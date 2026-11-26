# Helper functions for the PfVarX workflow, kept apart from the rules (Snakemake's
# modularisation guideline). Included by the Snakefile before any rule file; the functions
# read the globals defined there (R, P, CFG, PLAN, SAMPLES, LAYERS, COMPUTE, ...) and are
# called when rules are defined or when the DAG is built.

# Layer 5 locations and analysis lists, used by the eligibility and target functions.
REQUESTED = PLAN.get("analyses_requested", [])
INTERNAL = PLAN.get("analyses_internal", [])
DOWN = R / "downstream"
ELIG = R / "eligibility"


# --- Resources and the per-sample QC gates ---


PLANNED_MEM = {j["job"]: j["mem_mb"] for j in PLAN.get("resources", {}).get("jobs", [])}
# Per-rule thread and memory requests: config/registry/job_resources.yaml, the file the
# preflight resource plan also reads (pfvarx.planning.resources).
from pfvarx.planning import resources as _job_resources

RULE_REQUESTS = _job_resources.rule_requests(ROOT)
ZIPPER_HEAP = _job_resources.load(ROOT)["zipper_heap"]


def mem(rule, host=False, job=None):
    """Attempt-scaled memory of *rule* within the per-job ceiling (retries only on memory-scalable rules).

    ``host=True`` marks a bwa-mem2 job: with Human DNA removal (wgs only) it loads the
    combined parasite+host index, so the request grows by ``layer1.bwa_mem2_host_extra_mem_mb``.
    ``job`` names the job class whose request the frozen plan may have raised to the memory
    this computer measured in earlier runs, so Snakemake never starts more of them at once
    than ``max_mem_mb`` holds.
    """
    base = RULE_REQUESTS[rule]["mem_mb"]
    cap = COMPUTE["per_job"]["max_mem_mb"]
    factor = COMPUTE["retry_memory_factor"]
    if host and HOST:
        base += P["layer1.bwa_mem2_host_extra_mem_mb"]
    if job is not None:
        base = max(base, PLANNED_MEM.get(job, 0))
    return lambda wc, attempt: int(min(base * factor ** (attempt - 1), cap))


def host_filter(counts):
    """Shell pipe stage after bwa-mem2 that removes host read pairs (empty without a host FASTA).

    *counts* is a path pattern with {sample} and optionally {lane}; the stage is a params
    function so the wildcards are filled in.
    """
    def stage(wc):
        if not HOST:
            return ""
        path = counts.format(sample=wc.sample, lane=getattr(wc, "lane", ""))
        return f"| python -m pfvarx.references.host --counts '{path}' "
    return stage


def threads(rule):
    """Threads requested for *rule*, capped at compute.yaml per_job.max_threads."""
    return min(RULE_REQUESTS[rule]["threads"], COMPUTE["per_job"]["max_threads"])


def zipper_heap_mb(wc, resources):
    """Java heap of fgbio ZipperBams, which shares its job with bwa-mem2 (zipper_heap in job_resources.yaml)."""
    return max(ZIPPER_HEAP["heap_min_mb"], resources.mem_mb // ZIPPER_HEAP["heap_divisor"])


def gate(sample):
    """Layer 1 -> 2 gate decision for one sample (forces the per-sample checkpoint)."""
    out = checkpoints.l1_gate.get(sample=sample).output[0]
    return json.loads(Path(out).read_text())


def passing(stage):
    """Samples allowed into *stage* after their own QC gate (failure isolation)."""
    out = []
    for s in sorted(SAMPLES):
        if s in EXCLUDED:
            continue
        g = gate(s)
        if g["permitted"].get(stage, False):
            out.append(s)
    return out


def final_targets(wildcards):
    targets = []
    if 1 in LAYERS and PLAN["input_type"] == "FASTQ":
        targets += [str(R / "alignments" / f"{s}.bam") for s in passing("L1_OUTPUT")]
    if 2 in LAYERS:
        targets += [str(R / "variants" / "samples" / f"{s}.vcf.gz") for s in passing("L2")]
    if 3 not in LAYERS:
        return targets
    # The cohort status is always written. When no sample is permitted into Layer 3 (every
    # sample failed QC or input validation), the cohort, annotation and analyses are skipped
    # with an explicit outcome instead of running the merge on an empty file list.
    targets.append(str(R / "cohort" / "cohort_status.json"))
    if not passing("L3"):
        return targets
    targets += [str(R / "cohort" / "cohort.vcf.gz"), str(R / "cohort" / "missingness.tsv"),
                str(R / "cohort" / "cohort_qc.json")]
    if 4 in LAYERS:
        targets += [str(R / "annotation" / "cohort.ann.vcf.gz"), str(R / "annotation" / "annotation_qc.json")]
    if 5 in LAYERS:
        targets += da_targets()
    return targets


# --- Layer 1: what each library is (resolved at validation, frozen in the plan) ---


def library(sample):
    """Resolved library facts of one sample: technique, route, panel, protocol, state, count unit."""
    return SAMPLES[sample].get("assay") or {}


def umi(sample, name):
    """A UMI extraction or consensus setting: the sample's library protocol, else the layer-1 default."""
    setting = (library(sample).get("umi") or {}).get(name)
    if setting is not None and setting.get("value") is not None:
        return setting["value"]
    if name == "read_structure":  # runs frozen before library protocols (config schema 3)
        return CFG.get("assay", {}).get("read_structure")
    return P["layer1.fgbio_" + name]


def count_unit(sample):
    return library(sample).get("count_unit") or ("consensus_molecules" if DESIGN == "targeted_umi" else "reads")


def entry_steps(sample):
    return library(sample).get("entry_steps") or []


# --- Layer 1: FASTQ QC (l1_fastq_qc.smk) ---


def lane(wc):
    return next(x for x in SAMPLES[wc.sample]["lanes"] if x["lane"] == wc.lane)


def adapter_args(sample):
    """fastp adapter options for the sample's library kit (resolved at validation).

    Known kits are passed explicitly: fastp's detection can choose its poly-T FlowCell1
    sequence on AT-rich reads and trim real T runs. Paired-end overlap trimming stays on.
    """
    adapters = library(sample).get("adapters") or {}
    if adapters.get("kit", "AUTO") == "AUTO" or not adapters.get("read1"):
        return ["--detect_adapter_for_pe"]
    return [f"--adapter_sequence {adapters['read1']}", f"--adapter_sequence_r2 {adapters['read2']}"]


def fastp_args(sample):
    a = []
    if P["layer1.fastp_low_complexity_filter"]:
        a.append("--low_complexity_filter")
    if P["layer1.fastp_dedup"]:
        a.append("--dedup")
    a += adapter_args(sample)
    a += [f"--qualified_quality_phred {P['layer1.fastp_qualified_quality_phred']}",
          f"--unqualified_percent_limit {P['layer1.fastp_unqualified_percent_limit']}",
          f"--n_base_limit {P['layer1.fastp_n_base_limit']}",
          f"--length_required {P['layer1.fastp_length_required']}"]
    if P["layer1.fastp_cut_right"]:
        a.append(f"--cut_right --cut_window_size {P['layer1.fastp_cut_window_size']} "
                 f"--cut_mean_quality {P['layer1.fastp_cut_mean_quality']}")
    if CFG["assay"]["two_colour_chemistry"]:
        a.append("--trim_poly_g")
    else:
        a.append("--disable_trim_poly_g")
    if DESIGN == "targeted_umi":
        # UMIs are inline at the 5' end: no quality trimming or filtering before extraction.
        # Pairs too short to hold the fixed read-structure segments (adapter dimers) are removed
        # here, because fgbio FastqToBam stops on the first such read. The read structure is the
        # sample's library protocol.
        from pfvarx.inputs.validate import runtime_read_structures

        _, min_len, trim_safe = runtime_read_structures(umi(sample, "read_structure"))
        a = [x for x in a if not x.startswith(("--cut_right", "--length_required"))]
        a += ["--disable_quality_filtering", f"--length_required {min_len}"]
        if not trim_safe:
            # A fixed non-template segment ends the read: 3' trimming would corrupt it.
            a = [x for x in a if x not in ("--detect_adapter_for_pe", "--trim_poly_g", "--disable_trim_poly_g")
                 and not x.startswith("--adapter_sequence")]
            a += ["--disable_adapter_trimming", "--disable_trim_poly_g"]
    return " ".join(a)


# --- Layer 1: alignment and imported alignments (l1_align.smk, l1_umi.smk) ---


def rg(wc):
    s = SAMPLES[wc.sample]
    ln = next(x for x in s["lanes"] if x["lane"] == wc.lane)
    return (f"@RG\\tID:{wc.sample}.{wc.lane}\\tSM:{wc.sample}\\tLB:{ln['library_id']}"
            f"\\tPL:{CFG['assay']['platform']}\\tPU:{wc.sample}.{wc.lane}")


def umi_group_input(wc):
    """Reads to group by UMI: the aligned lanes (FASTQ entry) or the imported raw-UMI BAM."""
    if PLAN["input_type"] == "FASTQ":
        return [str(R / "work" / "umi" / "lanes" / f"{wc.sample}.{x['lane']}.mapped.bam")
                for x in SAMPLES[wc.sample]["lanes"]]
    return [str(R / "work" / "entry" / f"{wc.sample}.bam")]


def entry_trim_source(wc):
    """Imported alignment to primer-trim: new consensus reads (raw UMI import) or the import itself."""
    if "umi_consensus" in entry_steps(wc.sample):
        return str(R / "work" / "umi" / f"{wc.sample}.consensus.mapped.bam")
    return str(R / "work" / "entry" / f"{wc.sample}.bam")


def entry_final_source(wc):
    """Last product for an imported alignment, given the steps its documented state still needs."""
    steps = entry_steps(wc.sample)
    if "trim_primers" in steps:
        return str(R / "work" / "entry" / f"{wc.sample}.trimmed.bam")
    if "umi_consensus" in steps:
        return str(R / "work" / "umi" / f"{wc.sample}.consensus.mapped.bam")
    return str(R / "work" / "entry" / f"{wc.sample}.bam")


# --- Layer 1: QC metrics and gate (l1_metrics.smk) ---


def gate_inputs(wc):
    d = R / "qc" / "samples" / wc.sample
    files = {"flagstat": str(d / f"{wc.sample}.flagstat.json")}
    if HOST and PLAN["input_type"] == "FASTQ":
        files["host"] = [str(d / f"{wc.sample}.{x['lane']}.host.json") for x in SAMPLES[wc.sample]["lanes"]]
    if PLAN["input_type"] == "FASTQ":
        files["fastp"] = [str(d / f"{wc.sample}.{x['lane']}.fastp.json") for x in SAMPLES[wc.sample]["lanes"]]
    if DESIGN == "wgs":
        files["mosdepth"] = str(d / "core.mosdepth.region.dist.txt")
        if PLAN["input_type"] == "FASTQ":
            files["markdup"] = str(d / f"{wc.sample}.markdup.txt")
    else:
        files["targets"] = str(d / "targets.regions.bed.gz")
        files["ontarget"] = str(d / f"{wc.sample}.ontarget.txt")
    # UMI family sizes exist only where this run builds the consensus.
    if DESIGN == "targeted_umi" and (PLAN["input_type"] == "FASTQ" or "umi_consensus" in entry_steps(wc.sample)):
        files["umi"] = str(d / f"{wc.sample}.umi_family_sizes.txt")
    return files


# --- Layer 3: cohort (l3_cohort.smk) ---


def cohort_status_inputs(wc):
    if passing("L3"):
        return {"samples": str(R / "cohort" / "cohort_samples.txt")}
    return {}


# --- Layer 5: downstream analyses (l5_downstream.smk) ---


def common_inputs(wc=None):
    base = {"matrix": str(R / "cohort" / "evidence_matrix.tsv.gz"), "sites": str(R / "cohort" / "sites.tsv"),
            "missing": str(R / "cohort" / "missingness.tsv"), "samples": str(R / "cohort" / "cohort_samples.txt")}
    if 4 in LAYERS:
        base["ann"] = str(R / "annotation" / "cohort.ann.vcf.gz")
    return base


def eligibility_inputs(wc):
    d = {"gates": [str(R / "qc" / "samples" / s / "qc_decision.json") for s in sorted(SAMPLES) if s not in EXCLUDED],
         "samples": str(R / "cohort" / "cohort_samples.txt"),
         "cohort_qc": str(R / "cohort" / "cohort_qc.json")}
    if "DA001" in REQUESTED + INTERNAL:
        d["clonality"] = str(DOWN / "DA001" / "clonality.tsv")
    return d


def source_manifest(wc):
    src = CFG["analyses"].get("spatial_outcome_source")
    return [str(DOWN / src / "status_manifest.json")] if wc.da == "DA007" and src else []


def da_targets():
    """Phase two: expand only analyses whose cohort decision is EXECUTE."""
    out = [str(ELIG / "cohort_decisions.json")]
    if "DA001" in REQUESTED + INTERNAL:
        out.append(str(DOWN / "DA001" / "status_manifest.json"))
    decisions = json.loads(Path(checkpoints.da_eligibility.get().output.decisions).read_text())
    for da, d in decisions["analyses"].items():
        if da != "DA001" and d["decision"] == "EXECUTE":
            out.append(str(DOWN / da / "status_manifest.json"))
    return out

