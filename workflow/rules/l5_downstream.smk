# Layer 5: downstream analyses with two-phase eligibility (checkpoint), then dynamic expansion.

REQUESTED = PLAN.get("analyses_requested", [])
INTERNAL = PLAN.get("analyses_internal", [])
DOWN = R / "downstream"
ELIG = R / "eligibility"


def common_inputs(wc=None):
    base = {"matrix": str(R / "cohort" / "evidence_matrix.tsv.gz"), "sites": str(R / "cohort" / "sites.tsv"),
            "missing": str(R / "cohort" / "missingness.tsv"), "samples": str(R / "cohort" / "cohort_samples.txt")}
    if 4 in LAYERS:
        base["ann"] = str(R / "annotation" / "cohort.ann.vcf.gz")
    return base


rule da001_clonality:
    """Clonal-status evidence (DA001), run when requested or as an internal prerequisite."""
    input:
        unpack(common_inputs),
        gates=lambda wc: [str(R / "qc" / "samples" / s / "qc_decision.json") for s in sorted(SAMPLES) if s not in EXCLUDED],
    output:
        manifest=str(DOWN / "DA001" / "status_manifest.json"),
        clonality=str(DOWN / "DA001" / "clonality.tsv"),
    resources:
        mem_mb=mem(2000),
        runtime=RUNTIME,
    log:
        str(LOGS / "downstream" / "DA001.log"),
    benchmark:
        str(BENCH / "downstream" / "DA001.tsv")
    shell:
        "python -m pfvarx.downstream.run --da DA001 --run-dir {RUN:q} --root {ROOT:q} > {log:q} 2>&1"


def eligibility_inputs(wc):
    d = {"gates": [str(R / "qc" / "samples" / s / "qc_decision.json") for s in sorted(SAMPLES) if s not in EXCLUDED],
         "samples": str(R / "cohort" / "cohort_samples.txt")}
    if "DA001" in REQUESTED + INTERNAL:
        d["clonality"] = str(DOWN / "DA001" / "clonality.tsv")
    return d


checkpoint da_eligibility:
    input:
        unpack(eligibility_inputs),
    output:
        decisions=str(ELIG / "cohort_decisions.json"),
        matrix=str(ELIG / "sample_analysis_eligibility.tsv"),
    log:
        str(LOGS / "downstream" / "eligibility.log"),
    shell:
        "python -m pfvarx.downstream.run --eligibility --run-dir {RUN:q} --root {ROOT:q} > {log:q} 2>&1"


def source_manifest(wc):
    src = CFG["analyses"].get("spatial_outcome_source")
    return [str(DOWN / src / "status_manifest.json")] if wc.da == "DA007" and src else []


rule da_execute:
    input:
        unpack(common_inputs),
        decisions=str(ELIG / "cohort_decisions.json"),
        source=source_manifest,
    output:
        manifest=str(DOWN / "{da}" / "status_manifest.json"),
    wildcard_constraints:
        da="DA00[234567]",
    threads: threads(2)
    resources:
        mem_mb=mem(2000),
        runtime=RUNTIME,
    log:
        str(LOGS / "downstream" / "{da}.log"),
    benchmark:
        str(BENCH / "downstream" / "{da}.tsv")
    shell:
        "python -m pfvarx.downstream.run --da {wildcards.da} --run-dir {RUN:q} --root {ROOT:q} "
        "--threads {threads} > {log:q} 2>&1"


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
