"""Layer 6 harvester and renderer: results -> dashboard_data.json -> one offline HTML file.

The harvester only reads upstream records; every warning, threshold, reason code and
interpretation statement shown in the dashboard comes from those records or from the
fixed interpretation registry below. Missing, failed, skipped, not-requested and
not-applicable values are distinct states and are never rendered as zero.
"""

from __future__ import annotations

import csv
import datetime as dt
import gzip
import hashlib
import json
import math
import re
from collections import Counter
from importlib import resources
from pathlib import Path
from typing import Any

from pfvarx import __version__
from pfvarx.diagnostics.events import read_events
from pfvarx.io import read_json, read_tsv, write_json

DA_NAMES = {"DA002": "Drug-resistance markers", "DA003": "Genetic diversity", "DA001": "COI and clonal status",
            "DA004": "Relatedness and IBD", "DA005": "Population structure", "DA006": "Phylogenetics",
            "DA007": "Spatial-temporal", "DA008": "Copy-number markers",
            "DA009": "pfhrp2/3 joint deletion"}

# Fixed, analysis-specific interpretation consequences of QC warnings (Layer 6 spec, section 6).
CONSEQUENCES = {
    "DA002": "Low depth or callability reduces the number of callable codons; affected markers are UNKNOWN, "
             "never wild type. Mixed-allele fractions near the reporting floor are least reliable.",
    "DA003": "Missingness and low depth reduce the sites shared by all samples and can bias diversity downward.",
    "DA001": "Few informative loci or low minor-clone sensitivity can make a polyclonal infection look monoclonal; "
             "an uncertain classification is reported as INDETERMINATE.",
    "DA004": "Marker density, genotype uncertainty and infection complexity widen or bias relatedness estimates.",
    "DA005": "Missingness and the dominant-allele representation of polyclonal samples can shift PCA positions and FST.",
    "DA006": "Only MONOCLONAL samples are eligible; missing data and few informative sites weaken support values.",
    "DA007": "Missing locations or dates remove samples; spatial sampling bias is not modelled.",
    "DA008": "Low or uneven depth widens the copy-number interval; a locus that cannot be judged is UNKNOWN, "
             "never single copy or present. pfhrp2/3 lie in subtelomeric repeats where depth is least reliable.",
    "DA009": "A gene that cannot be judged leaves the sample INDETERMINATE, never 'no deletion'. On a targeted "
             "panel the call is amplicon presence, not copy number: a mutated primer site or a failed PCR looks "
             "like a deletion, which the control amplicons reduce but do not remove.",
}

VIEW_REGISTRY = {
    "DA002.prevalence": {"estimand": "marker prevalence over infections (callable denominator)",
                         "views": ["forest", "bar", "heatmap", "table", "map"], "uncertainty": "Wilson interval (simple random sampling)",
                         "fields": ["marker_id", "group", "callable", "mutant_or_mixed", "prevalence_infections"],
                         "limitations": "molecular association, not treatment outcome"},
    "DA002.calls": {"estimand": "per-sample marker state", "views": ["table", "heatmap"],
                    "uncertainty": "codon read counts and Wilson interval of the alt fraction",
                    "limitations": "UNKNOWN = uncallable; never wild type"},
    "DA001.clonality": {"estimand": "clonal status and COI posterior (MOIRE)", "views": ["coi", "histogram", "box", "stacked", "table"],
                        "uncertainty": "posterior interval of COI; P(COI>1)", "limitations": "MOIRE settings provisional"},
    "DA003.diversity": {"estimand": "expected heterozygosity and nucleotide diversity per callable base",
                        "views": ["comparison", "histogram", "box", "table"], "uncertainty": "bootstrap over samples",
                        "limitations": "dominant-allele representation"},
    "DA004.pairs": {"estimand": "pairwise IBD / relatedness", "views": ["table", "heatmap", "network"],
                    "uncertainty": "Dcifer CI; hmmIBD point estimate", "limitations": "network threshold is a display setting"},
    "DA005.pca": {"estimand": "PCA coordinates of analysed samples", "views": ["scatter", "components", "table", "fst"],
                  "uncertainty": "none (descriptive)", "limitations": "filtering highlights; it does not refit PCA"},
    "DA006.tree": {"estimand": "maximum-likelihood tree of MONOCLONAL samples", "views": ["tree", "distance", "table"],
                   "uncertainty": "UFBoot / SH-aLRT support", "limitations": "no tree when skipped"},
    "DA008.calls": {"estimand": "per-sample copy-number state of each locus (depth ratio to the core genome)",
                    "views": ["table", "bar"], "uncertainty": "bootstrap interval of the depth ratio",
                    "limitations": "DNA-level estimate; intermediate ratios in mixed infections; UNKNOWN never counted"},
    "DA008.prevalence": {"estimand": "share of judged infections with each amplification or deletion",
                         "views": ["bar", "table"], "uncertainty": "Wilson interval (simple random sampling)",
                         "limitations": "provisional thresholds"},
    "DA009.joint": {"estimand": "joint pfhrp2/pfhrp3 state of each infection (single, dual or no deletion)",
                    "views": ["bar", "table"], "uncertainty": "none (per-sample states)",
                    "limitations": "INDETERMINATE whenever either gene is UNKNOWN; evidence class differs by design"},
    "DA009.prevalence": {"estimand": "share of judged infections with pfhrp2, pfhrp3, dual or any deletion",
                         "views": ["bar", "table"], "uncertainty": "Wilson interval (simple random sampling)",
                         "limitations": "pooled across evidence classes; the stratified table is published beside it"},
    "DA009.cooccurrence": {"estimand": "whether the two deletions co-occur more than independence predicts",
                           "views": ["table"], "uncertainty": "Fisher exact test p-value",
                           "limitations": "describes this cohort; not evidence of selection"},
    "DA007.summary": {"estimand": "prevalence by country and month", "views": ["map", "heatmap", "line", "table"],
                      "uncertainty": "counts shown; GAM CI where fitted", "limitations": "aggregated; preferential sampling"},
}


def _read_tsv_safe(path: Path) -> list[dict[str, str]]:
    try:
        return read_tsv(path) if path.is_file() else []
    except ValueError:
        return []


def _benchmarks(run: Path) -> list[dict[str, Any]]:
    out = []
    for f in sorted((run / "benchmarks").rglob("*.tsv")):
        rows = list(csv.DictReader(f.open(), delimiter="\t"))
        for r in rows:
            out.append({"rule": f.parent.name, "job": f.stem, "wall_s": _num(r.get("s")),
                        "max_rss_mb": _num(r.get("max_rss")), "cpu_time_s": _num(r.get("cpu_time")),
                        "io_in_mb": _num(r.get("io_in")), "io_out_mb": _num(r.get("io_out")),
                        "source": "Snakemake benchmark directive (measured)"})
    return out


# Benchmark folder -> (step shown in the report, pipeline layer, packages that do the work).
def _timeline(run: Path, status: dict[str, Any]) -> dict[str, Any]:
    """When the run started and ended, from its id and its final status."""
    try:
        started = dt.datetime.strptime(run.name[:16], "%Y%m%dT%H%M%SZ").replace(tzinfo=dt.UTC)
    except ValueError:
        return {"started_utc": None, "finished_utc": None, "elapsed_min": None}
    ended = status.get("updated_utc") if status.get("status") not in (None, "RUNNING", "PENDING") else None
    elapsed = (dt.datetime.fromisoformat(ended) - started).total_seconds() / 60 if ended else None
    return {"started_utc": started.isoformat(), "finished_utc": ended,
            "elapsed_min": None if elapsed is None else round(elapsed, 1)}


def _planned_vs_actual(run: Path, est: dict[str, Any] | None, measured: list[dict[str, Any]]) -> dict[str, Any] | None:
    """The estimate frozen with the plan beside what the run measured: time, peak RAM and storage."""
    if not est:
        return None
    from pfvarx.planning.history import RULE_TO_JOB, STORAGE_FOLDERS, _folder_mb

    actual: dict[str, dict[str, float]] = {}
    for b in measured:
        job = RULE_TO_JOB.get(b["rule"])
        if job is None or b.get("wall_s") is None:
            continue
        a = actual.setdefault(job, {"wall_s": 0.0, "rss_mb": 0.0})
        a["wall_s"] += b["wall_s"]
        a["rss_mb"] = max(a["rss_mb"], b.get("max_rss_mb") or 0.0)
    steps = [{"step": s["name"], "jobs": s["jobs"], "expected_min": s["minutes"],
              "actual_min": round(actual[s["job"]]["wall_s"] / 60, 1) if s["job"] in actual else None,
              "expected_peak_ram_mb": s["peak_ram_mb"],
              "actual_peak_ram_mb": round(actual[s["job"]]["rss_mb"]) if s["job"] in actual else None,
              "basis": s["basis"]} for s in est.get("steps", [])]
    planned = est.get("storage_mb") or {}
    storage: list[dict[str, Any]] = [{"folder": f, "expected_mb": planned.get(f),
                                      "actual_mb": round(_folder_mb(run / "results" / f), 1)} for f in STORAGE_FOLDERS
                                     if (run / "results" / f).is_dir() or planned.get(f)]
    actual_total = round(sum(float(r["actual_mb"]) for r in storage), 1)
    return {"steps": steps, "storage": storage, "expected_wall_min": est.get("wall_minutes"),
            "expected_peak_ram_mb": est.get("peak_ram_mb"), "expected_storage_mb": est.get("storage_total_mb"),
            "actual_storage_mb": actual_total,
            "basis": est.get("basis"), "confidence": est.get("confidence")}


STEP_INFO: dict[str, tuple[str, str, tuple[str, ...]]] = {
    "fastp": ("Read QC and trimming", "L1", ("fastp",)),
    "bwa_mem2": ("Read alignment", "L1", ("bwa-mem2", "samtools")),
    "merge_lanes": ("Merge lanes", "L1", ("samtools",)),
    "markdup": ("Mark duplicates", "L1", ("picard", "samtools")),
    "insert_metrics_preclip": ("Insert metrics before primer clipping", "L1", ("samtools",)),
    "trim_primers": ("Primer trimming", "L1", ("fgbio", "samtools")),
    "alignment_entry": ("Aligned-input entry", "L1", ("samtools",)),
    "entry_finalize": ("Aligned-input finalise", "L1", ("samtools",)),
    "umi_fastq_to_bam": ("UMI: reads to BAM", "L1", ("fgbio",)),
    "umi_align_zipper": ("UMI: align and merge tags", "L1", ("bwa-mem2", "fgbio", "samtools")),
    "umi_group": ("UMI: group reads", "L1", ("fgbio",)),
    "umi_consensus": ("UMI: consensus reads", "L1", ("fgbio",)),
    "umi_realign": ("UMI: realign consensus", "L1", ("bwa-mem2", "fgbio", "samtools")),
    "umi_preclip_and_trim": ("UMI: clip and trim primers", "L1", ("fgbio", "samtools")),
    "flagstat": ("Alignment statistics", "L1", ("samtools",)),
    "mosdepth_core": ("Core-genome depth", "L1", ("mosdepth",)),
    "mosdepth": ("Target depth", "L1", ("mosdepth",)),
    "l1_gate": ("QC gate decision", "L1", ("python",)),
    "pfvarx_call": ("Variant calling (PfVarX caller)", "L2", ("python", "pysam", "numpy", "scipy")),
    "normalise_sample": ("Normalise sample VCF", "L2", ("bcftools",)),
    "cohort_merge": ("Cohort merge", "L3", ("bcftools",)),
    "cohort_status": ("Cohort call status", "L3", ("python", "pysam")),
    "cohort_evidence": ("Cohort evidence", "L3", ("python", "pysam")),
    "snpeff_build": ("Annotation database build", "L4", ("snpeff",)),
    "snpeff_annotate": ("Variant annotation", "L4", ("snpeff",)),
    "da_eligibility": ("Analysis eligibility", "L5", ("python",)),
}
ANALYSIS_TOOLS: dict[str, tuple[str, tuple[str, ...]]] = {
    "DA001": ("Complexity of infection", ("r-base", "r-moire")),
    "DA002": ("Drug-resistance markers", ("python",)),
    "DA003": ("Genetic diversity", ("scikit-allel",)),
    "DA004": ("Relatedness and IBD", ("hmmibd", "r-dcifer")),
    "DA005": ("Population structure", ("scikit-allel",)),
    "DA006": ("Phylogenetics", ("iqtree",)),
    "DA007": ("Spatial and temporal summary", ("r-base", "r-mgcv")),
    "DA008": ("Copy-number markers", ("python", "pysam")),
    "DA009": ("pfhrp2/pfhrp3 joint deletion", ("python", "pysam")),
}


def _label_steps(measured: list[dict[str, Any]], tools_installed: dict[str, Any]) -> None:
    """Add the step name, layer and tools (with the run's recorded versions) to each measurement."""
    def with_version(package: str) -> str:
        if package == "python":
            return "Python (PfVarX)"
        recorded = next((tools_installed[f"{name} (package)"] for name in (package, f"{package}-slim", f"{package}-minimal")
                         if f"{name} (package)" in tools_installed), None)
        return f"{package} {str(recorded).split(' ')[0]}" if recorded else package

    for row in measured:
        step, layer, packages = STEP_INFO.get(row["rule"], (row["rule"].replace("_", " "), None, ()))
        if row["rule"] == "downstream":
            name, packages = ANALYSIS_TOOLS.get(row["job"], (row["job"], ()))
            step, layer = f"{row['job']} {name}", "L5"
        row.update({"step": step, "layer": layer, "tools": ", ".join(with_version(p) for p in packages) or None})


def _num(v: Any) -> float | None:
    try:
        value = None if v in (None, "", "-", "NA") else float(v)
        return value if value is not None and math.isfinite(value) else None
    except (ValueError, TypeError):
        return None


def _redact(value: Any) -> Any:
    """Redact embedded report/log text without loading environment variables."""
    if isinstance(value, dict):
        return {k: ("<REDACTED>" if re.search(r"(?i)(password|passwd|secret|token|api[_-]?key|credential|authorization|access[_-]?key)", str(k))
                    else _redact(v)) for k, v in value.items()}
    if isinstance(value, list):
        return [_redact(v) for v in value]
    if not isinstance(value, str):
        return value
    value = re.sub(r"(?is)-----BEGIN [^-]*PRIVATE KEY-----.*?(?:-----END [^-]*PRIVATE KEY-----|$)", "<REDACTED PRIVATE KEY>", value)
    if "\n" in value:
        value = re.sub(r"(?im)^\s*[A-Za-z0-9+/=]{64,}\s*$", "<REDACTED ENCODED BLOCK>", value)
    value = re.sub(r"(?i)\bBearer\s+[A-Za-z0-9._~+/-]+=*", "Bearer <REDACTED>", value)
    value = re.sub(r"(?i)((?:password|passwd|secret|token|api[_-]?key|credential|authorization|access[_-]?key)[\w-]*[\"']?\s*(?:[=:]\s*|\s+))(?:\"[^\"]*\"|'[^']*'|[^\s,;]+)", r"\1<REDACTED>", value)
    return re.sub(r"(https?://)[^\s/@:]+:[^\s/@]+@", r"\1<REDACTED>@", value)


def _log_content(run: Path) -> list[dict[str, Any]]:
    out = []
    remaining = 512 * 1024
    paths = sorted((p for p in (run / "logs").rglob("*") if p.is_file() and p.suffix in {".log", ".txt", ".jsonl"}),
                   key=lambda p: p.as_posix())
    for path in paths[:60]:
        if not path.resolve().is_relative_to(run.resolve()) or remaining <= 0:
            continue
        try:
            size = path.stat().st_size
            limit = min(64 * 1024, remaining)
            offset = max(0, size - limit)
            with path.open("rb") as handle:
                handle.seek(offset)
                content = handle.read(limit)
            if offset:
                content = content.split(b"\n", 1)[-1]  # avoid embedding a partial credential line
            text = content.decode("utf-8", errors="replace")
            out.append({"path": str(path.relative_to(run)), "text": _redact(text),
                        "size_bytes": size, "truncated": offset > 0,
                        "scope": "redacted tail" if offset else "redacted full content"})
            remaining -= len(content)
        except OSError:
            continue
    return out


def _optional_json(path: Path, warnings: list[str], default: Any = None) -> Any:
    if not path.is_file():
        return default
    try:
        value = read_json(path)
        expected = list if isinstance(default, list) else dict
        if not isinstance(value, expected):
            warnings.append(f"Invalid optional report source {path.name}: expected {expected.__name__}")
            return default
        return value
    except (OSError, ValueError) as exc:
        warnings.append(f"Cannot read optional report source {path.name}: {type(exc).__name__}")
        return default


def _distance_matrix(path: Path, warnings: list[str]) -> list[dict[str, Any]]:
    if not path.is_file():
        return []
    try:
        lines = path.read_text().strip().splitlines()
        n = int(lines[0].strip())
        rows = [line.split() for line in lines[1:]]
        if len(rows) != n or any(len(row) != n + 1 for row in rows):
            raise ValueError("invalid IQ-TREE distance matrix dimensions")
        return [{"sample_1": row[0], "sample_2": rows[j][0], "distance": _num(value)}
                for row in rows for j, value in enumerate(row[1:])]
    except (OSError, ValueError, IndexError):
        warnings.append("Cannot read optional IQ-TREE tree.mldist distance matrix")
        return []


def _vcf_values(path: Path, warnings: list[str]) -> dict[tuple[str, int, str, str], dict[str, Any]]:
    """Read optional published fields from plain, gzip or BGZF VCF without re-calling."""
    out: dict[tuple[str, int, str, str], dict[str, Any]] = {}
    if not path.is_file():
        return out
    try:
        opener = gzip.open if path.suffix == ".gz" else open
        with opener(path, "rt") as handle:
            for line in handle:
                if line.startswith("#") or not line.strip():
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 8:
                    warnings.append(f"Malformed optional VCF row in {path.name}")
                    continue
                info = dict(value.split("=", 1) for value in fields[7].split(";") if "=" in value)
                alts = fields[4].split(",")
                call = dict(zip(fields[8].split(":"), fields[9].split(":"))) if len(fields) > 9 else {}
                gt = call.get("GT")
                genotype = gt if gt and any(part != "." for part in re.split(r"[|/]", gt)) else None
                for alt in alts:
                    entries = info.get("ANN", "").split(",") if info.get("ANN") else []
                    matched = [v for v in entries if len(alts) == 1 or v.split("|")[0] == alt]
                    consequences = sorted({v.split("|")[1] for v in matched if len(v.split("|")) > 1})
                    out[(fields[0], int(fields[1]), fields[3], alt)] = {
                        "quality": _num(fields[5]), "genotype": genotype,
                        "annotation": "; ".join(matched) or None,
                        "consequence": " & ".join(consequences) or None}
    except (OSError, ValueError, EOFError) as exc:
        warnings.append(f"Cannot read optional VCF fields from {path.name}: {type(exc).__name__}")
    return out


def _data_origin(root: Path, eff: dict[str, Any]) -> tuple[str, str]:
    catalog = eff["config"]["analyses"].get("resistance_catalog")
    if catalog:
        p = root / catalog
        if p.is_file() and "#status=TEST_FIXTURE" in p.read_text()[:4096]:
            return "synthetic", "Configured catalogue is explicitly TEST_FIXTURE"
    if str(eff["config"]["reference"].get("release_label", "")).upper().startswith("SYNTHETIC"):
        return "synthetic", "Reference release is explicitly marked SYNTHETIC"
    ref = root / eff["config"]["reference"]["fasta"]
    gff = eff["config"]["reference"].get("gff")
    source = root / gff if gff else ref.with_suffix(".gff")
    if source.is_file():
        with source.open() as handle:
            for _ in range(100):
                line = handle.readline()
                if not line:
                    break
                fields = line.split("\t")
                if len(fields) >= 3 and fields[1] == "synthetic":
                    return "synthetic", "Reference annotation explicitly identifies its source as synthetic"
    return "not_declared", "No explicit synthetic declaration detected; data origin is not inferred from sample names"


def harvest(root: Path, run: Path) -> dict[str, Any]:
    res = run / "results"
    warnings: list[str] = []
    eff = read_json(run / "frozen_config.json")
    plan = read_json(run / "frozen_plan.json")
    status = read_json(run / "status.json")
    snap = read_json(run / "snapshot.json")
    meta = {r["sample_id"]: r for r in _read_tsv_safe(root / eff["config"]["samples_tsv"])}
    samples, qc_matrix = [], []
    cohort_qc = (_optional_json(res / "cohort" / "cohort_qc.json", warnings) or {}).get("samples", {})
    all_ids = sorted({s["sample_id"] for s in plan["samples"]} | set(plan.get("samples_failed_validation", [])))
    excl = read_json(run / "runtime_exclusions.json")["samples"] if (run / "runtime_exclusions.json").is_file() else []
    for sid in all_ids:
        gate_p = res / "qc" / "samples" / sid / "qc_decision.json"
        gate = _optional_json(gate_p, warnings)
        if gate is not None and "verdict" not in gate:
            warnings.append(f"Invalid optional QC decision for {sid}: missing verdict")
            gate = None
        cq_p = res / "qc" / "samples" / sid / f"{sid}.caller_qc.json"
        cq = _optional_json(cq_p, warnings)
        if sid in plan.get("samples_failed_validation", []):
            processing = "BLOCKED_INPUT_VALIDATION"
        elif sid in excl:
            processing = "FAILED"
        elif gate is None:
            processing = "NOT_REACHED"
        elif gate["verdict"] == "FAIL":
            processing = "BLOCKED_QC"
        elif (res / "alignments" / f"{sid}.bam").is_file() and (
                plan["final_output"] == "FO001" or (res / "variants" / "samples" / f"{sid}.vcf.gz").is_file()):
            processing = "COMPLETED"
        elif status.get("status") in ("RUNNING", "PENDING", "READY"):
            processing = status["status"]
        else:
            processing = "INCOMPLETE"
        m = meta.get(sid, {})
        samples.append({"sample_id": sid, "processing": processing,
                        "qc_verdict": gate["verdict"] if gate else "NOT_EVALUATED",
                        "cohort_qc_verdict": (cohort_qc.get(sid) or {}).get("verdict"),
                        "label": gate.get("label") if gate else None,
                        "permitted": gate.get("permitted") if gate else {},
                        "metrics": gate.get("metrics") if gate else {},
                        "uncalibrated_gates": gate.get("uncalibrated_gates", []) if gate else [],
                        "caller": None if cq is None else {
                            "median_depth": cq.get("median_usable_depth_core"), "callable_fraction": cq.get("callable_fraction_core_caller"),
                            "status_counts": cq.get("status_counts"), "calibration": cq.get("calibration"),
                            "mixture_status": (cq.get("mixture") or {}).get("status"), "technical_k": (cq.get("mixture") or {}).get("k_selected"),
                            "technical_weights": ((cq.get("mixture") or {}).get("fit") or {}).get("weights"),
                            "error_model": (cq.get("error_model") or {}).get("snv_any")},
                        "population": m.get(eff["config"]["analyses"].get("population_column") or "population", ""),
                        "country": m.get("country", ""), "latitude": _num(m.get("latitude")),
                        "longitude": _num(m.get("longitude")), "collection_date": m.get("collection_date", "")})
        for d in [*(gate or {}).get("decisions", []), *(cohort_qc.get(sid) or {}).get("decisions", [])]:
            qc_matrix.append({"sample_id": sid, **{k: d.get(k) for k in (
                "gate_id", "gate", "stage", "metric", "unit", "observed", "fail_threshold", "warn_threshold",
                "direction", "status", "threshold_source", "calibrated", "policy_version", "mandatory", "reason")}})
    variants = []
    variant_files = []
    annotations = _vcf_values(res / "annotation" / "cohort.ann.vcf.gz", warnings)
    for sid in all_ids:
        ev = res / "variants" / "samples" / f"{sid}.evidence.tsv"
        if ev.is_file():
            variant_files.append(str(ev.relative_to(run)))
        vcf_values = _vcf_values(res / "variants" / "samples" / f"{sid}.vcf.gz", warnings)
        for r in _read_tsv_safe(ev):
            variant_key = (r["chrom"], int(r["pos"]), r["ref"], r["alt"])
            extra = vcf_values.get(variant_key, {})
            annotation = annotations.get(variant_key, {})
            variants.append({"sample_id": sid, "chrom": r["chrom"], "pos": int(r["pos"]), "ref": r["ref"], "alt": r["alt"],
                             "type": r["kind"], "status": r["status"], "wsaf": _num(r["wsaf"]), "wsaf_lo": _num(r["wsaf_lo"]),
                             "wsaf_hi": _num(r["wsaf_hi"]), "depth": int(r["n"]), "alt_reads": int(r["y"]),
                             "screen": r["screen"] if r["screen"] in ("inf", "-inf") else _num(r["screen"]),
                             "quality": extra.get("quality"), "genotype": extra.get("genotype"),
                             "annotation": annotation.get("annotation") or extra.get("annotation"),
                             "consequence": annotation.get("consequence") or extra.get("consequence"),
                             "source": str(ev.relative_to(run))})
    cap = int(eff["parameters"]["layer6.max_embedded_records"]["value"])
    counts = {"types": dict(Counter(v["type"] for v in variants)),
              "statuses": dict(Counter(v["status"] for v in variants)),
              "chromosomes": dict(Counter(v["chrom"] for v in variants)),
              "scope": "all evidence records before embedding cap", "sources": variant_files}
    variants_total = len(variants) if variant_files else None
    if variants_total is not None and variants_total > cap:  # bounded embedding: the omitted count is shown, never hidden
        variants.sort(key=lambda v: (v["status"] != "PASS", v["sample_id"], v["chrom"], v["pos"]))
        variants = variants[:cap]
    ann = _optional_json(res / "annotation" / "annotation_qc.json", warnings)
    downstream: dict[str, Any] = {}
    requested = plan.get("analyses_requested", [])
    internal = plan.get("analyses_internal", [])
    for da in sorted(set(requested) | set(internal)):
        d = res / "downstream" / da
        man = _optional_json(d / "status_manifest.json", warnings)
        entry: dict[str, Any] = {"name": DA_NAMES[da], "requested": da in requested, "internal": da in internal,
                                 "manifest": man, "terminal": man.get("terminal_state", "NOT_RECORDED") if man else "NOT_RUN_OR_FAILED",
                                 "consequence": CONSEQUENCES[da]}
        if entry["terminal"] not in {"COMPLETED", "COMPLETED_WITH_EXCLUSIONS"}:
            downstream[da] = entry
            continue
        if da == "DA002":
            entry["calls"] = _read_tsv_safe(d / "marker_calls.tsv")
            entry["prevalence"] = _read_tsv_safe(d / "prevalence.tsv")
        elif da == "DA003":
            entry["fws"] = _read_tsv_safe(d / "within_host_fws.tsv")
            entry["diversity"] = _read_tsv_safe(d / "population_diversity.tsv")
        elif da == "DA001":
            entry["clonality"] = _read_tsv_safe(d / "clonality.tsv")
        elif da == "DA004":
            entry["pairs"] = _read_tsv_safe(d / "pairwise_relatedness.tsv")
        elif da == "DA005":
            entry["pca"] = _read_tsv_safe(d / "pca.tsv")
            entry["fst"] = _read_tsv_safe(d / "fst.tsv")
            entry["pca_model"] = _optional_json(d / "pca_model.json", warnings)
        elif da == "DA006":
            tree = d / "tree.treefile"
            entry["newick"] = tree.read_text().strip() if tree.is_file() else None
            entry["gates"] = _optional_json(d / "runtime_gates.json", warnings)
            entry["distance_matrix"] = _distance_matrix(d / "tree.mldist", warnings)
        elif da == "DA008":
            entry["calls"] = _read_tsv_safe(d / "cnv_calls.tsv")
            entry["prevalence"] = _read_tsv_safe(d / "prevalence.tsv")
        elif da == "DA009":
            entry["gene_calls"] = _read_tsv_safe(d / "gene_calls.tsv")
            entry["joint"] = _read_tsv_safe(d / "joint_calls.tsv")
            entry["prevalence"] = _read_tsv_safe(d / "joint_prevalence.tsv")
            entry["prevalence_by_evidence"] = _read_tsv_safe(d / "joint_prevalence_by_evidence.tsv")
            entry["cooccurrence"] = _read_tsv_safe(d / "cooccurrence.tsv")
            entry["trend"] = _read_tsv_safe(d / "trend.tsv")
        elif da == "DA007":
            entry["summary"] = _read_tsv_safe(d / "country_month_summary.tsv")
            entry["grid"] = _optional_json(d / "grid_summary.json", warnings, [])
        downstream[da] = entry
    elig = _optional_json(res / "eligibility" / "sample_analysis_eligibility.json", warnings, [])
    events = read_events(run / "logs" / "events.jsonl")
    # Derived, structured diagnostics from upstream records (QC fails/warnings, exclusions, skipped analyses).
    for s in samples:
        for d in [q for q in qc_matrix if q["sample_id"] == s["sample_id"] and q["status"] in ("FAIL", "WARNING")]:
            events.append({"event_id": f"qc-{s['sample_id']}-{d['gate_id']}", "run_id": run.name, "sample_id": s["sample_id"],
                           "analysis_id": None, "layer": "L1", "stage": d["stage"], "status": d["status"],
                           "severity": "ERROR" if d["status"] == "FAIL" else "WARNING", "code": d["reason"],
                           "failure_class": "QC_FAILURE", "observed": d["observed"],
                           "expected": f"{d['direction']} fail={d['fail_threshold']} warn={d['warn_threshold']}",
                           "criterion": d["gate_id"], "cause_confirmed": None,
                           "causes_suspected": ["not established by the reported metric"],
                           "remedies": ["inspect the metric and its upstream tool report"],
                           "affected_downstream": [k for k, v in s["permitted"].items() if not v],
                           "log_paths": [], "retry_history": [], "tool": None, "rule": "l1_gate"})
    for r in elig:
        if r["eligibility"] not in ("ELIGIBLE",):
            events.append({"event_id": f"elig-{r['sample_id']}-{r['analysis_id']}-{r['scale']}", "run_id": run.name,
                           "sample_id": r["sample_id"], "analysis_id": r["analysis_id"], "layer": "L5",
                           "stage": f"{r['analysis_id']} eligibility", "status": r["eligibility"], "severity": "INFO",
                           "code": ",".join(r["reason_codes"]), "failure_class": r["failure_class"] or "SCIENTIFIC_EXCLUSION",
                           "observed": r["evidence"], "expected": r["rule_id"], "criterion": r["rule_id"],
                           "cause_confirmed": "scientific eligibility rule", "causes_suspected": [],
                           "remedies": ["the sample remains eligible for other analyses where their rules permit"],
                           "affected_downstream": [r["analysis_id"]], "log_paths": [], "retry_history": [],
                           "tool": None, "rule": "da_eligibility"})
    resources_ = {"plan": plan.get("resources"), "measured": _benchmarks(run), "passes": status.get("passes", [])}
    executed_rules = sorted({b["rule"] for b in resources_["measured"]})
    for row in resources_["measured"]:
        row["average_cpu_cores"] = row["cpu_time_s"] / row["wall_s"] if row["cpu_time_s"] is not None and row["wall_s"] else None
    # The pipeline version is the one that produced the run (recorded in its snapshot), not the
    # version installed when the report is rebuilt; the renderer version is kept separately.
    tools = {**snap.get("tools_installed", {}),
             **{f"{name} (package)": version for name, version in snap.get("packages", {}).items()}}
    _label_steps(resources_["measured"], tools)
    resources_["timeline"] = _timeline(run, status)
    resources_["planned_vs_actual"] = _planned_vs_actual(run, (plan.get("resources") or {}).get("estimate"),
                                                         resources_["measured"])
    software = {"pfvarx_version": snap.get("pfvarx_version", __version__),
                "tools_installed": tools,
                "rules_with_measured_execution": executed_rules, "git": snap.get("git"),
                "pixi_lock_sha256": snap["files"].get("pixi.lock"), "effective_config_hash": eff["hash"],
                "reference_sha256": plan["reference"].get("reference_sha256"),
                "qc_policy": eff["quality"]["policy_id"], "parameters": eff["parameters"],
                "frozen_files": len(snap["files"]), "warnings": snap.get("warnings", []),
                "configuration": eff["config"], "compute": eff["compute"],
                "environment": {k: snap.get(k) for k in ("python", "platform", "created_utc")}}
    validation = _optional_json(root / "metadata" / "validation_report.json", warnings)
    geo_path = root / "config" / "geo" / "ne_110m_countries.geojson"
    input_files = []
    for sample in plan["samples"]:
        for lane in sample.get("lanes", []):
            for key in ("fastq_1", "fastq_2"):
                if lane.get(key):
                    input_files.append({"sample_id": sample["sample_id"], "lane": lane.get("lane"),
                                        "role": key, "path": lane[key], "size_bytes": _file_size(root, lane[key])})
        if sample.get("bam"):
            input_files.append({"sample_id": sample["sample_id"], "role": "alignment",
                                "path": sample["bam"], "size_bytes": _file_size(root, sample["bam"])})
    reference = eff["config"]["reference"]
    keys = ("fasta", "gff", "core_mask_bed") + (("host_fasta",) if reference.get("host_removal") else ())
    for key in keys:
        value = reference.get(key)
        if value:
            input_files.append({"role": "reference." + key, "path": value, "size_bytes": _file_size(root, value)})
    origin, origin_basis = _data_origin(root, eff)
    data = {
        "schema": "pfvarx.dashboard/2", "data_kind": origin, "data_origin_basis": origin_basis, "generated_at": dt.datetime.now(dt.UTC).isoformat(),
        "project_id": eff["config"]["project_id"], "run_id": run.name, "purpose": eff["config"]["purpose"],
        "final_output": eff["config"]["final_output"],
        "design": plan.get("design") or eff["config"].get("assay_design") or eff["config"]["assay"].get("design"),
        "input_type": eff["config"]["input_type"], "run_status": status, "layers": plan["layers"],
        "panel": plan.get("panel"),
        "libraries": {s["sample_id"]: s.get("assay", {}) for s in plan.get("samples", [])},
        "analysis_selection": plan.get("analysis_selection", "explicit"),
        "analyses_not_applicable": plan.get("analyses_not_applicable", {}),
        "requested_analyses": requested, "internal_analyses": internal, "samples": samples, "qc_matrix": qc_matrix,
        "variants": variants, "variant_summary": counts, "variants_total": variants_total, "variants_embedded": len(variants),
        "variants_full_table": "results/variants/samples/<sample>.evidence.tsv", "annotation": ann, "downstream": downstream, "eligibility": elig, "diagnostics": events,
        "resources": resources_, "software": software, "report_warnings": warnings,
        "reference": eff["config"]["reference"], "log_content": _log_content(run),
        "inputs": {"files": input_files, "validation": None if validation is None else {
            "status": validation["status"], "samples": validation["samples"],
            "reports": [{**r, "size_bytes": _file_size(root, r["path"])} for r in validation["reports"]]}},
        "views": VIEW_REGISTRY, "calibration": "UNCALIBRATED (Phase 3 calibration not yet applied)",
        "geo": _optional_json(geo_path, warnings),
    }
    return _redact(data)


# Output files shown in the dashboard: (pattern on the path inside results/, what the file is).
OUTPUT_DESCRIPTIONS: list[tuple[str, str]] = [
    (r"^report/pfvarx_dashboard\.html$", "This dashboard (opens offline in any browser)"),
    (r"^report/dashboard_data\.json$", "All data shown in this dashboard, machine-readable"),
    (r"^report/clonality_report\.(tsv|txt)$", "Clonality summary per sample"),
    (r"^alignments/.+\.bam$", "Final alignment (BAM) used for variant calling"),
    (r"^alignments/.+\.(bai|csi)$", "Index of the alignment"),
    (r"^qc/centralized_quality_matrix\.json$", "Every QC gate for every sample"),
    (r"\.fastp\.html$", "fastp read-QC report (per-position quality, adapters)"),
    (r"\.fastp\.json$", "fastp read-QC numbers"),
    (r"\.flagstat\.json$", "Alignment statistics (samtools flagstat)"),
    (r"\.markdup\.txt$", "Duplicate-marking metrics (Picard)"),
    (r"mosdepth\.summary\.txt$", "Depth summary in the core genome (mosdepth)"),
    (r"mosdepth\..*dist\.txt$", "Depth distribution (mosdepth)"),
    (r"targets\.regions\.bed\.gz$", "Depth per target region (mosdepth; target gates)"),
    (r"qc_decision\.json$", "QC verdict of the sample and the reason for each gate"),
    (r"\.caller_qc\.json$", "Variant-caller model fit for the sample"),
    (r"\.entry_contract\.json$", "Checks done on an imported BAM/CRAM"),
    (r"^variants/samples/.+\.vcf\.gz$", "Variant calls of the sample (VCF)"),
    (r"^variants/samples/.+\.vcf\.gz\.tbi$", "Index of the VCF"),
    (r"\.evidence\.tsv$", "Every candidate site with its evidence and status (full table)"),
    (r"\.callable\.bed$", "Regions where the sample could be called"),
    (r"\.focal_loss\.bed$", "Regions with focal loss of callability"),
    (r"^cohort/cohort\.vcf\.gz$", "Cohort VCF: all samples merged"),
    (r"^cohort/cohort_qc\.json$", "Cohort QC (gate G5): each sample's completeness in the cohort, and site call rates"),
    (r"^cohort/", "Cohort-level table or status"),
    (r"^annotation/cohort\.ann\.vcf\.gz$", "Annotated cohort VCF (SnpEff)"),
    (r"^annotation/", "Annotation output"),
    (r"^eligibility/", "Which samples were eligible for each analysis"),
    (r"^downstream/DA00\d/status_manifest\.json$", "Run status and denominators of the analysis"),
    (r"^downstream/", "Downstream analysis result"),
]
OUTPUT_GROUPS = {"report": "Report", "alignments": "Alignments", "qc": "Quality control", "variants": "Variant calls",
                 "cohort": "Cohort", "annotation": "Annotation", "eligibility": "Eligibility",
                 "downstream": "Downstream analyses"}


def _describe_output(rel: str) -> str:
    for pattern, text in OUTPUT_DESCRIPTIONS:
        if re.search(pattern, rel):
            return text
    return ""


def _outputs(run: Path) -> dict[str, Any]:
    """Where the run's results are and every final output file, with links that work from the report folder.

    Temporary work files (results/work) are left out. Links are relative to results/report/, so they
    open when this dashboard is opened from the run folder on the computer that produced it.
    """
    res = run / "results"
    files: list[dict[str, Any]] = []
    planned = ["report/pfvarx_dashboard.html", "report/dashboard_data.json", "qc/centralized_quality_matrix.json"]
    seen = set()
    for path in sorted(p for p in res.rglob("*") if p.is_file() and not p.is_symlink()):
        rel = path.relative_to(res).as_posix()
        if rel.startswith("work/"):
            continue
        seen.add(rel)
        files.append(_output_entry(rel, path.stat().st_size))
    for rel in planned:  # written after this list is made
        if rel not in seen:
            files.append(_output_entry(rel, None))
    files.sort(key=lambda f: (list(OUTPUT_GROUPS).index(f["group_id"]) if f["group_id"] in OUTPUT_GROUPS else 99, f["path"]))
    folders = [{"label": OUTPUT_GROUPS.get(d.name, d.name), "path": f"results/{d.name}/", "href": f"../{d.name}/"}
               for d in sorted(res.iterdir(), key=lambda d: (list(OUTPUT_GROUPS).index(d.name) if d.name in OUTPUT_GROUPS else 99, d.name))
               if d.is_dir() and d.name != "work"]
    folders += [{"label": "Logs of every step", "path": "logs/", "href": "../../logs/"},
                {"label": "Time and memory of every step", "path": "benchmarks/", "href": "../../benchmarks/"}]
    run_files = [{"label": "Settings used for this run (frozen)", "path": "frozen_config.json", "href": "../../frozen_config.json"},
                 {"label": "Software and versions used", "path": "snapshot.json", "href": "../../snapshot.json"}]
    return {"run_dir": str(run.resolve()), "results_dir": str(res.resolve()), "results_href": "../",
            "folders": folders, "run_files": [f for f in run_files if (run / f["path"]).is_file()], "files": files}


def _output_entry(rel: str, size: int | None) -> dict[str, Any]:
    group = rel.split("/", 1)[0]
    return {"group_id": group, "group": OUTPUT_GROUPS.get(group, group), "path": "results/" + rel,
            "file": rel.rsplit("/", 1)[-1], "size_bytes": size, "description": _describe_output(rel),
            "href": "../" + rel}


def _file_size(root: Path, value: str) -> int | None:
    try:
        path = root / value
        return path.stat().st_size if path.is_file() else None
    except OSError:
        return None


def build_report(root: Path, run: Path) -> dict[str, Any]:
    if read_json(run / "frozen_config.json")["config"].get("schema_version") == "pfvarx.config/2":
        raise ValueError("Use the original PfVarX release to rebuild a dashboard with legacy analysis IDs")
    data = harvest(root, run)
    # Report rendering can be refreshed independently of the frozen scientific run.
    asset_names = ("build.py", "dashboard_template.html", "dashboard.css", "core.js", "exports.js", "tables.js", "charts.js", "views.js")
    data["software"]["report_renderer"] = {
        "sha256": hashlib.sha256(b"".join(resources.files("pfvarx.reporting").joinpath(name).read_bytes()
                                         for name in asset_names)).hexdigest(),
        "version": __version__,
        "scope": "Current reporting package; frozen upstream scientific provenance is preserved",
    }
    res = run / "results"
    from pfvarx.reporting.clonality import write_report as write_clonality

    data["clonality"] = write_clonality(res, [s["sample_id"] for s in data["samples"]
                                              if (res / "qc" / "samples" / s["sample_id"]).is_dir()], run.name)
    data["outputs"] = _outputs(run)
    write_json(res / "report" / "dashboard_data.json", data)
    write_json(res / "qc" / "centralized_quality_matrix.json",
               {"schema": "pfvarx.qc_matrix/1", "run_id": run.name, "rows": data["qc_matrix"],
                "samples": [{k: s[k] for k in ("sample_id", "processing", "qc_verdict", "label", "uncalibrated_gates")}
                            for s in data["samples"]]})
    template = resources.files("pfvarx.reporting").joinpath("dashboard_template.html").read_text(encoding="utf-8")
    payload = json.dumps(data, default=str, separators=(",", ":")).replace("</", "<\\/")
    css = resources.files("pfvarx.reporting").joinpath("dashboard.css").read_text(encoding="utf-8")
    scripts = "\n".join(resources.files("pfvarx.reporting").joinpath(name).read_text(encoding="utf-8")
                        for name in ("core.js", "exports.js", "tables.js", "charts.js", "views.js"))
    html = (template.replace("/*__PFVARX_STYLE__*/", css)
            .replace("/*__PFVARX_SCRIPT__*/", scripts.replace("</script", "<\\/script"))
            .replace("/*__PFVARX_DATA__*/null", payload))
    out = res / "report" / "pfvarx_dashboard.html"
    out.write_text(html, encoding="utf-8")
    return {"dashboard": str(out.relative_to(root)), "bytes": out.stat().st_size,
            "dashboard_data": str((res / "report" / "dashboard_data.json").relative_to(root)),
            "qc_matrix": str((res / "qc" / "centralized_quality_matrix.json").relative_to(root))}
