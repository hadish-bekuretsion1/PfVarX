"""Time, memory and storage learned from earlier runs on this computer (MEASURED_HISTORY).

After each completed run PfVarX records, from its Snakemake benchmark records, every job's
wall time, CPU time and peak memory together with the input size of the sample it processed,
and the size of each result folder. A new run's estimate scales those measurements by its own
input sizes. Nothing is invented: a step with no history on this computer is UNKNOWN, and
synthetic test data is never learned from (it says nothing about real genomes).

The history is kept per computer, outside the workspaces, so every study on the same machine
benefits: ``$PFVARX_HISTORY_DIR`` or ``~/.local/share/pfvarx/history/<host>.jsonl``.
Estimates are advice for planning and the progress page; they never change results.
"""

from __future__ import annotations

import datetime as dt
import json
import os
import socket
import statistics
from pathlib import Path
from typing import Any

# Benchmark folder (one per workflow step) -> job class used by the resource plan.
RULE_TO_JOB = {
    "fastp": "fastp", "bwa_mem2": "bwa_mem2_align", "merge_lanes": "metrics", "markdup": "markdup",
    "insert_metrics_preclip": "metrics", "trim_primers": "fgbio_clip", "alignment_entry": "alignment_entry",
    "entry_finalize": "alignment_entry", "umi_fastq_to_bam": "fgbio_umi", "umi_align_zipper": "fgbio_umi",
    "umi_group": "fgbio_umi", "umi_consensus": "fgbio_umi", "umi_realign": "fgbio_umi",
    "umi_preclip_and_trim": "fgbio_umi", "flagstat": "metrics", "mosdepth_core": "metrics", "mosdepth": "metrics",
    "l1_gate": "metrics", "pfvarx_call": "pfvarx_call", "normalise_sample": "pfvarx_call",
    "cohort_merge": "cohort_merge", "cohort_status": "cohort_merge", "cohort_evidence": "cohort_merge",
    "cohort_qc": "cohort_merge", "snpeff_build": "snpeff", "snpeff_annotate": "snpeff",
    "da_eligibility": "downstream", "downstream": "downstream",
}
# Job classes run once per sample; the others run once for the cohort.
PER_SAMPLE = {"fastp", "bwa_mem2_align", "markdup", "metrics", "fgbio_clip", "alignment_entry", "fgbio_umi",
              "pfvarx_call"}
JOB_NAMES = {
    "fastp": "Read QC and trimming (fastp)", "bwa_mem2_align": "Read alignment (bwa-mem2)",
    "markdup": "Mark duplicates (Picard)", "metrics": "Alignment QC (samtools, mosdepth)",
    "fgbio_clip": "Primer trimming (fgbio)", "fgbio_umi": "UMI consensus (fgbio)",
    "alignment_entry": "Imported BAM checks (samtools)", "pfvarx_call": "Variant calling (PfVarX caller, bcftools)",
    "cohort_merge": "Cohort merge, evidence and QC (bcftools)", "snpeff": "Annotation (SnpEff)",
    "downstream": "Downstream analyses", "dashboard": "Report dashboard",
}
STORAGE_FOLDERS = ("alignments", "qc", "variants", "cohort", "annotation", "downstream", "report")


def history_file() -> Path:
    base = Path(os.environ.get("PFVARX_HISTORY_DIR") or Path.home() / ".local/share/pfvarx/history")
    return base / f"{socket.gethostname()}.jsonl"


def sample_input_gb(root: Path, plan: dict[str, Any]) -> dict[str, float]:
    """Input size (GB) of each sample: its FASTQ files, or its BAM/CRAM."""
    sizes: dict[str, float] = {}
    for s in plan.get("samples", []):
        total = 0
        for lane in s.get("lanes", []):
            for f in (lane.get("fastq_1"), lane.get("fastq_2")):
                if f and (root / f).is_file():
                    total += (root / f).stat().st_size
        if s.get("bam") and (root / s["bam"]).is_file():
            total += (root / s["bam"]).stat().st_size
        sizes[s["sample_id"]] = total / 1e9
    return sizes


def _folder_mb(path: Path) -> float:
    return sum(p.stat().st_size for p in path.rglob("*") if p.is_file() and not p.is_symlink()) / 1e6 if path.is_dir() else 0.0


def record(root: Path, run: Path) -> dict[str, Any] | None:
    """Add a finished real-data run to this computer's history; synthetic runs are skipped."""
    from pfvarx.io import read_json
    from pfvarx.reporting.build import _benchmarks, _data_origin

    eff = read_json(run / "frozen_config.json")
    if _data_origin(root, eff)[0] == "synthetic":
        return None
    plan = read_json(run / "frozen_plan.json")
    sizes = sample_input_gb(root, plan)
    total_gb = sum(sizes.values())
    jobs = []
    for b in _benchmarks(run):
        job = RULE_TO_JOB.get(b["rule"])
        if job is None or b["wall_s"] is None:
            continue
        sample = next((s for s in sizes if b["job"] == s or b["job"].startswith(s + ".")), None)
        jobs.append({"job": job, "rule": b["rule"], "sample": sample,
                     "input_gb": sizes.get(sample or "", total_gb) if job in PER_SAMPLE else total_gb,
                     "wall_s": b["wall_s"], "cpu_s": b["cpu_time_s"], "rss_mb": b["max_rss_mb"]})
    storage = {name: _folder_mb(run / "results" / name) for name in STORAGE_FOLDERS}
    entry = {"schema": "pfvarx.history/1", "run_id": run.name, "recorded_utc": dt.datetime.now(dt.UTC).isoformat(),
             "design": plan.get("design"), "input_type": plan.get("input_type"), "final_output": plan.get("final_output"),
             "samples": len(sizes), "input_gb": round(total_gb, 4), "cpus": os.cpu_count(),
             "threads_per_job": eff["compute"]["per_job"]["max_threads"], "jobs": jobs, "storage_mb": storage}
    path = history_file()
    path.parent.mkdir(parents=True, exist_ok=True)
    known = {json.loads(line).get("run_id") for line in path.read_text().splitlines() if line.strip()} if path.is_file() else set()
    if run.name not in known:
        with path.open("a") as handle:
            handle.write(json.dumps(entry) + "\n")
    return entry


def load_history() -> list[dict[str, Any]]:
    path = history_file()
    if not path.is_file():
        return []
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]


def measured_peak_mb(job: str, design: str | None = None, input_type: str | None = None,
                     history: list[dict[str, Any]] | None = None) -> float | None:
    """Highest memory one job of this class used in earlier real-data runs on this computer
    (same design and input type when there are such runs)."""
    runs = history if history is not None else load_history()
    same = [r for r in runs if r.get("design") == design and r.get("input_type") == input_type] or runs
    peaks = [j["rss_mb"] for r in same for j in r["jobs"] if j["job"] == job and j.get("rss_mb")]
    return max(peaks) if peaks else None


def estimate(job_classes: list[str], sample_gb: dict[str, float], max_jobs: int, design: str | None = None,
             input_type: str | None = None, history: list[dict[str, Any]] | None = None,
             layers: list[int] | None = None) -> dict[str, Any]:
    """Expected time, peak memory and storage for each job class of a new run.

    Per-sample steps scale by each sample's input size (median seconds per input GB measured
    here before); cohort steps scale by the total input size. Wall-clock time assumes samples
    run ``max_jobs`` at a time. Every number names how many earlier jobs it rests on.
    """
    runs = history if history is not None else load_history()
    same = [r for r in runs if (design is None or r.get("design") == design)
            and (input_type is None or r.get("input_type") == input_type)] or runs
    n_samples = max(len(sample_gb), 1)
    total_gb = sum(sample_gb.values())
    steps: list[dict[str, Any]] = []
    unknown: list[str] = []
    sample_minutes: dict[str, float] = {s: 0.0 for s in sample_gb}
    cohort_minutes = 0.0
    for job in job_classes:
        if job == "dashboard":
            continue
        # One measurement per earlier (run, sample) for this job class: summed wall time of its rules.
        per_unit: dict[tuple[str, str | None], dict[str, float]] = {}
        # Prefer runs with the same design and input type; fall back to any run that ran this step.
        source = [r for r in same if any(j["job"] == job for j in r["jobs"])]
        fallback = not source
        source = source or [r for r in runs if any(j["job"] == job for j in r["jobs"])]
        for r in source:
            for j in r["jobs"]:
                if j["job"] != job:
                    continue
                key = (r["run_id"], j["sample"] if job in PER_SAMPLE else None)
                unit = per_unit.setdefault(key, {"wall_s": 0.0, "rss_mb": 0.0, "input_gb": j["input_gb"]})
                unit["wall_s"] += j["wall_s"] or 0.0
                unit["rss_mb"] = max(unit["rss_mb"], j["rss_mb"] or 0.0)
        rates = [u["wall_s"] / u["input_gb"] for u in per_unit.values() if u["input_gb"] and u["input_gb"] > 0.005]
        if not rates:
            unknown.append(job)
            steps.append({"job": job, "name": JOB_NAMES.get(job, job), "per_sample": job in PER_SAMPLE,
                          "jobs": n_samples if job in PER_SAMPLE else 1, "minutes": None, "peak_ram_mb": None,
                          "basis": "no earlier run of this step on this computer", "confidence": "UNKNOWN"})
            continue
        rate = statistics.median(rates)  # seconds per input GB
        peak = max(u["rss_mb"] for u in per_unit.values())
        if job in PER_SAMPLE:
            per = {s: rate * gb / 60 for s, gb in sample_gb.items()}
            for s, minutes in per.items():
                sample_minutes[s] += minutes
            minutes = sum(per.values())
        else:
            minutes = rate * total_gb / 60
            cohort_minutes += minutes
        steps.append({"job": job, "name": JOB_NAMES.get(job, job), "per_sample": job in PER_SAMPLE,
                      "jobs": n_samples if job in PER_SAMPLE else 1, "minutes": round(minutes, 1),
                      "minutes_per_sample": {s: round(rate * gb / 60, 1) for s, gb in sample_gb.items()} if job in PER_SAMPLE else None,
                      "peak_ram_mb": round(peak), "seconds_per_input_gb": round(rate, 1),
                      "basis": f"median of {len(rates)} earlier job(s) on this computer"
                               + (" (runs with another design or input type)" if fallback else ""),
                      "confidence": "MEASURED_HISTORY"})
    # Samples run max_jobs at a time: the longest chains decide the wall clock.
    lanes = [0.0] * max(1, min(max_jobs, n_samples))
    for minutes in sorted(sample_minutes.values(), reverse=True):
        lanes[lanes.index(min(lanes))] += minutes
    wall = max(lanes) + cohort_minutes
    storage = {}
    # Only the folders this final output writes: variants (layer 2), cohort (3), annotation (4), analyses (5).
    produced = {"variants": 2, "cohort": 3, "annotation": 4, "downstream": 5}
    for folder in STORAGE_FOLDERS:
        if layers is not None and folder in produced and produced[folder] not in layers:
            continue
        # Only runs that produced this folder; same design and input type first.
        made = [r for r in same if r.get("input_gb") and r["storage_mb"].get(folder)]
        made = made or [r for r in runs if r.get("input_gb") and r["storage_mb"].get(folder)]
        ratios = [r["storage_mb"][folder] / r["input_gb"] for r in made]
        storage[folder] = round(statistics.median(ratios) * total_gb) if ratios else None
    return {"schema": "pfvarx.estimate/1", "input_gb": round(total_gb, 3), "samples": len(sample_gb),
            "max_jobs": max_jobs, "steps": steps, "unknown_steps": unknown,
            "wall_minutes": None if unknown else round(wall, 1),
            "wall_minutes_known_steps": round(wall, 1),
            "peak_ram_mb": max((int(s["peak_ram_mb"] or 0) for s in steps), default=0) or None,
            "storage_mb": storage, "storage_total_mb": sum(v for v in storage.values() if v) or None,
            "history_runs": len(same),
            "confidence": "UNKNOWN" if unknown or not same else "MEASURED_HISTORY",
            "basis": f"{len(same)} earlier real-data run(s) on this computer" if same else
                     "no earlier real-data run on this computer yet; this run will be learned from"}


def _size(mb: float | None) -> str:
    if mb is None:
        return "unknown"
    return f"{mb / 1024:.1f} GB" if mb >= 1024 else f"{mb:.0f} MB"


def _minutes(minutes: float | None) -> str:
    if minutes is None:
        return "unknown"
    if minutes < 1:
        return "under 1 min"
    whole = int(round(minutes))
    return f"{whole // 60} h {whole % 60:02d} min" if whole >= 60 else f"{whole} min"


def estimate_table(est: dict[str, Any], final_output: str | None = None, input_type: str | None = None) -> str:
    """The estimate as a plain-text table: shown by 'workflow-check' in the terminal and in the GUI."""
    rows = [("Step (tools)", "Jobs", "Time", "Peak RAM", "Based on")]
    for s in est.get("steps", []):
        rows.append((s["name"], str(s["jobs"]), _minutes(s["minutes"]), _size(s["peak_ram_mb"]), s["basis"]))
    widths = [max(len(r[i]) for r in rows) for i in range(4)]
    lines = [f"Expected resources for {final_output or 'this final output'}"
             + (f" from {input_type} input" if input_type else "")
             + f": {est.get('samples', 0)} sample(s), {est.get('input_gb', 0):.2f} GB of input, "
             + f"{est.get('max_jobs', 1)} job(s) at a time", ""]
    for i, r in enumerate(rows):
        lines.append("  ".join(r[c].ljust(widths[c]) for c in range(4)) + "  " + r[4])
        if i == 0:
            lines.append("-" * (sum(widths) + 8 + len("Based on")))
    lines.append("")
    total = est.get("wall_minutes")
    if total is None and est.get("unknown_steps"):
        lines.append(f"Total time: unknown (no earlier run on this computer for: {', '.join(est['unknown_steps'])}); "
                     f"the steps with history need about {_minutes(est.get('wall_minutes_known_steps'))}")
    else:
        lines.append(f"Total time: about {_minutes(total)}")
    lines.append(f"Peak memory of one job: {_size(est.get('peak_ram_mb'))}")
    storage = est.get("storage_mb") or {}
    lines.append("Storage of the results: " + (", ".join(f"{k} {_size(v)}" for k, v in storage.items() if v)
                                               or "unknown") + f"; total {_size(est.get('storage_total_mb'))}")
    lines.append(f"Basis: {est.get('basis', '')}. Learned from real-data runs on this computer and scaled by "
                 "input size; a guide for planning, never used to change results.")
    return "\n".join(lines)
