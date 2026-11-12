"""Summarise measured per-job resources (Snakemake benchmark records) into a profile.

The profile (config/profiles/resources.yaml) records the hardware, dataset size and
executor it was measured on; the planner uses it only as MEASURED_PROFILE evidence for
the same job classes, and never generalises it to untested hardware without saying so.
"""

from __future__ import annotations

import os
import platform
from pathlib import Path
from typing import Any

import yaml

from pfvarx.io import read_json
from pfvarx.reporting.build import _benchmarks

RULE_TO_JOB = {"fastp": "fastp", "bwa_mem2": "bwa_mem2_align", "markdup": "markdup", "pfvarx_call": "pfvarx_call",
               "downstream": "downstream", "trim_primers": "metrics"}


def summarise(root: Path, run: Path) -> dict[str, Any]:
    rows = _benchmarks(run)
    jobs: dict[str, dict[str, Any]] = {}
    for r in rows:
        job = RULE_TO_JOB.get(r["rule"])
        if not job or r["max_rss_mb"] is None:
            continue
        cur = jobs.setdefault(job, {"peak_rss_mb": 0.0, "wall_min": 0.0, "n": 0})
        cur["peak_rss_mb"] = max(cur["peak_rss_mb"], r["max_rss_mb"])
        cur["wall_min"] = max(cur["wall_min"], (r["wall_s"] or 0) / 60)
        cur["n"] += 1
    plan = read_json(run / "frozen_plan.json")
    basis = {"run_id": run.name, "input_gb": plan.get("input_gb"), "samples": len(plan["samples"]),
             "cpu_count": os.cpu_count(), "platform": platform.platform(), "executor": "local"}
    for j in jobs.values():
        j["basis"] = f"max over {j['n']} jobs in run {run.name} ({basis['input_gb']} GB input, {basis['cpu_count']} CPUs)"
    prof = {"schema": "pfvarx.resource_profile/1", "measured_on": basis, "jobs": jobs}
    out = root / "config" / "profiles" / "resources.yaml"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(yaml.safe_dump(prof, sort_keys=True))
    return {"profile": str(out.relative_to(root)), **prof}
