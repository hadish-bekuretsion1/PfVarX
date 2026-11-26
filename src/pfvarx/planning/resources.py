"""Compute requests of the workflow rules and job classes, read from ``config/registry/job_resources.yaml``.

The Snakemake rules (``workflow/rules/common.smk``) and the preflight resource plan
(``planning/plan.py``) both read that file, so the requests are written in one place only.
A job class with rules requests the largest threads and memory among its rules.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml

RESOURCES_FILE = Path("config") / "registry" / "job_resources.yaml"


def load(root: Path) -> dict[str, Any]:
    """The whole resource file of the workspace (or pipeline checkout) at *root*."""
    data: dict[str, Any] = yaml.safe_load((root / RESOURCES_FILE).read_text())
    return data


def job_classes(root: Path) -> dict[str, dict[str, Any]]:
    """Job class -> {layer, threads, mem_mb, per_gb_mb, host_index}; requests are the maximum over its rules."""
    out: dict[str, dict[str, Any]] = {}
    for name, entry in load(root)["job_classes"].items():
        rules = entry.get("rules") or {}
        threads = max((r["threads"] for r in rules.values()), default=entry.get("threads", 1))
        mem = max((r["mem_mb"] for r in rules.values()), default=entry.get("mem_mb", 0))
        out[name] = {"layer": entry["layer"], "threads": int(threads), "mem_mb": int(mem),
                     "per_gb_mb": int(entry.get("per_gb_mb", 0)), "host_index": bool(entry.get("host_index", False)),
                     "rules": sorted(rules)}
    return out


def rule_requests(root: Path) -> dict[str, dict[str, int]]:
    """Rule name -> {threads, mem_mb} for every rule listed under a job class."""
    out: dict[str, dict[str, int]] = {}
    for entry in load(root)["job_classes"].values():
        for rule, req in (entry.get("rules") or {}).items():
            out[rule] = {"threads": int(req["threads"]), "mem_mb": int(req["mem_mb"])}
    return out


def planning(root: Path) -> dict[str, Any]:
    """Planning factors (default rule memory, measured-memory headroom, temporary-disk factor, ...)."""
    return dict(load(root)["planning"])
