"""Content-addressed snapshot of the run-defining files (configuration, policy, code)."""

from __future__ import annotations

import datetime as dt
import json
import platform
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

from pfvarx import __version__
from pfvarx.io import read_json, sha256

SCHEMA = "pfvarx.snapshot/2"
REQUIRED = ("pixi.toml", "pyproject.toml", "config/config.yaml", "config/quality.yaml",
            "config/compute.yaml", "workflow/Snakefile")


def run_defining_files(root: Path) -> list[Path]:
    files = {root / r for r in REQUIRED}
    for pattern, base in (("*.yaml", "config/registry"), ("*.tsv", "config/registry"),
                          ("*.yaml", "config/parameters"), ("*.py", "src/pfvarx"),
                          ("*.py", "workflow/scripts"), ("*.R", "workflow/scripts"),
                          ("*.smk", "workflow/rules")):
        files.update((root / base).rglob(pattern))
    for name in ("bootstrap.py", "pfvarx", "config/installation.json"):
        if (root / name).is_file():
            files.add(root / name)
    if (root / "pixi.lock").is_file():
        files.add(root / "pixi.lock")
    missing = [r for r in REQUIRED if not (root / r).is_file()]
    if missing:
        raise ValueError(f"required run-defining files missing: {missing}")
    return sorted((f for f in files if f.is_file()), key=lambda p: p.as_posix())


def git_identity(root: Path) -> dict[str, Any]:
    if not shutil.which("git"):
        return {"status": "UNKNOWN", "commit": None, "dirty": None}
    try:
        commit = subprocess.run(["git", "-C", str(root), "rev-parse", "HEAD"], check=True,
                                capture_output=True, text=True).stdout.strip()
        dirty = subprocess.run(["git", "-C", str(root), "status", "--porcelain"], check=True,
                               capture_output=True, text=True).stdout.strip()
        return {"status": "KNOWN", "commit": commit, "dirty": bool(dirty)}
    except subprocess.CalledProcessError:
        return {"status": "UNKNOWN", "commit": None, "dirty": None}


def tool_versions() -> dict[str, str]:
    probes = {"samtools": ["samtools", "--version"], "bcftools": ["bcftools", "--version"],
              "fastp": ["fastp", "--version"], "bwa-mem2": ["bwa-mem2", "version"],
              "snakemake": ["snakemake", "--version"], "mosdepth": ["mosdepth", "--version"]}
    out = {}
    for name, cmd in probes.items():
        if not shutil.which(cmd[0]):
            out[name] = "NOT_INSTALLED"
            continue
        r = subprocess.run(cmd, capture_output=True, text=True)
        text = (r.stdout or r.stderr).strip().splitlines()
        out[name] = text[-1] if name == "bwa-mem2" and text else (text[0] if text else "UNKNOWN")
    return out


# Conda packages whose versions decide scientific results. Read from the environment's own
# package records, so a tool that is never probed on the command line is still recorded.
SCIENTIFIC_PACKAGES = (
    "samtools", "bcftools", "htslib", "pysam", "fastp", "bwa-mem2", "mosdepth", "picard", "picard-slim", "fgbio", "fgbio-minimal",
    "snpeff", "snakemake-minimal", "iqtree", "hmmibd", "rtg-tools", "scikit-allel", "numpy", "scipy",
    "r-base", "r-moire", "r-dcifer", "r-mgcv",
)


def package_versions(prefix: Path | None = None) -> dict[str, str]:
    """Exact version and build of each installed scientific package (empty outside a conda env)."""
    meta = (prefix or Path(sys.prefix)) / "conda-meta"
    out: dict[str, str] = {}
    if not meta.is_dir():
        return out
    for record in meta.glob("*.json"):
        try:
            info = json.loads(record.read_text())
        except (OSError, ValueError):
            continue
        if info.get("name") in SCIENTIFIC_PACKAGES:
            out[info["name"]] = f"{info.get('version')} ({info.get('build')})"
    return dict(sorted(out.items()))


def build(root: Path) -> dict[str, Any]:
    manifest = {f.relative_to(root).as_posix(): sha256(f) for f in run_defining_files(root)}
    prepared = root / "prepared" / "manifest.json"
    warnings = []
    if "pixi.lock" not in manifest:
        warnings.append("PIXI_LOCK_MISSING")
    gid = git_identity(root)
    if gid["status"] == "UNKNOWN":
        warnings.append("GIT_REVISION_UNKNOWN")
    elif gid["dirty"]:
        warnings.append("GIT_WORKTREE_DIRTY")
    return {"schema": SCHEMA, "created_utc": dt.datetime.now(dt.UTC).isoformat(),
            "accelerators": accelerators(),
            "pfvarx_version": __version__, "python": platform.python_version(),
            "platform": platform.platform(), "git": gid, "files": manifest,
            "prepared": read_json(prepared) if prepared.is_file() else None,
            "tools_installed": tool_versions(), "packages": package_versions(), "warnings": warnings}


def accelerators() -> dict[str, Any]:
    """Graphics processors present on this computer, recorded for provenance only.

    PfVarX computes on CPUs: every pinned tool (fastp, bwa-mem2, samtools, fgbio, mosdepth,
    SnpEff, IQ-TREE, hmmIBD and the R packages) is a CPU program, and the native caller uses
    NumPy and SciPy. A run on a machine with a GPU uses its CPUs and gives identical results.
    This is recorded so a report says which machine produced it, and so a reader never has to
    wonder whether an accelerator changed a number. ``used_for_compute`` is always false.
    """
    found: list[dict[str, str]] = []
    detection = "none"
    if shutil.which("nvidia-smi"):
        detection = "nvidia-smi"
        try:
            result = subprocess.run(["nvidia-smi", "--query-gpu=name,memory.total,driver_version",
                                     "--format=csv,noheader"], capture_output=True, text=True, timeout=10)
            for line in (result.stdout or "").strip().splitlines():
                parts = [p.strip() for p in line.split(",")]
                if parts and parts[0]:
                    found.append({"vendor": "NVIDIA", "name": parts[0],
                                  "memory": parts[1] if len(parts) > 1 else "",
                                  "driver": parts[2] if len(parts) > 2 else ""})
        except (OSError, subprocess.SubprocessError):
            detection = "nvidia-smi-failed"
    elif sys.platform == "darwin" and shutil.which("system_profiler"):
        detection = "system_profiler"
        try:
            result = subprocess.run(["system_profiler", "-json", "SPDisplaysDataType"],
                                    capture_output=True, text=True, timeout=30)
            for item in json.loads(result.stdout or "{}").get("SPDisplaysDataType", []):
                name = item.get("sppci_model") or item.get("_name")
                if name:
                    found.append({"vendor": item.get("spdisplays_vendor", "Apple"), "name": str(name),
                                  "memory": str(item.get("spdisplays_vram", "")), "driver": ""})
        except (OSError, subprocess.SubprocessError, json.JSONDecodeError):
            detection = "system_profiler-failed"
    return {"detected": found, "detection": detection, "count": len(found), "used_for_compute": False,
            "note": "Recorded for provenance. PfVarX computes on CPUs only; no result depends on a GPU."}


def drift(root: Path, snapshot: dict[str, Any]) -> list[dict[str, str]]:
    issues = []
    for rel, digest in snapshot["files"].items():
        p = root / rel
        if not p.is_file():
            issues.append({"path": rel, "issue": "MISSING"})
        elif sha256(p) != digest:
            issues.append({"path": rel, "issue": "CONTENT_CHANGED"})
    return issues
