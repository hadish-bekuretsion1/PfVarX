"""Objective and resource planning (L0.7). Read-only: it never creates output artefacts.

Resource requests per job class are FRAMEWORK_POLICY starting requests. Runtime
estimates are reported only from measured profiles (``config/profiles/resources.yaml``,
written by the benchmarking workspace); otherwise the estimate is UNKNOWN.
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path
from typing import Any

import yaml

from pfvarx.configuration.models import Bundle
from pfvarx.constants import Design, FinalOutput, InputType, Reason

# Starting requests, summarising the rules in workflow/rules (tests/unit/test_resource_plan.py
# keeps the two equal); memory scales with input size where a predictor is known.
JOB_CLASSES: dict[str, dict[str, Any]] = {
    "fastp": {"layer": 1, "threads": 4, "mem_mb": 2000, "per_gb_mb": 0},
    "bwa_mem2_index": {"layer": 0, "threads": 1, "mem_mb": 4000, "per_gb_mb": 0},
    "bwa_mem2_align": {"layer": 1, "threads": 4, "mem_mb": 6000, "per_gb_mb": 0},
    "markdup": {"layer": 1, "threads": 1, "mem_mb": 4000, "per_gb_mb": 500},
    "fgbio_umi": {"layer": 1, "threads": 4, "mem_mb": 6000, "per_gb_mb": 500},
    "alignment_entry": {"layer": 1, "threads": 1, "mem_mb": 4000, "per_gb_mb": 0},
    "fgbio_clip": {"layer": 1, "threads": 1, "mem_mb": 3000, "per_gb_mb": 0},
    "metrics": {"layer": 1, "threads": 2, "mem_mb": 1000, "per_gb_mb": 0},
    "pfvarx_call": {"layer": 2, "threads": 8, "mem_mb": 5000, "per_gb_mb": 0},
    "cohort_merge": {"layer": 3, "threads": 1, "mem_mb": 2000, "per_gb_mb": 0},
    "snpeff": {"layer": 4, "threads": 1, "mem_mb": 4000, "per_gb_mb": 0},
    "downstream": {"layer": 5, "threads": 2, "mem_mb": 2000, "per_gb_mb": 0},
    "dashboard": {"layer": 6, "threads": 1, "mem_mb": 1000, "per_gb_mb": 0},
}
# Job classes whose bwa-mem2 step loads the combined parasite+host index (Human DNA removal, wgs only).
HOST_INDEX_JOBS = ("bwa_mem2_align",)


def layers_for(cfg_final: FinalOutput, input_type: InputType) -> list[int]:
    start = {InputType.FASTQ: 1, InputType.BAM: 2, InputType.CRAM: 2, InputType.VCF: 3}[input_type]
    stop = {FinalOutput.FO001: 1, FinalOutput.FO002: 2, FinalOutput.FO003: 3,
            FinalOutput.FO004: 5}[cfg_final]
    layers = [0] + [layer for layer in range(start, stop + 1)]
    if cfg_final == FinalOutput.FO004:
        layers = sorted(set(layers) | {3, 4})
    return layers + [6]


def jobs_for(layers: list[int], design: Design, input_type: InputType = InputType.FASTQ) -> list[str]:
    jobs = ["dashboard"]
    if input_type in (InputType.BAM, InputType.CRAM):
        # Imported alignments may still need primer trimming or UMI consensus (by documented state).
        jobs += ["alignment_entry", "metrics"]
        jobs += ["fgbio_umi"] if design == Design.TARGETED_UMI else []
        jobs += ["fgbio_clip"] if design.targeted else []
    elif 1 in layers:
        jobs += ["fastp", "bwa_mem2_align", "metrics"]
        jobs += ["markdup"] if design == Design.WGS else []
        jobs += ["fgbio_umi"] if design == Design.TARGETED_UMI else []
        jobs += ["fgbio_clip"] if design == Design.TARGETED_NO_UMI else []
    if 2 in layers:
        jobs.append("pfvarx_call")
    if 3 in layers:
        jobs.append("cohort_merge")
    if 4 in layers:
        jobs.append("snpeff")
    if 5 in layers:
        jobs.append("downstream")
    return jobs


def measured_profiles(root: Path) -> dict[str, Any]:
    p = root / "config" / "profiles" / "resources.yaml"
    return yaml.safe_load(p.read_text()) if p.is_file() else {}


def resource_plan(root: Path, bundle: Bundle, input_gb: float, n_samples: int) -> dict[str, Any]:
    comp = bundle.compute
    layers = layers_for(bundle.config.final_output, bundle.config.input_type)
    if bundle.config.analyses.selection == "auto" and bundle.config.reference.gff is None and 4 in layers:
        layers.remove(4)
    jobs = jobs_for(layers, bundle.config.assay_design, bundle.config.input_type)
    profiles = measured_profiles(root).get("jobs", {})
    host_extra = 0
    if (bundle.config.reference.host_fasta_used is not None and bundle.config.input_type == InputType.FASTQ
            and bundle.config.assay_design == Design.WGS):
        from pfvarx.configuration.loader import load_parameters

        host_extra = int(load_parameters(root)["layer1.bwa_mem2_host_extra_mem_mb"]["value"])
    from pfvarx.planning.history import measured_peak_mb

    rows, blocked = [], []
    warnings = []
    for name in jobs:
        jc = JOB_CLASSES[name]
        threads = min(jc["threads"], comp.per_job.max_threads)
        mem = jc["mem_mb"] + int(jc["per_gb_mb"] * input_gb / max(n_samples, 1))
        mem += host_extra if name in HOST_INDEX_JOBS else 0
        prof = profiles.get(name)
        basis = "MEASURED_PROFILE" if prof else "FRAMEWORK_POLICY_REQUEST"
        # Memory a job of this class really used here before, with 20 % headroom (MEASURED_HISTORY).
        # Kept within the per-job ceiling; a tighter ceiling is a warning, not a block.
        peak = measured_peak_mb(name, bundle.config.assay_design.value, bundle.config.input_type.value)
        if peak and int(peak * 1.2) > mem:
            needed = int(peak * 1.2)
            mem, basis = min(needed, comp.per_job.max_mem_mb), "MEASURED_HISTORY"
            if needed > comp.per_job.max_mem_mb:
                warnings.append({"reason": "RESOURCE_MEMORY_TIGHT", "job": name,
                                 "detail": f"{name} used up to {int(peak)} MB in an earlier run on this computer; "
                                           f"per_job.max_mem_mb is {comp.per_job.max_mem_mb}. Raise it if the "
                                           "machine has the memory, or the job may run out of memory"})
        if prof:
            mem = max(mem, int(prof["peak_rss_mb"] * 1.2))
            runtime = {"value_min": prof["wall_min"], "confidence": "MEASURED_PROFILE",
                       "basis": prof.get("basis", "benchmark workspace")}
        else:
            runtime = {"value_min": None, "confidence": "UNKNOWN",
                       "basis": Reason.RESOURCE_ESTIMATE_UNAVAILABLE.value}
        row = {"job": name, "layer": jc["layer"], "threads": threads, "mem_mb": mem,
               "runtime": runtime, "memory_basis": basis}
        if mem > comp.per_job.max_mem_mb:
            row["blocked"] = Reason.RESOURCE_MEMORY_EXCEEDED.value
            blocked.append({"job": name, "reason": Reason.RESOURCE_MEMORY_EXCEEDED.value,
                            "requested_mb": mem, "limit_mb": comp.per_job.max_mem_mb,
                            "remedies": ["raise per_job.max_mem_mb and max_mem_mb if the machine has the memory",
                                         "split inputs into smaller samples/lanes"]
                            + (["host-read removal loads the parasite+host bwa-mem2 index "
                                f"(+{host_extra} MB, layer1.bwa_mem2_host_extra_mem_mb): run on a machine with more "
                                "memory, or set reference.host_removal to false"] if host_extra and name in HOST_INDEX_JOBS
                               else [])})
        rows.append(row)
    host_cpus = len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else (os.cpu_count() or 1)
    host_mem = _host_mem_mb()
    free = _available_mem_mb()
    at_once = min(comp.max_jobs, max(n_samples, 1))
    biggest = max((r["mem_mb"] for r in rows), default=0)
    if free is not None and biggest and biggest * at_once > free:
        fits = max(1, free // biggest)
        warnings.append({"reason": "RESOURCE_MEMORY_IN_USE",
                         "detail": f"{at_once} job(s) at once may need {biggest * at_once} MB, but only {free} MB of "
                                   f"memory is free now (other programs use the rest). Close other programs or set "
                                   f"max_jobs to {fits} in config/compute.yaml, or the system may stop the run "
                                   "when memory runs out"})
    if comp.max_cores > host_cpus:
        warnings.append({"reason": "RESOURCE_CPU_LIMIT", "detail": f"max_cores {comp.max_cores} > host CPUs {host_cpus}; "
                         "jobs will oversubscribe the CPUs (slower, not unsafe)"})
    if host_mem is not None and comp.max_mem_mb > host_mem:
        blocked.append({"job": "all", "reason": Reason.RESOURCE_MEMORY_EXCEEDED.value, "requested_mb": comp.max_mem_mb,
                        "limit_mb": host_mem, "basis": "compute.yaml max_mem_mb exceeds the detected memory limit",
                        "remedies": [f"set max_mem_mb <= {host_mem}"]})
    tmp_need = int(input_gb * 3 * 1024)
    if tmp_need > comp.tmp_disk_mb:
        blocked.append({"job": "all", "reason": Reason.RESOURCE_DISK_EXHAUSTED.value,
                        "requested_mb": tmp_need, "limit_mb": comp.tmp_disk_mb,
                        "basis": "3x input size (FRAMEWORK_POLICY heuristic, not measured)",
                        "remedies": ["increase tmp_disk_mb", "free disk space"]})
    total = None
    if all(r["runtime"]["value_min"] is not None for r in rows):
        total = sum(r["runtime"]["value_min"] for r in rows) * n_samples / max(comp.max_jobs, 1)
    return {"layers": layers, "jobs": rows, "blocked": blocked, "feasible": not blocked, "warnings": warnings,
            "host": {"cpus": host_cpus, "mem_mb": host_mem},
            "estimated_total_runtime_min": total,
            "estimate_confidence": "MEASURED_PROFILE" if total is not None else "UNKNOWN"}


def _available_mem_mb() -> int | None:
    """Memory free for new work now (Linux MemAvailable), or None where it cannot be read."""
    try:
        for line in Path("/proc/meminfo").read_text().splitlines():
            if line.startswith("MemAvailable:"):
                return int(line.split()[1]) // 1024
    except OSError:
        return None
    return None


def _host_mem_mb() -> int | None:
    limits: list[int] = []
    if sys.platform == "darwin":
        try:
            result = subprocess.run(["sysctl", "-n", "hw.memsize"], check=True,
                                    capture_output=True, text=True, timeout=5)
            size = int(result.stdout.strip()) // (1024 * 1024)
            if size > 0:
                limits.append(size)
        except (OSError, ValueError, subprocess.SubprocessError):
            pass
        allocation = os.environ.get("SLURM_MEM_PER_NODE", "")
        if os.environ.get("SLURM_JOB_ID") and allocation.isdigit() and int(allocation) > 0:
            limits.append(int(allocation))
        return min(limits) if limits else None
    try:
        for line in Path("/proc/meminfo").read_text().splitlines():
            if line.startswith("MemTotal:"):
                limits.append(int(line.split()[1]) // 1024)
                break
    except (OSError, ValueError):
        pass
    mount = Path("/sys/fs/cgroup")
    try:
        group = next(line[3:] for line in Path("/proc/self/cgroup").read_text().splitlines()
                     if line.startswith("0::"))
        current = (mount / group.lstrip("/")).resolve()
        while current.is_relative_to(mount):
            try:
                value = (current / "memory.max").read_text().strip()
                if value.isdigit():
                    limits.append(int(value) // (1024 * 1024))
            except OSError:
                pass
            if current == mount:
                break
            current = current.parent
    except (OSError, StopIteration):
        pass
    # SLURM's per-node memory request is expressed in MiB.
    allocation = os.environ.get("SLURM_MEM_PER_NODE", "")
    if os.environ.get("SLURM_JOB_ID") and allocation.isdigit() and int(allocation) > 0:
        limits.append(int(allocation))
    return min(limits) if limits else None

