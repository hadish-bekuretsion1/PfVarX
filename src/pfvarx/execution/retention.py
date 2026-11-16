"""Conservative cleanup of completed-run intermediates only; preview by default."""

from __future__ import annotations

import fcntl
import stat
from contextlib import contextmanager
from pathlib import Path
from typing import Any

from pfvarx.io import read_json, write_json


@contextmanager
def workspace_lock(root: Path):
    lock = root / ".pfvarx.lock"
    if lock.is_symlink():
        raise ValueError("Workspace lock must not be a symlink")
    with lock.open("a") as handle:
        try:
            fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise ValueError("Workspace has active execution/cleanup; wait for it to finish") from exc
        try:
            yield
        finally:
            fcntl.flock(handle, fcntl.LOCK_UN)


def cleanup_run(root: Path, run: Path, apply: bool = False) -> dict[str, Any]:
    root = root.resolve()
    run = run.absolute()
    if run.is_symlink() or not run.resolve().is_relative_to(root / "runs"):
        raise ValueError("Invalid run cleanup boundary")
    if read_json(run / "status.json")["status"] != "COMPLETED":
        raise ValueError("Cleanup requires a COMPLETED run; failed/interrupted evidence is retained")
    work = run / "results/work"
    if (run / "results").is_symlink() or work.is_symlink():
        raise ValueError("Intermediate directory symlink refused")
    cfg = read_json(run / "frozen_config.json")["config"]
    plan = read_json(run / "frozen_plan.json")
    inputs = [root / str(v) for v in cfg.get("reference", {}).values() if isinstance(v, str)]
    for sample in plan["samples"]:
        if sample.get("bam"):
            inputs.append(root / sample["bam"])
        for lane in sample.get("lanes", []):
            inputs.extend(root / lane[k] for k in ("fastq_1", "fastq_2") if lane.get(k))
    for value in [
        cfg.get("samples_tsv"),
        cfg.get("assay", {}).get("targets_bed"),  # runs frozen before config schema 4
        cfg.get("assay", {}).get("primers_tsv"),
        cfg.get("assay", {}).get("library_catalog"),
        *[(plan.get("panel") or {}).get(k) for k in ("targets_bed", "primers_tsv", "probes_bed")],
        cfg.get("analyses", {}).get("resistance_catalog"),
    ]:
        if value:
            inputs.append(root / value)
    inputs += [root / "config" / name for name in ("config.yaml", "quality.yaml", "compute.yaml")]
    protected = {p.resolve() for p in inputs}
    eligible: list[dict[str, Any]] = []
    skipped: list[str] = []
    if work.exists():
        for p in sorted(work.rglob("*")):
            if p.is_symlink():
                skipped.append(str(p.relative_to(run)))
                continue
            if not p.is_file():
                continue
            if any(part.is_symlink() for part in [p, *list(p.parents)[: len(p.relative_to(work).parts)]]):
                skipped.append(str(p.relative_to(run)))
                continue
            if p.resolve() in protected or not stat.S_ISREG(p.lstat().st_mode):
                skipped.append(str(p.relative_to(run)))
                continue
            eligible.append({"path": str(p.relative_to(run)), "bytes": p.stat().st_size})
    result: dict[str, Any] = {
        "schema": "pfvarx.cleanup/1",
        "run_id": run.name,
        "apply": apply,
        "files": eligible,
        "bytes_eligible": sum(p["bytes"] for p in eligible),
        "skipped": skipped,
        "deleted": [],
    }
    if apply:
        for item in eligible:
            p = run / item["path"]
            if p.is_symlink() or not p.resolve().is_relative_to(work.resolve()) or p.resolve() in protected:
                raise ValueError("Intermediate changed during cleanup; stopped")
            p.unlink()
            result["deleted"].append(item["path"])
        if work.is_dir():
            for p in sorted(work.rglob("*"), key=lambda p: len(p.parts), reverse=True):
                if p.is_dir() and not p.is_symlink() and not any(p.iterdir()):
                    p.rmdir()
        write_json(run / "cleanup.json", result)
    return result
