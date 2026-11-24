"""Atomic run creation with NEW, RESUME and EXTEND modes (L0.6).

* NEW     - fresh run directory, frozen configuration and snapshot.
* RESUME  - continue an interrupted run; refused unless the effective configuration
            hash and every run-defining file hash are identical.
* EXTEND  - new run that may reuse per-sample artefacts of a previous run, but only
            when the reuse key (reference, QC policy, layer 1-2 parameters and the
            code that produced them) is identical. Anything else is recomputed.
"""

from __future__ import annotations

import datetime as dt
import json
import os
import re
import shutil
import socket
import uuid
from pathlib import Path
from typing import Any

from pfvarx.configuration.loader import canonical_json, sha256_text
from pfvarx.io import read_json, sha256, write_json
from pfvarx.provenance import snapshot as snap

RUN_ID = re.compile(r"^\d{8}T\d{6}Z-[A-Za-z0-9_-]+-[0-9a-f]{8}$")
REUSE_CODE = ("src/pfvarx/caller/", "src/pfvarx/caller_adapter.py", "src/pfvarx/genomic_evidence.py", "workflow/rules/l1_", "workflow/rules/l2_", "workflow/scripts/l1_",
              "workflow/scripts/l2_")


class ReuseError(ValueError):
    pass


def reuse_key(eff: dict[str, Any], snapshot: dict[str, Any]) -> str:
    code = {k: v for k, v in snapshot["files"].items() if k.startswith(REUSE_CODE)}
    params = {k: v for k, v in eff["parameters"].items() if k.startswith(("layer0.", "layer1.", "layer2."))}
    key = {"reference": (snapshot.get("prepared") or {}).get("reference_sha256"),
           "assay": eff["config"]["assay"], "policy": eff["quality"]["policy_id"],
           "policy_file": snapshot["files"].get(f"config/registry/qc_policy_{eff['quality']['policy_id']}.yaml"),
           "params": params, "code": code, "caller": eff["config"].get("caller")}
    return sha256_text(canonical_json(key))


def _new_id(project: str) -> str:
    return (dt.datetime.now(dt.UTC).strftime("%Y%m%dT%H%M%SZ")
            + f"-{project}-{uuid.uuid4().hex[:8]}")


def create_run(root: Path, eff: dict[str, Any], plan: dict[str, Any], mode: str = "NEW",
               previous: str | None = None) -> Path:
    root = root.resolve()
    runs = root / "runs"
    runs.mkdir(exist_ok=True)
    current = snap.build(root)
    if mode == "RESUME":
        if previous is None:
            raise ReuseError("RESUME requires --run-id")
        run = runs / previous
        if not RUN_ID.fullmatch(previous) or not run.is_dir():
            raise ReuseError(f"run {previous} not found")
        frozen = read_json(run / "frozen_config.json")
        old = read_json(run / "snapshot.json")
        if frozen["hash"] != eff["hash"]:
            raise ReuseError("RUN_REUSE_INCOMPATIBLE: effective configuration changed; use NEW")
        changed = [k for k, v in old["files"].items() if current["files"].get(k) != v]
        if changed:
            raise ReuseError(f"RUN_REUSE_INCOMPATIBLE: run-defining files changed: {changed[:10]}")
        _attempt(run, "RESUME")
        return run
    run_id = _new_id(eff["config"]["project_id"])
    staging = runs / f".staging-{run_id}"
    staging.mkdir()
    try:
        for rel in current["files"]:
            dst = staging / "frozen" / rel
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(root / rel, dst)
            if sha256(dst) != current["files"][rel]:
                raise ValueError(f"file changed while freezing: {rel}")
        write_json(staging / "snapshot.json", current)
        write_json(staging / "frozen_config.json", eff)
        write_json(staging / "frozen_plan.json", plan)
        write_json(staging / "status.json", {"run_id": run_id, "mode": mode, "status": "PENDING",
                                             "previous": previous,
                                             "reuse_key": reuse_key(eff, current)})
        (staging / "results").mkdir()
        (staging / "logs").mkdir()
        os.rename(staging, runs / run_id)
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    run = runs / run_id
    _attempt(run, mode)
    if mode == "EXTEND":
        if previous is None:
            raise ReuseError("EXTEND requires --run-id of the run to extend")
        link_reusable(run, runs / previous)
    return run


def _attempt(run: Path, mode: str) -> None:
    attempts = run / "attempts.jsonl"
    with attempts.open("a") as h:
        h.write(json.dumps({"mode": mode, "at": dt.datetime.now(dt.UTC).isoformat()}) + "\n")


def link_reusable(run: Path, prev: Path) -> list[str]:
    """Hard-link per-sample Layer 1-2 artefacts from *prev* when the reuse keys match."""
    new_key = read_json(run / "status.json")["reuse_key"]
    old_key = read_json(prev / "status.json")["reuse_key"]
    if new_key != old_key:
        raise ReuseError("RUN_REUSE_INCOMPATIBLE: reference, policy, parameters or caller code differ")
    reused = []
    for sub in ("alignments", "qc/samples", "variants/samples"):
        src = prev / "results" / sub
        if not src.is_dir():
            continue
        for f in src.rglob("*"):
            if f.is_file():
                dst = run / "results" / f.relative_to(prev / "results")
                dst.parent.mkdir(parents=True, exist_ok=True)
                os.link(f, dst)
                reused.append(str(dst.relative_to(run)))
    write_json(run / "reused_artifacts.json", {"from_run": prev.name, "files": reused})
    return reused


def set_status(run: Path, status: str, **extra: Any) -> None:
    data = read_json(run / "status.json")
    data.update({"status": status, **extra,
                 "updated_utc": dt.datetime.now(dt.UTC).isoformat()})
    write_json(run / "status.json", data)


def process_identity() -> dict[str, Any]:
    """Who is running a run: recorded with RUNNING so a stopped run can be recognised later."""
    return {"pid": os.getpid(), "host": socket.gethostname()}


def _alive(pid: int) -> bool:
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True  # exists, owned by someone else
    cmdline = Path(f"/proc/{pid}/cmdline")
    if cmdline.is_file():  # Linux: the number may have been reused by another program
        return b"pfvarx" in cmdline.read_bytes()
    return True


def run_status(run: Path) -> dict[str, Any]:
    """status.json, with RUNNING replaced by INTERRUPTED when the run's process has ended.

    A run stopped by Ctrl+C, a kill or the out-of-memory killer, or by a restart of the
    computer, never writes its final status itself. The file is not changed here.
    """
    data = read_json(run / "status.json") if (run / "status.json").is_file() else {}
    if (data.get("status") == "RUNNING" and data.get("pid") and data.get("host") == socket.gethostname()
            and not _alive(int(data["pid"]))):
        data = {**data, "status": "INTERRUPTED", "recorded_status": "RUNNING",
                "interruption": "the run's process is no longer running (stopped with Ctrl+C, killed, out of "
                                "memory, or the computer restarted); resume it with "
                                f"'./pfvarx pipeline --mode RESUME --run-id {run.name}'"}
    return data


def verify(root: Path, run_id: str) -> dict[str, Any]:
    run = root / "runs" / run_id
    if not RUN_ID.fullmatch(run_id) or not run.is_dir():
        raise ValueError(f"run {run_id} not found")
    snapshot = read_json(run / "snapshot.json")
    tampered = [rel for rel, d in snapshot["files"].items()
                if not (run / "frozen" / rel).is_file() or sha256(run / "frozen" / rel) != d]
    return {"run_id": run_id, "frozen_integrity": "PASS" if not tampered else "FAIL",
            "tampered": tampered, "workspace_drift": snap.drift(root, snapshot),
            "status": read_json(run / "status.json")}


def latest(root: Path) -> Path | None:
    runs = sorted(p for p in (root / "runs").glob("*") if p.is_dir() and RUN_ID.fullmatch(p.name))
    return runs[-1] if runs else None
