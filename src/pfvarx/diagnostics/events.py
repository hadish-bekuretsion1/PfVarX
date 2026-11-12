"""Structured diagnostic events (append-only JSONL) with redaction.

A record never asserts a root cause from an exit code alone: ``cause_confirmed`` is
set only by the caller when evidence establishes it; otherwise causes are listed as
``causes_suspected`` (hypotheses).
"""

from __future__ import annotations

import datetime as dt
import json
import os
import re
import uuid
from pathlib import Path
from typing import Any

from pfvarx.constants import FailureClass, Severity

SECRET = re.compile(r"(?i)(token|password|secret|api[_-]?key)\s*[=:]\s*\S+")
REMEDIES = {
    "RESOURCE_MEMORY_EXCEEDED": ["increase the validated memory allocation within compute.yaml limits",
                                 "reduce max_jobs", "use a smaller window_size_bp profile"],
    "RESOURCE_TIME_LIMIT_EXCEEDED": ["increase per_job.max_runtime_min", "reduce input size per job"],
    "RESOURCE_DISK_EXHAUSTED": ["free temporary disk", "point TMPDIR to a larger volume"],
    "TOOL_NOT_FOUND": ["run the command through 'pixi run' so the locked environment is active"],
    "FASTQ_TRUNCATED": ["re-transfer the file and verify its checksum"],
    "FASTQ_MATE_NAME_MISMATCH": ["check that R1/R2 files belong to the same lane and are in the same order"],
    "REFERENCE_MISMATCH": ["use the exact reference FASTA the alignment was made against"],
    "CLONALITY_POLYCLONAL": ["exclude from conventional phylogenetics; the sample remains eligible for "
                             "polyclonal-compatible analyses"],
    "MONOCLONAL_COUNT_BELOW_MINIMUM": ["add monoclonal samples; no tree is produced below the minimum"],
    "QC_GATE_FAIL": ["inspect the metric in the QC matrix; re-sequence or re-prepare the library if confirmed"],
}


def redact(text: str) -> str:
    return SECRET.sub(lambda m: f"{m.group(1)}=<REDACTED>", text)


def event(*, run_id: str, layer: str, stage: str, status: str, severity: Severity, code: str,
          failure_class: FailureClass | None, sample_id: str | None = None,
          analysis_id: str | None = None, rule: str | None = None, tool: str | None = None,
          observed: Any = None, expected: Any = None, criterion: str | None = None,
          inputs: list[str] | None = None, outputs: list[str] | None = None,
          log_paths: list[str] | None = None, exit_code: int | None = None,
          cause_confirmed: str | None = None, causes_suspected: list[str] | None = None,
          affected_downstream: list[str] | None = None, interpretation_impact: str | None = None,
          retry_history: list[dict[str, Any]] | None = None, project_id: str | None = None
          ) -> dict[str, Any]:
    return {
        "event_id": uuid.uuid4().hex, "timestamp": dt.datetime.now(dt.UTC).isoformat(),
        "project_id": project_id, "run_id": run_id, "sample_id": sample_id,
        "analysis_id": analysis_id, "layer": layer, "stage": stage, "rule": rule, "tool": tool,
        "status": status, "severity": severity.value, "code": code,
        "failure_class": failure_class.value if failure_class else None,
        "observed": observed, "expected": expected, "criterion": criterion,
        "inputs": inputs or [], "outputs": outputs or [], "log_paths": log_paths or [],
        "exit_code": exit_code, "cause_confirmed": cause_confirmed,
        "causes_suspected": causes_suspected or ([] if cause_confirmed else
                                                 ["not established from available evidence"]),
        "remedies": REMEDIES.get(code, []), "affected_downstream": affected_downstream or [],
        "interpretation_impact": interpretation_impact, "retry_history": retry_history or [],
    }


def append(path: Path, record: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    line = redact(json.dumps(record, default=str, sort_keys=True))
    with path.open("a", encoding="utf-8") as handle:
        handle.write(line + "\n")
        handle.flush()
        os.fsync(handle.fileno())


def read_events(path: Path) -> list[dict[str, Any]]:
    if not path.is_file():
        return []
    return [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines() if line]
