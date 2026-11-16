"""Evaluate versioned QC gates. A gate without a threshold is NOT_EVALUATED, never PASS."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml

from pfvarx.constants import QC, Design, ParameterSource, Reason


def load_policy(root: Path, policy_id: str) -> dict[str, Any]:
    path = root / "config" / "registry" / f"qc_policy_{policy_id}.yaml"
    if not path.is_file():
        raise ValueError(f"QC policy {policy_id!r} is not registered ({path} missing)")
    policy = yaml.safe_load(path.read_text())
    if policy.get("policy_id") != policy_id:
        raise ValueError("QC policy file id does not match its name")
    return policy


def evaluate_gate(gate: dict[str, Any], metrics: dict[str, Any], design: str,
                  policy_version: str) -> dict[str, Any]:
    rec = {"gate_id": gate["id"], "gate": gate["gate"], "stage": gate["stage"],
           "metric": gate["metric"], "unit": gate["unit"], "observed": metrics.get(gate["metric"]),
           "fail_threshold": gate.get("fail"), "warn_threshold": gate.get("warn"),
           "direction": gate["direction"], "threshold_source": gate["threshold_source"],
           "policy_version": policy_version, "mandatory": gate["mandatory"],
           "warning_permitted_for": gate.get("warning_permitted_for", []),
           "calibrated": gate["threshold_source"] in (ParameterSource.EMPIRICALLY_CALIBRATED,
                                                       ParameterSource.ASSAY_VALIDATED,
                                                       ParameterSource.FRAMEWORK_POLICY,
                                                       ParameterSource.LITERATURE_SUPPORTED),
           "reason": None}
    if Design.parse(design) not in [Design.parse(d) for d in gate["designs"]]:
        rec["status"] = QC.NOT_APPLICABLE
        return rec
    obs = rec["observed"]
    if gate["direction"] == "descriptive":
        rec["status"] = QC.NOT_APPLICABLE if obs is None else QC.PASS
        rec["reason"] = "DESCRIPTIVE_ONLY"
        return rec
    if obs is None:
        rec["status"] = QC.NOT_EVALUATED
        rec["reason"] = "METRIC_UNAVAILABLE"
        return rec
    if gate["threshold_source"] == ParameterSource.UNSET_BLOCKING or gate.get("fail") is None:
        rec["status"] = QC.NOT_EVALUATED
        rec["reason"] = Reason.QC_GATE_UNCALIBRATED.value
        return rec
    value = float(obs)
    lo = gate["direction"] == "min"
    fail = float(gate["fail"])
    warn = None if gate.get("warn") is None else float(gate["warn"])
    if (lo and value < fail) or (not lo and value > fail):
        rec["status"], rec["reason"] = QC.FAIL, Reason.QC_GATE_FAIL.value
    elif warn is not None and ((lo and value < warn) or (not lo and value > warn)):
        rec["status"], rec["reason"] = QC.WARNING, Reason.QC_GATE_WARNING.value
    else:
        rec["status"] = QC.PASS
    return rec


def evaluate_sample(policy: dict[str, Any], metrics: dict[str, Any], design: str,
                    stages: set[str] | None = None) -> dict[str, Any]:
    gates = [g for g in policy["gates"] if stages is None or g["gate"] in stages]
    decisions = [evaluate_gate(g, metrics, design, policy["version"]) for g in gates]
    mandatory_fail = [d for d in decisions if d["status"] == QC.FAIL and d["mandatory"]]
    warnings = [d for d in decisions if d["status"] == QC.WARNING
                or (d["status"] == QC.FAIL and not d["mandatory"])]
    if mandatory_fail:
        verdict = QC.FAIL
    elif warnings:
        verdict = QC.WARNING
    elif any(d["status"] == QC.PASS for d in decisions):
        verdict = QC.PASS
    else:
        verdict = QC.NOT_EVALUATED
    permitted: set[str] | None = None
    for d in warnings:
        allowed = set(d["warning_permitted_for"])
        permitted = allowed if permitted is None else permitted & allowed
    uncalibrated = [d["gate_id"] for d in decisions
                    if d["status"] in (QC.PASS, QC.WARNING, QC.FAIL) and not d["calibrated"]]
    return {"verdict": verdict.value, "decisions": decisions,
            "warning_permitted_for": sorted(permitted) if permitted is not None else None,
            "uncalibrated_gates": uncalibrated,
            "label": ("Analysed with warning — interpret affected results with caution"
                      if verdict == QC.WARNING else None)}


def stage_permitted(qc: dict[str, Any], stage: str) -> bool:
    """May a sample with this QC record enter *stage* (e.g. 'L2', 'DA002')?"""
    if qc["verdict"] == QC.FAIL:
        return False
    if qc["verdict"] == QC.WARNING:
        allowed = qc.get("warning_permitted_for") or []
        return stage in allowed or stage.split(".")[0] in allowed
    return True
