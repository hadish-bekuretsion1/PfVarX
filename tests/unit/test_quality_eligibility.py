from pathlib import Path

from pfvarx.eligibility.engine import SampleState, apply_strict, assess, internal_prerequisites, load_registry
from pfvarx.quality.policy import evaluate_gate, evaluate_sample, load_policy, stage_permitted

REPO = Path(__file__).resolve().parents[2]
PARAMS = {"da006_min_monoclonal": 4, "da005_min_group_size": 3, "da003_min_group_size": 3}
CTX = {"run_id": "r", "project_id": "p", "policy_version": "pfvarx_qc_v1", "population_column": "population",
       "run_policy": "CONTINUE_INDEPENDENT"}


def test_unset_threshold_is_not_evaluated_never_pass() -> None:
    gate = {"id": "x", "gate": "G1", "stage": "s", "metric": "m", "unit": "u", "direction": "min", "fail": None,
            "warn": None, "threshold_source": "UNSET_BLOCKING", "designs": ["WGS"], "mandatory": True}
    assert evaluate_gate(gate, {"m": 1.0}, "WGS", "1")["status"] == "NOT_EVALUATED"
    assert evaluate_gate(gate, {"m": 1.0}, "TARGETED_UMI", "1")["status"] == "NOT_APPLICABLE"


def test_fail_blocks_and_warning_permission() -> None:
    policy = load_policy(REPO, "pfvarx_qc_v1")
    bad = evaluate_sample(policy, {"mapped_fraction": 0.2}, "WGS")
    assert bad["verdict"] == "FAIL" and not stage_permitted(bad, "L2")
    warn = evaluate_sample(policy, {"mapped_fraction": 0.9, "median_depth_core": 15, "duplicate_fraction": 0.1,
                                    "callable_fraction_core": 0.95, "reads_retained_fraction": 0.95}, "WGS")
    assert warn["verdict"] == "WARNING" and warn["label"].startswith("Analysed with warning")
    assert stage_permitted(warn, "L2") and stage_permitted(warn, "DA002") and not stage_permitted(warn, "DA006")


def _states(clon: list[str | None]) -> dict[str, SampleState]:
    ok = {"verdict": "PASS", "permitted": {}, "uncalibrated_gates": []}
    return {f"S{i}": SampleState(f"S{i}", "WGS", ok, c, {"population": "A" if i % 2 else "B"})
            for i, c in enumerate(clon)}


def _da006(clon: list[str | None], policy: str = "CONTINUE_INDEPENDENT") -> dict:
    reg = load_registry(REPO)
    res = assess(["DA006", "DA002"], {}, reg, _states(clon), PARAMS, {**CTX, "run_policy": policy})
    return {c["analysis_id"]: c for c in apply_strict(res["cohort"], policy)} | {"_matrix": res["matrix"]}


def test_da006_eligible() -> None:
    assert _da006(["MONOCLONAL"] * 5)["DA006"]["decision"] == "EXECUTE"


def test_da006_partially_eligible_excludes_polyclonal_and_indeterminate() -> None:
    r = _da006(["MONOCLONAL"] * 4 + ["POLYCLONAL", "INDETERMINATE"])
    assert r["DA006"]["decision"] == "EXECUTE" and set(r["DA006"]["excluded"]) == {"S4", "S5"}
    reasons = {m["sample_id"]: m["reason_codes"] for m in r["_matrix"] if m["analysis_id"] == "DA006"}
    assert reasons["S4"] == ["CLONALITY_POLYCLONAL"] and reasons["S5"] == ["CLONALITY_INDETERMINATE"]


def test_da006_all_polyclonal_or_indeterminate_is_skipped_not_failed() -> None:
    r = _da006(["POLYCLONAL"] * 3 + ["INDETERMINATE"] * 3)
    assert r["DA006"]["decision"] == "SKIPPED_INELIGIBLE"
    assert r["DA002"]["decision"] == "EXECUTE"  # independent analysis continues


def test_da006_waits_for_clonality_evidence() -> None:
    assert _da006([None] * 5)["DA006"]["decision"] == "PENDING_PREREQUISITE"


def test_strict_run_stops_unstarted_analyses() -> None:
    r = _da006(["POLYCLONAL"] * 6, policy="STRICT")
    assert r["DA006"]["decision"] == "SKIPPED_INELIGIBLE"
    assert r["DA002"]["decision"] == "BLOCKED_MISSING_PREREQUISITE" and "STRICT_RUN_STOP" in r["DA002"]["reason_codes"]


def test_internal_prerequisite_inserted() -> None:
    reg = load_registry(REPO)
    assert internal_prerequisites(["DA006"], reg, None) == ["DA001"]
    assert internal_prerequisites(["DA007"], reg, "DA002") == ["DA002"]
    assert internal_prerequisites(["DA002"], reg, None) == []


def test_failed_sample_isolated() -> None:
    reg = load_registry(REPO)
    states = _states(["MONOCLONAL"] * 5)
    states["S0"].processing = "BLOCKED_UPSTREAM"
    res = assess(["DA002"], {}, reg, states, PARAMS, CTX)
    rec = {m["sample_id"]: m for m in res["matrix"]}
    assert rec["S0"]["eligibility"] == "INELIGIBLE" and rec["S1"]["eligibility"] == "ELIGIBLE"
    assert res["cohort"][0]["decision"] == "EXECUTE" and "S0" in res["cohort"][0]["excluded"]
