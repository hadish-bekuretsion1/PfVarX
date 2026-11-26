"""Scientific eligibility engine (L0.5): sample-by-analysis matrix and cohort decisions.

Eligibility is a scientific decision recorded with evidence and reason codes. It is
never inferred from a software exception, and no user setting can bypass it.
"""

from __future__ import annotations

import datetime as dt
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

import yaml

from pfvarx.configuration.loader import frozen_param
from pfvarx.constants import (
    QC,
    Clonality,
    Design,
    Eligibility,
    FailureClass,
    Reason,
    Scale,
    Terminal,
)
from pfvarx.quality.policy import stage_permitted

RULE_VERSION = "1.1.0"  # 1.1.0: panel territory and callable-marker evidence


def load_registry(root: Path) -> dict[str, Any]:
    data = yaml.safe_load((root / "config" / "registry" / "downstream_analyses.yaml").read_text())
    return data


@dataclass
class SampleState:
    sample_id: str
    design: str
    qc: dict[str, Any]
    clonality: str | None = None  # MONOCLONAL | POLYCLONAL | INDETERMINATE | None (not assessed)
    metadata: dict[str, str] = field(default_factory=dict)
    processing: str = "COMPLETED"
    cohort_qc: dict[str, Any] | None = None  # G5 decision (Layer 3); None when not evaluated


@dataclass
class EligibilityRecord:
    run_id: str
    project_id: str
    sample_id: str
    analysis_id: str
    scale: str
    eligibility: str
    reason_codes: list[str]
    failure_class: str | None
    evidence: dict[str, Any]
    rule_id: str
    rule_version: str
    policy_version: str
    dependencies: list[str]
    caveats: list[str]
    decided_at: str


def _now() -> str:
    return dt.datetime.now(dt.UTC).isoformat()


def sample_eligibility(da: str, contract: dict[str, Any], s: SampleState, scale: str,
                       ctx: dict[str, Any]) -> EligibilityRecord:
    reasons: list[str] = []
    failure: str | None = None
    caveats: list[str] = []
    evidence: dict[str, Any] = {"design": s.design, "qc_verdict": s.qc.get("verdict"),
                                "clonality": s.clonality}
    status = Eligibility.ELIGIBLE
    if s.processing != "COMPLETED":
        status, failure = Eligibility.INELIGIBLE, FailureClass.QC_FAILURE
        reasons.append(Reason.QC_UPSTREAM_FAILED.value)
    elif not stage_permitted(s.qc, da):
        status, failure = Eligibility.INELIGIBLE, FailureClass.QC_FAILURE
        reasons.append(Reason.QC_GATE_FAIL.value if s.qc.get("verdict") == QC.FAIL
                       else Reason.QC_GATE_WARNING.value)
    elif s.cohort_qc is not None and not stage_permitted(s.cohort_qc, da):
        status, failure = Eligibility.INELIGIBLE, FailureClass.QC_FAILURE
        reasons.append(Reason.COHORT_QC_FAIL.value if s.cohort_qc.get("verdict") == QC.FAIL
                       else Reason.COHORT_QC_WARNING.value)
    if s.cohort_qc is not None:
        evidence["cohort_qc_verdict"] = s.cohort_qc.get("verdict")
        evidence["cohort_missing_fraction"] = (s.cohort_qc.get("metrics") or {}).get("cohort_missing_fraction")
    if Design.parse(s.design) not in [Design.parse(d) for d in contract["designs"]]:
        status, failure = Eligibility.UNSUPPORTED, FailureClass.SCIENTIFIC_EXCLUSION
        reasons.append(Reason.DESIGN_NOT_SUPPORTED.value)
    needs = contract.get("needs") or {}
    if needs.get("callable_markers"):
        # Callable evidence: the sample must cover whole catalogue codons in its callable regions.
        found = (ctx.get("callable_markers") or {}).get(s.sample_id)
        evidence["callable_markers"] = found
        if status == Eligibility.ELIGIBLE and found is not None and found < int(needs["callable_markers"]):
            status, failure = Eligibility.INELIGIBLE, FailureClass.SCIENTIFIC_EXCLUSION
            reasons.append(Reason.NO_CALLABLE_MARKERS.value)
    if status == Eligibility.ELIGIBLE and "DA001" in contract.get("prerequisites", []):
        if s.clonality is None:
            status = Eligibility.CONDITIONALLY_ELIGIBLE
            reasons.append(Reason.CLONALITY_EVIDENCE_MISSING.value)
        elif s.clonality not in contract["clonality"]:
            status, failure = Eligibility.INELIGIBLE, FailureClass.SCIENTIFIC_EXCLUSION
            reasons.append(Reason.CLONALITY_POLYCLONAL.value if s.clonality == Clonality.POLYCLONAL
                           else Reason.CLONALITY_INDETERMINATE.value)
    for col in contract.get("required_metadata", []):
        key = ctx.get("population_column") if col == "population" else col
        if status == Eligibility.ELIGIBLE and (not key or not s.metadata.get(key)):
            status, failure = Eligibility.INELIGIBLE, FailureClass.SCIENTIFIC_EXCLUSION
            reasons.append(Reason.GROUP_ASSIGNMENT_MISSING.value if col == "population"
                           else Reason.METADATA_DATE_INVALID.value if col == "collection_date"
                           else Reason.METADATA_GEOGRAPHY_INVALID.value)
    if s.qc.get("verdict") == QC.WARNING and status == Eligibility.ELIGIBLE:
        caveats.append(s.qc.get("label") or "analysed with warning")
    if s.qc.get("uncalibrated_gates"):
        caveats.append("QC thresholds are UNCALIBRATED_PROVISIONAL: "
                       + ",".join(s.qc["uncalibrated_gates"]))
    return EligibilityRecord(
        run_id=ctx["run_id"], project_id=ctx["project_id"], sample_id=s.sample_id,
        analysis_id=da, scale=scale, eligibility=status.value, reason_codes=reasons,
        failure_class=failure, evidence=evidence, rule_id=f"{da}.sample", rule_version=RULE_VERSION,
        policy_version=ctx["policy_version"], dependencies=list(contract.get("prerequisites", [])),
        caveats=caveats, decided_at=_now())


def cohort_decision(da: str, contract: dict[str, Any], records: list[EligibilityRecord],
                    scale: str, samples: dict[str, SampleState], params: dict[str, Any],
                    ctx: dict[str, Any]) -> dict[str, Any]:
    eligible = [r.sample_id for r in records if r.eligibility == Eligibility.ELIGIBLE]
    conditional = [r.sample_id for r in records if r.eligibility == Eligibility.CONDITIONALLY_ELIGIBLE]
    excluded = [r.sample_id for r in records if r.eligibility not in
                (Eligibility.ELIGIBLE, Eligibility.CONDITIONALLY_ELIGIBLE)]
    reasons: list[str] = []
    terminal: str | None = None
    if scale not in contract["scales"]:
        return {"analysis_id": da, "scale": scale, "decision": Terminal.SKIPPED_INELIGIBLE.value,
                "reason_codes": [Reason.SCALE_NOT_SUPPORTED.value], "eligible": [], "excluded": excluded,
                "conditional": conditional}
    if conditional:
        decision = "PENDING_PREREQUISITE"
    else:
        decision = "EXECUTE"
    # What was sequenced: the panel (or core genome) must be able to answer the question.
    territory = ctx.get("territory") or {}
    needs = contract.get("needs") or {}
    if territory and needs.get("catalogue_markers_in_territory") is not None and \
            territory.get("catalogue_markers_in_territory", 0) < int(needs["catalogue_markers_in_territory"]):
        terminal = Terminal.SKIPPED_INELIGIBLE.value
        reasons.append(Reason.PANEL_LACKS_MARKERS.value)
    loci_param = needs.get("independent_loci_param")
    if territory and loci_param and territory.get("independent_loci") is not None and \
            territory["independent_loci"] < int(params[loci_param]):
        terminal = Terminal.SKIPPED_INELIGIBLE.value
        reasons.append(Reason.TOO_FEW_INDEPENDENT_LOCI.value)
    if da == "DA006" and not conditional:
        minimum = int(params["da006_min_monoclonal"])
        if len(eligible) < minimum:
            terminal = Terminal.SKIPPED_INELIGIBLE.value
            reasons.append(Reason.MONOCLONAL_COUNT_BELOW_MINIMUM.value)
    if da == "DA004" and not conditional and len(eligible) < 2:
        terminal = Terminal.SKIPPED_INELIGIBLE.value
        reasons.append(Reason.NO_ELIGIBLE_PAIR.value)
    if da in ("DA003", "DA005") and scale in (Scale.POPULATION, Scale.BETWEEN_POPULATIONS) and not conditional:
        col = ctx.get("population_column")
        groups: dict[str, int] = {}
        for sid in eligible:
            g = samples[sid].metadata.get(col or "", "")
            if g:
                groups[g] = groups.get(g, 0) + 1
        minimum = int(params["da005_min_group_size" if da == "DA005" else "da003_min_group_size"])
        good = {g: n for g, n in groups.items() if n >= minimum}
        need = 2 if scale == Scale.BETWEEN_POPULATIONS or da == "DA005" else 1
        if len(good) < need:
            terminal = Terminal.SKIPPED_INELIGIBLE.value
            reasons.append(Reason.GROUP_TOO_SMALL.value)
    if da == "DA005" and not conditional and len(eligible) < int(frozen_param(params, "da005_min_samples")):
        terminal = Terminal.SKIPPED_INELIGIBLE.value
        reasons.append(Reason.GROUP_TOO_SMALL.value)
    if not eligible and not conditional and terminal is None:
        terminal = Terminal.SKIPPED_INELIGIBLE.value
        reasons.append("NO_ELIGIBLE_SAMPLES")
    if terminal:
        decision = terminal
    return {"analysis_id": da, "scale": scale, "decision": decision, "reason_codes": reasons,
            "eligible": eligible, "excluded": excluded, "conditional": conditional,
            "territory": {k: v for k, v in territory.items() if not k.startswith("_")} or None,
            "rule_id": f"{da}.cohort", "rule_version": RULE_VERSION}


def assess(da_list: list[str], scales: dict[str, list[str]], registry: dict[str, Any],
           samples: dict[str, SampleState], params: dict[str, Any], ctx: dict[str, Any]
           ) -> dict[str, Any]:
    matrix: list[dict[str, Any]] = []
    cohort: list[dict[str, Any]] = []
    for da in da_list:
        contract = registry["analyses"][da]
        for scale in scales.get(da) or [contract["scales"][0]]:
            recs = [sample_eligibility(da, contract, s, scale, ctx) for s in samples.values()]
            matrix += [asdict(r) for r in recs]
            cohort.append(cohort_decision(da, contract, recs, scale, samples, params, ctx))
    return {"schema": "pfvarx.eligibility/1", "matrix": matrix, "cohort": cohort,
            "run_policy": ctx.get("run_policy", "CONTINUE_INDEPENDENT")}


def internal_prerequisites(requested: list[str], registry: dict[str, Any],
                           spatial_source: str | None) -> list[str]:
    """Analyses that must run internally (e.g. DA001 for DA006) even if not requested."""
    need: list[str] = []
    for da in requested:
        for pre in registry["analyses"][da].get("prerequisites", []):
            pre = spatial_source if pre == "SOURCE" else pre
            if pre and pre not in requested and pre not in need:
                need.append(pre)
    return need


def apply_strict(cohort: list[dict[str, Any]], run_policy: str) -> list[dict[str, Any]]:
    """Under STRICT, any SKIPPED_INELIGIBLE stops analyses that have not yet executed."""
    if run_policy != "STRICT":
        return cohort
    if not any(c["decision"] == Terminal.SKIPPED_INELIGIBLE for c in cohort):
        return cohort
    out = []
    for c in cohort:
        c = dict(c)
        if c["decision"] in ("EXECUTE", "PENDING_PREREQUISITE"):
            c["decision"] = Terminal.BLOCKED_MISSING_PREREQUISITE.value
            c["reason_codes"] = c["reason_codes"] + [Reason.STRICT_RUN_STOP.value]
        out.append(c)
    return out
