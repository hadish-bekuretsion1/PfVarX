"""Layer 5 dispatcher: the eligibility checkpoint and execution of one analysis.

    python -m pfvarx.downstream.run --eligibility --run-dir RUN --root ROOT
    python -m pfvarx.downstream.run --da DA006 --run-dir RUN --root ROOT
    python -m pfvarx.downstream.run --cohort-status [--empty] --run-dir RUN --root ROOT

A scientific exclusion ends in a terminal status manifest, never in a non-zero exit;
a non-zero exit is reserved for genuine execution failures.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

from pfvarx.constants import Reason, Severity, Terminal
from pfvarx.diagnostics.events import append, event
from pfvarx.downstream.common import Context, load_context, write_manifest
from pfvarx.eligibility.engine import SampleState, apply_strict, assess, load_registry
from pfvarx.io import read_tsv, write_json, write_tsv

COHORT_STATUS_SCHEMA = "pfvarx.cohort_status/1"


def write_cohort_status(ctx: Context, empty: bool) -> None:
    """Record the Layer 3 outcome. Without an admitted sample, skip Layers 3-5 explicitly.

    An empty cohort is a scientific outcome (every sample was blocked by input validation,
    its QC gate or a run-time failure), not an execution failure: every requested analysis
    receives a SKIPPED_INELIGIBLE manifest so reports and denominators stay complete.
    """
    registered = sorted({s["sample_id"] for s in ctx.plan["samples"]}
                        | set(ctx.plan.get("samples_failed_validation", [])))
    blocked: dict[str, str] = {}
    for sid in registered:
        if sid in ctx.plan.get("samples_failed_validation", []):
            blocked[sid] = "INPUT_VALIDATION_FAILED"
            continue
        try:
            gate = ctx.gate(sid)
        except FileNotFoundError:
            blocked[sid] = "QC_NOT_EVALUATED"
            continue
        if not gate.get("permitted", {}).get("L3"):
            blocked[sid] = Reason.QC_GATE_FAIL.value if gate.get("verdict") == "FAIL" else Reason.QC_GATE_WARNING.value
    members = [] if empty else ctx.samples
    status = "SKIPPED_NO_ELIGIBLE_SAMPLES" if empty else "COMPLETED"
    skipped: list[str] = []
    layers = ctx.plan["layers"]
    if empty and 5 in layers:
        registry = load_registry(ctx.root)
        das = list(dict.fromkeys([*ctx.plan.get("analyses_requested", []), *ctx.plan.get("analyses_internal", [])]))
        for da in das:
            write_manifest(ctx, da, Terminal.SKIPPED_INELIGIBLE.value,
                           reason_codes=[Reason.NO_ELIGIBLE_SAMPLES.value, Reason.NO_COHORT.value],
                           requested=[], eligible=[], excluded={}, analysed=[], outputs=[],
                           caveats=["No sample passed input validation and QC into the cohort; "
                                    "no scientific output was produced for this analysis."],
                           method=registry["analyses"][da])
            skipped.append(da)
    write_json(ctx.results / "cohort" / "cohort_status.json", {
        "schema": COHORT_STATUS_SCHEMA, "status": status, "run_id": ctx.run.name,
        "cohort_samples": members, "blocked_samples": blocked,
        "layers_skipped": [layer for layer in (3, 4, 5) if layer in layers] if empty else [],
        "analyses_skipped": skipped,
        "recommended_action": ("Review the per-sample QC decisions in results/qc/samples/; no sample "
                               "qualified for cohort analysis." if empty else None),
    })


def clonality_map(ctx: Context) -> dict[str, str] | None:
    path = ctx.results / "downstream" / "DA001" / "clonality.tsv"
    if not path.is_file():
        return None
    return {r["sample_id"]: r["clonality"] for r in read_tsv(path)}


def run_eligibility(ctx: Context) -> None:
    registry = load_registry(ctx.root)
    requested = list(ctx.plan.get("analyses_requested", []))
    internal = list(ctx.plan.get("analyses_internal", []))
    clon = clonality_map(ctx)
    cohort_qc_path = ctx.results / "cohort" / "cohort_qc.json"
    cohort_qc = json.loads(cohort_qc_path.read_text())["samples"] if cohort_qc_path.is_file() else {}
    states: dict[str, SampleState] = {}
    for sid in sorted({s["sample_id"] for s in ctx.plan["samples"]}):
        try:
            gate = ctx.gate(sid)
        except FileNotFoundError:
            gate = {"verdict": "NOT_EVALUATED", "permitted": {}}
        processing = "COMPLETED" if sid in ctx.samples else "BLOCKED_UPSTREAM"
        states[sid] = SampleState(sample_id=sid, design=ctx.design.value, qc=gate,
                                  clonality=(clon or {}).get(sid) if clon is not None else None,
                                  metadata=ctx.metadata.get(sid, {}), processing=processing,
                                  cohort_qc=cohort_qc.get(sid))
    scales = ctx.eff["config"]["analyses"]["scales"]
    from pfvarx.eligibility.territory import callable_markers, describe

    territory = describe(ctx)
    run_ctx = {"run_id": ctx.run.name, "project_id": ctx.eff["config"]["project_id"],
               "territory": territory,
               "callable_markers": {sid: callable_markers(ctx, sid, territory) for sid in ctx.samples},
               "policy_version": ctx.eff["quality"]["policy_id"],
               "population_column": ctx.eff["config"]["analyses"]["population_column"],
               "run_policy": ctx.eff["config"]["run_policy"]}
    das = [d for d in requested + internal if d != "DA001"]
    res = assess(das, scales, registry, states, ctx.params, run_ctx)
    if 4 not in ctx.plan["layers"]:
        for decision in res["cohort"]:
            if decision["analysis_id"] in {"DA002", "DA007", "DA008", "DA009"}:
                # DA002, DA008 and DA009 place their markers and genes with the GFF; DA007 needs DA002.
                reason = "SOURCE_OUTCOME_INCOMPLETE" if decision["analysis_id"] == "DA007" else "ANNOTATION_NOT_CONFIGURED"
                decision.update(decision=Terminal.BLOCKED_MISSING_PREREQUISITE.value,
                                reason_codes=[*decision["reason_codes"], reason],
                                excluded=sorted(states), eligible=[], conditional=[])
        for row in res["matrix"]:
            if row["analysis_id"] in {"DA002", "DA007", "DA008", "DA009"}:
                row["eligibility"] = "INELIGIBLE"
                row["failure_class"] = "SCIENTIFIC_EXCLUSION"
                row["reason_codes"].append(
                    "SOURCE_OUTCOME_INCOMPLETE" if row["analysis_id"] == "DA007" else "ANNOTATION_NOT_CONFIGURED"
                )
    cohort = apply_strict(res["cohort"], run_ctx["run_policy"])
    analyses: dict[str, Any] = {}
    for c in cohort:
        a = analyses.setdefault(c["analysis_id"], {"scales": {}, "decision": None})
        a["scales"][c["scale"]] = c
    for da, a in analyses.items():
        decs = [c["decision"] for c in a["scales"].values()]
        a["decision"] = "EXECUTE" if "EXECUTE" in decs else (
            "PENDING_PREREQUISITE" if "PENDING_PREREQUISITE" in decs else decs[0])
        if a["decision"] == "PENDING_PREREQUISITE":
            # Clonality evidence should exist here (DA001 ran before this checkpoint); if it
            # does not, the prerequisite is missing.
            a["decision"] = Terminal.BLOCKED_MISSING_PREREQUISITE.value
        if a["decision"] not in ("EXECUTE",):
            first = next(iter(a["scales"].values()))
            excluded: dict[str, list[str]] = {}
            for r in res["matrix"]:
                if r["analysis_id"] == da and r["eligibility"] != "ELIGIBLE":
                    excluded.setdefault(",".join(r["reason_codes"]) or r["eligibility"], []).append(r["sample_id"])
            write_manifest(ctx, da, a["decision"], reason_codes=sorted({x for c in a["scales"].values()
                                                                          for x in c["reason_codes"]} | {x for r in res["matrix"]
                                                                            if r["analysis_id"] == da for x in r["reason_codes"]}),
                           requested=sorted(states), eligible=first["eligible"], excluded=excluded,
                           analysed=[], outputs=[], caveats=["No scientific output was produced for this analysis."],
                           method=registry["analyses"][da], scale=",".join(a["scales"]))
    if "DA001" in requested + internal:
        analyses["DA001"] = {"decision": "EXECUTED_AS_PREREQUISITE" if "DA001" not in requested else "EXECUTE",
                             "scales": {}, "reported": "DA001" in requested}
    out = ctx.results / "eligibility"
    cols = ["run_id", "project_id", "sample_id", "analysis_id", "scale", "eligibility", "reason_codes",
            "failure_class", "rule_id", "rule_version", "policy_version", "caveats", "decided_at"]
    rows = [{**r, "reason_codes": ",".join(r["reason_codes"]), "caveats": " | ".join(r["caveats"])}
            for r in res["matrix"]]
    write_tsv(out / "sample_analysis_eligibility.tsv", rows, cols)
    write_json(out / "sample_analysis_eligibility.json", res["matrix"])
    write_json(out / "cohort_decisions.json", {"schema": "pfvarx.cohort_decisions/1", "analyses": analyses,
                                               "run_policy": run_ctx["run_policy"],
                                               "internal_prerequisites": internal,
                                               "cohort_scale": cohort_scale(ctx, das)})


def analysis_module(da: str) -> str:
    """Module implementing an analysis."""
    modules = {"DA002": "resistance", "DA003": "diversity", "DA001": "coi", "DA004": "relatedness",
               "DA005": "structure", "DA006": "phylogeny", "DA007": "spatial", "DA008": "cnv",
               "DA009": "hrp"}
    return modules[da]


def needs_matrix(da: str) -> bool:
    """Whether an analysis reads the cohort evidence matrix (``USES_COHORT_MATRIX``).

    A few hundred samples over tens of thousands of sites is gigabytes of evidence, and the
    analyses that work from the alignments (DA002, DA007, DA008, DA009) never read it. An
    analysis that does not declare itself gets the matrix, so a new one cannot lose its data
    by omission.
    """
    import importlib

    module = importlib.import_module("pfvarx.downstream." + analysis_module(da))
    return bool(getattr(module, "USES_COHORT_MATRIX", True))


def cohort_scale(ctx: Context, das: list[str]) -> dict[str, Any]:
    """How much memory this cohort's evidence will need, and a warning while there is still time.

    The analyses that read the evidence matrix hold one cell per sample per cohort site. The
    cost per cell is measured (``planning.cohort_matrix_bytes_per_cell``), and the sample and
    site counts are known here, so a cohort too large for ``per_job.max_mem_mb`` can be said
    before the analyses start rather than discovered as an out-of-memory kill. It is a warning,
    not a block: the estimate is approximate and an analysis may still succeed.
    """
    from pfvarx.downstream.common import cohort_sites, matrix_memory_mb
    from pfvarx.planning import resources

    samples, sites = len(ctx.samples), cohort_sites(ctx.run)
    try:
        per_cell = int(resources.planning(ctx.root)["cohort_matrix_bytes_per_cell"])
    except (OSError, KeyError):  # a workspace frozen before the factor existed
        per_cell = 32
    needed = matrix_memory_mb(samples, sites, per_cell)
    limit = int(ctx.eff["compute"]["per_job"]["max_mem_mb"])
    reading = sorted(da for da in das if needs_matrix(da))
    scale: dict[str, Any] = {"samples": samples, "cohort_sites": sites, "cells": samples * sites,
                             "bytes_per_cell": per_cell, "evidence_matrix_mb": needed,
                             "per_job_max_mem_mb": limit, "analyses_reading_the_matrix": reading,
                             "warning": None}
    # Half the job's memory: an analysis needs working memory beyond the matrix it reads.
    if reading and needed > limit / 2:
        scale["warning"] = {
            "reason": "RESOURCE_MEMORY_TIGHT",
            "detail": (f"{samples} samples x {sites} cohort sites is about {needed} MB of evidence, and "
                       f"per_job.max_mem_mb is {limit} MB. {', '.join(reading)} read it and need working "
                       "memory besides. Raise per_job.max_mem_mb and max_mem_mb if the machine has it, "
                       "request fewer analyses, or analyse a smaller cohort."),
        }
        append(ctx.run / "logs" / "events.jsonl", event(
            run_id=ctx.run.name, layer="L5", stage="eligibility", status="WARNING",
            severity=Severity.WARNING, code="RESOURCE_MEMORY_TIGHT", failure_class=None,
            observed=f"{needed} MB", expected=f"<= {limit // 2} MB",
            criterion="cohort evidence matrix within half of per_job.max_mem_mb",
            interpretation_impact="none: this is a resource estimate, not a scientific decision",
            project_id=ctx.eff["config"]["project_id"]))
    return scale


def run_da(ctx: Context, da: str, threads: int) -> None:
    import importlib

    try:
        importlib.import_module("pfvarx.downstream." + analysis_module(da)).run(ctx, threads=threads)
    except (Exception, SystemExit) as exc:
        if ctx.eff["config"]["run_policy"] != "CONTINUE_INDEPENDENT":
            raise
        import traceback

        traceback.print_exc()
        if da == "DA001":
            # Release the eligibility checkpoint with no assessed clonality.
            write_tsv(ctx.out(da) / "clonality.tsv", [], ["sample_id", "clonality"])
        write_manifest(
            ctx, da, Terminal.FAILED_EXECUTION.value,
            reason_codes=["METHOD_EXECUTION_FAILED"], requested=ctx.samples,
            eligible=[], excluded={}, analysed=[], outputs=[],
            caveats=[f"{type(exc).__name__}: {exc}"],
            method=load_registry(ctx.root)["analyses"][da],
        )



def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--run-dir", required=True)
    ap.add_argument("--root", required=True)
    ap.add_argument("--da")
    ap.add_argument("--eligibility", action="store_true")
    ap.add_argument("--cohort-status", action="store_true")
    ap.add_argument("--empty", action="store_true", help="with --cohort-status: no sample entered Layer 3")
    ap.add_argument("--threads", type=int, default=1)
    a = ap.parse_args(argv)
    with_matrix = bool(a.da) and needs_matrix(a.da)
    ctx = load_context(Path(a.root), Path(a.run_dir), with_matrix=with_matrix)
    if a.cohort_status:
        write_cohort_status(ctx, a.empty)
    elif a.eligibility:
        run_eligibility(ctx)
    elif a.da:
        run_da(ctx, a.da, a.threads)
    else:
        ap.error("need --eligibility or --da")
    sys.stdout.write(json.dumps({"status": "ok"}) + "\n")


if __name__ == "__main__":
    main()
