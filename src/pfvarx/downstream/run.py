"""Layer 5 dispatcher: the eligibility checkpoint and execution of one analysis.

    python -m pfvarx.downstream.run --eligibility --run-dir RUN --root ROOT
    python -m pfvarx.downstream.run --da DA006 --run-dir RUN --root ROOT

A scientific exclusion ends in a terminal status manifest, never in a non-zero exit;
a non-zero exit is reserved for genuine execution failures.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

from pfvarx.constants import Terminal
from pfvarx.downstream.common import Context, load_context, write_manifest
from pfvarx.eligibility.engine import SampleState, apply_strict, assess, load_registry
from pfvarx.io import read_tsv, write_json, write_tsv


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
    states: dict[str, SampleState] = {}
    for sid in sorted({s["sample_id"] for s in ctx.plan["samples"]}):
        try:
            gate = ctx.gate(sid)
        except FileNotFoundError:
            gate = {"verdict": "NOT_EVALUATED", "permitted": {}}
        processing = "COMPLETED" if sid in ctx.samples else "BLOCKED_UPSTREAM"
        states[sid] = SampleState(sample_id=sid, design=ctx.eff["config"]["assay"]["design"], qc=gate,
                                  clonality=(clon or {}).get(sid) if clon is not None else None,
                                  metadata=ctx.metadata.get(sid, {}), processing=processing)
    scales = ctx.eff["config"]["analyses"]["scales"]
    run_ctx = {"run_id": ctx.run.name, "project_id": ctx.eff["config"]["project_id"],
               "policy_version": ctx.eff["quality"]["policy_id"],
               "population_column": ctx.eff["config"]["analyses"]["population_column"],
               "run_policy": ctx.eff["config"]["run_policy"]}
    das = [d for d in requested + internal if d != "DA001"]
    res = assess(das, scales, registry, states, ctx.params, run_ctx)
    if 4 not in ctx.plan["layers"]:
        for decision in res["cohort"]:
            if decision["analysis_id"] in {"DA002", "DA007"}:
                reason = "ANNOTATION_NOT_CONFIGURED" if decision["analysis_id"] == "DA002" else "SOURCE_OUTCOME_INCOMPLETE"
                decision.update(decision=Terminal.BLOCKED_MISSING_PREREQUISITE.value,
                                reason_codes=[*decision["reason_codes"], reason],
                                excluded=sorted(states), eligible=[], conditional=[])
        for row in res["matrix"]:
            if row["analysis_id"] in {"DA002", "DA007"}:
                row["eligibility"] = "INELIGIBLE"
                row["failure_class"] = "SCIENTIFIC_EXCLUSION"
                row["reason_codes"].append(
                    "ANNOTATION_NOT_CONFIGURED" if row["analysis_id"] == "DA002" else "SOURCE_OUTCOME_INCOMPLETE"
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
                                               "internal_prerequisites": internal})


def run_da(ctx: Context, da: str, threads: int) -> None:
    import importlib

    modules = {"DA002": "resistance", "DA003": "diversity", "DA001": "coi", "DA004": "relatedness",
               "DA005": "structure", "DA006": "phylogeny", "DA007": "spatial"}
    try:
        importlib.import_module("pfvarx.downstream." + modules[da]).run(ctx, threads=threads)
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
    ap.add_argument("--threads", type=int, default=1)
    a = ap.parse_args(argv)
    ctx = load_context(Path(a.root), Path(a.run_dir), with_matrix=not a.eligibility)
    if a.eligibility:
        run_eligibility(ctx)
    elif a.da:
        run_da(ctx, a.da, a.threads)
    else:
        ap.error("need --eligibility or --da")
    sys.stdout.write(json.dumps({"status": "ok"}) + "\n")


if __name__ == "__main__":
    main()
