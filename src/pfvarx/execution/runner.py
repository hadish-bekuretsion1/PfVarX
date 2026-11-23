"""Build the frozen plan, launch Snakemake, isolate failed samples and finalise the run (L0.8).

Execution is resource-bounded by Snakemake (--cores, --jobs, --resources mem_mb). A
first pass runs with --keep-going. If per-sample jobs failed, those samples are written
to ``runtime_exclusions.json`` (with their diagnostics) and a second pass re-evaluates
cohort membership and denominators without them. Reporting always runs, even when every
downstream analysis was skipped.
"""

from __future__ import annotations

import os
import subprocess
import tempfile
from pathlib import Path
from typing import Any, TypedDict

from pfvarx.caller_adapter import check_backend
from pfvarx.configuration.loader import effective_config, load_bundle
from pfvarx.constants import Design, FailureClass, Severity
from pfvarx.diagnostics.events import append, event
from pfvarx.eligibility.engine import internal_prerequisites, load_registry
from pfvarx.inputs.registry import validate_all
from pfvarx.io import read_json, write_json
from pfvarx.planning.plan import resource_plan
from pfvarx.references.prepare import verify_prepared
from pfvarx.runs import registry as runs


def build_plan(root: Path) -> dict[str, Any]:
    from pfvarx.installation import require_tools

    require_tools(root)  # readiness, dry run and pipeline need the tools of the saved configuration
    bundle = load_bundle(root)
    cfg = bundle.config
    if cfg.final_output.value != "FO001":
        check_backend(cfg.caller)
    validation = validate_all(root, bundle)
    prep = verify_prepared(root, bundle)
    samples: dict[str, dict[str, Any]] = {}
    for row in validation["registry"]:
        s = samples.setdefault(row["sample_id"], {"sample_id": row["sample_id"], "lanes": [],
                                                  "bam": row.get("bam", ""), "vcf": row.get("vcf", ""),
                                                  "assay": validation["libraries"].get(row["sample_id"], {})})
        s["lanes"].append({"lane": row["lane"], "fastq_1": row.get("fastq_1", ""), "fastq_2": row.get("fastq_2", ""),
                           "library_id": row["library_id"], "flowcell_type": row.get("flowcell_type", "")})
    failed_inputs = sorted(s for s, st in validation["samples"].items() if st == "FAIL")
    registry = load_registry(root)
    requested = list(cfg.analyses.requested)
    internal = internal_prerequisites(requested, registry, cfg.analyses.spatial_outcome_source)
    gb = 0.0
    for s in samples.values():
        for ln in s["lanes"]:
            for f in (ln["fastq_1"], ln["fastq_2"]):
                if f and (root / f).is_file():
                    gb += (root / f).stat().st_size / 1e9
        if s["bam"] and (root / s["bam"]).is_file():
            gb += (root / s["bam"]).stat().st_size / 1e9
    resources = resource_plan(root, bundle, gb, max(1, len(samples)))
    manifest = root / "prepared" / "manifest.json"
    ref = read_json(manifest) if manifest.is_file() else {}
    return {
        "schema": "pfvarx.plan/1",
        "input_type": cfg.input_type.value, "design": cfg.assay_design.value,
        "panel": validation["panel"],
        "final_output": cfg.final_output.value, "layers": resources["layers"],
        "samples": [s for k, s in sorted(samples.items()) if k not in failed_inputs],
        "samples_failed_validation": failed_inputs,
        "validation_status": validation["status"],
        "reference": {"combined_fasta": ref.get("combined_fasta"), "parasite_fasta": ref.get("parasite_fasta"),
                      "reference_sha256": ref.get("reference_sha256"), "host": ref.get("host")},
        "preparation_current": prep["status"] == "PASS",
        "prepared_verification": prep,
        "analysis_selection": cfg.analyses.selection,
        "analyses_not_applicable": {
            da: "DESIGN_NOT_SUPPORTED" for da, contract in registry["analyses"].items()
            if cfg.analyses.selection == "auto"
            and cfg.assay_design.value not in [Design.parse(d).value for d in contract["designs"]]
        },
        "analyses_requested": requested, "analyses_internal": internal,
        "resources": resources, "input_gb": round(gb, 4),
    }


def snakemake_cmd(root: Path, run: Path, eff: dict[str, Any], dry: bool) -> list[str]:
    comp = eff["compute"]
    cmd = ["snakemake", "-s", str(root / "workflow" / "Snakefile"), "--directory", str(run),
           "--config", f"run_dir={run}", f"root={root}", "--cores", str(comp["max_cores"]),
           "--jobs", str(comp["max_jobs"]), "--resources", f"mem_mb={comp['max_mem_mb']}",
           "--default-resources", f"mem_mb={min(1000, comp['per_job']['max_mem_mb'])}",
           f"runtime={comp['per_job']['max_runtime_min']}",
           *([f"tmpdir={comp['tmp_dir']}"] if comp.get("tmp_dir") else []),
           # Greedy scheduling: PfVarX jobs are few per sample with fixed limits, so Snakemake's ILP
           # optimiser adds nothing and can stall 10 s per decision before falling back to greedy.
           "--scheduler", "greedy",
           "--rerun-incomplete", "--printshellcmds", "--show-failed-logs"]
    if eff["config"].get("retention", {}).get("intermediates") == "keep":
        cmd += ["--notemp"]
    if dry:
        cmd += ["-n", "--quiet", "rules"]
    else:
        cmd += ["--keep-going"]
    return cmd


def dry_run(root: Path) -> dict[str, Any]:
    """Pipeline preflight: plan + Snakemake -n in a temporary directory outside the workspace."""
    plan = build_plan(root)
    eff = effective_config(root)
    if plan["validation_status"] == "FAIL" or not plan["samples"] or not plan["resources"]["feasible"]:
        return {"status": "NOT_READY", "reason": "Fix validation/resources and require admitted samples", "plan": plan}
    if not plan["preparation_current"]:
        return {"status": "NOT_READY", "reason": "run 'pixi run prepare' first", "plan": plan}
    with tempfile.TemporaryDirectory(prefix="pfvarx-dryrun-") as tmp:
        run = Path(tmp)
        write_json(run / "frozen_config.json", eff)
        write_json(run / "frozen_plan.json", plan)
        r = subprocess.run(snakemake_cmd(root, run, eff, dry=True), capture_output=True, text=True)
    return {"status": "OK" if r.returncode == 0 else "FAILED", "snakemake_exit": r.returncode,
            "jobs_before_checkpoints": r.stdout[-4000:], "stderr": r.stderr[-2000:],
            "note": ("snakemake -n checks the DAG up to the per-sample QC checkpoints and the downstream "
                     "eligibility checkpoint; jobs after them are decided at run time (deferred gates)."),
            "resources": plan["resources"], "outputs_created": []}


def failed_samples(run: Path, plan: dict[str, Any]) -> dict[str, str]:
    res = run / "results"
    out = {}
    for s in plan["samples"]:
        sid = s["sample_id"]
        gate = res / "qc" / "samples" / sid / "qc_decision.json"
        if not gate.is_file():
            out[sid] = "L1"
            continue
        g = read_json(gate)
        if 2 in plan["layers"] and g["permitted"].get("L2") and not (res / "variants" / "samples" / f"{sid}.vcf.gz").is_file():
            out[sid] = "L2"
    return out


ExecutionPass = TypedDict("ExecutionPass", {"pass": int, "exit_code": int,
                                            "failed_samples": dict[str, str], "log": str})

def execute(root: Path, mode: str = "NEW", previous: str | None = None) -> dict[str, Any]:
    from pfvarx.execution.retention import workspace_lock
    with workspace_lock(root):
        return _execute(root, mode, previous)


def _execute(root: Path, mode: str = "NEW", previous: str | None = None) -> dict[str, Any]:
    eff = effective_config(root)
    plan = build_plan(root)
    if plan["validation_status"] == "FAIL":
        raise ValueError("validation FAILED for global inputs; run 'pixi run validate' and fix the reported files")
    if not plan["preparation_current"]:
        raise ValueError("prepared reference is missing or stale; run 'pixi run prepare'")
    if not plan["resources"]["feasible"]:
        raise ValueError(f"resource plan infeasible: {plan['resources']['blocked']}")
    if not plan["samples"]:
        raise ValueError("No admitted samples; fix the registered inputs before execution")
    if not (root / "pixi.lock").is_file():
        raise ValueError("pixi.lock missing: a reproducible environment is required for 'pipeline'")
    run = runs.create_run(root, eff, plan, mode=mode, previous=previous)
    events = run / "logs" / "events.jsonl"
    for sid in plan["samples_failed_validation"]:
        append(events, event(run_id=run.name, layer="L0", stage="input_validation", status="BLOCKED",
                             severity=Severity.ERROR, code="INPUT_VALIDATION_FAILED",
                             failure_class=FailureClass.INPUT_FAILURE, sample_id=sid,
                             cause_confirmed="input validation failed (see metadata/validation_report.json)",
                             project_id=eff["config"]["project_id"]))
    if mode == "RESUME":
        for manifest in (run / "results/downstream").glob("DA*/status_manifest.json"):
            if read_json(manifest).get("terminal_state") == "FAILED_EXECUTION":
                history = run / "logs" / "downstream" / "failed_manifests"
                history.mkdir(parents=True, exist_ok=True)
                # Retain the old decision and force Snakemake to retry the failed rule.
                with tempfile.NamedTemporaryFile(dir=history, prefix=manifest.parent.name + "-", suffix=".json",
                                                 delete=False) as archived:
                    archived.write(manifest.read_bytes())
                manifest.unlink()
        exclusions = run / "runtime_exclusions.json"
        if exclusions.is_file():
            # A resumed run retries samples that failed at run time (for example after a transient
            # tool or resource failure). The earlier exclusion record is kept with the logs.
            history = run / "logs" / "runtime_exclusions"
            history.mkdir(parents=True, exist_ok=True)
            with tempfile.NamedTemporaryFile(dir=history, prefix="exclusions-", suffix=".json",
                                             delete=False) as archived:
                archived.write(exclusions.read_bytes())
            exclusions.unlink()
        # Without this the finished target rule is treated as up to date and nothing is retried.
        (run / "workflow.complete").unlink(missing_ok=True)
    runs.set_status(run, "RUNNING")
    env = dict(os.environ)
    env.pop("PFVARX_JSON_OUT", None)  # the CLI result file belongs to this command, not to workflow jobs
    passes: list[ExecutionPass] = []
    log_dir = run / "logs"
    prior_passes = [
        int(path.stem.removeprefix("snakemake.pass"))
        for path in log_dir.glob("snakemake.pass*.log")
        if path.stem.removeprefix("snakemake.pass").isdigit()
    ]
    offset = max(prior_passes, default=0)
    for attempt in (1, 2):
        log = log_dir / f"snakemake.pass{offset + attempt}.log"
        with log.open("x") as h:
            r = subprocess.run(snakemake_cmd(root, run, eff, dry=False), stdout=h, stderr=subprocess.STDOUT, env=env)
        failed = failed_samples(run, plan)
        passes.append({"pass": offset + attempt, "exit_code": r.returncode, "failed_samples": failed, "log": str(log)})
        if r.returncode != 0 and len(failed) == len(plan["samples"]):
            # Nothing completed for any sample: a workflow-level error (e.g. DAG construction),
            # not a per-sample failure. Do not re-run; report it as an execution failure.
            append(events, event(run_id=run.name, layer="L0", stage="workflow execution", status="FAILED",
                                 severity=Severity.CRITICAL, code="TOOL_NONZERO_EXIT",
                                 failure_class=FailureClass.SOFTWARE_FAILURE, rule="snakemake",
                                 log_paths=[str(log.relative_to(run))], exit_code=r.returncode,
                                 interpretation_impact="no sample completed; see the Snakemake log",
                                 project_id=eff["config"]["project_id"]))
            break
        if r.returncode == 0 or not failed or attempt == 2:
            break
        prior = read_json(run / "runtime_exclusions.json")["samples"] if (run / "runtime_exclusions.json").is_file() else []
        write_json(run / "runtime_exclusions.json", {"samples": sorted(set(prior) | set(failed)), "stages": failed})
        for sid, stage in failed.items():
            append(events, event(run_id=run.name, layer=stage, stage=f"{stage} processing", status="FAILED",
                                 severity=Severity.ERROR, code="TOOL_NONZERO_EXIT",
                                 failure_class=FailureClass.SOFTWARE_FAILURE, sample_id=sid,
                                 log_paths=[str(p.relative_to(run)) for p in (run / "logs").rglob(f"*{sid}*")],
                                 affected_downstream=["L3", "L4", "L5"],
                                 interpretation_impact="sample removed from cohort; denominators re-evaluated",
                                 project_id=eff["config"]["project_id"]))
    from pfvarx.reporting.build import build_report
    downstream_failures = sorted(
        manifest.parent.name for manifest in (run / "results/downstream").glob("DA*/status_manifest.json")
        if read_json(manifest).get("terminal_state") == "FAILED_EXECUTION"
    )
    any_failed = any(p["failed_samples"] for p in passes) or bool(downstream_failures)
    if passes[-1]["exit_code"] != 0:
        final = "FAILED_EXECUTION"
    elif passes[-1]["exit_code"] == 0 and not any_failed:
        final = "COMPLETED"
    else:
        final = "COMPLETED_WITH_FAILURES"
    runs.set_status(run, final, passes=passes, downstream_failures=downstream_failures)
    report = build_report(root, run)
    link = root / "results"
    if link.is_symlink() or not link.exists():
        if link.is_symlink():
            link.unlink()
        link.symlink_to(run.relative_to(root) / "results")
    cleanup = None
    if final == "COMPLETED" and eff["config"].get("retention", {}).get("intermediates", "successful") == "successful":
        from pfvarx.execution.retention import cleanup_run
        cleanup = cleanup_run(root, run, apply=True)
    return {"run_id": run.name, "status": final, "passes": passes, "downstream_failures": downstream_failures, "report": report, "cleanup": cleanup}

