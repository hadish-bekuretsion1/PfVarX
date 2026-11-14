"""PfVarX command-line interface (Layer 0). Every command prints JSON.

Read-only commands: check-layout, validate (writes only metadata/validation_report.json),
ready, dry-run, workflow-check, run-verify, status. Mutating commands: config-init, init,
add-input, fetch-input, prepare, pipeline, workflow-run, report, fixture, benchmark-*.
"""

from __future__ import annotations

import json
import os
import subprocess
from pathlib import Path
from typing import Any

import typer

from pfvarx import __version__
from pfvarx.configuration.loader import ConfigError, effective_config, load_bundle

app = typer.Typer(help="PfVarX: assay-aware SNP/small-indel analysis of P. falciparum", no_args_is_help=True,
                  add_completion=False)


def _root() -> Path:
    root = Path.cwd().resolve()
    if not (root / "pixi.toml").is_file() or not (root / "workflow" / "Snakefile").is_file():
        raise typer.BadParameter("run PfVarX from the repository/workspace root (pixi.toml and workflow/ present)")
    return root


def emit(data: Any, code: int = 0) -> None:
    """Print the command result as JSON on stdout (pixi's own banner goes to stderr).

    If PFVARX_JSON_OUT is set, the same JSON is also written to that file, so scripts can
    read the result without parsing a terminal stream (e.g. one captured with 2>&1).
    """
    text = json.dumps(data, indent=2, default=str)
    typer.echo(text)
    target = os.environ.get("PFVARX_JSON_OUT")
    if target:
        path = Path(target)
        path.parent.mkdir(parents=True, exist_ok=True)
        tmp = path.with_name(path.name + ".tmp")
        tmp.write_text(text + "\n")
        tmp.replace(path)  # atomic: a reader never sees a half-written file
    if code:
        raise typer.Exit(code)


def guard(fn: Any) -> Any:
    try:
        return fn()
    except ConfigError as exc:
        emit({"status": "CONFIG_ERROR", "problems": exc.problems}, 2)
    except subprocess.CalledProcessError as exc:
        emit({"status": "ERROR", "error": f"Required tool exited {exc.returncode}; inspect stage/tool logs and verify installed environment"}, 2)
    except (ValueError, OSError) as exc:
        emit({"status": "ERROR", "error": str(exc)}, 2)


@app.command()
def version() -> None:
    """Print the package version."""
    typer.echo(__version__)


@app.command()
def configure(gui: bool = typer.Option(False, "--gui"), port: int = typer.Option(0, min=0, max=65535),
              no_browser: bool = typer.Option(False, "--no-browser")) -> None:
    """Open the local configuration form; YAML remains available for headless use."""
    from pfvarx.configuration.web import serve

    if not gui:
        raise typer.BadParameter("use --gui, or edit the YAML created by config-init")
    guard(lambda: serve(_root(), port, not no_browser))


@app.command("config-init")
def config_init(purpose: str = typer.Option(..., "--purpose", help="analysis | benchmarking")) -> None:
    """Copy the purpose-specific YAML templates into config/ without overwriting."""
    from pfvarx.initialization.layout import config_init as ci

    root = _root()
    emit({"created": [str(p.relative_to(root)) for p in guard(lambda: ci(root, purpose))]})


@app.command()
def init() -> None:
    """Create purpose- and assay-specific input folders only."""
    from pfvarx.initialization.layout import check_release_label, init_inputs

    root = _root()

    def create() -> Any:
        cfg = load_bundle(root).config
        check_release_label(cfg)  # the GUI applies the same gate before it creates folders
        return init_inputs(root, cfg)
    emit(guard(create))


@app.command("check-layout")
def check_layout() -> None:
    """Read-only check of required folders and files."""
    from pfvarx.initialization.layout import check_layout as cl

    root = _root()
    res = guard(lambda: cl(root, load_bundle(root).config))
    emit(res, 0 if res["status"] == "PASS" else 1)


@app.command()
def validate() -> None:
    """Validate every input read-only; write metadata/validation_report.json."""
    from pfvarx.execution.preflight import validate_inputs

    root = _root()
    rep = guard(lambda: validate_inputs(root))
    brief = {"status": rep["status"], "samples": rep["samples"],
             "failures": [f for r in rep["reports"] for f in r["findings"] if f["level"] == "FAIL"],
             "warnings": [f for r in rep["reports"] for f in r["findings"] if f["level"] == "WARNING"],
             "report": "metadata/validation_report.json"}
    emit(brief, 0 if rep["status"] in ("PASS", "PARTIAL") else 1)


@app.command()
def prepare(dry_run: bool = typer.Option(False, "--dry-run")) -> None:
    """Build only the missing or stale reference derivatives under prepared/."""
    from pfvarx.execution.preflight import prepare_inputs

    root = _root()
    if dry_run:
        emit(guard(lambda: prepare_inputs(root, True)))
        return
    from pfvarx.execution.retention import workspace_lock
    def prepare_locked():
        with workspace_lock(root):
            return prepare_inputs(root, dry_run)
    emit(guard(prepare_locked))


@app.command()
def ready() -> None:
    """Read-only readiness: validation, preparation, eligibility pre-plan and resource feasibility."""
    from pfvarx.execution.runner import build_plan

    root = _root()
    plan = guard(lambda: build_plan(root))
    ok = (plan["validation_status"] in ("PASS", "PARTIAL") and plan["preparation_current"]
          and plan["resources"]["feasible"] and bool(plan["samples"]))
    emit({"ready": ok, "validation": plan["validation_status"], "preparation_current": plan["preparation_current"],
          "samples": [s["sample_id"] for s in plan["samples"]], "samples_failed_validation": plan["samples_failed_validation"],
          "layers": plan["layers"], "analyses_requested": plan["analyses_requested"],
          "analyses_internal_prerequisites": plan["analyses_internal"], "resources": plan["resources"],
          "deferred_gates": ["per-sample QC (Layer 1 checkpoint)", "clonal status (DA001) where required",
                             "downstream eligibility checkpoint", "DA006 alignment informativeness"]},
         0 if ok else 1)


@app.command("dry-run")
def dry_run() -> None:
    """Scientific preflight plus 'snakemake -n' in a temporary directory; creates no outputs."""
    from pfvarx.execution.runner import dry_run as dr

    res = guard(lambda: dr(_root()))
    emit(res, 0 if res["status"] == "OK" else 1)


@app.command()
def pipeline(mode: str = typer.Option("NEW", "--mode", help="NEW | RESUME | EXTEND"),
             run_id: str | None = typer.Option(None, "--run-id")) -> None:
    """Execute the single configured final-output mode."""
    from pfvarx.execution.runner import execute

    root = _root()
    res = guard(lambda: execute(root, mode.upper(), run_id))
    emit(res, 0 if res["status"] == "COMPLETED" else 3)


@app.command()
def report(run_id: str | None = typer.Option(None, "--run-id")) -> None:
    """Rebuild the offline dashboard from stored results (no analysis is re-run)."""
    from pfvarx.reporting.build import build_report
    from pfvarx.runs.registry import latest

    root = _root()
    run = root / "runs" / run_id if run_id else latest(root)
    if run is None:
        emit({"status": "ERROR", "error": "no runs"}, 2)
    emit(guard(lambda: build_report(root, run)))  # type: ignore[arg-type]


@app.command("run-verify")
def run_verify(run_id: str = typer.Option(..., "--run-id")) -> None:
    """Verify frozen-file integrity and report workspace drift (read-only)."""
    from pfvarx.runs.registry import verify

    res = guard(lambda: verify(_root(), run_id))
    emit(res, 0 if res["frozen_integrity"] == "PASS" else 1)


@app.command()
def status() -> None:
    """Show the effective configuration hash and the latest run status."""
    from pfvarx.io import read_json
    from pfvarx.runs.registry import latest

    root = _root()
    eff = guard(lambda: effective_config(root))
    run = latest(root)
    emit({"effective_config_hash": eff["hash"], "latest_run": run.name if run else None,
          "latest_status": read_json(run / "status.json") if run else None})


@app.command()
def fixture(design: str = typer.Option("WGS", "--design"), purpose: str = typer.Option("analysis", "--purpose"),
            depth: int = typer.Option(40, "--depth"), seed: int = typer.Option(7, "--seed"),
            n_mono: int = typer.Option(6, "--n-mono"), n_poly: int = typer.Option(3, "--n-poly"),
            low_depth_sample: bool = typer.Option(False, "--low-depth-sample"),
            variant_scale: int = typer.Option(1, "--variant-scale", min=1,
                                              help="multiply simulated variant sites per contig")) -> None:
    """Write a synthetic SOFTWARE-TEST dataset into inputs/ (not biological data)."""
    from pfvarx.fixtures.synthetic import build

    root = _root()
    if any((root / "inputs" / "fastq").glob("*.fastq.gz")):
        emit({"status": "ERROR", "error": "inputs/fastq already contains FASTQ; refusing to overwrite"}, 2)
    emit(build(root, design=design, purpose=purpose, depth=depth, seed=seed, n_mono=n_mono, n_poly=n_poly,
               low_depth_sample=low_depth_sample, variant_scale=variant_scale))


@app.command("benchmark-partition")
def benchmark_partition(holdout_fraction: float = typer.Option(0.25, "--holdout-fraction"),
                        seed: int = typer.Option(20261106, "--seed")) -> None:
    """Lock the benchmarking / held-out partition of truth samples (once, before tuning)."""
    from pfvarx.benchmarking.partitions import lock

    root = _root()
    bundle = guard(lambda: load_bundle(root))
    if bundle.config.purpose.value != "benchmarking":
        emit({"status": "ERROR", "error": "benchmark-partition is only available in a benchmarking workspace"}, 2)
    emit(guard(lambda: lock(root, holdout_fraction, seed)))


@app.command("benchmark-evaluate")
def benchmark_evaluate(run_id: str | None = typer.Option(None, "--run-id"),
                       out: str = typer.Option("benchmarks/evaluation", "--out"),
                       holdout: bool = typer.Option(False, "--holdout", help="score the held-out partition (once)")) -> None:
    """Score a completed benchmarking run against truth (benchmarking workspace only)."""
    from pfvarx.benchmarking.evaluate import evaluate
    from pfvarx.runs.registry import latest

    root = _root()
    bundle = guard(lambda: load_bundle(root))
    if bundle.config.purpose.value != "benchmarking":
        emit({"status": "ERROR", "error": "benchmark-evaluate is only available in a benchmarking workspace"}, 2)
    run = root / "runs" / run_id if run_id else latest(root)
    emit(guard(lambda: evaluate(root, run, root / out, holdout)))  # type: ignore[arg-type]


@app.command("benchmark-resources")
def benchmark_resources(run_id: str | None = typer.Option(None, "--run-id")) -> None:
    """Summarise measured per-job resources of a run and write config/profiles/resources.yaml."""
    from pfvarx.benchmarking.resources import summarise
    from pfvarx.runs.registry import latest

    root = _root()
    run = root / "runs" / run_id if run_id else latest(root)
    emit(guard(lambda: summarise(root, run)))  # type: ignore[arg-type]


@app.command("workspace-init")
def workspace_init(config_only: bool = typer.Option(False, "--config-only")) -> None:
    """Initialize only profile/config-selected input and metadata folders; never overwrite configuration."""
    from pfvarx.initialization.workspace import initialize
    emit(guard(lambda: initialize(_root(), config_only=config_only)))


@app.command("discover-inputs")
def discover_inputs(directory: Path = typer.Option(..., "--directory")) -> None:
    """Generate a review-required candidate sample sheet from unambiguous R1/R2 filenames."""
    from pfvarx.inputs.organize import discover
    emit(guard(lambda: discover(_root(), directory)))


@app.command("register-inputs")
def register_inputs(manifest: Path = typer.Option(..., "--manifest"),
                    mode: str = typer.Option("reference", "--mode", help="reference | link; never copy reads")) -> None:
    """Register a reviewed TSV, preserving originals and refusing nonempty sheet overwrite."""
    from pfvarx.inputs.organize import register
    emit(guard(lambda: register(_root(), manifest, mode)))


@app.command("fetch-input")
def fetch_input(url: str = typer.Option(..., "--url"), sha256: str | None = typer.Option(None, "--sha256"),
                role: str = typer.Option("raw", "--role", help="raw (reads or alignments) | reference | targets | metadata | truth"),
                name: str = typer.Option(..., "--name")) -> None:
    """Download an HTTPS resource into its input folder (same checks as the GUI download)."""
    from pfvarx.execution.retention import workspace_lock
    from pfvarx.inputs.transfer import download

    root = _root()

    def run() -> Any:
        with workspace_lock(root):
            return download(root, url, sha256, role, name)
    emit(guard(run))


@app.command("add-input")
def add_input(files: list[Path] = typer.Option(..., "--file", help="file on this computer; repeat for several"),
              role: str = typer.Option("raw", "--role", help="raw (reads or alignments) | reference | targets | metadata | truth"),
              name: str | None = typer.Option(None, "--name", help="destination filename (single --file only)")) -> None:
    """Copy files from disk into their input folder (same checks as the GUI upload)."""
    from pfvarx.execution.retention import workspace_lock
    from pfvarx.inputs.transfer import import_file

    root = _root()
    if name is not None and len(files) != 1:
        raise typer.BadParameter("--name needs exactly one --file")

    def run() -> Any:
        saved = []
        with workspace_lock(root):
            for source in files:
                saved.append(import_file(root, source, role, name))
        return {"saved": saved}
    emit(guard(run))


@app.command("verify-prepared")
def prepared_verify() -> None:
    """Verify registered inputs, source/derived hashes, indexes and sequence dictionary before execution."""
    from pfvarx.execution.preflight import verify_inputs
    result = guard(lambda: verify_inputs(_root()))
    emit(result, 0 if result["status"] == "PASS" else 1)


@app.command("orchestrate")
def orchestrate_command(through: str = typer.Option("execute", "--through")) -> None:
    """Initialize -> validate -> prepare -> verify -> preflight -> execute/report/retain, stopping on failure."""
    from pfvarx.execution.preflight import orchestrate
    result = guard(lambda: orchestrate(_root(), through))
    emit(result, 0 if result["status"] in ("READY_FOR_CONFIGURATION", "VALIDATED", "PREPARED", "VERIFIED", "COMPLETED") else 3)


@app.command("cleanup")
def cleanup_command(run_id: str = typer.Option(..., "--run-id"), apply: bool = typer.Option(False, "--apply")) -> None:
    """Preview completed-run intermediate cleanup; --apply unlinks only safe results/work files."""
    from pfvarx.execution.retention import cleanup_run, workspace_lock
    from pfvarx.runs.registry import RUN_ID
    root = _root()
    if not RUN_ID.fullmatch(run_id):
        emit({"status": "ERROR", "error": "Invalid run ID"}, 2)
    def clean():
        with workspace_lock(root):
            return cleanup_run(root, root / "runs" / run_id, apply)
    emit(guard(clean))


@app.command("input-directories")
def input_directories() -> None:
    """Create input directories after editing the configuration YAML."""
    init()


@app.command("input-files")
def input_files() -> None:
    """List configured input folders for copying or registering files."""
    from pfvarx.initialization.layout import required_dirs

    root = _root()
    bundle = guard(lambda: load_bundle(root))
    emit({"directories": [{"path": path + "/", "exists": (root / path).is_dir()}
                          for path in required_dirs(bundle.config)]})


@app.command("validate-input-files")
def validate_input_files() -> None:
    """Validate files registered in the saved study configuration."""
    validate()


def _run_checks(stages: tuple[str, ...]) -> None:
    from pfvarx.configuration.web import run_preflight

    result = guard(lambda: run_preflight(_root(), stages))
    emit(result, 0 if result["passed"] else 1)


@app.command("prepare-reference")
def prepare_reference() -> None:
    """Prepare references, verify artifacts and check readiness."""
    _run_checks(("prepare", "verify-prepared", "ready"))


@app.command("workflow-readiness-checks")
def workflow_readiness_checks() -> None:
    """Validate, prepare, verify and check readiness in order."""
    from pfvarx.configuration.web import PREFLIGHT_STAGES

    _run_checks(PREFLIGHT_STAGES)


@app.command("workflow-check")
def workflow_check() -> None:
    """Check the workflow plan without creating outputs (the GUI 'Check workflow' button)."""
    dry_run()


@app.command("workflow-run")
def workflow_run() -> None:
    """Check the workflow, then run it as a new analysis (the GUI 'Run workflow' button)."""
    from pfvarx.execution.runner import dry_run as dr
    from pfvarx.execution.runner import execute

    root = _root()
    check = guard(lambda: dr(root))
    if check["status"] != "OK":
        emit({"status": "NOT_STARTED", "reason": "workflow check failed; nothing was run", "workflow_check": check}, 1)
    result = guard(lambda: execute(root, "NEW", None))
    emit({**result, "workflow_check": {"status": "OK"}}, 0 if result["status"] == "COMPLETED" else 3)


if __name__ == "__main__":
    app()
