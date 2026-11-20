"""Install the tools the saved configuration needs (GUI: Install required tools).

A workspace bootstrapped on demand (``python3 bootstrap.py`` without ``--final-output``)
starts with the small ``control`` environment: configuration, GUI and validation, no
pipeline tools. ``./pfvarx install-tools`` then reads the saved configuration (final
output, design, input type and analyses with their prerequisites), chooses the matching
locked environment with the same rules as ``bootstrap.py --final-output``, installs it
with ``pixi install --locked`` and records it in ``config/installation.json``; from then
on ``./pfvarx`` and the GUI's steps run in it. Pressing it again after the configuration
changed installs only the newly needed environment; an unchanged configuration installs
nothing.

Steps that run tools (Inputs Preparation, readiness, dry run, pipeline) call
:func:`require_tools` and stop with a plain instruction when the tools for the saved
configuration are not installed yet. Workspaces installed up front
(``--final-output``/``--profile``) are unaffected.
"""

from __future__ import annotations

import os
import runpy
import shutil
import subprocess
from pathlib import Path
from typing import Any

from pfvarx.io import read_json, write_json

ON_DEMAND = "on-demand"
CONTROL = "control"
INSTALLATION = "config/installation.json"
LOG = "metadata/install_tools.log"


def metadata(root: Path) -> dict[str, Any]:
    path = root / INSTALLATION
    return read_json(path) if path.is_file() else {}


def on_demand(root: Path) -> bool:
    return metadata(root).get("profile") == ON_DEMAND


def environment_python(root: Path, environment: str) -> Path:
    return root / ".pixi" / "envs" / environment / "bin" / "python"


def required(root: Path) -> dict[str, Any]:
    """What the saved configuration needs: the environment and the settings that decide it."""
    from pfvarx.configuration.loader import load_bundle
    from pfvarx.eligibility.engine import internal_prerequisites, load_registry

    cfg = load_bundle(root).config
    bootstrap = runpy.run_path(str(root / "bootstrap.py"))
    if cfg.final_output.value == "FO004" and cfg.analyses.selection == "auto":
        analyses = sorted(bootstrap["METHODS"])  # eligibility decides per run; install every method's tools
    else:
        requested = list(cfg.analyses.requested)
        analyses = sorted(set(requested) | set(internal_prerequisites(requested, load_registry(root),
                                                                      cfg.analyses.spatial_outcome_source)))
    environment = bootstrap["environment_for"](cfg.assay_design.value, analyses, False, cfg.input_type.value,
                                               cfg.final_output.value)
    return {"environment": environment, "final_output": cfg.final_output.value, "design": cfg.assay_design.value,
            "input_type": cfg.input_type.value, "analyses": analyses}


def status(root: Path) -> dict[str, Any]:
    """Whether the tools for the saved configuration are installed."""
    meta = metadata(root)
    if meta.get("profile") != ON_DEMAND:
        return {"mode": "up-front", "current": True, "installed_environment": meta.get("environment"),
                "detail": "Tools were installed with the workspace (bootstrap --final-output or --profile)."}
    need = required(root)
    installed = meta.get("environment", CONTROL)
    current = installed == need["environment"] and environment_python(root, installed).is_file()
    return {"mode": ON_DEMAND, "current": current, "installed_environment": installed,
            "required_environment": need["environment"], "required_for": need,
            "detail": (f"Installed for this configuration (environment {installed})." if current else
                       f"Not installed yet for this configuration: {need['final_output']}, {need['design']}, "
                       f"{need['input_type']} needs environment {need['environment']}. Click Install required "
                       "tools on the Input directories tab, or run ./pfvarx install-tools.")}


def require_tools(root: Path) -> None:
    """Stop a tool-running step when the tools for the saved configuration are not installed."""
    if not on_demand(root):
        return
    state = status(root)
    if not state["current"]:
        raise ValueError("TOOLS_NOT_INSTALLED: " + state["detail"])


def install(root: Path) -> dict[str, Any]:
    """Install the environment the saved configuration needs (nothing when it is installed)."""
    meta = metadata(root)
    if meta.get("profile") != ON_DEMAND:
        return {"status": "NOT_NEEDED", **status(root)}
    need = required(root)
    environment = need["environment"]
    if meta.get("environment") == environment and environment_python(root, environment).is_file():
        return {"status": "ALREADY_INSTALLED", "environment": environment, "required_for": need}
    pixi = shutil.which("pixi")
    if pixi is None:
        raise ValueError("Pixi is not on PATH; install the reviewed Pixi version, then run ./pfvarx install-tools")
    log = root / LOG
    log.parent.mkdir(parents=True, exist_ok=True)
    env = {k: v for k, v in os.environ.items() if not k.startswith("PIXI_") and k != "CONDA_PREFIX"}
    with log.open("w") as handle:
        result = subprocess.run([pixi, "install", "--locked", "-e", environment], cwd=root, env=env,
                                stdout=handle, stderr=subprocess.STDOUT)
    if result.returncode or not environment_python(root, environment).is_file():
        raise ValueError(f"Installing environment {environment} failed (exit {result.returncode}); see {LOG}. "
                         "It needs internet access to the conda channels, or a filled Pixi cache.")
    meta.update(environment=environment, tools_installed_for=need)
    meta.setdefault("environments_installed", [])
    if environment not in meta["environments_installed"]:
        meta["environments_installed"].append(environment)
    write_json(root / INSTALLATION, meta)
    return {"status": "INSTALLED", "environment": environment, "required_for": need, "log": LOG}


def stage_environment(root: Path, base: dict[str, str]) -> dict[str, str]:
    """Process environment for a pipeline stage: the installed tools environment's ``bin`` first on PATH.

    The GUI keeps running in the environment it started in (on demand: ``control``), so the
    stages it starts after Install required tools must find samtools, bwa-mem2 and the rest there.
    """
    env = dict(base)
    environment = metadata(root).get("environment")
    if environment:
        prefix = root / ".pixi" / "envs" / environment
        if (prefix / "bin").is_dir():
            env["PATH"] = str(prefix / "bin") + os.pathsep + env.get("PATH", "")
            env["CONDA_PREFIX"] = str(prefix)
    return env


def stage_python(root: Path, fallback: str) -> str:
    """The Python that runs pipeline stages: the installed tools environment's, else *fallback*."""
    environment = metadata(root).get("environment")
    if environment:
        candidate = environment_python(root, environment)
        if candidate.is_file():
            return str(candidate)
    return fallback
