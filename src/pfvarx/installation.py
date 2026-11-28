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

import json
import os
import runpy
import shutil
import subprocess
from pathlib import Path
from typing import Any

from pfvarx.io import read_json, write_json

# The tools and packages a configuration can need: (conda package, name shown, what PfVarX uses
# it for). Listed in pipeline order; only those in the chosen environment are shown.
TOOLS = [
    ("python", "Python", "runs PfVarX itself (configuration, validation, QC gates, reports)"),
    ("pysam", "pysam", "reads FASTA, BAM and VCF files in validation and the caller"),
    ("snakemake-minimal", "Snakemake", "runs the workflow steps in order"),
    ("fastp", "fastp", "read trimming and read QC (FASTQ input)"),
    ("bwa-mem2", "BWA-MEM2", "aligns reads to the reference (and to the human genome with Human DNA removal)"),
    ("samtools", "SAMtools", "sorts, merges, indexes and counts alignments; FASTA index"),
    ("picard-slim", "Picard", "marks duplicates (WGS) and writes the sequence dictionary"),
    ("mosdepth", "mosdepth", "depth and callable coverage for the QC gates"),
    ("fgbio-minimal", "fgbio", "targeted designs: primer trimming, UMI grouping and consensus reads"),
    ("numpy", "NumPy", "variant caller arithmetic (FO002 and later)"),
    ("scipy", "SciPy", "variant caller statistics and mixture fitting (FO002 and later)"),
    ("numba", "Numba", "compiles the caller's pair-HMM indel rescoring for speed (FO002 and later)"),
    ("bcftools", "BCFtools", "normalises per-sample VCFs and merges the cohort (FO002 and later)"),
    ("snpeff", "SnpEff", "annotates the cohort VCF with gene effects (FO004)"),
    ("r-base", "R", "runs the R analysis scripts (DA001, DA004, DA007)"),
    ("r-moire", "MOIRE", "complexity of infection and clonality (DA001)"),
    ("r-dcifer", "Dcifer", "relatedness between infections (DA004)"),
    ("hmmibd", "hmmIBD", "identity by descent between infections (DA004)"),
    ("scikit-allel", "scikit-allel", "population structure (DA005)"),
    ("iqtree", "IQ-TREE", "phylogeny (DA006)"),
    ("r-mgcv", "mgcv", "spatial and temporal models (DA007)"),
    ("rtg-tools", "RTG Tools", "benchmark comparison against truth VCFs (benchmarking)"),
]
_LOCK_CACHE: dict[tuple[str, float, str], dict[str, dict[str, str]]] = {}

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


def required(root: Path, cfg: Any = None) -> dict[str, Any]:
    """What a configuration needs (the saved one by default): the environment and the settings that decide it."""
    from pfvarx.configuration.loader import load_bundle
    from pfvarx.eligibility.engine import internal_prerequisites, load_registry

    cfg = cfg if cfg is not None else load_bundle(root).config
    bootstrap = runpy.run_path(str(root / "bootstrap.py"))
    if cfg.final_output.value == "FO004" and cfg.analyses.selection == "auto":
        analyses = sorted(bootstrap["METHODS"])  # eligibility decides per run; install every method's tools
    elif cfg.final_output.value == "FO004":
        requested = list(cfg.analyses.requested)
        analyses = sorted(set(requested) | set(internal_prerequisites(requested, load_registry(root),
                                                                      cfg.analyses.spatial_outcome_source)))
    else:
        analyses = []
    rtg = cfg.purpose.value == "benchmarking" and bool(metadata(root).get("rtg"))
    environment = bootstrap["environment_for"](cfg.assay_design.value, analyses, rtg, cfg.input_type.value,
                                               cfg.final_output.value)
    return {"environment": environment, "final_output": cfg.final_output.value, "design": cfg.assay_design.value,
            "input_type": cfg.input_type.value, "analyses": analyses}


def host_platform(root: Path) -> str:
    """Use the same stdlib-only platform detection as bootstrap."""
    return runpy.run_path(str(root / "bootstrap.py"))["host_platform"]()


def locked_versions(root: Path) -> dict[str, dict[str, str]]:
    """Locked versions for this host only; never combine builds for different platforms."""
    import yaml

    lock = root / "pixi.lock"
    target = host_platform(root)
    key = (str(lock), lock.stat().st_mtime, target)
    if key not in _LOCK_CACHE:
        loader = getattr(yaml, "CSafeLoader", yaml.SafeLoader)
        data = yaml.load(lock.read_text(), Loader=loader)  # noqa: S506 - a safe loader
        versions: dict[str, str] = {}
        for package in data.get("packages", []):
            url = package.get("conda") or ""
            if url:
                name, version, _build = url.rsplit("/", 1)[-1].rsplit("-", 2)
                versions[url] = version
                versions.setdefault(name, version)
        environments: dict[str, dict[str, str]] = {}
        for environment, record in data.get("environments", {}).items():
            if target not in record.get("packages", {}):
                raise ValueError(f"Environment {environment} has no lock for {target}; use a source release "
                                 "with this platform locked. Do not copy an installed .pixi directory between hosts.")
            names: dict[str, str] = {}
            for item in record.get("packages", {}).get(target, []):
                url = item.get("conda") or ""
                if url:
                    names[url.rsplit("/", 1)[-1].rsplit("-", 2)[0]] = versions.get(url, "")
            environments[environment] = names
        _LOCK_CACHE.clear()
        _LOCK_CACHE[key] = environments
    return _LOCK_CACHE[key]


def tools_for(root: Path, environment: str, need: dict[str, Any] | None = None) -> list[dict[str, Any]]:
    """The tools and packages *environment* installs, with locked version and what each is for."""
    prefix = root / ".pixi" / "envs" / environment
    packages = locked_versions(root).get(environment, {})
    versions = {}
    selected = locked_versions(root)[need["environment"]] if need is not None else packages
    tracked = {name for name, _, _ in TOOLS if name in selected}
    if need is not None:
        if need["input_type"] != "FASTQ":
            tracked.discard("fastp")
            if need["design"] != "targeted_umi":
                tracked.discard("bwa-mem2")
    for record in (prefix / "conda-meta").glob("*.json"):
        # Conda record names are <package>-<version>-<build>.json.
        package = record.stem.rsplit("-", 2)[0]
        if package not in packages or package not in tracked:
            continue
        try:
            info = json.loads(record.read_bytes())
        except (OSError, ValueError):
            continue
        if isinstance(info, dict) and info.get("name") == package and isinstance(info.get("version"), str):
            versions[package] = info["version"]
    return [{"package": package, "name": name, "version": packages.get(package, selected[package]), "used_for": purpose,
             "installed_version": versions.get(package), "installed": package in packages and versions.get(package) == packages[package]}
            for package, name, purpose in TOOLS if package in tracked]


def _complete(root: Path, environment: str, tools: list[dict[str, Any]]) -> bool:
    return bool(tools) and environment_python(root, environment).is_file() and all(t["installed"] for t in tools)



def plan(root: Path, values: dict[str, Any] | None = None) -> dict[str, Any]:
    """The tools a configuration needs: the saved one, or the GUI form's current values (not saved)."""
    cfg = None
    if values is not None:
        from pfvarx.configuration.loader import validate_bundle

        cfg = validate_bundle(root, values).config
    meta = metadata(root)
    # A freshly exported up-front workspace may still lack its editable configuration.
    # In that case its explicitly installed profile is the only available plan.
    upfront = meta.get("profile") != ON_DEMAND
    need = None if upfront and cfg is None and not (root / "config/config.yaml").exists() else required(root, cfg)
    if upfront:
        environment = meta.get("environment", "default")
    else:
        assert need is not None
        environment = need["environment"]
    tools = tools_for(root, environment, need)
    return {"mode": "up-front" if upfront else ON_DEMAND, "environment": environment, "required_for": need,
            "tools": tools, "installed": _complete(root, environment, tools)}


def status(root: Path) -> dict[str, Any]:
    """Verify packages required by the saved configuration in its execution environment."""
    selection = plan(root)
    environment = selection["environment"]
    installed = metadata(root).get("environment", "default")
    current = installed == environment and selection["installed"]
    if current:
        detail = f"Required packages installed for this configuration in {environment}."
    elif selection["mode"] == ON_DEMAND:
        detail = (f"Required packages for {environment} are not ready. Click Install required tools "
                  "on the Input directories tab, or run ./pfvarx install-tools.")
    else:
        detail = f"Package verification incomplete. Run pixi install --locked -e {environment} in this workspace."
    return {"mode": selection["mode"], "current": current, "installed_environment": installed,
            "required_environment": environment, "required_for": selection["required_for"],
            "tools": selection["tools"], "detail": detail}


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
    if meta.get("environment") == environment and _complete(root, environment, tools_for(root, environment, need)):
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
    tools = tools_for(root, environment, need)
    if not _complete(root, environment, tools):
        missing = ", ".join(t["name"] for t in tools if not t["installed"])
        raise ValueError(f"Pixi finished but package verification failed for {environment}: {missing}. See {LOG}.")
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

