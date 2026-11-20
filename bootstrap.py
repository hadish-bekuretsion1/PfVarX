"""Export/install a minimal configured workspace; Python 3.11+, no Python packages needed.

Works from a full source extraction or a sparse Git clone. No code runs merely by cloning.
Two ways to install the tools, both following the final output:

* On demand (the default, no --final-output or --profile): only the small "control"
  environment is installed (configuration, GUI, validation). After saving the configuration,
  "Install required tools" (./pfvarx install-tools) installs exactly what it needs.
* Up front (--final-output FO00N, or --profile): the tools are installed now, for machines
  that have no internet access later (for example HPC compute nodes). FO001 installs
  alignment and QC tools only; FO002 adds the variant caller and bcftools; FO003 the cohort
  step; FO004 the annotation and only the tools of the requested downstream analyses.
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parent
METHODS = {
    "DA002": "resistance",
    "DA003": "diversity",
    "DA001": "coi",
    "DA004": "relatedness",
    "DA005": "structure",
    "DA006": "phylogeny",
    "DA007": "spatial",
}
R_SCRIPTS = {"DA001": "da001_moire.R", "DA004": "da004_dcifer.R", "DA007": "da007_gam.R"}
ROOT_FILES = {
    "README.md",
    "CHANGELOG.md",
    "LICENSE",
    ".gitignore",
    "pixi.toml",
    "pixi.lock",
    "pyproject.toml",
    "bootstrap.py",
    "pfvarx",
}


DESIGNS = ("wgs", "targeted_umi", "targeted_no_umi")
FINAL_OUTPUTS = ("FO001", "FO002", "FO003", "FO004")
# The profile each final output installs (benchmark and development are chosen by --profile).
FINAL_OUTPUT_PROFILE = {"FO001": "sample", "FO002": "sample", "FO003": "cohort", "FO004": "downstream"}


ON_DEMAND = "on-demand"


def profile_for(final_output: str | None, profile: str | None, default: str = "sample") -> str:
    """The installation profile from --final-output and/or --profile (they must agree)."""
    if final_output is None:
        return profile or default
    if final_output not in FINAL_OUTPUTS:
        raise ValueError("Unknown final output; use FO001, FO002, FO003 or FO004")
    expected = FINAL_OUTPUT_PROFILE[final_output]
    if profile not in (None, expected):
        raise ValueError(f"--final-output {final_output} installs the {expected} profile; omit --profile {profile}")
    return expected


def design_of(value: str) -> str:
    """The assay design (config ``assay_design``); the upper-case spellings of earlier releases are accepted."""
    text = value.strip().lower()
    text = {"targeted_non_umi": "targeted_no_umi"}.get(text, text)
    if text not in DESIGNS:
        raise ValueError("Unknown design; use wgs, targeted_umi or targeted_no_umi")
    return text


def environment_for(design: str, analyses: list[str], rtg: bool = False, input_type: str = "FASTQ",
                    final_output: str | None = None) -> str:
    targeted = design_of(design) != "wgs"
    if final_output == "FO001":
        # Alignment and QC only: no caller, bcftools, numpy or scipy (pixi feature "calling").
        return "alignment-qc-targeted" if targeted else "alignment-qc"
    features = set()
    if targeted:
        # fgbio trims primers and builds UMI consensus for FASTQ and for imported BAMs whose
        # documented state still needs it; consensus reads are realigned with bwa-mem2.
        features.add("targeted")
    if analyses:
        features.add("annotation")
    if any(d in analyses for d in ("DA001", "DA004", "DA005", "DA006")):
        features.add("coi")
    if "DA004" in analyses:
        features.add("related")
    if "DA005" in analyses:
        features.add("structure")
    if "DA006" in analyses:
        features.add("phylogeny")
    if "DA007" in analyses:
        features.add("spatial")
    if rtg:
        features.add("bench")
    if not features:
        return "core" if input_type == "FASTQ" else "aligned"
    if input_type == "FASTQ" or targeted:
        features.add("reads")
    return "selected-" + "-".join(sorted(features))


def source_names(source: Path) -> list[str]:
    names = {
        p.relative_to(source).as_posix()
        for p in source.rglob("*")
        if p.is_file()
        and not any(x in {".git", ".pixi", "__pycache__"} for x in p.relative_to(source).parts)
        and p.relative_to(source).parts[0] not in {"inputs", "runs", "results", "work", "prepared", "metadata"}
    }
    if (source / ".git").exists():
        r = subprocess.run(
            ["git", "-C", str(source), "ls-tree", "-r", "--name-only", "HEAD"], capture_output=True, text=True
        )
        if r.returncode == 0:
            names.update(r.stdout.splitlines())
    return sorted(names)


def export_workspace(
    source: Path,
    destination: Path,
    profile: str,
    design: str,
    analyses: list[str],
    rtg: bool = False,
    input_type: str = "FASTQ",
    final_output: str | None = None,
) -> dict:
    if input_type not in {"FASTQ", "BAM", "CRAM"}:
        raise ValueError("Select FASTQ, BAM or CRAM; VCF entry is not implemented")
    if final_output is not None:
        profile = profile_for(final_output, profile)
        if final_output == "FO001" and input_type != "FASTQ":
            raise ValueError("FO001 is FASTQ-only: BAM/CRAM enter at Layer 2; choose FO002 or later")
    if profile not in {"sample", "cohort", "downstream", "benchmark", "development", ON_DEMAND}:
        raise ValueError("Unknown profile")
    if profile == ON_DEMAND and (analyses or rtg or final_output):
        raise ValueError("An on-demand workspace takes its analyses and final output from the saved configuration")
    design = design_of(design)
    if set(analyses) - set(METHODS):
        raise ValueError("Unknown analysis; choose DA001–DA007")
    if profile in {"sample", "cohort"} and analyses:
        raise ValueError("Select downstream or benchmark profile for analyses")
    if profile == "downstream" and not analyses:
        raise ValueError("Choose --analyses or --auto-analyses for the downstream profile")
    if rtg and profile not in {"benchmark", "development"}:
        raise ValueError("--rtg requires benchmarking")
    if destination.exists():
        raise ValueError(f"Workspace exists: {destination}; choose a new directory")
    selected = set(METHODS) if profile == ON_DEMAND else set(analyses)
    if selected & {"DA004", "DA005", "DA006"}:
        selected.add("DA001")
    if "DA007" in selected:
        selected.add("DA002")
    modules = {METHODS[d] for d in selected}
    if selected & {"DA003", "DA001", "DA004", "DA005"}:
        modules.add("coi")  # diversity imports the shared FWS helper; no MOIRE tool unless requested/prerequisite
    if selected & {"DA001", "DA004", "DA005"}:
        modules.add("diversity")  # shared dominant-allele helper
    if profile == "development":
        env = "default"
    elif profile == ON_DEMAND:
        env = "control"  # the tools follow with ./pfvarx install-tools for the saved configuration
    else:
        env = environment_for(design, sorted(selected), rtg, input_type, final_output)
    mode = final_output or ("FO004" if selected else ("FO003" if profile == "cohort" else "FO002"))
    metadata = {
        "schema": "pfvarx.installation/2",
        "profile": profile,
        "design": design,
        "input_type": input_type,
        "environment": env,
        "analyses_requested": analyses,
        "analyses_available": sorted(selected),
        "max_final_output": mode,
        "purpose": "benchmarking" if profile == "benchmark" else "analysis",
        "rtg": rtg,
    }
    metadata["source_revision"] = None
    metadata["source_dirty"] = None
    if (source / ".git").exists():
        revision = subprocess.run(["git", "-C", str(source), "rev-parse", "HEAD"], capture_output=True, text=True)
        if revision.returncode == 0:
            metadata["source_revision"] = revision.stdout.strip()
            dirty = subprocess.run(
                ["git", "-C", str(source), "status", "--porcelain", "--untracked-files=no"],
                capture_output=True,
                text=True,
            )
            metadata["source_dirty"] = bool(dirty.stdout.strip())
    names = []
    for name in source_names(source):
        parts = Path(name).parts
        keep = name in ROOT_FILES or name.startswith(
            (
                "src/pfvarx/",
                "config/templates/",
                "config/parameters/",
                "config/registry/",
                "config/geo/",
                "docs/",
                "workflow/",
            )
        )
        if profile == "development":
            keep |= name.startswith(("tests/", "scripts/", ".github/", "containers/"))
        if not keep:
            continue
        if profile not in ("development", ON_DEMAND):  # on demand: all code, tools installed later
            if name.startswith("src/pfvarx/downstream/") and (
                not selected or Path(name).stem not in modules | {"__init__", "common", "run"}
            ):
                continue
            if name.startswith("src/pfvarx/annotation/") and not selected:
                continue
            if mode == "FO001" and (name.startswith("src/pfvarx/caller/") or name == "workflow/rules/l2_caller.smk"):
                continue  # no per-sample calling in an FO001 installation
            if name.startswith("src/pfvarx/cohort/") and mode not in {"FO003", "FO004"}:
                continue
            if name.startswith("src/pfvarx/benchmarking/") and profile != "benchmark":
                continue
            if name.startswith("config/geo/") and not selected & {"DA002", "DA007"}:
                continue
            if name.startswith("workflow/rules/l3_") and mode not in {"FO003", "FO004"}:
                continue
            if name.startswith(("workflow/rules/l4_", "workflow/rules/l5_")) and not selected:
                continue
            if name == "workflow/rules/l1_umi.smk" and design != "targeted_umi":
                continue
            if name.startswith("workflow/scripts/da") and Path(name).name not in {
                R_SCRIPTS[d] for d in selected if d in R_SCRIPTS
            }:
                continue
        if len(parts) > 1 and parts[0] == "src" and Path(name).suffix not in {".py", ".html", ".css", ".js", ".json"}:
            continue
        names.append(name)
    # Resolve all source content before exposing a partially written destination.
    content = {}
    for name in names:
        p = source / name
        if p.is_symlink():
            raise ValueError("Source symlink refused: " + name)
        if p.is_file():
            content[name] = p.read_bytes()
        else:
            content[name] = subprocess.run(
                ["git", "-C", str(source), "show", "HEAD:" + name], check=True, capture_output=True
            ).stdout
    if "src/pfvarx/cli.py" not in content:
        raise ValueError("Incomplete source: missing CLI")
    try:
        destination.mkdir(parents=True)
        for name, data in content.items():
            target = destination / name
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(data)
        (destination / "pfvarx").chmod(0o755)
        (destination / "config/installation.json").write_text(json.dumps(metadata, indent=2) + "\n")
    except BaseException:
        shutil.rmtree(destination, ignore_errors=True)
        raise
    return metadata


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--workspace", type=Path, default=ROOT.parent / "pfvarx-analysis",
                   help="New workspace path (default: pfvarx-analysis beside the source checkout)")
    p.add_argument("--final-output", choices=list(FINAL_OUTPUTS), default=None,
                   help="install now only what this final output needs (FO004 also needs --analyses or "
                        "--auto-analyses); without it, tools are installed on demand by ./pfvarx install-tools")
    p.add_argument(
        "--profile", choices=["sample", "cohort", "downstream", "benchmark", "development"], default=None,
        help="installation profile (default sample, or the one --final-output implies)"
    )
    p.add_argument("--design", type=design_of, default="wgs", help="wgs | targeted_umi | targeted_no_umi")
    selection = p.add_mutually_exclusive_group()
    selection.add_argument("--analyses", nargs="*", default=[])
    selection.add_argument("--auto-analyses", action="store_true")
    p.add_argument("--input-type", choices=["FASTQ", "BAM", "CRAM"], default="FASTQ")
    p.add_argument("--rtg", action="store_true")
    p.add_argument("--no-install", action="store_true")
    a = p.parse_args()
    try:
        destination = a.workspace.expanduser().resolve()
        if a.auto_analyses:
            a.analyses = sorted(METHODS)  # eligibility is decided per run from the panel and callable evidence
        profile = profile_for(a.final_output, a.profile, default=ON_DEMAND)
        meta = export_workspace(ROOT, destination, profile, a.design, a.analyses, a.rtg, a.input_type, a.final_output)
        if a.auto_analyses:
            meta["auto_analyses"] = True
            (destination / "config/installation.json").write_text(json.dumps(meta, indent=2) + "\n")
        if not a.no_install:
            subprocess.run(["pixi", "install", "--locked", "-e", meta["environment"]], cwd=destination, check=True)
            subprocess.run(
                ["pixi", "run", "--locked", "-e", meta["environment"], "pfvarx", "workspace-init", "--config-only"], cwd=destination, check=True
            )
        print(
            json.dumps(
                {
                    "workspace": str(destination),
                    "installation": meta,
                    "configuration_ready": not a.no_install,
                    "next": (
                        f"Run pixi install --locked -e {meta['environment']}, then ./pfvarx workspace-init --config-only"
                        if a.no_install
                        else ("Configure with ./pfvarx configure --gui, save, then click Install required tools "
                              "(or edit YAML, run ./pfvarx init, then ./pfvarx install-tools)"
                              if profile == ON_DEMAND else
                              "Configure with ./pfvarx configure --gui, or edit YAML and run ./pfvarx init")
                    ),
                },
                indent=2,
            )
        )
    except (ValueError, OSError, subprocess.CalledProcessError) as exc:
        p.exit(2, str(exc) + "\n")


if __name__ == "__main__":
    main()
