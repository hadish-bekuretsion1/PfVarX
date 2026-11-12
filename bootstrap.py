"""Export/install a minimal configured workspace; Python 3.11+, no Python packages needed.

Works from a full source extraction or a sparse Git clone. No code runs merely by cloning.
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


def environment_for(design: str, analyses: list[str], rtg: bool = False, input_type: str = "FASTQ") -> str:
    features = set()
    if design != "WGS" and input_type == "FASTQ":
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
    if input_type == "FASTQ":
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
) -> dict:
    if input_type not in {"FASTQ", "BAM", "CRAM"}:
        raise ValueError("Select FASTQ, BAM or CRAM; VCF entry is not implemented")
    if profile not in {"sample", "cohort", "downstream", "benchmark", "development"}:
        raise ValueError("Unknown profile")
    if design not in {"WGS", "TARGETED_NON_UMI", "TARGETED_UMI"}:
        raise ValueError("Unknown design")
    if set(analyses) - set(METHODS):
        raise ValueError("Unknown analysis; choose DA001–DA007")
    if profile in {"sample", "cohort"} and analyses:
        raise ValueError("Select downstream or benchmark profile for analyses")
    if profile == "downstream" and not analyses:
        raise ValueError("downstream requires --analyses")
    if design != "WGS" and set(analyses) & {"DA005", "DA006"}:
        raise ValueError("DA005/DA006 require WGS")
    if rtg and profile not in {"benchmark", "development"}:
        raise ValueError("--rtg requires benchmarking")
    if destination.exists():
        raise ValueError(f"Workspace exists: {destination}; choose a new directory")
    selected = set(analyses)
    if selected & {"DA004", "DA005", "DA006"}:
        selected.add("DA001")
    if "DA007" in selected:
        selected.add("DA002")
    modules = {METHODS[d] for d in selected}
    if selected & {"DA003", "DA001", "DA004", "DA005"}:
        modules.add("coi")  # diversity imports the shared FWS helper; no MOIRE tool unless requested/prerequisite
    if selected & {"DA001", "DA004", "DA005"}:
        modules.add("diversity")  # shared dominant-allele helper
    env = "default" if profile == "development" else environment_for(design, sorted(selected), rtg, input_type)
    mode = "FO004" if selected else ("FO003" if profile == "cohort" else "FO002")
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
        if profile != "development":
            if name.startswith("src/pfvarx/downstream/") and (
                not selected or Path(name).stem not in modules | {"__init__", "common", "run"}
            ):
                continue
            if name.startswith("src/pfvarx/annotation/") and not selected:
                continue
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
            if name == "workflow/rules/l1_umi.smk" and (design != "TARGETED_UMI" or input_type != "FASTQ"):
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
    p.add_argument(
        "--profile", choices=["sample", "cohort", "downstream", "benchmark", "development"], default="sample"
    )
    p.add_argument("--design", choices=["WGS", "TARGETED_NON_UMI", "TARGETED_UMI"], default="WGS")
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
            a.analyses = [da for da in sorted(METHODS) if a.design == "WGS" or da not in {"DA005", "DA006"}]
        meta = export_workspace(ROOT, destination, a.profile, a.design, a.analyses, a.rtg, a.input_type)
        if a.auto_analyses:
            meta["auto_analyses"] = True
            (destination / "config/installation.json").write_text(json.dumps(meta, indent=2) + "\n")
        if not a.no_install:
            subprocess.run(["pixi", "install", "--locked", "-e", meta["environment"]], cwd=destination, check=True)
            subprocess.run(
                ["pixi", "run", "-e", meta["environment"], "pfvarx", "workspace-init", "--config-only"], cwd=destination, check=True
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
                        else "Configure with ./pfvarx configure --gui, or edit YAML and run ./pfvarx init"
                    ),
                },
                indent=2,
            )
        )
    except (ValueError, OSError, subprocess.CalledProcessError) as exc:
        p.exit(2, str(exc) + "\n")


if __name__ == "__main__":
    main()
