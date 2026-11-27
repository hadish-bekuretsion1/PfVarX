"""Build a deterministic software-only source archive from a PfVarX checkout."""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ROOT_FILES = ("README.md", "CHANGELOG.md", "LICENSE", "bootstrap.py", "pfvarx", "pyproject.toml", "pixi.toml", "pixi.lock", ".gitignore")
SOP_HTML = ("docs/sop/PfVarX_SOP_Manual.html", "docs/README.html", "docs/configuration-guide.html",
            "docs/PfVarX_Sample_Analysis_Guide.html", "docs/PfVarX_Benchmarking_Guide.html")
REQUIRED = (*ROOT_FILES, *SOP_HTML, ".github/workflows/ci.yaml", "src/pfvarx/cli.py", "workflow/Snakefile",
            "config/templates/analysis.yaml", "config/templates/benchmarking.yaml",
            "config/templates/quality.yaml", "config/templates/compute.yaml",
            "docs/analysis/README.md", "docs/benchmarking/README.md", "scripts/package.py", "scripts/verify.py")
SUFFIXES = {"src": {".py", ".html", ".css", ".js", ".json"}, "workflow": {".smk", ".py", ".R"},
            "tests": {".py", ".txt"}, "scripts": {".py"}, "docs": {".md"},
            "config/templates": {".yaml"}, "config/parameters": {".yaml"},
            "config/registry": {".yaml", ".tsv", ".bed"}, "config/geo": {".geojson", ".md"}}
CACHES = {"__pycache__", ".pytest_cache", ".mypy_cache", ".ruff_cache", ".pixi", ".snakemake", ".git"}
# Development-only material (guides, notes, experiment scripts) is never part of a release.
DEVELOPMENT_ONLY = {"development"}
RUNTIME = {"inputs", "prepared", "metadata", "runs", "results", "benchmarks", "work", "logs", "dist", "build"}


def source_files(root: Path) -> list[Path]:
    """Select software-owned paths; runtime paths and user config are never exported."""
    missing = [name for name in REQUIRED if not (root / name).is_file()]
    if missing:
        raise ValueError(f"Incomplete software tree; missing: {missing}")
    paths = {root / name for name in (*ROOT_FILES, *SOP_HTML, ".github/workflows/ci.yaml", "workflow/Snakefile")}
    for directory, suffixes in SUFFIXES.items():
        for path in (root / directory).rglob("*"):
            if path.suffix in suffixes and not any(p in CACHES for p in path.relative_to(root).parts):
                if path.is_file() or path.is_symlink():
                    paths.add(path)
    leaked = sorted(p.relative_to(root).as_posix() for p in paths if p.relative_to(root).parts[0] in DEVELOPMENT_ONLY)
    if leaked:
        raise ValueError(f"Development-only files selected for release: {leaked}")
    for path in paths:
        # Reject external references even if their filename looks like source code.
        if path.is_symlink() or not path.resolve().is_relative_to(root.resolve()):
            raise ValueError(f"Source symlinks are not distributable: {path.relative_to(root)}")
    return sorted(paths, key=lambda p: p.relative_to(root).as_posix())


def check_pristine(root: Path) -> dict:
    selected = {p.relative_to(root).as_posix() for p in source_files(root)}
    present = {p.relative_to(root).as_posix() for p in root.rglob("*") if p.is_file() or p.is_symlink()}
    extras = sorted(present - selected)
    bad_dirs = sorted(p.relative_to(root).as_posix() for p in root.rglob("*") if p.is_dir() and
                      (p.name in CACHES or (p.parent == root and p.name in RUNTIME)))
    if extras or bad_dirs:
        raise ValueError(f"Tree is not pristine; extra files={extras}, runtime/cache directories={bad_dirs}")
    return {"status": "PASS", "software_files": len(selected), "runtime_files": 0}


def pipeline_version(root: Path) -> str:
    """The authoritative release version, read from the package without importing it."""
    match = re.search(r'^__version__\s*=\s*"([^"]+)"', (root / "src/pfvarx/__init__.py").read_text(), re.MULTILINE)
    if match is None:
        raise ValueError("Missing __version__ in src/pfvarx/__init__.py")
    return match[1]


def build_archive(root: Path, output: Path) -> dict:
    files = source_files(root)
    version = pipeline_version(root)
    output.mkdir(parents=True, exist_ok=True)
    archive = output / f"pfvarx-{version}-source.zip"
    temporary = archive.with_suffix(".zip.tmp")
    with zipfile.ZipFile(temporary, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as bundle:
        for path in files:
            name = "pfvarx/" + path.relative_to(root).as_posix()
            info = zipfile.ZipInfo(name, date_time=(2020, 1, 1, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.create_system = 3
            info.external_attr = (0o100755 if path.name == "pfvarx" else 0o100644) << 16
            bundle.writestr(info, path.read_bytes())
    temporary.replace(archive)
    digest = hashlib.sha256(archive.read_bytes()).hexdigest()
    checksum = archive.with_suffix(".zip.sha256")
    checksum.write_text(f"{digest}  {archive.name}\n")
    return {"archive": str(archive), "sha256": digest, "checksum": str(checksum), "files": len(files)}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "dist")
    parser.add_argument("--check", action="store_true", help="assert that this checkout contains only source files")
    args = parser.parse_args()
    result = check_pristine(ROOT) if args.check else build_archive(ROOT, args.output.resolve())
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
