"""Input-file transfer rules shared by the configuration GUI and the CLI.

Both interfaces place a file in ``inputs/<role>/`` through the same checks: the folder must
belong to the saved configuration and already exist, the filename must be plain and carry an
extension accepted for that folder, and an existing destination is never overwritten.
"""

from __future__ import annotations

import os
import re
import shutil
import tempfile
from pathlib import Path
from typing import Any

from pfvarx.configuration.loader import load_bundle
from pfvarx.initialization.layout import required_dirs
from pfvarx.inputs.formats import FOLDER_FORMATS, check_extension
from pfvarx.inputs.organize import fetch
from pfvarx.io import sha256, within

ALIGNMENT_FORMATS = {"BAM": (".bam", ".bai", ".csi"), "CRAM": (".cram", ".crai", ".csi")}


def folder_extensions(root: Path) -> dict[str, tuple[str, ...]]:
    """Accepted extensions per input folder for the saved configuration."""
    cfg = load_bundle(root).config
    extensions = dict(FOLDER_FORMATS)
    if cfg.input_type.value in ALIGNMENT_FORMATS:
        extensions["alignments"] = ALIGNMENT_FORMATS[cfg.input_type.value]
    return extensions


def destination_for(root: Path, name: str, role: str = "raw") -> Path:
    """Checked destination of one input file; ``raw`` means the reads or alignments folder."""
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", name):
        raise ValueError("Use a plain filename with letters, digits, dots, underscores or hyphens")
    cfg = load_bundle(root).config
    if role == "raw":
        role = "fastq" if cfg.input_type.value == "FASTQ" else "alignments"
    if f"inputs/{role}" not in required_dirs(cfg) or role == "partitions":
        raise ValueError("Input folder is unavailable for the saved configuration")
    check_extension(name, folder_extensions(root)[role], role)
    folder = within(root, f"inputs/{role}")
    if not folder.is_dir():
        raise ValueError("Save the configuration and create input directories first")
    destination = folder / name
    if destination.exists() or destination.is_symlink():
        raise ValueError("Destination exists; refusing overwrite")
    return destination


def download(root: Path, url: str, expected: str | None, role: str, name: str) -> dict[str, Any]:
    """Download one HTTPS resource to its checked destination."""
    destination = destination_for(root, name, role)
    return fetch(root, url, expected, destination.parent.name, destination.name)


def import_file(root: Path, source: Path, role: str = "raw", name: str | None = None) -> dict[str, Any]:
    """Copy one file from disk into its input folder (the CLI counterpart of a GUI upload)."""
    source = source.expanduser()
    if not source.is_file():
        raise ValueError(f"Input file not found: {source}")
    if source.stat().st_size == 0:
        raise ValueError(f"Input file is empty: {source}")
    destination = destination_for(root, name or source.name, role)
    with tempfile.NamedTemporaryFile(dir=destination.parent, prefix=".import-") as staged:
        with source.open("rb") as handle:
            shutil.copyfileobj(handle, staged, 1024 * 1024)
        staged.flush()
        os.link(staged.name, destination)  # atomic and refuses a racing existing destination
    destination.chmod(0o644)
    return {"path": str(destination.relative_to(root)), "bytes": destination.stat().st_size,
            "sha256": sha256(destination), "source": str(source.resolve())}
