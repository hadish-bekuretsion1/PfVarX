"""Recognise input files by their content and structure, never by their name.

Every route that accepts an input file (GUI upload and download, ``add-input``,
``fetch-input``, registration, folder settings in ``config.yaml`` and validation) asks this
module what a file *is*. The extension does not decide anything, with one exception that the
tools force on us: fastp reads a FASTQ as gzip only when its name ends in ``.gz`` (a
gzip-compressed FASTQ without it yields zero reads, silently), so compressed content and the
``.gz`` suffix must agree for reads.
"""

from __future__ import annotations

import gzip
import zlib
from pathlib import Path

# Content kinds recognised here, with the plain-language name used in messages.
KINDS = {
    "FASTA": "FASTA sequence",
    "FASTQ": "FASTQ reads",
    "BAM": "BAM alignment",
    "CRAM": "CRAM alignment",
    "ALIGNMENT_INDEX": "alignment index (.bai/.csi/.crai)",
    "BED": "BED intervals",
    "GFF": "GFF3 annotation",
    "VCF": "VCF variants",
    "SAMPLE_SHEET": "sample sheet (TSV with sample_id)",
    "PRIMERS": "primer table (chrom, left_start, left_end, right_start, right_end)",
    "CATALOGUE": "resistance marker catalogue (TSV with marker_id)",
    "LIBRARY_CATALOG": "library catalog (YAML with schema_version pfvarx.library_catalog/1)",
    "TSV": "tab-separated table",
    "JSON": "JSON",
    "EMPTY": "empty file",
    "UNKNOWN": "unrecognised content",
}

# What each configuration setting must contain.
SETTING_KINDS = {
    "reference.fasta": ("FASTA",),
    "reference.core_mask_bed": ("BED",),
    "reference.gff": ("GFF",),
    "reference.host_fasta": ("FASTA",),
    "assay.library_catalog": ("LIBRARY_CATALOG",),
    "samples_tsv": ("SAMPLE_SHEET",),
    "analyses.resistance_catalog": ("CATALOGUE",),
}

# What each input folder may hold (GUI upload/download, add-input, fetch-input).
FOLDER_KINDS = {
    "reference": ("FASTA", "BED", "GFF"),
    "targets": ("BED", "PRIMERS"),
    "metadata": ("SAMPLE_SHEET", "CATALOGUE", "LIBRARY_CATALOG", "TSV"),
    "fastq": ("FASTQ",),
    "bam": ("BAM", "CRAM", "ALIGNMENT_INDEX"),  # BAM/CRAM input; results/alignments/ holds the pipeline's own
    "truth": ("VCF", "BED", "JSON"),
}
ALIGNMENT_KINDS = {"BAM": ("BAM", "ALIGNMENT_INDEX"), "CRAM": ("CRAM", "ALIGNMENT_INDEX")}

PRIMER_HEADER = ["chrom", "left_start", "left_end", "right_start", "right_end"]
# Inputs the pipeline reads or copies as plain text (FASTA is copied byte for byte before indexing).
PLAIN_TEXT_ONLY = {"FASTA", "BED", "GFF", "PRIMERS", "SAMPLE_SHEET", "CATALOGUE", "LIBRARY_CATALOG"}
_SAMPLE = 65536


def _head(path: Path) -> tuple[bytes, bool]:
    """Up to 64 KiB of decompressed content, and whether the file is gzip/BGZF compressed."""
    with path.open("rb") as handle:
        raw = handle.read(_SAMPLE)
    if raw[:2] != b"\x1f\x8b":
        return raw, False
    out = b""
    decompressor = zlib.decompressobj(31)
    data = raw
    try:
        while data and len(out) < _SAMPLE:  # BGZF is a series of gzip members
            out += decompressor.decompress(data, _SAMPLE - len(out))
            data = decompressor.unused_data if decompressor.eof else decompressor.unconsumed_tail
            if decompressor.eof:
                decompressor = zlib.decompressobj(31)
    except zlib.error:
        pass
    return out, True


def is_compressed(path: Path) -> bool:
    with path.open("rb") as handle:
        return handle.read(2) == b"\x1f\x8b"


def _lines(data: bytes) -> list[str]:
    text = data.decode("ascii", errors="replace")
    lines = text.splitlines()
    if len(data) >= _SAMPLE and lines:
        lines = lines[:-1]  # the last line may be cut by the sample size
    return lines


def _int(value: str) -> bool:
    return value.strip().isdigit()


def _is_bed(rows: list[list[str]]) -> bool:
    return bool(rows) and all(len(r) >= 3 and _int(r[1]) and _int(r[2]) and int(r[1]) < int(r[2])
                              for r in rows[:200])


def _is_gff(rows: list[list[str]]) -> bool:
    return bool(rows) and all(len(r) == 9 and _int(r[3]) and _int(r[4]) and r[6] in "+-.?" for r in rows[:200])


def detect(path: Path) -> str:
    """The content kind of *path* (a key of ``KINDS``), from its bytes alone."""
    if path.stat().st_size == 0:
        return "EMPTY"
    with path.open("rb") as handle:
        magic = handle.read(4)
    if magic == b"CRAM":
        return "CRAM"
    if magic == b"BAI\x01":
        return "ALIGNMENT_INDEX"
    data, compressed = _head(path)
    if data[:4] == b"BAM\x01":
        return "BAM"
    if data[:4] == b"CSI\x01":
        return "ALIGNMENT_INDEX"
    lines = [line for line in _lines(data) if line.strip()]
    if not lines:
        return "UNKNOWN"
    first = lines[0]
    if first.startswith(">"):
        return "FASTA"
    if first.startswith("@") and len(lines) >= 3 and lines[2].startswith("+"):
        return "FASTQ"
    if first.startswith("##fileformat=VCF"):
        return "VCF"
    if first.lstrip().startswith(("{", "[")):
        return "JSON"
    if first.startswith("##gff-version"):
        return "GFF"
    if any(line.split("#", 1)[0].strip().replace('"', "").replace("'", "")
           == "schema_version: pfvarx.library_catalog/1" for line in lines[:50]):
        return "LIBRARY_CATALOG"
    body = [line for line in lines if not line.startswith(("#", "track", "browser"))]
    rows = [line.split("\t") for line in body]
    if compressed and rows and all(len(r) == 6 and all(_int(c) for c in r) for r in rows[:50]):
        return "ALIGNMENT_INDEX"  # CRAM index: gzip-compressed rows of six integers
    if _is_gff(rows):
        return "GFF"
    if _is_bed(rows):
        return "BED"
    if rows:
        header = [c.strip() for c in rows[0]]
        if header[:5] == PRIMER_HEADER:
            return "PRIMERS"
        if "marker_id" in header:
            return "CATALOGUE"
        if "sample_id" in header:
            return "SAMPLE_SHEET"
        if len(header) > 1:
            return "TSV"
    return "UNKNOWN"


def describe(kinds: tuple[str, ...]) -> str:
    return " or ".join(KINDS[k] for k in kinds)


def content_problem(path: Path, accepted: tuple[str, ...], name: str | None = None) -> str | None:
    """Why *path* cannot be used where *accepted* content is required, or None if it can.

    *name* is the file name the content will have (for a staged upload or download); it
    matters only for the FASTQ .gz rule that fastp imposes."""
    shown = name or path.name
    kind = detect(path)
    if kind not in accepted:
        return f"{shown} contains {KINDS[kind]}; this input needs {describe(accepted)}"
    if kind in PLAIN_TEXT_ONLY and is_compressed(path):
        return (f"{shown} is gzip-compressed {KINDS[kind]}; PfVarX reads this input as plain text. "
                "Decompress it (gunzip) and use the uncompressed file")
    if kind == "FASTQ":
        gz_name = shown.lower().endswith(".gz")
        if is_compressed(path) and not gz_name:
            return (f"{shown} is gzip-compressed but its name does not end in .gz; fastp reads it as plain "
                    "text and finds no reads. Rename it to end in .gz")
        if not is_compressed(path) and gz_name:
            return f"{shown} ends in .gz but is not compressed; rename it without .gz or compress it with gzip"
    return None


def require(path: Path, accepted: tuple[str, ...], name: str | None = None) -> None:
    problem = content_problem(path, accepted, name)
    if problem:
        raise ValueError(problem)


def open_text(path: Path):  # noqa: ANN201 - text handle for plain or gzip content
    """Open a text input, decompressing when the *content* is gzip, whatever its name."""
    if is_compressed(path):
        return gzip.open(path, "rt", encoding="ascii", errors="replace")
    return path.open("rt", encoding="ascii", errors="replace")


def find_in_folder(folder: Path, accepted: tuple[str, ...]) -> list[Path]:
    """Files in *folder* whose content is one of *accepted*, by name order."""
    if not folder.is_dir():
        return []
    found = []
    for path in sorted(folder.iterdir()):
        if path.is_file() and not path.name.startswith("."):
            try:
                if detect(path) in accepted:
                    found.append(path)
            except OSError:
                continue
    return found
