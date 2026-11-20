"""Optional Human DNA removal by competitive alignment (L0.4 and L1.2), wgs design only.

Targeted designs (targeted_no_umi, targeted_umi) are ineligible: their reads come from
parasite-specific amplicons or probes, and the configuration refuses host_removal for them.

When ``reference.host_fasta`` is set, ``prepare`` writes ``combined.fasta`` as the parasite
contigs followed by the host contigs, each host name prefixed with :data:`HOST_PREFIX`. Every
read pair is aligned once with bwa-mem2 to that combined reference, so a read that resembles
both genomes goes to the better match. :func:`filter_stream` then removes every read pair
whose primary alignment (either mate) is on a host contig, drops secondary and supplementary
records on host contigs from the pairs it keeps, removes host hits from their ``SA``/``XA``
tags and writes a header with the parasite contigs only. The parasite contigs come first in
the combined reference, so their reference IDs are unchanged and the output matches
``parasite.fasta`` exactly.

The bwa-mem2 index of a human-sized reference needs far more memory to build (tens of GB)
than most workstations have, so PfVarX never builds it: the user builds it once from
``prepared/reference/combined.fasta`` on a large machine and copies the index files back.
:func:`index_problems` checks that a supplied index was built from that exact FASTA.
"""

from __future__ import annotations

import argparse
import json
import sys
from collections.abc import Iterable, Iterator
from pathlib import Path
from typing import Any

HOST_PREFIX = "host_"
HOST_INDEX_COMMAND = "bwa-mem2 index combined.fasta"


def _headers(path: Path) -> Iterator[str]:
    with path.open("rb") as handle:
        for line in handle:
            if line.startswith(b">"):
                yield line[1:].split()[0].decode()


def write_combined(parasite: Path, host: Path, out: Path) -> dict[str, Any]:
    """Write parasite contigs, then host contigs renamed ``host_<name>``; returns counts."""
    parasite_names = list(_headers(parasite))
    clash = [n for n in parasite_names if n.startswith(HOST_PREFIX)]
    if clash:
        raise ValueError(f"parasite contig names may not start with {HOST_PREFIX!r}: {clash[:3]}")
    host_contigs = 0
    with out.open("wb") as dest:
        data = parasite.read_bytes().rstrip(b"\n")
        dest.write(data + b"\n")
        with host.open("rb") as source:
            for line in source:
                if line.startswith(b">"):
                    host_contigs += 1
                    line = b">" + HOST_PREFIX.encode() + line[1:].lstrip()
                dest.write(line if line.endswith(b"\n") else line + b"\n")
    if host_contigs == 0:
        raise ValueError(f"{host}: no FASTA contigs")
    return {"parasite_contigs": len(parasite_names), "host_contigs": host_contigs, "prefix": HOST_PREFIX}


def _fai(path: Path) -> list[tuple[str, int]]:
    rows = []
    for line in path.read_text().splitlines():
        cells = line.split("\t")
        rows.append((cells[0], int(cells[1])))
    return rows


def _ann(path: Path) -> list[tuple[str, int]]:
    """Contig names and lengths recorded in a bwa(-mem2) ``.ann`` file."""
    lines = path.read_text().splitlines()
    count = int(lines[0].split()[1])
    out = []
    for i in range(count):
        name = lines[1 + 2 * i].split()[1]
        length = int(lines[2 + 2 * i].split()[1])
        out.append((name, length))
    return out


def index_problems(combined: Path, suffixes: Iterable[str]) -> list[str]:
    """Why a supplied bwa-mem2 index cannot be used with ``combined`` (empty when it can)."""
    missing = [s for s in suffixes if not Path(f"{combined}{s}").is_file() or Path(f"{combined}{s}").stat().st_size == 0]
    if missing:
        return [f"bwa-mem2 index files missing: {', '.join('combined.fasta' + s for s in missing)}"]
    fai = Path(f"{combined}.fai")
    try:
        if _ann(Path(f"{combined}.ann")) != _fai(fai):
            return ["the bwa-mem2 index was built from a different FASTA than prepared/reference/combined.fasta "
                    "(contig names or lengths differ)"]
    except (OSError, ValueError, IndexError) as exc:
        return [f"cannot read combined.fasta.ann or combined.fasta.fai: {exc}"]
    return []


def index_instructions(prep: Path) -> list[str]:
    """Numbered steps for building the host index on another machine."""
    return [
        f"Copy {prep / 'combined.fasta'} to a machine with enough memory (building a human-sized bwa-mem2 "
        "index is commonly reported to need about 64 GB of RAM).",
        f"There, in the folder holding the copy, run: {HOST_INDEX_COMMAND}  (bwa-mem2 2.2.1, as pinned here).",
        "Copy the five index files back next to the original: combined.fasta.0123, .amb, .ann, "
        f".bwt.2bit.64 and .pac, into {prep}/.",
        "Run ./pfvarx prepare again: it checks that the index matches combined.fasta and seals it.",
    ]


def _clean_hits(value: str) -> str:
    """Drop host entries from an SA (``;``-separated) or XA tag."""
    return ";".join(h for h in value.split(";") if h and not h.startswith(HOST_PREFIX)) + (";" if value.endswith(";") else "")


def filter_stream(source: Iterable[Any], write: Any, n_parasite: int) -> dict[str, int]:
    """Filter name-grouped alignments (bwa-mem2 output order); returns template counts."""
    counts = {"templates": 0, "host_templates": 0, "kept_templates": 0, "host_secondary_records_dropped": 0}
    group: list[Any] = []

    def flush() -> None:
        if not group:
            return
        counts["templates"] += 1
        primaries = [r for r in group if not r.is_secondary and not r.is_supplementary]
        if any(not r.is_unmapped and r.reference_id >= n_parasite for r in primaries):
            counts["host_templates"] += 1
            return
        counts["kept_templates"] += 1
        for r in group:
            if not r.is_unmapped and r.reference_id >= n_parasite:
                counts["host_secondary_records_dropped"] += 1
                continue
            for tag in ("SA", "XA"):
                if r.has_tag(tag):
                    cleaned = _clean_hits(str(r.get_tag(tag)))
                    r.set_tag(tag, cleaned if cleaned.strip(";") else None)
            write(r)

    for record in source:
        if group and record.query_name != group[0].query_name:
            flush()
            group = []
        group.append(record)
    flush()
    return counts


def parasite_header(header: dict[str, Any]) -> tuple[dict[str, Any], int]:
    """The input header with host @SQ lines removed; host contigs must follow every parasite contig."""
    sq = header.get("SQ", [])
    n = sum(1 for s in sq if not s["SN"].startswith(HOST_PREFIX))
    if any(s["SN"].startswith(HOST_PREFIX) for s in sq[:n]):
        raise ValueError("host contigs must follow the parasite contigs in the combined reference")
    out = dict(header)
    out["SQ"] = sq[:n]
    out.setdefault("PG", []).append({"ID": "pfvarx.host_filter", "PN": "pfvarx host filter",
                                     "CL": "remove read pairs aligned to host contigs (" + HOST_PREFIX + "*)"})
    return out, n


def main(argv: list[str] | None = None) -> int:
    """``python -m pfvarx.references.host --counts FILE``: SAM/BAM on stdin, uncompressed BAM on stdout."""
    import pysam

    parser = argparse.ArgumentParser(description=(__doc__ or "").splitlines()[0])
    parser.add_argument("--counts", required=True, type=Path, help="JSON file for template counts")
    parser.add_argument("--input", default="-", help="name-grouped SAM/BAM from bwa-mem2 (default stdin)")
    parser.add_argument("--output", default="-", help="output BAM (default stdout)")
    args = parser.parse_args(argv)
    source = pysam.AlignmentFile(args.input, "r")
    header, n = parasite_header(source.header.to_dict())
    with pysam.AlignmentFile(args.output, "wbu", header=header) as sink:
        counts = filter_stream(source, sink.write, n)
    counts["host_fraction"] = counts["host_templates"] / counts["templates"] if counts["templates"] else 0.0  # type: ignore[assignment]
    args.counts.write_text(json.dumps(counts, indent=2) + "\n")
    print(f"host filter: {json.dumps(counts)}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
