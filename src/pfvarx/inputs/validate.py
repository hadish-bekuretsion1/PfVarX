"""Read-only input validation (L0.4). Nothing here writes or modifies a source file.

Each check returns ``Finding`` records with a stable reason code, the file and,
where meaningful, the record number, so that a malformed input is localised.
"""

from __future__ import annotations

import gzip
import hashlib
import re
from collections.abc import Iterator
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import IO, Any

from pfvarx.constants import Design, InputType, Reason
from pfvarx.inputs.content import content_problem, open_text
from pfvarx.io import read_tsv

PHRED33_MIN, PHRED33_MAX = 33, 33 + 93


@dataclass
class Finding:
    level: str  # FAIL | WARNING
    reason: str
    path: str
    detail: str
    record: int | None = None


@dataclass
class FileReport:
    path: str
    kind: str
    status: str = "PASS"
    findings: list[Finding] = field(default_factory=list)
    facts: dict[str, Any] = field(default_factory=dict)

    def fail(self, reason: Reason, detail: str, record: int | None = None) -> None:
        self.findings.append(Finding("FAIL", reason.value, self.path, detail, record))
        self.status = "FAIL"

    def warn(self, reason: Reason, detail: str) -> None:
        self.findings.append(Finding("WARNING", reason.value, self.path, detail))
        if self.status == "PASS":
            self.status = "WARNING"

    def as_dict(self) -> dict[str, Any]:
        d = asdict(self)
        return d


def _open(path: Path) -> IO[str]:
    """Text handle; gzip is recognised from the content, not from the file name."""
    return open_text(path)  # type: ignore[no-any-return]


def _fastq_records(handle: IO[str]) -> Iterator[tuple[str, str, str]]:
    while True:
        header = handle.readline()
        if not header:
            return
        seq = handle.readline()
        plus = handle.readline()
        qual = handle.readline()
        if not qual:
            raise EOFError("truncated record")
        if not header.startswith("@") or not plus.startswith("+"):
            raise ValueError("record does not follow @/+ structure")
        seq, qual = seq.rstrip("\n"), qual.rstrip("\n")
        if len(seq) != len(qual):
            raise ValueError("sequence and quality lengths differ")
        yield header[1:].split()[0], seq, qual


def _mate_key(name: str) -> str:
    return re.sub(r"/[12]$", "", name)


def validate_fastq_pair(r1: Path, r2: Path, encoding_records: int) -> FileReport:
    rep = FileReport(path=f"{r1} | {r2}", kind="FASTQ_PAIR")
    for p in (r1, r2):
        if not p.is_file():
            rep.fail(Reason.FILE_MISSING, f"{p} does not exist")
        elif problem := content_problem(p, ("FASTQ",)):
            rep.fail(Reason.CONTENT_MISMATCH, problem)
    if rep.status == "FAIL":
        return rep
    n = 0
    qmin, qmax = 255, 0
    try:
        with _open(r1) as h1, _open(r2) as h2:
            it1, it2 = _fastq_records(h1), _fastq_records(h2)
            while True:
                a = next(it1, None)
                b = next(it2, None)
                if a is None and b is None:
                    break
                n += 1
                if a is None or b is None:
                    rep.fail(Reason.FASTQ_MATE_COUNT_MISMATCH,
                             f"{'R1' if a is None else 'R2'} ended first", n)
                    break
                if _mate_key(a[0]) != _mate_key(b[0]):
                    rep.fail(Reason.FASTQ_MATE_NAME_MISMATCH, f"{a[0]} vs {b[0]}", n)
                    break
                if n <= encoding_records:
                    for q in (a[2], b[2]):
                        codes = [ord(c) for c in q]
                        if codes:  # an empty read (already trimmed away) carries no quality codes
                            qmin, qmax = min(qmin, min(codes)), max(qmax, max(codes))
    except EOFError:
        rep.fail(Reason.FASTQ_TRUNCATED, "file ends inside a record", n + 1)
    except (ValueError, gzip.BadGzipFile, OSError) as exc:
        rep.fail(Reason.FASTQ_MALFORMED_RECORD, str(exc), n + 1)
    if n and qmin <= qmax and (qmin < PHRED33_MIN or qmax > PHRED33_MAX):
        rep.fail(Reason.FASTQ_QUALITY_ENCODING, f"quality codes {qmin}-{qmax} outside Phred+33")
    if n == 0 and rep.status != "FAIL":
        rep.fail(Reason.FASTQ_MALFORMED_RECORD, "no records")
    rep.facts = {"read_pairs": n, "quality_code_range": [qmin, qmax]}
    return rep


def fasta_contigs(path: Path) -> dict[str, dict[str, Any]]:
    """Return {name: {length, md5}} using the SAM M5 convention (upper-case, no whitespace)."""
    contigs: dict[str, dict[str, Any]] = {}
    name: str | None = None
    md5 = hashlib.md5()
    length = 0
    with _open(path) as handle:
        for line in handle:
            line = line.rstrip("\n\r")
            if line.startswith(">"):
                if name is not None:
                    contigs[name] = {"length": length, "md5": md5.hexdigest()}
                name = line[1:].split()[0]
                if name in contigs:
                    raise ValueError(f"duplicate contig {name}")
                md5, length = hashlib.md5(), 0
            elif line:
                if name is None:
                    raise ValueError("sequence before first header")
                seq = re.sub(r"\s", "", line).upper()
                if re.search(r"[^ACGTNRYKMSWBDHV*-]", seq):
                    raise ValueError(f"invalid characters in {name}")
                md5.update(seq.encode())
                length += len(seq)
    if name is not None:
        contigs[name] = {"length": length, "md5": md5.hexdigest()}
    if not contigs:
        raise ValueError("no contigs")
    return contigs


def validate_fasta(path: Path) -> FileReport:
    rep = FileReport(path=str(path), kind="FASTA")
    if not path.is_file():
        rep.fail(Reason.FILE_MISSING, "reference FASTA not found")
        return rep
    try:
        rep.facts["contigs"] = fasta_contigs(path)
    except ValueError as exc:
        rep.fail(Reason.REFERENCE_MALFORMED, str(exc))
    return rep


def validate_bed(path: Path, contigs: dict[str, dict[str, Any]], kind: str = "BED") -> FileReport:
    rep = FileReport(path=str(path), kind=kind)
    if not path.is_file():
        rep.fail(Reason.FILE_MISSING, "BED not found")
        return rep
    n = 0
    total = 0
    names: list[str] = []
    intervals: list[tuple[str, int, int]] = []
    for i, line in enumerate(path.read_text().splitlines(), start=1):
        if not line.strip() or line.startswith(("#", "track", "browser")):
            continue
        cells = line.split("\t")
        if len(cells) > 3 and len(names) < 200:
            names.append(cells[3])
        try:
            chrom, start, end = cells[0], int(cells[1]), int(cells[2])
        except (IndexError, ValueError):
            rep.fail(Reason.BED_MALFORMED, "need chrom<TAB>start<TAB>end integers", i)
            continue
        if chrom not in contigs:
            rep.fail(Reason.REFERENCE_MISMATCH, f"contig {chrom} not in reference", i)
        elif not 0 <= start < end <= contigs[chrom]["length"]:
            rep.fail(Reason.BED_MALFORMED, "0-based half-open interval out of range", i)
        n += 1
        total += max(0, end - start)
        intervals.append((chrom, start, end))
    if n == 0 and rep.status == "PASS":
        rep.fail(Reason.BED_MALFORMED, "no intervals")
    rep.facts = {"intervals": n, "bases": total}
    if kind == "CORE_MASK" and n and rep.status != "FAIL":
        concerns = core_mask_concerns(intervals, names, contigs)
        if concerns:
            rep.warn(Reason.CORE_MASK_SUSPECT, "; ".join(concerns))
    return rep


GENE_LIKE = re.compile(r"(PF3D7_\d{7}|PF3D7_[A-Z]{3}\d+|exon_|-CDS|\.\d+-p\d|gene|mRNA)", re.IGNORECASE)


def core_mask_concerns(intervals: list[tuple[str, int, int]], names: list[str],
                       contigs: dict[str, dict[str, Any]]) -> list[str]:
    """Signs that a 'core mask' is really a gene or feature list. A region mask is a few dozen
    non-overlapping intervals covering most of the genome (about 89% of 3D7 for the standard
    core of Miles et al. 2016); a gene BED has thousands of named, often overlapping intervals."""
    concerns = []
    if names and sum(bool(GENE_LIKE.search(x)) for x in names) >= max(1, len(names) // 2):
        concerns.append(f"column 4 holds gene or feature names (for example {names[0]!r}): this looks like a gene "
                        "BED, not a core-genome mask")
    if len(intervals) > 1000:
        concerns.append(f"{len(intervals)} intervals: a core-genome mask normally has fewer than 200")
    last: dict[str, int] = {}
    for chrom, start, end in sorted(intervals):
        if start < last.get(chrom, -1):
            concerns.append("intervals overlap: a region mask has non-overlapping intervals")
            break
        last[chrom] = max(last.get(chrom, -1), end)
    genome = sum(int(c["length"]) for c in contigs.values())
    covered = sum(end - start for chrom, start, end in merge_intervals(intervals))
    if genome and not 0.70 <= covered / genome <= 0.97:
        concerns.append(f"covers {covered / genome:.0%} of the reference; the standard 3D7 core covers about 89%")
    return concerns


def merge_intervals(intervals: list[tuple[str, int, int]]) -> list[tuple[str, int, int]]:
    merged: list[tuple[str, int, int]] = []
    for chrom, start, end in sorted(intervals):
        if merged and merged[-1][0] == chrom and start <= merged[-1][2]:
            merged[-1] = (chrom, merged[-1][1], max(merged[-1][2], end))
        else:
            merged.append((chrom, start, end))
    return merged


PRIMER_HEADER = ["chrom", "left_start", "left_end", "right_start", "right_end"]


def validate_gff(path: Path, contigs: dict[str, dict[str, Any]]) -> FileReport:
    """GFF3 structure: nine tab-separated columns, reference contigs, coordinates inside them,
    and the gene/CDS features that annotation and resistance codons need."""
    rep = FileReport(path=str(path), kind="GFF")
    if not path.is_file():
        rep.fail(Reason.FILE_MISSING, "GFF not found")
        return rep
    features: dict[str, int] = {}
    unknown: set[str] = set()
    with _open(path) as handle:
        for i, line in enumerate(handle, start=1):
            if line.startswith("##FASTA"):
                break
            if not line.strip() or line.startswith("#"):
                continue
            cells = line.rstrip("\n").split("\t")
            if len(cells) != 9 or not cells[3].isdigit() or not cells[4].isdigit():
                rep.fail(Reason.GFF_MALFORMED, "need nine tab-separated columns with integer start and end", i)
                break
            seqid, start, end = cells[0], int(cells[3]), int(cells[4])
            if seqid not in contigs:
                unknown.add(seqid)
            elif not 1 <= start <= end <= contigs[seqid]["length"]:
                rep.fail(Reason.GFF_MALFORMED, f"{seqid}:{start}-{end} is outside the contig (GFF is one-based)", i)
                break
            features[cells[2]] = features.get(cells[2], 0) + 1
    if unknown:
        rep.fail(Reason.REFERENCE_MISMATCH, "contigs not in the reference FASTA: " + ", ".join(sorted(unknown)[:5]))
    # PlasmoDB types genes as protein_coding_gene, ncRNA_gene or pseudogene; any *gene type counts.
    has_genes = any(kind.lower().endswith("gene") for kind in features)
    if rep.status != "FAIL" and not (has_genes and features.get("CDS")):
        rep.fail(Reason.GFF_MALFORMED, "no gene and CDS features; annotation and resistance codons need both")
    rep.facts = {"features": dict(sorted(features.items()))}
    return rep


CATALOG_COLUMNS = ("marker_id", "gene_id", "gene_name", "chrom", "codon", "codon_start", "strand", "ref_aa", "alt_aa")


def validate_catalog(path: Path) -> FileReport:
    """Resistance catalogue structure: '#key=value' metadata (with status) and the marker columns."""
    rep = FileReport(path=str(path), kind="CATALOGUE")
    if not path.is_file():
        rep.fail(Reason.FILE_MISSING, "catalogue not found")
        return rep
    meta: dict[str, str] = {}
    header: list[str] | None = None
    markers = 0
    for line in path.read_text(errors="replace").splitlines():
        if line.startswith("#") and "=" in line:
            key, value = line[1:].split("=", 1)
            meta[key.strip()] = value.strip()
        elif line.strip() and header is None:
            header = line.split("\t")
        elif line.strip():
            markers += 1
    missing = [c for c in CATALOG_COLUMNS if c not in (header or [])]
    if missing:
        rep.fail(Reason.CATALOG_MALFORMED, "missing columns: " + ", ".join(missing))
    elif not markers:
        rep.fail(Reason.CATALOG_MALFORMED, "no marker rows")
    if "status" not in meta:
        rep.warn(Reason.CATALOG_MALFORMED, "no '#status=' line; DA002 runs only for MAINTAINER_VERIFIED or TEST_FIXTURE")
    rep.facts = {"markers": markers, "status": meta.get("status")}
    return rep


def validate_primers(path: Path, contigs: dict[str, dict[str, Any]]) -> FileReport:
    """fgbio TrimPrimers table: header required, 1-based inclusive coordinates."""
    rep = FileReport(path=str(path), kind="PRIMERS")
    if not path.is_file():
        rep.fail(Reason.FILE_MISSING, "primer table not found")
        return rep
    lines = [ln for ln in path.read_text().splitlines() if ln.strip()]
    if not lines or lines[0].split("\t")[:5] != PRIMER_HEADER:
        rep.fail(Reason.PRIMER_TABLE_MALFORMED, f"first line must be {PRIMER_HEADER}", 1)
        return rep
    for i, line in enumerate(lines[1:], start=2):
        c = line.split("\t")
        try:
            chrom, ls, le, rs, re_ = c[0], int(c[1]), int(c[2]), int(c[3]), int(c[4])
        except (IndexError, ValueError):
            rep.fail(Reason.PRIMER_TABLE_MALFORMED, "non-integer coordinate", i)
            continue
        if chrom not in contigs:
            rep.fail(Reason.REFERENCE_MISMATCH, f"contig {chrom} not in reference", i)
            continue
        if not (1 <= ls <= le < rs <= re_ <= contigs[chrom]["length"]):
            rep.fail(Reason.PRIMER_TABLE_MALFORMED,
                     "need 1 <= left_start <= left_end < right_start <= right_end <= contig length", i)
    rep.facts = {"amplicons": len(lines) - 1}
    return rep


def validate_read_structure(structure: str, read_lengths: tuple[int, int] | None) -> FileReport:
    rep = FileReport(path="assay.read_structure", kind="UMI_READ_STRUCTURE")
    parts = structure.split()
    pattern = re.compile(r"^(\d+|\+)[TBMSC](?:(\d+|\+)[TBMSC])*$")
    if len(parts) != 2 or not all(pattern.match(p) for p in parts):
        rep.fail(Reason.UMI_READ_STRUCTURE, "need two fgbio read structures, e.g. '8M+T 8M+T'")
        return rep
    if not any("M" in p for p in parts):
        rep.fail(Reason.UMI_READ_STRUCTURE, "no molecular-barcode (M) segment")
    if read_lengths:
        for p, length in zip(parts, read_lengths):
            segments = re.findall(r"(\d+)([TBMSC])", p)
            fixed = sum(int(n) for n, _ in segments)
            # Only non-template segments (UMI, barcode, skip) consume bases before the template.
            non_template = sum(int(n) for n, kind in segments if kind != "T")
            if "+" not in p and fixed != length:
                rep.fail(Reason.UMI_READ_STRUCTURE, f"structure {p} totals {fixed}, reads are {length}")
            if "T" not in p or non_template >= length:
                rep.fail(Reason.UMI_READ_STRUCTURE, f"structure {p} leaves no template bases")
    return rep


SEGMENT = re.compile(r"(\d+|\+)([TBMSC])")


def runtime_read_structures(structure: str) -> tuple[list[str], int, bool]:
    """Read structures as passed to fgbio after 3' trimming, the shortest usable read, and trim safety.

    Adapter and poly-G trimming shorten reads at the 3' end, so a fixed-length final template
    segment ('8M143T') no longer matches every read and fgbio stops on the first shorter one.
    Such a segment is relaxed to '+T' (the remaining bases are the template). The second value
    is the smallest read length that still contains every fixed segment (plus one base for '+').
    The third value is False when a read ends in a fixed non-template segment: 3' trimming would
    then remove or shift that segment, so it has to be switched off for this library.
    """
    relaxed: list[str] = []
    minimum = 0
    trim_safe = True
    for part in structure.split():
        segments = SEGMENT.findall(part)
        if segments and segments[-1][1] == "T":
            segments[-1] = ("+", "T")
        elif segments and segments[-1][0] != "+":
            trim_safe = False
        relaxed.append("".join(n + kind for n, kind in segments))
        fixed = sum(int(n) for n, _ in segments if n != "+")
        minimum = max(minimum, fixed + (1 if any(n == "+" for n, _ in segments) else 0))
    return relaxed, minimum, trim_safe


def first_read_length(path: Path) -> int:
    with _open(path) as handle:
        for _, seq, _q in _fastq_records(handle):
            return len(seq)
    return 0


def validate_alignment(path: Path, contigs: dict[str, dict[str, Any]], reference: Path) -> FileReport:
    """Use a temporary sidecar for unindexed CRAM references; never index the original."""
    import tempfile

    import pysam

    if path.suffix.lower() == ".cram" and reference.is_file() and not Path(str(reference) + ".fai").is_file():
        with tempfile.TemporaryDirectory(prefix="pfvarx-cram-validation-") as tmp:
            staged = Path(tmp) / "reference.fasta"
            staged.symlink_to(reference.resolve())
            pysam.faidx(str(staged))
            return _validate_alignment(path, contigs, staged)
    return _validate_alignment(path, contigs, reference)


def _validate_alignment(path: Path, contigs: dict[str, dict[str, Any]], reference: Path) -> FileReport:
    """Validate a user-supplied BAM/CRAM for Layer 2 entry (L1.9)."""
    import pysam

    rep = FileReport(path=str(path), kind="ALIGNMENT")
    if not path.is_file():
        rep.fail(Reason.FILE_MISSING, "alignment not found")
        return rep
    try:
        aln = pysam.AlignmentFile(str(path), reference_filename=str(reference))
    except (ValueError, OSError) as exc:
        rep.fail(Reason.REFERENCE_MALFORMED, f"cannot open alignment: {exc}")
        return rep
    with aln:
        header = aln.header.to_dict()
        hd = header.get("HD", {})
        if hd.get("SO") != "coordinate":
            rep.fail(Reason.ALIGNMENT_UNSORTED, f"@HD SO is {hd.get('SO')!r}, need 'coordinate'")
        sq = {s["SN"]: s for s in header.get("SQ", [])}
        for name, info in contigs.items():
            if name not in sq:
                rep.fail(Reason.REFERENCE_MISMATCH, f"contig {name} absent from header")
            elif sq[name]["LN"] != info["length"]:
                rep.fail(Reason.REFERENCE_MISMATCH, f"contig {name} length differs")
            elif "M5" in sq[name] and sq[name]["M5"] != info["md5"]:
                rep.fail(Reason.REFERENCE_MISMATCH, f"contig {name} M5 checksum differs")
        extra = sorted(set(sq) - set(contigs))
        if extra:
            rep.warn(Reason.REFERENCE_MISMATCH, f"extra header contigs (e.g. host decoy): {extra[:5]}")
        rgs = header.get("RG", [])
        if not rgs or any(not {"ID", "SM", "LB", "PL"} <= set(rg) for rg in rgs):
            rep.fail(Reason.ALIGNMENT_READ_GROUP, "every @RG needs ID, SM, LB and PL")
        samples = sorted({rg.get("SM") for rg in rgs if rg.get("SM")})
        programs = [pg.get("PN", pg.get("ID", "")) for pg in header.get("PG", [])]
        rep.facts = {"samples": samples, "programs": programs, "sort_order": hd.get("SO")}
        try:
            aln.check_index()
        except ValueError:
            rep.fail(Reason.ALIGNMENT_NO_INDEX, "no .bai/.csi/.crai index")
        if not any("MarkDuplicates" in p or "markdup" in p.lower() for p in programs):
            rep.warn(Reason.ALIGNMENT_PROCESSING_STATE,
                     "no duplicate-marking program in @PG; WGS entry will repeat duplicate marking")
    return rep


REQUIRED_SAMPLE_COLUMNS = {"sample_id", "lane", "library_id"}


def validate_samples(path: Path, root: Path, input_type: InputType, design: Design,
                     required_metadata: set[str]) -> tuple[FileReport, list[dict[str, str]]]:
    rep = FileReport(path=str(path), kind="SAMPLE_METADATA")
    if not path.is_file():
        rep.fail(Reason.FILE_MISSING, "samples.tsv not found")
        return rep, []
    try:
        rows = read_tsv(path)
    except ValueError as exc:
        rep.fail(Reason.METADATA_MISSING_FIELD, str(exc))
        return rep, []
    if not rows:
        rep.fail(Reason.METADATA_MISSING_FIELD, "no samples registered")
        return rep, []
    cols = set(rows[0])
    need = set(REQUIRED_SAMPLE_COLUMNS)
    need |= {"fastq_1", "fastq_2"} if input_type == InputType.FASTQ else set()
    need |= {"bam"} if input_type in (InputType.BAM, InputType.CRAM) else set()
    need |= {"vcf"} if input_type == InputType.VCF else set()
    for col in sorted(need - cols):
        rep.fail(Reason.METADATA_MISSING_FIELD, f"column {col} missing")
    seen: set[tuple[str, str]] = set()
    flowcells: dict[str, str] = {}
    for i, row in enumerate(rows, start=2):
        key = (row.get("sample_id", ""), row.get("lane", ""))
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]{0,63}", key[0]):
            rep.fail(Reason.METADATA_MISSING_FIELD, f"invalid sample_id {key[0]!r} (letters, digits, _ and - only)", i)
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]{0,31}", key[1]):
            rep.fail(Reason.METADATA_MISSING_FIELD, f"invalid lane {key[1]!r} (letters, digits, _ and - only)", i)
        if key in seen:
            rep.fail(Reason.METADATA_DUPLICATE_SAMPLE, f"duplicate sample/lane {key}", i)
        seen.add(key)
        for col in need - REQUIRED_SAMPLE_COLUMNS:
            value = row.get(col, "")
            if not value or not (root / value).is_file():
                rep.fail(Reason.FILE_MISSING, f"{col}={value!r} not found", i)
        for col in sorted(required_metadata):
            if not row.get(col):
                rep.warn(Reason.METADATA_MISSING_FIELD, f"{col} empty for {key[0]} (affects analyses needing it)")
        if input_type == InputType.FASTQ and design == Design.WGS:
            # Duplicates are marked once on the merged sample with one optical-duplicate
            # distance, so every lane of a sample must come from the same flow-cell type.
            first = flowcells.setdefault(key[0], row.get("flowcell_type", ""))
            if row.get("flowcell_type", "") != first:
                rep.fail(Reason.METADATA_MISSING_FIELD,
                         f"sample {key[0]} mixes flowcell_type {first!r} and {row.get('flowcell_type', '')!r}; "
                         "register lanes from different flow-cell types as separate samples", i)
    rep.facts = {"samples": len({r['sample_id'] for r in rows}), "lanes": len(rows),
                 "design": design.value}
    return rep, rows
