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
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="ascii", errors="replace")
    return path.open("rt", encoding="ascii", errors="replace")


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
                        qmin, qmax = min(qmin, min(codes)), max(qmax, max(codes))
    except EOFError:
        rep.fail(Reason.FASTQ_TRUNCATED, "file ends inside a record", n + 1)
    except (ValueError, gzip.BadGzipFile, OSError) as exc:
        rep.fail(Reason.FASTQ_MALFORMED_RECORD, str(exc), n + 1)
    if n and (qmin < PHRED33_MIN or qmax > PHRED33_MAX):
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
    for i, line in enumerate(path.read_text().splitlines(), start=1):
        if not line.strip() or line.startswith(("#", "track", "browser")):
            continue
        cells = line.split("\t")
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
    if n == 0 and rep.status == "PASS":
        rep.fail(Reason.BED_MALFORMED, "no intervals")
    rep.facts = {"intervals": n, "bases": total}
    return rep


PRIMER_HEADER = ["chrom", "left_start", "left_end", "right_start", "right_end"]


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
            fixed = sum(int(n) for n, _ in re.findall(r"(\d+)([TBMSC])", p))
            if "+" not in p and fixed != length:
                rep.fail(Reason.UMI_READ_STRUCTURE, f"structure {p} totals {fixed}, reads are {length}")
            if fixed >= length:
                rep.fail(Reason.UMI_READ_STRUCTURE, f"structure {p} leaves no template bases")
    return rep


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
    rep.facts = {"samples": len({r['sample_id'] for r in rows}), "lanes": len(rows),
                 "design": design.value}
    return rep, rows
