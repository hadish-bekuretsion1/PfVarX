"""Inputs are recognised by content and structure; file names and extensions decide nothing."""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from pfvarx.configuration.loader import resolve_input_defaults
from pfvarx.inputs import validate as v
from pfvarx.inputs.content import content_problem, detect, find_in_folder

CONTIGS = {"chr1": {"length": 1000}}


def _write(path: Path, text: str, compress: bool = False) -> Path:
    data = text.encode()
    path.write_bytes(gzip.compress(data) if compress else data)
    return path


@pytest.mark.parametrize(("text", "compress", "kind"), [
    (">chr1\nACGT\n", False, "FASTA"),
    ("@r1\nACGT\n+\nIIII\n", True, "FASTQ"),
    ("chr1\t0\t100\n", False, "BED"),
    ("##gff-version 3\nchr1\tsrc\tgene\t1\t90\t.\t+\t.\tID=g\n", False, "GFF"),
    ("chr1\tsrc\tCDS\t1\t90\t.\t-\t0\tID=c\n", False, "GFF"),
    ("##fileformat=VCFv4.2\n#CHROM\tPOS\n", False, "VCF"),
    ("sample_id\tlane\tlibrary_id\nS1\tL1\tLIB\n", False, "SAMPLE_SHEET"),
    ("chrom\tleft_start\tleft_end\tright_start\tright_end\nchr1\t1\t20\t80\t100\n", False, "PRIMERS"),
    ("#status=TEST_FIXTURE\nmarker_id\tgene_id\nm1\tg1\n", False, "CATALOGUE"),
    ('{"a": 1}', False, "JSON"),
    ("just some notes", False, "UNKNOWN"),
])
def test_content_kind_ignores_the_name(tmp_path: Path, text: str, compress: bool, kind: str) -> None:
    assert detect(_write(tmp_path / "misleading.name.bam", text, compress)) == kind


def test_binary_alignments_and_indexes(tmp_path: Path) -> None:
    assert detect(_write(tmp_path / "a.txt", "CRAM\x03\x00rest")) == "CRAM"
    (tmp_path / "b.txt").write_bytes(gzip.compress(b"BAM\x01" + b"\x00" * 40))
    assert detect(tmp_path / "b.txt") == "BAM"
    (tmp_path / "c.txt").write_bytes(b"BAI\x01" + b"\x00" * 8)
    assert detect(tmp_path / "c.txt") == "ALIGNMENT_INDEX"
    (tmp_path / "empty.fasta").write_bytes(b"")
    assert detect(tmp_path / "empty.fasta") == "EMPTY"


def test_wrong_content_and_compression_rules(tmp_path: Path) -> None:
    gff = _write(tmp_path / "core_genome.bed", "##gff-version 3\nchr1\ts\tgene\t1\t9\t.\t+\t.\tID=g\n")
    assert content_problem(gff, ("BED",)) == "core_genome.bed contains GFF3 annotation; this input needs BED intervals"
    reads = "@r1\nACGT\n+\nIIII\n"
    assert "does not end in .gz" in content_problem(_write(tmp_path / "S_R1.fastq", reads, True), ("FASTQ",))
    assert "is not compressed" in content_problem(_write(tmp_path / "S_R2.fastq.gz", reads), ("FASTQ",))
    assert content_problem(_write(tmp_path / "reads.gz", reads, True), ("FASTQ",)) is None
    assert "Decompress it" in content_problem(_write(tmp_path / "ref.fa.gz", ">c\nACGT\n", True), ("FASTA",))


def test_fastq_pair_validation_checks_content(tmp_path: Path) -> None:
    r1 = _write(tmp_path / "a.fastq.gz", "@r1/1\nACGT\n+\nIIII\n", True)
    r2 = _write(tmp_path / "b.fastq.gz", "chr1\t0\t10\n", True)
    rep = v.validate_fastq_pair(r1, r2, 100)
    assert rep.status == "FAIL" and rep.findings[0].reason == "CONTENT_MISMATCH"
    r2 = _write(tmp_path / "b.reads", "@r1/2\nACGT\n+\nIIII\n")  # plain FASTQ, any name
    assert v.validate_fastq_pair(r1, r2, 100).status == "PASS"


def test_gff_structure(tmp_path: Path) -> None:
    good = "##gff-version 3\nchr1\ts\tgene\t1\t90\t.\t+\t.\tID=g\nchr1\ts\tCDS\t1\t90\t.\t+\t0\tParent=g\n"
    assert v.validate_gff(_write(tmp_path / "a", good), CONTIGS).status == "PASS"
    plasmodb = good.replace("\tgene\t", "\tprotein_coding_gene\t")  # PlasmoDB gene type
    assert v.validate_gff(_write(tmp_path / "p", plasmodb), CONTIGS).status == "PASS"
    outside = good.replace("\t90\t.\t+\t0", "\t5000\t.\t+\t0")
    rep = v.validate_gff(_write(tmp_path / "b", outside), CONTIGS)
    assert rep.findings[0].reason == "GFF_MALFORMED" and "outside the contig" in rep.findings[0].detail
    other = good.replace("chr1", "chrX")
    assert v.validate_gff(_write(tmp_path / "c", other), CONTIGS).findings[0].reason == "REFERENCE_MISMATCH"
    genes_only = "##gff-version 3\nchr1\ts\tgene\t1\t90\t.\t+\t.\tID=g\n"
    assert "no gene and CDS" in v.validate_gff(_write(tmp_path / "d", genes_only), CONTIGS).findings[0].detail


def test_catalogue_structure(tmp_path: Path) -> None:
    header = "marker_id\tgene_id\tgene_name\tchrom\tcodon\tcodon_start\tstrand\tref_aa\talt_aa\n"
    good = "#status=TEST_FIXTURE\n" + header + "m1\tg1\tk13\tchr1\t1\t1\t+\tCys\tTyr\n"
    assert v.validate_catalog(_write(tmp_path / "a", good)).status == "PASS"
    rep = v.validate_catalog(_write(tmp_path / "b", "#status=X\nmarker_id\tgene_id\nm1\tg1\n"))
    assert rep.status == "FAIL" and "missing columns" in rep.findings[0].detail
    assert v.validate_catalog(_write(tmp_path / "c", header + "m1\tg1\tk13\tchr1\t1\t1\t+\tC\tY\n")).status == "WARNING"


def test_folder_setting_picks_the_file_by_content(tmp_path: Path) -> None:
    ref = tmp_path / "inputs/reference"
    ref.mkdir(parents=True)
    _write(ref / "genome_v68.txt", ">chr1\nACGT\n")
    _write(ref / "regions.list", "chr1\t0\t3\n")
    _write(ref / "annotation.data", "##gff-version 3\nchr1\ts\tgene\t1\t3\t.\t+\t.\tID=g\n")
    project = {"reference": {"fasta": "inputs/reference/", "core_mask_bed": "inputs/reference/",
                             "gff": "inputs/reference"}, "assay": {}, "samples_tsv": "inputs/metadata/"}
    resolved = resolve_input_defaults(project, tmp_path)
    assert resolved["reference"] == {"fasta": "inputs/reference/genome_v68.txt",
                                     "core_mask_bed": "inputs/reference/regions.list",
                                     "gff": "inputs/reference/annotation.data"}
    assert resolved["samples_tsv"] == "inputs/metadata/samples.tsv"  # nothing found: the usual name is reported
    _write(ref / "second.bed", "chr1\t0\t2\n")  # two BEDs: ambiguous, so validation must ask
    assert resolve_input_defaults(project, tmp_path)["reference"]["core_mask_bed"] == "inputs/reference/core_genome.bed"
    assert [p.name for p in find_in_folder(ref, ("BED",))] == ["regions.list", "second.bed"]
