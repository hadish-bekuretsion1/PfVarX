import gzip
from pathlib import Path

from pfvarx.inputs import validate as v


def _fq(path: Path, names: list[str], qual: str = "I", truncate: bool = False) -> None:
    with gzip.open(path, "wt") as h:
        for n in names:
            h.write(f"@{n}\nACGT\n+\n{qual * 4}\n")
        if truncate:
            h.write("@cut\nACGT\n+\n")


def test_valid_pair(tmp_path: Path) -> None:
    _fq(tmp_path / "a.fq.gz", ["r1/1", "r2/1"])
    _fq(tmp_path / "b.fq.gz", ["r1/2", "r2/2"])
    rep = v.validate_fastq_pair(tmp_path / "a.fq.gz", tmp_path / "b.fq.gz", 100)
    assert rep.status == "PASS" and rep.facts["read_pairs"] == 2


def test_truncated_mismatched_and_encoding(tmp_path: Path) -> None:
    _fq(tmp_path / "a.fq.gz", ["r1", "r2"], truncate=True)
    _fq(tmp_path / "b.fq.gz", ["r1", "r2", "r3"])
    rep = v.validate_fastq_pair(tmp_path / "a.fq.gz", tmp_path / "b.fq.gz", 100)
    assert rep.status == "FAIL" and rep.findings[0].reason == "FASTQ_TRUNCATED"
    _fq(tmp_path / "c.fq.gz", ["x1", "x2"])
    _fq(tmp_path / "d.fq.gz", ["x1", "y2"])
    rep = v.validate_fastq_pair(tmp_path / "c.fq.gz", tmp_path / "d.fq.gz", 100)
    assert rep.findings[0].reason == "FASTQ_MATE_NAME_MISMATCH" and rep.findings[0].record == 2
    _fq(tmp_path / "e.fq.gz", ["z"], qual=" ")
    _fq(tmp_path / "f.fq.gz", ["z"], qual=" ")
    rep = v.validate_fastq_pair(tmp_path / "e.fq.gz", tmp_path / "f.fq.gz", 100)
    assert any(f.reason == "FASTQ_QUALITY_ENCODING" for f in rep.findings)


def test_fasta_bed_primers(tmp_path: Path) -> None:
    (tmp_path / "r.fa").write_text(">c1\nACGTACGTAC\n>c2\nAAAA\n")
    fa = v.validate_fasta(tmp_path / "r.fa")
    assert fa.status == "PASS" and fa.facts["contigs"]["c1"]["length"] == 10
    contigs = fa.facts["contigs"]
    (tmp_path / "ok.bed").write_text("c1\t0\t10\n")
    (tmp_path / "bad.bed").write_text("c1\t5\t11\nc9\t0\t1\n")
    assert v.validate_bed(tmp_path / "ok.bed", contigs).status == "PASS"
    bad = v.validate_bed(tmp_path / "bad.bed", contigs)
    assert {f.reason for f in bad.findings} == {"BED_MALFORMED", "REFERENCE_MISMATCH"}
    (tmp_path / "p.tsv").write_text("chrom\tleft_start\tleft_end\tright_start\tright_end\nc1\t1\t3\t8\t10\n")
    assert v.validate_primers(tmp_path / "p.tsv", contigs).status == "PASS"
    (tmp_path / "p0.tsv").write_text("chrom\tleft_start\tleft_end\tright_start\tright_end\nc1\t0\t3\t8\t10\n")
    assert v.validate_primers(tmp_path / "p0.tsv", contigs).status == "FAIL"  # 0 is not 1-based


def test_read_structure() -> None:
    assert v.validate_read_structure("8M+T 8M+T", (108, 108)).status == "PASS"
    assert v.validate_read_structure("100T 100T", (100, 100)).status == "FAIL"  # no UMI segment
    assert v.validate_read_structure("8M92T 8M92T", (108, 108)).status == "FAIL"  # length mismatch


def test_fixed_length_read_structure_is_accepted_and_relaxed_for_trimmed_reads() -> None:
    # A fully fixed structure that matches the read length is valid (it used to be rejected).
    assert v.validate_read_structure("8M100T 8M100T", (108, 108)).status == "PASS"
    assert v.validate_read_structure("8M143T 8M143T", (151, 151)).status == "PASS"
    assert v.validate_read_structure("108M 108M", (108, 108)).status == "FAIL"  # no template bases
    # fgbio receives '+T' because adapter trimming shortens reads at the 3' end.
    assert v.runtime_read_structures("8M143T 8M143T") == (["8M+T", "8M+T"], 9, True)
    assert v.runtime_read_structures("3S8M+T +T") == (["3S8M+T", "+T"], 12, True)
    # A fixed non-template segment at the 3' end forbids trimming and keeps the exact length.
    assert v.runtime_read_structures("8M135T8M 8M135T8M") == (["8M135T8M", "8M135T8M"], 151, False)


def test_wgs_sample_must_not_mix_flowcell_types(tmp_path: Path) -> None:
    """Duplicates are marked on the merged sample with one optical distance (PFX-12)."""
    from pfvarx.constants import Design, InputType

    for name in ("a1", "a2", "b1", "b2"):
        (tmp_path / f"{name}.fq.gz").touch()
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(
        "sample_id\tlane\tfastq_1\tfastq_2\tlibrary_id\tflowcell_type\n"
        "S1\tL1\ta1.fq.gz\ta2.fq.gz\tlib\tpatterned\n"
        "S1\tL2\ta1.fq.gz\ta2.fq.gz\tlib\tunpatterned\n"
        "S2\tL1\tb1.fq.gz\tb2.fq.gz\tlib\tpatterned\n"
        "S2\tL2\tb1.fq.gz\tb2.fq.gz\tlib\tpatterned\n")
    rep, _ = v.validate_samples(sheet, tmp_path, InputType.FASTQ, Design.WGS, set())
    failures = [f for f in rep.findings if f.level == "FAIL"]
    assert len(failures) == 1 and failures[0].record == 3 and "mixes flowcell_type" in failures[0].detail
    targeted, _ = v.validate_samples(sheet, tmp_path, InputType.FASTQ, Design.TARGETED_NO_UMI, set())
    assert not [f for f in targeted.findings if f.level == "FAIL"]
