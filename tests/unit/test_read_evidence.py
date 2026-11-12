"""Missing optional BAM fields must not crash or supply allele evidence."""

from pathlib import Path

import pysam
import pytest

from pfvarx.caller.call import Candidate, Params, linkage_blocks, rescore_indels
from pfvarx.cohort.evidence import CohortSite, allele_counts
from pfvarx.configuration.loader import load_parameters
from pfvarx.downstream.resistance import codon_reads


@pytest.fixture()
def reads(tmp_path: Path) -> Path:
    path = tmp_path / "reads.bam"
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 100}]}
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        sequence_with_codon = "A" * 10 + "ATG" + "A" * 17
        for name, sequence, quality in [("valid", sequence_with_codon, "I" * 30),
                                         ("no_quality", sequence_with_codon, None),
                                         ("no_sequence", None, None)]:
            read = pysam.AlignedSegment(out.header)
            read.query_name = name
            read.query_sequence = sequence
            read.flag = 0
            read.reference_id = 0
            read.reference_start = 10
            read.mapping_quality = 60
            read.cigarstring = "30M"
            read.query_qualities = pysam.qualitystring_to_array(quality) if quality else None
            out.write(read)
    pysam.index(str(path))
    return path


def test_cohort_excludes_reads_without_sequence_or_quality(reads: Path) -> None:
    site: CohortSite = {"chrom": "chr1", "pos": 11, "ref": "A", "alt": "C", "type": "SNV",
                        "n_pass": 0, "statuses": {}}
    with pysam.AlignmentFile(str(reads)) as bam:
        assert allele_counts(bam, site, min_mapq=20, min_bq=13) == (1, 0, 1)


def test_codon_evidence_excludes_reads_without_sequence_or_quality(reads: Path) -> None:
    with pysam.AlignmentFile(str(reads)) as bam:
        assert codon_reads(bam, "chr1", [20, 21, 22], "+", min_mapq=20, min_bq=13) == ["ATG"]


@pytest.fixture()
def params() -> Params:
    root = Path(__file__).resolve().parents[2]
    return Params.from_frozen({"parameters": load_parameters(root), "config": {"random_seed": 20261106}})


def test_indel_rescoring_excludes_reads_without_quality(reads: Path, params: Params) -> None:
    candidate = Candidate("chr1", 20, "A", "AT", "INS", y=1, n=1)
    with pysam.AlignmentFile(str(reads)) as bam:
        rescore_indels(bam, {"chr1": "A" * 100}, [candidate], params)
    assert candidate.pairhmm_reads == 1


def test_linkage_handles_missing_sequence(reads: Path, params: Params) -> None:
    candidates = [Candidate("chr1", 20, "C", "A", "SNV", y=2, n=2, status="PASS"),
                  Candidate("chr1", 21, "C", "T", "SNV", y=2, n=2, status="PASS")]
    with pysam.AlignmentFile(str(reads)) as bam:
        linkage_blocks(bam, candidates, params)
    assert candidates[0].lhid == candidates[1].lhid
    assert candidates[0].lhid != "."
