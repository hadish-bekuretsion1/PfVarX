"""DA008 copy-number markers: depth ratios, calls, junction read pairs and the catalogue."""

from __future__ import annotations

import random
from pathlib import Path

import pysam
import pytest

from pfvarx.downstream import cnv

LENGTH = 20000
AMP = (5000, 7000)    # 0-based gene with twice the coverage (a tandem duplication)
DEL = (12000, 13000)  # gene with no reads (a deletion)
PARAMS = {"da008_amplification_ratio": 1.5, "da008_amplification_lower": 1.25, "da008_deletion_ratio": 0.1,
          "da008_present_ratio": 0.3, "da008_flank_min_ratio": 0.3, "da008_min_core_depth": 20,
          "da008_flank_bp": 2000, "da008_window_bp": 200, "da008_bootstrap": 200}


def _pair(header: pysam.AlignmentHeader, name: str, left: int, right: int, left_rev: bool, right_rev: bool,
          seq: str) -> list[pysam.AlignedSegment]:
    reads = []
    for pos, mate, rev, mate_rev, first in ((left, right, left_rev, right_rev, True), (right, left, right_rev, left_rev, False)):
        a = pysam.AlignedSegment(header)
        a.query_name, a.reference_id, a.reference_start = name, 0, pos
        a.query_sequence, a.cigarstring, a.mapping_quality = seq, "100M", 60
        a.query_qualities = pysam.qualitystring_to_array("I" * 100)
        a.is_paired, a.is_read1, a.is_read2 = True, first, not first
        a.is_reverse, a.mate_is_reverse = rev, mate_rev
        a.next_reference_id, a.next_reference_start = 0, mate
        a.template_length = (right + 100 - left) * (1 if first else -1)
        reads.append(a)
    return reads


@pytest.fixture(scope="module")
def fasta(tmp_path_factory: pytest.TempPathFactory) -> Path:
    path = tmp_path_factory.mktemp("ref") / "ref.fasta"
    rng = random.Random(2)
    path.write_text(">ctg\n" + "".join(rng.choices("ACGT", weights=[40, 10, 10, 40], k=LENGTH)) + "\n")
    pysam.faidx(str(path))
    return path


@pytest.fixture(scope="module")
def bam(tmp_path_factory: pytest.TempPathFactory) -> Path:
    d = tmp_path_factory.mktemp("cnv")
    header = pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6", "SO": "unsorted"},
                                              "SQ": [{"SN": "ctg", "LN": LENGTH}]})
    seq = "".join(random.Random(1).choices("ACGT", k=100))
    reads: list[pysam.AlignedSegment] = []
    n = 0
    for start in range(0, LENGTH - 400, 10):  # fragments of 400 bp every 10 bp: about 20x per base
        if DEL[0] - 400 < start < DEL[1]:
            continue  # nothing overlaps the deleted gene
        copies = 2 if AMP[0] <= start and start + 400 <= AMP[1] else 1
        for _ in range(copies):
            n += 1
            reads += _pair(header, f"f{n}", start, start + 300, False, True, seq)
    for k in range(6):  # duplication junction: outward-facing pairs across the gene ends
        reads += _pair(header, f"j{k}", AMP[0] + 20 + 10 * k, AMP[1] - 150 - 10 * k, True, False, seq)
    for k in range(4):  # deletion: pairs spanning the deleted gene
        reads += _pair(header, f"s{k}", DEL[0] - 150 - 10 * k, DEL[1] + 50 + 10 * k, False, True, seq)
    raw = d / "raw.bam"
    with pysam.AlignmentFile(str(raw), "wb", header=header) as out:
        for r in reads:
            out.write(r)
    path = d / "s.bam"
    pysam.sort("-o", str(path), str(raw))
    pysam.index(str(path))
    return path


def test_depth_ratios_and_calls(bam: Path) -> None:
    with pysam.AlignmentFile(str(bam)) as b:
        core = cnv.statistics.median(cnv.window_medians(cnv.depth_profile(b, "ctg", 1000, 4000, 20), 200))
        assert 15 <= core <= 25
        amp = cnv.measure(b, ("ctg", *AMP), core, PARAMS, 20, seed=1)
        gone = cnv.measure(b, ("ctg", *DEL), core, PARAMS, 20, seed=2)
        normal = cnv.measure(b, ("ctg", 15000, 16500), core, PARAMS, 20, seed=3)
    assert 1.8 <= amp["depth_ratio"] <= 2.2 and amp["ratio_lo"] >= 1.25
    assert amp["outward_pairs"] >= 6
    assert cnv.call_locus("AMPLIFICATION", amp["depth_ratio"], amp["ratio_lo"], amp["flank_ratio"], PARAMS) == ("AMPLIFIED", "")
    assert gone["depth_ratio"] == 0 and gone["flank_ratio"] >= 0.3 and gone["spanning_pairs"] >= 4
    assert cnv.call_locus("DELETION", gone["depth_ratio"], gone["ratio_lo"], gone["flank_ratio"], PARAMS) == ("DELETED", "")
    assert cnv.call_locus("AMPLIFICATION", normal["depth_ratio"], normal["ratio_lo"], normal["flank_ratio"], PARAMS)[0] == "SINGLE_COPY"
    assert cnv.call_locus("DELETION", normal["depth_ratio"], normal["ratio_lo"], normal["flank_ratio"], PARAMS)[0] == "PRESENT"


def test_unjudgeable_loci_are_unknown_never_present() -> None:
    assert cnv.call_locus("DELETION", 0.0, 0.0, 0.05, PARAMS) == ("UNKNOWN", "FLANKS_NOT_COVERED")
    assert cnv.call_locus("DELETION", 0.2, 0.1, 1.0, PARAMS) == ("UNKNOWN", "AMBIGUOUS_DEPTH")
    assert cnv.call_locus("AMPLIFICATION", 0.1, 0.05, 1.0, PARAMS) == ("UNKNOWN", "GENE_DEPTH_LOW")
    assert cnv.call_locus("AMPLIFICATION", None, None, None, PARAMS) == ("UNKNOWN", "NO_DEPTH")
    # A minor clone carrying a duplication: ratio between one and two copies, interval not clear of 1.25.
    assert cnv.call_locus("AMPLIFICATION", 1.3, 1.1, 1.0, PARAMS)[0] == "SINGLE_COPY"


def test_gene_location_from_gene_feature_or_cds(tmp_path: Path) -> None:
    gff = tmp_path / "g.gff"
    gff.write_text("##gff-version 3\n"
                   "Pf3D7_05_v3\tVEuPathDB\tprotein_coding_gene\t955955\t963095\t.\t+\t.\tID=PF3D7_0523000;Name=MDR1\n"
                   "ctg\tsynthetic\tCDS\t101\t400\t.\t+\t0\tID=PFX_A.1-CDS;Parent=PFX_A.1\n")
    assert cnv.gene_location(gff, "PF3D7_0523000") == ("Pf3D7_05_v3", 955954, 963095)
    assert cnv.gene_location(gff, "PFX_A") == ("ctg", 100, 400)
    assert cnv.gene_location(gff, "PF3D7_0831800") is None


def test_prevalence_counts_only_judged_infections() -> None:
    calls = [{"sample_id": s, "locus_id": "hrp2", "state": st, "population": "A"}
             for s, st in (("a", "DELETED"), ("b", "PRESENT"), ("c", "UNKNOWN"))]
    rows = cnv.prevalence(calls, [{"locus_id": "hrp2", "gene_name": "pfhrp2", "event": "DELETION"}])
    overall = next(r for r in rows if r["group"] == "ALL")
    assert (overall["judged"], overall["unknown"], overall["with_event"], overall["prevalence"]) == (2, 1, 1, "0.5000")


def test_shipped_catalogue_is_draft_and_complete() -> None:
    repo = Path(__file__).resolve().parents[2]
    from pfvarx.downstream.resistance import load_catalog

    meta, loci = load_catalog(repo / cnv.CATALOG)
    assert meta["status"] == "DRAFT_REQUIRES_VERIFICATION"
    assert {lc["gene_id"] for lc in loci} == {"PF3D7_0523000", "PF3D7_1408000", "PF3D7_1408100", "PF3D7_1224000",
                                              "PF3D7_0831800", "PF3D7_1372200"}
    assert all(set(cnv.CATALOG_COLUMNS) <= set(lc) and lc["event"] in cnv.EVENTS for lc in loci)


def test_unverified_catalogues_warn_before_the_run(tmp_path: Path) -> None:
    from pfvarx.execution.runner import catalogue_warnings

    repo = Path(__file__).resolve().parents[2]
    (tmp_path / "config/registry").mkdir(parents=True)
    for name in ("cnv_catalog_pfvarx_v1.tsv", "resistance_catalog_pfvarx_v1.tsv"):
        (tmp_path / "config/registry" / name).write_text((repo / "config/registry" / name).read_text())
    warnings = catalogue_warnings(tmp_path, ["DA002", "DA008"], None)
    assert [w["analysis"] for w in warnings] == ["DA002", "DA008"]
    assert all(w["reason"] == "CATALOG_NOT_VERIFIED" and "EXTEND" in w["detail"] for w in warnings)
    assert catalogue_warnings(tmp_path, ["DA003"], None) == []


def test_gc_corrected_ratios_match_the_simulated_copy_numbers(bam: Path, fasta: Path) -> None:
    with pysam.AlignmentFile(str(bam)) as b, pysam.FastaFile(str(fasta)) as fa:
        model = cnv.GcModel(b, fa, [("ctg", 1000, 4000), ("ctg", 14500, 19000)], 200, 300, 20, seed=1)
        assert len(model.points) == 300 and 15 <= (model.expected(0.2) or 0) <= 25
        amp = cnv.measure(b, ("ctg", *AMP), 20.0, PARAMS, 20, 1, fa, model)
        gone = cnv.measure(b, ("ctg", *DEL), 20.0, PARAMS, 20, 2, fa, model)
    assert 1.8 <= amp["depth_ratio"] <= 2.2 and amp["ratio_lo"] >= 1.25
    assert gone["depth_ratio"] == 0 and gone["flank_ratio"] >= 0.3
    assert cnv.gc_fraction("ACGTNN") == 0.5 and cnv.gc_fraction("NNN") == 0.0
