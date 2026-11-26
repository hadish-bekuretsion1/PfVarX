"""Layer 3 microhaplotypes: fragment haplotypes, floors, locus resolution and the locus basis."""

from __future__ import annotations

import gzip
from pathlib import Path

import pysam
import pytest

from pfvarx.cohort import microhap
from pfvarx.constants import Design
from pfvarx.downstream import loci as loci_module

LENGTH = 4000
POSITIONS = [999, 1049, 1099]  # 0-based variable positions inside one locus
PARAMS = {"microhap_min_depth": 10, "microhap_min_haplotype_reads": 2,
          "microhap_min_haplotype_fraction": 0.05, "microhap_max_variable_positions": 32,
          "microhap_min_loci_for_basis": 10}


def _read(header: pysam.AlignmentHeader, name: str, start: int, seq: str, quals: str,
          first: bool = True, mate: int | None = None) -> pysam.AlignedSegment:
    a = pysam.AlignedSegment(header)
    a.query_name, a.reference_id, a.reference_start = name, 0, start
    a.query_sequence, a.cigarstring, a.mapping_quality = seq, f"{len(seq)}M", 60
    a.query_qualities = pysam.qualitystring_to_array(quals)
    if mate is not None:
        a.is_paired, a.is_read1, a.is_read2 = True, first, not first
        a.next_reference_id, a.next_reference_start = 0, mate
    return a


def _bam(tmp_path: Path, reads: list[tuple[str, int, str, str, int | None, bool]]) -> Path:
    header = pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6", "SO": "unsorted"},
                                              "SQ": [{"SN": "ctg", "LN": LENGTH}]})
    raw, out = tmp_path / "raw.bam", tmp_path / "s.bam"
    with pysam.AlignmentFile(str(raw), "wb", header=header) as handle:
        for name, start, seq, quals, mate, first in reads:
            handle.write(_read(header, name, start, seq, quals, first, mate))
    pysam.sort("-o", str(out), str(raw))
    pysam.index(str(out))
    return out


def _spanning(name: str, bases: tuple[str, str, str], start: int = 950, quals: str | None = None
              ) -> tuple[str, int, str, str, int | None, bool]:
    """One read covering all three variable positions, carrying *bases* at them."""
    seq = list("A" * 200)
    for position, base in zip(POSITIONS, bases):
        seq[position - start] = base
    text = "".join(seq)
    return (name, start, text, quals or "I" * len(text), None, True)


def test_haplotypes_are_counted_per_fragment(tmp_path: Path) -> None:
    reads = [_spanning(f"h1_{i}", ("A", "C", "G")) for i in range(7)]
    reads += [_spanning(f"h2_{i}", ("A", "T", "G")) for i in range(3)]
    bam = _bam(tmp_path, reads)
    with pysam.AlignmentFile(str(bam)) as handle:
        counts, incomplete = microhap.fragment_haplotypes(handle, "ctg", POSITIONS, 20, 13)
    assert dict(counts) == {"ACG": 7, "ATG": 3} and incomplete == 0
    state, reason, kept = microhap.call_sample(counts, PARAMS)
    assert (state, reason) == ("CALLED", "")
    assert [(h, n) for h, n, _ in kept] == [("ACG", 7), ("ATG", 3)]
    assert kept[0][2] == pytest.approx(0.7)


def test_overlapping_mates_count_once_and_disagreement_is_unresolved(tmp_path: Path) -> None:
    # Both mates cover the locus (as amplicon pairs do): one fragment, not two.
    agree = [_spanning("pair_ok", ("A", "C", "G"), start=950)[0:4] + (1000, True),
             _spanning("pair_ok", ("A", "C", "G"), start=960)[0:4] + (950, False)]
    conflict = [_spanning("pair_bad", ("A", "C", "G"), start=950)[0:4] + (1000, True),
                _spanning("pair_bad", ("A", "T", "G"), start=960)[0:4] + (950, False)]
    bam = _bam(tmp_path, [(*a[:4], a[4], a[5]) for a in agree + conflict])
    with pysam.AlignmentFile(str(bam)) as handle:
        counts, incomplete = microhap.fragment_haplotypes(handle, "ctg", POSITIONS, 20, 13)
    assert dict(counts) == {"ACG": 1}, "the agreeing pair must count once"
    assert incomplete == 1, "mates that disagree leave the position unresolved"


def test_low_quality_base_makes_a_fragment_incomplete(tmp_path: Path) -> None:
    good = [_spanning(f"g{i}", ("A", "C", "G")) for i in range(10)]
    poor = _spanning("poor", ("A", "C", "G"))
    quals = list(poor[3])
    quals[POSITIONS[1] - 950] = "#"  # Phred 2 at one variable position
    bam = _bam(tmp_path, [*good, (poor[0], poor[1], poor[2], "".join(quals), None, True)])
    with pysam.AlignmentFile(str(bam)) as handle:
        counts, incomplete = microhap.fragment_haplotypes(handle, "ctg", POSITIONS, 20, 13)
    assert dict(counts) == {"ACG": 10} and incomplete == 1


def test_floors_separate_haplotypes_from_sequencing_error() -> None:
    from collections import Counter

    # 1 fragment in 100 is below both the read and the fraction floor.
    state, reason, kept = microhap.call_sample(Counter({"ACG": 99, "ATG": 1}), PARAMS)
    assert state == "CALLED" and [h for h, _, _ in kept] == ["ACG"]
    # 4 in 100 clears the read floor but not the 5% fraction floor.
    _, _, kept = microhap.call_sample(Counter({"ACG": 96, "ATG": 4}), PARAMS)
    assert [h for h, _, _ in kept] == ["ACG"]
    # 6 in 100 clears both.
    _, _, kept = microhap.call_sample(Counter({"ACG": 94, "ATG": 6}), PARAMS)
    assert [h for h, _, _ in kept] == ["ACG", "ATG"]


def test_too_few_fragments_is_uncallable_not_no_haplotype() -> None:
    from collections import Counter

    state, reason, kept = microhap.call_sample(Counter({"ACG": 9}), PARAMS)
    assert (state, reason, kept) == ("UNCALLABLE", "LOW_FRAGMENT_DEPTH", [])
    state, reason, _ = microhap.call_sample(Counter(), PARAMS)
    assert (state, reason) == ("UNCALLABLE", "LOW_FRAGMENT_DEPTH")


def test_locus_resolution_from_coordinates_and_its_refusals() -> None:
    lengths = {"ctg": LENGTH}
    place, reason = microhap.resolve_locus({"locus_id": "m1", "chrom": "ctg", "start": "1000", "end": "1200"},
                                           None, lengths)
    assert place == ("ctg", 999, 1200) and reason == ""  # 1-based inclusive -> 0-based half-open
    assert microhap.resolve_locus({"chrom": "other", "start": "1", "end": "9"}, None, lengths)[1] == \
        "CONTIG_NOT_IN_REFERENCE"
    assert microhap.resolve_locus({"chrom": "ctg", "start": "1", "end": "99999"}, None, lengths)[1] == \
        "COORDINATES_OUT_OF_RANGE"
    assert microhap.resolve_locus({"chrom": "ctg", "start": "x", "end": "9"}, None, lengths)[1] == \
        "COORDINATES_MALFORMED"
    assert microhap.resolve_locus({"locus_id": "m2"}, None, lengths)[1] == "NO_COORDINATES_OR_GENE"
    assert microhap.resolve_locus({"gene_id": "PF3D7_0100100"}, None, lengths)[1] == "GENE_LOCUS_NEEDS_GFF"


def test_a_locus_must_lie_in_what_the_run_sequenced() -> None:
    panel = {"ctg": [(900, 1200)]}
    assert microhap.territory_status(("ctg", 999, 1100), panel, Design.TARGETED_UMI) == (True, "")
    assert microhap.territory_status(("ctg", 500, 1100), panel, Design.TARGETED_UMI) == \
        (False, "OUTSIDE_PANEL_TARGETS")
    assert microhap.territory_status(("ctg", 500, 1100), {"ctg": [(0, 400)]}, Design.WGS) == \
        (False, "OUTSIDE_CORE_MASK")


def test_variable_positions_are_the_cohorts_pass_snvs_only() -> None:
    sites = [{"chrom": "ctg", "pos": "1000", "ref": "A", "alt": "C", "n_pass": "2"},
             {"chrom": "ctg", "pos": "1050", "ref": "A", "alt": "T", "n_pass": "1"},
             {"chrom": "ctg", "pos": "1060", "ref": "AT", "alt": "A", "n_pass": "3"},   # indel: excluded
             {"chrom": "ctg", "pos": "1070", "ref": "A", "alt": "G", "n_pass": "0"},    # no PASS sample
             {"chrom": "ctg", "pos": "3000", "ref": "A", "alt": "G", "n_pass": "2"},    # outside the locus
             {"chrom": "other", "pos": "1000", "ref": "A", "alt": "G", "n_pass": "2"}]  # other contig
    positions, capped = microhap.variable_positions(sites, ("ctg", 900, 1200), 32)
    assert positions == [999, 1049] and capped == ""
    few, capped = microhap.variable_positions(sites, ("ctg", 900, 1200), 1)
    assert few == [999] and capped == "VARIABLE_POSITIONS_CAPPED"


def test_a_locus_with_no_variable_site_is_recorded_not_dropped() -> None:
    positions, _ = microhap.variable_positions([], ("ctg", 900, 1200), 32)
    assert positions == []


class _Ctx:
    """The little of Context that the locus-basis helper reads."""

    def __init__(self, results: Path, params: dict[str, object]) -> None:
        self.results = results
        self.params = params


def _artifact(results: Path, rows: list[tuple[str, str, str, str]]) -> None:
    (results / "cohort").mkdir(parents=True, exist_ok=True)
    with gzip.open(results / "cohort" / "microhaplotypes.tsv.gz", "wt") as handle:
        handle.write("locus_id\tsample_id\thaplotype\tfragments\tfraction\tstate\treason\tincomplete_fragments\n")
        for locus, sample, haplotype, state in rows:
            handle.write(f"{locus}\t{sample}\t{haplotype}\t20\t0.5\t{state}\t\t0\n")


def test_microhaplotypes_become_the_locus_basis_when_there_are_enough(tmp_path: Path) -> None:
    rows = [(f"m{i}", s, hap, "CALLED") for i in range(10) for s, hap in (("S1", "ACG"), ("S2", "ATG"))]
    rows.append(("m0", "S1", "ATG", "CALLED"))      # a second haplotype in one sample
    rows.append(("m11", "S1", "", "UNCALLABLE"))    # not counted
    _artifact(tmp_path, rows)
    ctx = _Ctx(tmp_path, PARAMS)
    found = loci_module.allele_set(ctx, ["S1", "S2"], 0.05)  # type: ignore[arg-type]
    assert found.basis == loci_module.MICROHAPLOTYPE and found.loci == 10
    assert found.loci_per_sample == {"S1": 10, "S2": 10}
    assert {(r["locus"], r["allele"]) for r in found.rows if r["sample_id"] == "S1" and r["locus"] == "m0"} == \
        {("m0", "ACG"), ("m0", "ATG")}
    assert found.provenance["locus_basis"] == "MICROHAPLOTYPE"
    assert found.provenance["haplotypes"] >= 2


def test_too_few_microhaplotype_loci_keeps_the_biallelic_set_and_says_why(tmp_path: Path,
                                                                         monkeypatch: pytest.MonkeyPatch) -> None:
    _artifact(tmp_path, [("m0", "S1", "ACG", "CALLED"), ("m1", "S1", "ACG", "CALLED")])
    ctx = _Ctx(tmp_path, PARAMS)
    monkeypatch.setattr(loci_module, "biallelic_rows",
                        lambda c, eligible, floor: ([{"sample_id": "S1", "locus": "s1", "allele": "R"}],
                                                    {"S1": 40}, 40))
    found = loci_module.allele_set(ctx, ["S1"], 0.05)  # type: ignore[arg-type]
    assert found.basis == loci_module.BIALLELIC_SNV and found.loci == 40
    assert "only 2 microhaplotype loci" in found.reason
    assert found.provenance["microhaplotype_loci_available"] == 2


def test_no_artifact_falls_back_without_error(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    ctx = _Ctx(tmp_path, PARAMS)
    monkeypatch.setattr(loci_module, "biallelic_rows", lambda c, eligible, floor: ([], {"S1": 0}, 0))
    found = loci_module.allele_set(ctx, ["S1"], 0.05)  # type: ignore[arg-type]
    assert found.basis == loci_module.BIALLELIC_SNV and "no microhaplotype locus" in found.reason


def test_shipped_panel_has_the_documented_format_and_no_invented_loci() -> None:
    from pfvarx.downstream.resistance import load_catalog

    root = Path(__file__).resolve().parents[2]
    meta, rows = load_catalog(root / microhap.PANEL)
    assert rows == [], "the shipped panel must not contain invented locus coordinates"
    assert meta["status"] == "DRAFT_REQUIRES_VERIFICATION"
    header = next(line for line in (root / microhap.PANEL).read_text().splitlines()
                  if not line.startswith("#"))
    assert tuple(header.split("\t")) == microhap.PANEL_COLUMNS


def test_build_produces_locus_and_haplotype_rows(tmp_path: Path) -> None:
    """The whole path for one locus: resolve, check the territory, count fragments, call samples."""
    reads = [_spanning(f"a{i}", ("A", "C", "G")) for i in range(12)]
    reads += [_spanning(f"b{i}", ("A", "T", "G")) for i in range(8)]
    (tmp_path / "rich").mkdir()
    (tmp_path / "thin").mkdir()
    rich = _bam(tmp_path / "rich", reads)
    thin = _bam(tmp_path / "thin", [_spanning(f"c{i}", ("A", "C", "G")) for i in range(3)])
    sites = [{"chrom": "ctg", "pos": str(p + 1), "ref": "A", "alt": "C", "n_pass": "2"} for p in POSITIONS]
    panel = [{"locus_id": "m1", "chrom": "ctg", "start": "950", "end": "1150", "gene_id": "",
              "name": "amplicon 1", "source": "fixture"},
             {"locus_id": "m2", "chrom": "ctg", "start": "3000", "end": "3200", "gene_id": "",
              "name": "outside the panel", "source": "fixture"}]
    haplotypes, loci = microhap.build(
        sites=sites, panel=panel, samples=["RICH", "THIN"],
        bams={"RICH": rich, "THIN": thin}, callable_beds={}, lengths={"ctg": LENGTH},
        territory={"ctg": [(900, 1200)]}, design=Design.TARGETED_UMI, gff=None,
        params=PARAMS, min_mapq=20, min_bq=13)
    used = {row["locus_id"]: row for row in loci}
    assert used["m1"]["status"] == "USED" and used["m1"]["variable_positions"] == 3
    assert used["m1"]["positions"] == "1000,1050,1100"  # reported 1-based
    assert used["m1"]["samples_called"] == 1 and used["m1"]["samples_uncallable"] == 1
    assert used["m1"]["haplotypes"] == 2
    # A locus outside the panel's targets is recorded with its reason, not quietly dropped.
    assert used["m2"]["status"] == "EXCLUDED" and used["m2"]["reason"] == "OUTSIDE_PANEL_TARGETS"
    rich_rows = {r["haplotype"]: r for r in haplotypes if r["sample_id"] == "RICH" and r["state"] == "CALLED"}
    assert set(rich_rows) == {"ACG", "ATG"}
    assert rich_rows["ACG"]["fragments"] == 12 and rich_rows["ACG"]["fraction"] == pytest.approx(0.6)
    (thin_row,) = [r for r in haplotypes if r["sample_id"] == "THIN"]
    assert thin_row["state"] == "UNCALLABLE" and thin_row["reason"] == "LOW_FRAGMENT_DEPTH"


def test_build_skips_a_sample_where_the_positions_are_not_callable(tmp_path: Path) -> None:
    bam = _bam(tmp_path, [_spanning(f"a{i}", ("A", "C", "G")) for i in range(12)])
    bed = tmp_path / "callable.bed"
    bed.write_text("ctg\t0\t500\n")  # the locus lies outside the sample's callable region
    sites = [{"chrom": "ctg", "pos": str(p + 1), "ref": "A", "alt": "C", "n_pass": "2"} for p in POSITIONS]
    panel = [{"locus_id": "m1", "chrom": "ctg", "start": "950", "end": "1150", "gene_id": "", "name": "",
              "source": "fixture"}]
    haplotypes, loci = microhap.build(sites=sites, panel=panel, samples=["S1"], bams={"S1": bam},
                                      callable_beds={"S1": bed}, lengths={"ctg": LENGTH},
                                      territory={"ctg": [(900, 1200)]}, design=Design.WGS, gff=None,
                                      params=PARAMS, min_mapq=20, min_bq=13)
    assert [r["reason"] for r in haplotypes] == ["POSITION_NOT_CALLABLE"]
    assert loci[0]["samples_called"] == 0 and loci[0]["samples_uncallable"] == 1
