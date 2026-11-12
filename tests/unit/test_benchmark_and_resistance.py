from pathlib import Path

import numpy as np
import pytest

from pfvarx.benchmarking.evaluate import wild_cluster_interval
from pfvarx.downstream.resistance import CODE, resolve_codon


def test_wild_cluster_enumerates_webb_for_few_clusters() -> None:
    r = wild_cluster_interval(np.array([90.0, 80, 95, 85]), np.array([100.0, 100, 100, 100]))
    assert r["weights"] == "Webb six-point" and r["vectors"] == "enumerated"
    assert r["lo"] < r["estimate"] < r["hi"] and r["estimate"] == pytest.approx(0.875)


def test_wild_cluster_rademacher_for_many_clusters() -> None:
    r = wild_cluster_interval(np.full(12, 50.0), np.full(12, 100.0))
    assert r["weights"] == "Rademacher"


def test_codon_resolution_minus_strand(tmp_path: Path) -> None:
    (tmp_path / "g.gff").write_text("c1\tx\tCDS\t101\t106\t.\t-\t0\tID=c;Parent=G1.1\n"
                                    "c1\tx\tCDS\t201\t203\t.\t-\t0\tID=c;Parent=G1.1\n")
    chrom, pos, strand = resolve_codon(tmp_path / "g.gff", "G1", 1)  # type: ignore[misc]
    assert strand == "-" and pos == [202, 201, 200]
    chrom, pos, strand = resolve_codon(tmp_path / "g.gff", "G1", 2)  # type: ignore[misc]
    assert pos == [105, 104, 103]


def test_genetic_code_complete() -> None:
    assert len(CODE) == 64 and CODE["TGT"] == "Cys" and CODE["TAT"] == "Tyr"
