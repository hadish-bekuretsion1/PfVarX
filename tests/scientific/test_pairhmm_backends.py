"""The compiled pair-HMM kernel gives the same likelihoods and read decisions as the numpy reference."""

import os
import subprocess
import sys

import numpy as np
import pytest

from pfvarx.caller import pairhmm

pytest.importorskip("numba")

GAPS = [(1e-3, 0.1), (1e-2, 0.5), (0.3, 0.9)]


def _cases(count: int, seed: int) -> list[tuple[str, list[int], str, str]]:
    """AT-rich haplotypes (as in P. falciparum) with homopolymer indels, mismatches and N bases."""
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(count):
        length = int(rng.integers(30, 280))
        hap = "".join(rng.choice(list("AAAATTTTCGN"), size=length, p=[.2, .2, .1, .1, .1, .1, .05, .05, .04, .04, .02]))
        pos = int(rng.integers(5, length - 5))
        alt = hap[:pos] + "T" * int(rng.integers(1, 4)) + hap[pos:] if rng.random() < .5 else hap[:pos] + hap[pos + 2:]
        source = alt if rng.random() < .5 else hap
        start = int(rng.integers(0, max(1, len(source) // 4)))
        read = "".join(b if rng.random() > .03 else rng.choice(list("ACGTN")) for b in source[start:start + int(rng.integers(10, 151))])
        quals = [int(q) for q in rng.integers(2, 42, size=len(read))]
        out.append((read, quals, hap, alt))
    return out


@pytest.mark.parametrize(("gap_open", "gap_extend"), GAPS)
def test_compiled_kernel_matches_reference(gap_open: float, gap_extend: float) -> None:
    assert pairhmm.BACKEND == "numba"
    for read, quals, hap, alt in _cases(300, seed=int(gap_open * 1e4)):
        for target in (hap, alt):
            expected = pairhmm.log10_likelihood_reference(read, quals, target, gap_open, gap_extend)
            got = pairhmm.log10_likelihood(read, quals, target, gap_open, gap_extend)
            assert got == pytest.approx(expected, rel=1e-12, abs=1e-12)


def test_read_decisions_are_identical() -> None:
    """support() is what the caller uses: +1/-1/0 must agree read by read."""
    def reference(read: str, quals: list[int], ref: str, alt: str) -> int:
        lr = (pairhmm.log10_likelihood_reference(read, quals, alt, 1e-3, 0.1)
              - pairhmm.log10_likelihood_reference(read, quals, ref, 1e-3, 0.1))
        return 1 if lr >= 1.0 else -1 if lr <= -1.0 else 0

    cases = _cases(1000, seed=7)
    decisions = [pairhmm.support(read, quals, hap, alt, 1e-3, 0.1) for read, quals, hap, alt in cases]
    assert decisions == [reference(*case) for case in cases]
    assert {-1, 0, 1} <= set(decisions)  # the test reaches every outcome


def test_validation_is_unchanged() -> None:
    with pytest.raises(ValueError, match="empty"):
        pairhmm.log10_likelihood("", [], "ACGT", 1e-3, 0.1)
    with pytest.raises(ValueError, match="quality length"):
        pairhmm.log10_likelihood("ACG", [30], "ACGT", 1e-3, 0.1)
    with pytest.raises(ValueError, match="gap"):
        pairhmm.log10_likelihood("ACG", [30] * 3, "ACGT", 0.6, 0.1)


def test_reference_backend_can_be_forced() -> None:
    env = dict(os.environ, PFVARX_PAIRHMM_BACKEND="reference")
    out = subprocess.run([sys.executable, "-c", "from pfvarx.caller import pairhmm; print(pairhmm.BACKEND)"],
                         env=env, capture_output=True, text=True, check=True)
    assert out.stdout.strip() == "reference"
