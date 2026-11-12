import math

import numpy as np
import pytest
from scipy import stats

from pfvarx.caller import mixture, pairhmm
from pfvarx.caller.stats import betabinom_logpmf, betabinom_sf_log10, fit_betabinom, profile_interval


def test_betabinom_matches_scipy() -> None:
    mu, phi = 0.01, 50.0
    ref = stats.betabinom(40, mu * phi, (1 - mu) * phi)
    for k in (0, 1, 5, 40):
        assert betabinom_logpmf(np.array([k]), np.array([40]), mu, phi)[0] == pytest.approx(ref.logpmf(k), rel=1e-9)
    assert betabinom_sf_log10(5, 40, mu, phi) == pytest.approx(-math.log10(ref.sf(4)), rel=1e-6)


def test_fit_recovers_parameters_and_trims_tail() -> None:
    rng = np.random.default_rng(1)
    n = np.full(20000, 60)
    p = rng.beta(0.002 * 300, 0.998 * 300, size=n.size)
    y = rng.binomial(n, p)
    y[:20] = 30  # unmasked true variants in the "known-reference" set
    f = fit_betabinom(y, n, 0.025, 0.005, 200)
    assert f["mu"] == pytest.approx(0.002, rel=0.15) and f["phi_source"] == "SAMPLE_FIT"


def test_profile_interval_contains_truth() -> None:
    est, lo, hi = profile_interval(30, 100, 0.001, 300, 0.95)
    assert lo < 0.3 < hi and est == pytest.approx(0.3, abs=0.01)


def test_pairhmm_prefers_true_haplotype_and_supports_indel() -> None:
    hap = "ACGTTAGCATTTTTGACCAGTACGATCGA"
    alt = hap[:14] + "T" + hap[14:]  # one-base insertion in the T tract
    read_alt = alt[3:25]
    q = [35] * len(read_alt)
    assert pairhmm.log10_likelihood(read_alt, q, alt, 1e-3, 0.1) > pairhmm.log10_likelihood(read_alt, q, hap, 1e-3, 0.1)
    assert pairhmm.support(read_alt, q, hap, alt, 1e-3, 0.1) == 1
    read_ref = hap[3:25]
    assert pairhmm.support(read_ref, q, hap, alt, 1e-3, 0.1) == -1


def _sim(weights: list[float], sites: int, seed: int) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    k = len(weights)
    a = rng.integers(0, 2, size=(sites, k))
    a[a.sum(1) == 0, 0] = 1
    pi = a @ np.array(weights)
    n = rng.integers(40, 80, size=sites)
    return rng.binomial(n, pi * (1 - 3e-4) + (1 - pi) * 3e-4), n


def test_em_recovers_two_clone_weights() -> None:
    y, n = _sim([0.7, 0.3], 150, 3)
    f = mixture.fit(y, n, 2, 3e-4, 500.0, seed=1, restarts=10, agreement=1e-4, tol=1e-8, max_iter=500)
    assert f.weights[0] == pytest.approx(0.7, abs=0.05) and f.stable


def test_k_selection_monoclonal_and_mixture() -> None:
    y1, n1 = _sim([1.0], 80, 5)
    r1 = mixture.select_k(y1, n1, 3e-4, 500.0, 3, 0.05, 39, 11, 4, 1e-4, 1e-8, 300)
    assert r1["k_selected"] == 1 and "NOT biological COI" in r1["label"]
    y2, n2 = _sim([0.65, 0.35], 120, 7)
    r2 = mixture.select_k(y2, n2, 3e-4, 500.0, 3, 0.05, 39, 13, 4, 1e-4, 1e-8, 300)
    assert r2["k_selected"] >= 2


def test_negligible_component_collapses_to_smaller_k() -> None:
    y, n = _sim([1.0], 80, 5)
    y = y.copy()
    y[:2] = (n[:2] * 0.03).astype(int)  # two outlier sites that a tiny extra component could absorb
    r = mixture.select_k(y, n, 3e-4, 500.0, 2, 0.05, 39, 17, 4, 1e-4, 1e-8, 300, min_weight=0.05)
    assert r["k_selected"] == 1
    assert all(c["min_weight"] < 0.05 for c in r["boundary_collapses"])
