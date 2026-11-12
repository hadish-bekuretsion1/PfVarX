"""Numerical primitives for PfVarX-Call (L2.4): beta-binomial likelihoods and fitting.

Parameterisation follows Equations 3.2-3.5: epsilon ~ Beta(mu*phi, (1-mu)*phi),
rho = 1/(phi + 1). All functions work in the log domain and are vectorised with numpy.
"""

from __future__ import annotations

import math
from typing import NotRequired, TypedDict

import numpy as np
from scipy import optimize, special, stats


class ErrorFit(TypedDict):
    mu: float
    phi: float
    n_sites: int
    fitted: bool
    trimmed: int
    reason: NotRequired[str]
    phi_source: NotRequired[str]
    nonzero_sites: NotRequired[int]
    neg_log_lik: NotRequired[float]
    phi_bounds: NotRequired[list[float]]


def betabinom_logpmf(y: np.ndarray, n: np.ndarray, mu: np.ndarray | float, phi: float) -> np.ndarray:
    y = np.asarray(y, dtype=float)
    n = np.asarray(n, dtype=float)
    mu = np.clip(np.asarray(mu, dtype=float), 1e-12, 1 - 1e-12)
    a, b = mu * phi, (1.0 - mu) * phi
    return (special.gammaln(n + 1) - special.gammaln(y + 1) - special.gammaln(n - y + 1)
            + special.betaln(y + a, n - y + b) - special.betaln(a, b))


def betabinom_sf_log10(y: int, n: int, mu: float, phi: float) -> float:
    """-log10 P(Y >= y | n, mu, phi) (Equation 3.10). Exact summation of the upper tail."""
    if y <= 0:
        return 0.0
    if y > n:
        return math.inf
    k = np.arange(y, n + 1)
    lp = betabinom_logpmf(k, np.full_like(k, n), mu, phi)
    return float(-special.logsumexp(lp) / math.log(10))


def fit_betabinom(y: np.ndarray, n: np.ndarray, trim_fraction: float,
                  init_mu: float, init_phi: float, phi_bounds: tuple[float, float] = (10.0, 1e7),
                  min_nonzero_for_phi: int = 50) -> ErrorFit:
    """ML fit of (mu, phi) on known-reference positions after trimming the upper tail.

    Trimming is rank-based among positions with y > 0 (zero-count positions are never
    trimmed): the top *trim_fraction* of them by y/n are removed, so that unmasked true
    low-frequency variants do not inflate the error model. phi is bounded below at 10
    (intra-sample correlation rho <= 0.09) because near-all-zero data cannot identify
    stronger overdispersion; the bound is FRAMEWORK_POLICY and recorded.
    """
    y = np.asarray(y, dtype=float)
    n = np.asarray(n, dtype=float)
    keep = n > 0
    y, n = y[keep], n[keep]
    if len(y) < 50 or y.sum() == 0:
        return {"mu": init_mu, "phi": init_phi, "n_sites": int(len(y)), "fitted": False,
                "trimmed": 0, "reason": "too few sites or no error events"}
    nz = np.where(y > 0)[0]
    drop = int(np.ceil(trim_fraction * len(nz)))
    sel = np.ones(len(y), dtype=bool)
    if drop:
        order = nz[np.argsort(-(y[nz] / n[nz]), kind="stable")]
        sel[order[:drop]] = False
    y, n = y[sel], n[sel]

    if len(nz) - drop < min_nonzero_for_phi:
        # Too few error events to identify overdispersion: shrink phi fully to the mode
        # default and estimate mu only (proposal: shrinkage toward phi(m)).
        res1 = optimize.minimize_scalar(
            lambda t: -float(np.sum(betabinom_logpmf(y, n, special.expit(t), init_phi))),
            bounds=(-18.0, 0.0), method="bounded")
        return {"mu": float(special.expit(res1.x)), "phi": float(init_phi), "n_sites": int(len(y)),
                "fitted": True, "trimmed": int(drop), "phi_source": "MODE_DEFAULT_SHRINKAGE",
                "nonzero_sites": int(len(nz) - drop)}

    def nll(theta: np.ndarray) -> float:
        mu = special.expit(theta[0])
        phi = math.exp(theta[1])
        return -float(np.sum(betabinom_logpmf(y, n, mu, phi)))

    start = np.array([special.logit(max(min(init_mu, 0.5), 1e-6)),
                      math.log(min(max(init_phi, phi_bounds[0]), phi_bounds[1]))])
    res = optimize.minimize(nll, start, method="L-BFGS-B",
                            bounds=[(-18.0, 0.0), (math.log(phi_bounds[0]), math.log(phi_bounds[1]))])
    return {"mu": float(special.expit(res.x[0])), "phi": float(math.exp(res.x[1])),
            "n_sites": int(len(y)), "fitted": bool(res.success), "trimmed": int(drop),
            "neg_log_lik": float(res.fun), "phi_bounds": list(phi_bounds),
            "phi_source": "SAMPLE_FIT", "nonzero_sites": int(len(nz) - drop)}


def profile_interval(y: int, n: int, eps: float, phi: float, level: float) -> tuple[float, float, float]:
    """Profile-likelihood MLE and interval for the DNA fraction pi (Equations 3.18-3.19)."""
    grid = np.linspace(0.0, 1.0, 2001)
    mu = grid * (1 - eps) + (1 - grid) * eps
    ll = betabinom_logpmf(np.full_like(grid, y), np.full_like(grid, n), mu, phi)
    best = int(np.argmax(ll))
    crit = stats.chi2.ppf(level, df=1) / 2.0
    inside = grid[ll >= ll[best] - crit]
    return float(grid[best]), float(inside.min()), float(inside.max())


def fisher_strand_neglog10(ref_f: int, ref_r: int, alt_f: int, alt_r: int) -> float:
    if alt_f + alt_r == 0 or ref_f + ref_r == 0:
        return 0.0
    p = stats.fisher_exact([[ref_f, ref_r], [alt_f, alt_r]]).pvalue
    return float(-math.log10(max(p, 1e-300)))
