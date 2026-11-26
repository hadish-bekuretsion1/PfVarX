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
    variant_like_excluded: NotRequired[int]
    truncation: NotRequired[list[float]]


def betabinom_logpmf(y: np.ndarray, n: np.ndarray, mu: np.ndarray | float, phi: float) -> np.ndarray:
    y = np.asarray(y, dtype=float)
    n = np.asarray(n, dtype=float)
    mu = np.clip(np.asarray(mu, dtype=float), 1e-12, 1 - 1e-12)
    a, b = mu * phi, (1.0 - mu) * phi
    return (special.gammaln(n + 1) - special.gammaln(y + 1) - special.gammaln(n - y + 1)
            + special.betaln(y + a, n - y + b) - special.betaln(a, b))


def _sum_logpmf(y: np.ndarray, n: np.ndarray, mu: float, phi: float,
                index: tuple[np.ndarray, np.ndarray, np.ndarray] | None = None) -> float:
    """Sum of beta-binomial log-probabilities over sites.

    With *index* (unique y, unique n, inverse map from ``np.unique``), the log-probability is
    computed once per distinct (y, n) pair and expanded back to every site in the original
    order before summing: the same numbers are summed in the same order, so the result is
    identical to the direct computation, at the cost of only the distinct pairs."""
    if index is None:
        return float(np.sum(betabinom_logpmf(y, n, mu, phi)))
    uy, un, inverse = index
    return float(np.sum(betabinom_logpmf(uy, un, mu, phi)[inverse]))


def _unique_pairs(y: np.ndarray, n: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    pairs, inverse = np.unique(np.stack([y, n], axis=1), axis=0, return_inverse=True)
    return pairs[:, 0], pairs[:, 1], inverse.reshape(-1)


def betabinom_sf_log10(y: int, n: int, mu: float, phi: float) -> float:
    """-log10 P(Y >= y | n, mu, phi) (Equation 3.10). Exact summation of the upper tail."""
    if y <= 0:
        return 0.0
    if y > n:
        return math.inf
    k = np.arange(y, n + 1)
    lp = betabinom_logpmf(k, np.full_like(k, n), mu, phi)
    return float(-special.logsumexp(lp) / math.log(10))


def _truncation_ymax(n: np.ndarray, exclude: tuple[int, float]) -> np.ndarray:
    """Largest count kept at depth n when sites with >= min_reads and >= fraction are excluded."""
    min_reads, fraction = exclude
    return np.maximum(min_reads - 1, np.ceil(fraction * n) - 1).astype(int)


def _log_kept_mass(un: np.ndarray, ymax: np.ndarray, mu: float, phi: float) -> np.ndarray:
    """log P(Y <= ymax(n) | n) for each distinct depth (the truncation normaliser)."""
    out = np.empty(len(un))
    for i, (depth, top) in enumerate(zip(un, ymax, strict=True)):
        k = np.arange(0, min(int(top), int(depth)) + 1)
        out[i] = special.logsumexp(betabinom_logpmf(k, np.full(len(k), depth), mu, phi))
    return out


def fit_betabinom(y: np.ndarray, n: np.ndarray, trim_fraction: float,
                  init_mu: float, init_phi: float, *, phi_bounds: tuple[float, float],
                  min_sites: int, min_nonzero_for_phi: int,
                  exclude: tuple[int, float] | None = None) -> ErrorFit:
    """ML fit of (mu, phi) on known-reference positions after trimming the upper tail.

    Trimming is rank-based among positions with y > 0 (zero-count positions are never
    trimmed): the top *trim_fraction* of them by y/n are removed, so that unmasked true
    low-frequency variants do not inflate the error model. phi is kept within *phi_bounds*
    (layer2.overdispersion_phi_bounds; the lower bound because near-all-zero data cannot
    identify stronger overdispersion); the bounds are recorded with the fit. Fewer than
    *min_sites* sites (layer2.error_fit_min_sites) keep the starting values.
    """
    y = np.asarray(y, dtype=float)
    n = np.asarray(n, dtype=float)
    keep = n > 0
    y, n = y[keep], n[keep]
    excluded = 0
    if exclude is not None:
        # Sites whose count looks like a real allele (>= min_reads and >= fraction of reads)
        # are left out; the likelihood of the rest is truncated to match (see below). In a
        # mixed infection this keeps minor-clone variants from inflating the error model.
        variant_like = y > _truncation_ymax(n, exclude)
        excluded = int(variant_like.sum())
        y, n = y[~variant_like], n[~variant_like]
    if len(y) < min_sites or y.sum() == 0:
        return {"mu": init_mu, "phi": init_phi, "n_sites": int(len(y)), "fitted": False,
                "trimmed": 0, "reason": "too few sites or no error events"}
    nz = np.where(y > 0)[0]
    drop = int(np.ceil(trim_fraction * len(nz)))
    sel = np.ones(len(y), dtype=bool)
    if drop:
        order = nz[np.argsort(-(y[nz] / n[nz]), kind="stable")]
        sel[order[:drop]] = False
    y, n = y[sel], n[sel]
    index = _unique_pairs(y, n)
    depths, depth_inverse = np.unique(n, return_inverse=True)
    ymax = _truncation_ymax(depths, exclude) if exclude is not None else depths.astype(int)

    def truncated(total: float, mu: float, phi: float) -> float:
        """Subtract each site's log P(Y <= ymax(n)): the truncated (censored-subset) likelihood."""
        if exclude is None:
            return total
        return total - float(np.sum(_log_kept_mass(depths, ymax, mu, phi)[depth_inverse]))

    if len(nz) - drop < min_nonzero_for_phi:
        # Too few error events to identify overdispersion: shrink phi fully to the mode
        # default and estimate mu only (proposal: shrinkage toward phi(m)).
        res1 = optimize.minimize_scalar(
            lambda t: -truncated(_sum_logpmf(y, n, special.expit(t), init_phi, index), special.expit(t), init_phi),
            bounds=(-18.0, 0.0), method="bounded")
        shrunk: ErrorFit = {"mu": float(special.expit(res1.x)), "phi": float(init_phi), "n_sites": int(len(y)),
                            "fitted": True, "trimmed": int(drop), "phi_source": "MODE_DEFAULT_SHRINKAGE",
                            "nonzero_sites": int(len(nz) - drop)}
        if exclude is not None:
            shrunk["variant_like_excluded"] = excluded
            shrunk["truncation"] = [float(exclude[0]), float(exclude[1])]
        return shrunk

    def nll(theta: np.ndarray) -> float:
        mu = special.expit(theta[0])
        phi = math.exp(theta[1])
        return -truncated(_sum_logpmf(y, n, mu, phi, index), mu, phi)

    start = np.array([special.logit(max(min(init_mu, 0.5), 1e-6)),
                      math.log(min(max(init_phi, phi_bounds[0]), phi_bounds[1]))])
    res = optimize.minimize(nll, start, method="L-BFGS-B",
                            bounds=[(-18.0, 0.0), (math.log(phi_bounds[0]), math.log(phi_bounds[1]))])
    result: ErrorFit = {"mu": float(special.expit(res.x[0])), "phi": float(math.exp(res.x[1])),
                        "n_sites": int(len(y)), "fitted": bool(res.success), "trimmed": int(drop),
                        "neg_log_lik": float(res.fun), "phi_bounds": list(phi_bounds),
                        "phi_source": "SAMPLE_FIT", "nonzero_sites": int(len(nz) - drop)}
    if exclude is not None:
        result["variant_like_excluded"] = excluded
        result["truncation"] = [float(exclude[0]), float(exclude[1])]
    return result


def profile_interval(y: int, n: int, eps: float, phi: float, level: float,
                     grid_points: int) -> tuple[float, float, float]:
    """Profile-likelihood MLE and interval for the DNA fraction pi (Equations 3.18-3.19), on
    *grid_points* values in [0, 1] (layer2.wsaf_profile_grid_points)."""
    grid = np.linspace(0.0, 1.0, grid_points)
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
