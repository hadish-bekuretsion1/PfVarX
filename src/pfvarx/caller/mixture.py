"""Haploid clone-mixture model (L2.6; Equations 3.17-3.26).

Pass one of the two-pass design: global number of technical components K and their
DNA proportions w are estimated from mixture-informative biallelic SNV sites with the
per-locus beta-binomial emission. K is selected by a parametric bootstrap of the
likelihood-ratio statistic. The result is a *technical mixture model of one sample's
reads*; it is never reported as biological complexity of infection.
"""

from __future__ import annotations

import itertools
import math
from dataclasses import dataclass, field
from typing import Any

import numpy as np
from scipy import optimize, special

from pfvarx.caller.stats import betabinom_logpmf


@dataclass
class MixtureFit:
    k: int
    weights: list[float]
    config_prior: list[float]
    log_lik: float
    converged: bool
    iterations: int
    starts_agreeing: int
    stable: bool
    seeds: list[int] = field(default_factory=list)


def _configs(k: int) -> np.ndarray:
    return np.array(list(itertools.product((0, 1), repeat=k)), dtype=float)  # (2^K, K)


def _site_loglik(y: np.ndarray, n: np.ndarray, w: np.ndarray, eps: float, phi: float,
                 configs: np.ndarray) -> np.ndarray:
    pi = configs @ w  # (C,)
    mu = pi * (1 - eps) + (1 - pi) * eps
    return betabinom_logpmf(y[:, None], n[:, None], mu[None, :], phi)  # (L, C)


def _em(y: np.ndarray, n: np.ndarray, k: int, eps: float, phi: float, rng: np.random.Generator,
        tol: float, max_iter: int) -> tuple[np.ndarray, np.ndarray, float, bool, int]:
    configs = _configs(k)
    c = len(configs)
    w = rng.dirichlet(np.ones(k))
    prior = np.full(c, 1.0 / c)
    prev = -math.inf
    ll = -math.inf
    converged = False
    it = 0
    for it in range(1, max_iter + 1):
        site = _site_loglik(y, n, w, eps, phi, configs) + np.log(prior)[None, :]
        norm = special.logsumexp(site, axis=1)
        ll = float(norm.sum())
        gamma = np.exp(site - norm[:, None])  # E-step (Equation 3.21)
        prior = np.clip(gamma.mean(axis=0), 1e-12, None)
        prior /= prior.sum()

        def q(theta: np.ndarray) -> tuple[float, np.ndarray]:
            """Negative M-step objective (Equation 3.22) and analytic gradient (softmax simplex)."""
            ww = special.softmax(np.concatenate([[0.0], theta]))
            pi = configs @ ww
            mu = np.clip(pi * (1 - eps) + (1 - pi) * eps, 1e-12, 1 - 1e-12)
            a_, b_ = mu * phi, (1 - mu) * phi
            ll_ = betabinom_logpmf(y[:, None], n[:, None], mu[None, :], phi)
            dmu = phi * (special.digamma(y[:, None] + a_[None, :]) - special.digamma(a_)[None, :]
                         - special.digamma(n[:, None] - y[:, None] + b_[None, :])
                         + special.digamma(b_)[None, :])
            gmu = (gamma * dmu).sum(axis=0) * (1 - 2 * eps)  # (C,)
            gw = configs.T @ gmu  # dQ/dw_k
            jac = ww[1:] * (gw[1:] - gw @ ww)
            return -float(np.sum(gamma * ll_)), -jac

        theta0 = np.log(np.clip(w[1:], 1e-9, None)) - math.log(max(w[0], 1e-9))
        if k > 1:
            res = optimize.minimize(q, theta0, jac=True, method="L-BFGS-B",
                                    bounds=[(-30.0, 30.0)] * (k - 1), options={"maxiter": 50})
            w = special.softmax(np.concatenate([[0.0], res.x]))
        else:
            w = np.array([1.0])
        if abs(ll - prev) <= tol * max(1.0, abs(ll)):
            converged = True
            break
        prev = ll
    order = np.argsort(-w)
    w = w[order]
    # re-index configuration prior to the new component order
    configs = _configs(k)
    reordered = configs[:, order]
    mapping = [int(np.where((configs == row).all(axis=1))[0][0]) for row in reordered]
    new_prior = np.zeros_like(prior)
    for old_idx, new_idx in enumerate(mapping):
        new_prior[new_idx] = prior[old_idx]
    return w, new_prior, ll, converged, it


def fit(y: np.ndarray, n: np.ndarray, k: int, eps: float, phi: float, seed: int, restarts: int,
        agreement: float, tol: float, max_iter: int, _doubled: bool = False) -> MixtureFit:
    """Best of *restarts* EM runs; accepted as stable when >= 2 starts agree within *agreement*.

    If unstable, the number of starts is doubled once; remaining instability is reported.
    """
    best = None
    lls = []
    seeds = [seed + i for i in range(restarts)]
    for s in seeds:
        w, prior, ll, conv, it = _em(y, n, k, eps, phi, np.random.default_rng(s), tol, max_iter)
        lls.append(ll)
        if best is None or ll > best[2]:
            best = (w, prior, ll, conv, it)
    assert best is not None
    agree = sum(1 for v in lls if best[2] - v <= agreement)
    stable = agree >= 2 or k == 1 or restarts == 1
    if not stable and not _doubled:
        return fit(y, n, k, eps, phi, seed + 10_000, restarts * 2, agreement, tol, max_iter, True)
    return MixtureFit(k, [float(v) for v in best[0]], [float(v) for v in best[1]], float(best[2]),
                      bool(best[3]), int(best[4]), agree, stable, seeds)


def simulate(n: np.ndarray, f: MixtureFit, eps: float, phi: float, rng: np.random.Generator) -> np.ndarray:
    configs = _configs(f.k)
    idx = rng.choice(len(configs), size=len(n), p=np.asarray(f.config_prior))
    pi = configs[idx] @ np.asarray(f.weights)
    mu = pi * (1 - eps) + (1 - pi) * eps
    p = rng.beta(np.maximum(mu * phi, 1e-9), np.maximum((1 - mu) * phi, 1e-9))
    return rng.binomial(n.astype(int), p)


def select_k(y: np.ndarray, n: np.ndarray, eps: float, phi: float, k_max: int, alpha: float,
             replicates: int, seed: int, restarts: int, agreement: float, tol: float,
             max_iter: int, min_weight: float = 0.0) -> dict[str, Any]:
    """Sequential parametric-bootstrap LRT (Equations 3.24-3.26).

    Smallest K whose test of K versus K+1 is not rejected at alpha. Replicate fits use
    fewer restarts (max(3, restarts//4)); both counts are recorded.
    """
    rng = np.random.default_rng(seed)
    fits = {1: fit(y, n, 1, eps, phi, seed, 1, agreement, tol, max_iter)}
    tests: list[dict[str, int | float | str | bool | None]] = []
    chosen = None
    for k in range(1, k_max + 1):
        if k + 1 > k_max:
            chosen = k
            tests.append({"k": k, "p_value": None, "note": "K_MAX reached; larger K not tested"})
            break
        fits[k + 1] = fit(y, n, k + 1, eps, phi, seed + 100 * k, restarts, agreement, tol, max_iter)
        lr = 2 * (fits[k + 1].log_lik - fits[k].log_lik)
        boot_restarts = max(2, restarts // 10)
        boot_iter = min(max_iter, 200)
        # Besag-Clifford sequential stopping: once h exceedances are seen, p >= alpha is
        # certain and the remaining replicates cannot change the decision.
        h = math.ceil(alpha * (replicates + 1))
        exceed = 0
        done = 0
        for b in range(replicates):
            yb = simulate(n, fits[k], eps, phi, rng)
            s = seed + 1_000_000 + 1000 * k + b
            l0 = fit(yb, n, k, eps, phi, s, boot_restarts, agreement, tol, boot_iter).log_lik
            l1 = fit(yb, n, k + 1, eps, phi, s + 7, boot_restarts, agreement, tol, boot_iter).log_lik
            exceed += int(2 * (l1 - l0) >= lr)
            done += 1
            if exceed >= h:
                break
        sequential_stop = done < replicates
        p = (1 + exceed) / (replicates + 1) if not sequential_stop else max(alpha, (1 + exceed) / (done + 1))
        tests.append({"k": k, "lr": lr, "p_value": p, "replicates_planned": replicates,
                      "replicates_run": done, "sequential_stop": sequential_stop,
                      "replicate_restarts": boot_restarts, "replicate_max_iter": boot_iter})
        if p >= alpha:
            chosen = k
            break
    assert chosen is not None
    # Identifiability boundary: a component with negligible weight is not a supported clone;
    # fall back to the largest smaller K whose fit has no such component, and record it.
    collapsed = []
    while chosen > 1 and min(fits[chosen].weights) < min_weight:
        collapsed.append({"k": chosen, "min_weight": min(fits[chosen].weights)})
        chosen -= 1
    f = fits[chosen]
    return {"k_selected": chosen, "fit": f.__dict__, "tests": tests, "boundary_collapses": collapsed,
            "min_component_weight": min_weight,
            "fits": {str(k): v.__dict__ for k, v in fits.items()},
            "label": "technical mixture components; NOT biological COI"}
