"""Vectorised pair hidden Markov model (L2.5; Equations 3.11-3.16).

Semi-global: the read is aligned end to end, the haplotype locally (uniform start
over haplotype positions, free end), the convention used for read-to-haplotype
likelihoods in local-reassembly callers. Computation is in probability space with
per-row scaling; the recurrence along the haplotype for the deletion state is a
first-order linear filter solved with ``scipy.signal.lfilter``.
"""

from __future__ import annotations

import math

import numpy as np
from scipy.signal import lfilter


def _encode(seq: str) -> np.ndarray:
    return np.frombuffer(seq.upper().encode(), dtype=np.uint8).copy()


def log10_likelihood(read: str, quals: np.ndarray | list[int], hap: str, gap_open: float,
                     gap_extend: float, insert_emission: float = 0.25) -> float:
    """log10 Pr(read | haplotype)."""
    if not read or not hap:
        raise ValueError("empty read or haplotype")
    if len(quals) != len(read):
        raise ValueError("quality length differs from read length")
    if not (0 < gap_open < 0.5 and 0 < gap_extend < 1):
        raise ValueError("invalid gap parameters")
    r = _encode(read)
    h = _encode(hap)
    q = np.asarray(quals, dtype=float)
    err = 10.0 ** (-q / 10.0)
    n, m = len(r), len(h)
    delta, eps = gap_open, gap_extend
    mm = 1.0 - 2.0 * delta
    gm = 1.0 - eps
    # Row 0: deletions may start anywhere (local start on the haplotype).
    M = np.zeros(m + 1)
    Ins = np.zeros(m + 1)
    D = np.full(m + 1, 1.0 / m)
    D[m] = 0.0
    log_scale = 0.0
    for i in range(1, n + 1):
        match = (h == r[i - 1])
        e = np.where(match, 1.0 - err[i - 1], err[i - 1] / 3.0)
        e[h == 78] = 0.25  # 'N' in haplotype
        if r[i - 1] == 78:  # 'N' in read
            e[:] = 0.25
        newM = np.zeros(m + 1)
        newM[1:] = e * (mm * M[:-1] + gm * (Ins[:-1] + D[:-1]))
        newI = insert_emission * (delta * M + eps * Ins)
        # D[i, j] = delta * M[i, j-1] + eps * D[i, j-1]  -> linear recurrence along j
        drive = np.zeros(m + 1)
        drive[1:] = delta * newM[:-1]
        newD = lfilter([1.0], [1.0, -eps], drive)
        scale = max(newM.max(), newI.max(), newD.max(), 1e-300)
        M, Ins, D = newM / scale, newI / scale, newD / scale
        log_scale += math.log(scale)
    total = M[1:].sum() + Ins[1:].sum()
    return (math.log(max(total, 1e-300)) + log_scale) / math.log(10)


def support(read: str, quals: list[int], ref_hap: str, alt_hap: str, gap_open: float,
            gap_extend: float, min_log10_ratio: float = 1.0) -> int:
    """+1 read supports alt, -1 supports ref, 0 uninformative."""
    lr = (log10_likelihood(read, quals, alt_hap, gap_open, gap_extend)
          - log10_likelihood(read, quals, ref_hap, gap_open, gap_extend))
    if lr >= min_log10_ratio:
        return 1
    if lr <= -min_log10_ratio:
        return -1
    return 0
