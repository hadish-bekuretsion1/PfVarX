"""Pair hidden Markov model (L2.5; Equations 3.11-3.16).

Semi-global: the read is aligned end to end, the haplotype locally (uniform start
over haplotype positions, free end), the convention used for read-to-haplotype
likelihoods in local-reassembly callers. Computation is in probability space with
per-row scaling.

Two implementations of the same recurrences:

* ``log10_likelihood_reference``: vectorised numpy, the deletion-state recurrence along
  the haplotype solved with ``scipy.signal.lfilter``. It is the reference implementation
  and is kept unchanged.
* a Numba-compiled kernel with the same arithmetic written as loops, used when Numba is
  installed. It only removes per-row Python and numpy call overhead; the model, its
  parameters and the decision rule in ``support`` are the same.

``log10_likelihood`` uses the compiled kernel when available, otherwise the reference.
Set ``PFVARX_PAIRHMM_BACKEND=reference`` to force the reference implementation (for
checking that both give the same calls). ``BACKEND`` names the implementation in use.
"""

from __future__ import annotations

import math
import os
import warnings
from collections.abc import Callable
from typing import Any

import numpy as np
from scipy.signal import lfilter


def _encode(seq: str) -> np.ndarray:
    return np.frombuffer(seq.upper().encode(), dtype=np.uint8).copy()


def _check(read: str, quals: np.ndarray | list[int], hap: str, gap_open: float, gap_extend: float) -> None:
    if not read or not hap:
        raise ValueError("empty read or haplotype")
    if len(quals) != len(read):
        raise ValueError("quality length differs from read length")
    if not (0 < gap_open < 0.5 and 0 < gap_extend < 1):
        raise ValueError("invalid gap parameters")


def log10_likelihood_reference(read: str, quals: np.ndarray | list[int], hap: str, gap_open: float,
                               gap_extend: float, insert_emission: float = 0.25) -> float:
    """log10 Pr(read | haplotype): the vectorised numpy reference implementation."""
    _check(read, quals, hap, gap_open, gap_extend)
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


def _kernel(r: np.ndarray, h: np.ndarray, err: np.ndarray, delta: float, eps: float,
            insert_emission: float) -> float:
    """The reference recurrences as explicit loops (compiled by Numba when it is installed)."""
    n, m = r.size, h.size
    mm = 1.0 - 2.0 * delta
    gm = 1.0 - eps
    M = np.zeros(m + 1)
    Ins = np.zeros(m + 1)
    D = np.full(m + 1, 1.0 / m)
    D[m] = 0.0
    newM = np.empty(m + 1)
    newI = np.empty(m + 1)
    newD = np.empty(m + 1)
    log_scale = 0.0
    for i in range(1, n + 1):
        ri = r[i - 1]
        e_match = 1.0 - err[i - 1]
        e_mismatch = err[i - 1] / 3.0
        newM[0] = 0.0
        for j in range(1, m + 1):
            hj = h[j - 1]
            if ri == 78 or hj == 78:  # 'N' in read or haplotype
                e = 0.25
            elif hj == ri:
                e = e_match
            else:
                e = e_mismatch
            newM[j] = e * (mm * M[j - 1] + gm * (Ins[j - 1] + D[j - 1]))
        for j in range(m + 1):
            newI[j] = insert_emission * (delta * M[j] + eps * Ins[j])
        # D[i, j] = delta * M[i, j-1] + eps * D[i, j-1]
        newD[0] = 0.0
        for j in range(1, m + 1):
            newD[j] = delta * newM[j - 1] + eps * newD[j - 1]
        scale = max(newM.max(), newI.max(), newD.max(), 1e-300)
        for j in range(m + 1):
            M[j] = newM[j] / scale
            Ins[j] = newI[j] / scale
            D[j] = newD[j] / scale
        log_scale += math.log(scale)
    total = M[1:].sum() + Ins[1:].sum()
    return (math.log(max(total, 1e-300)) + log_scale) / math.log(10.0)


def _compiled() -> Callable[..., float] | None:
    """The Numba-compiled kernel, or None without Numba (or with the reference forced)."""
    if os.environ.get("PFVARX_PAIRHMM_BACKEND", "").lower() == "reference":
        return None
    try:
        import numba
    except ImportError:
        return None
    try:
        # Compiled code is cached next to this file (or in NUMBA_CACHE_DIR) for later runs.
        # fastmath stays off: it would let the compiler reorder floating-point arithmetic.
        jitted: Any = numba.njit(cache=True)(_kernel)
    except RuntimeError:  # no writable cache location: compile in every process instead
        warnings.warn("Numba cannot cache the pair-HMM kernel; it is compiled in each process",
                      RuntimeWarning, stacklevel=2)
        jitted = numba.njit(cache=False)(_kernel)
    return jitted  # type: ignore[no-any-return]


_KERNEL = _compiled()
BACKEND = "numba" if _KERNEL is not None else "reference"


def log10_likelihood(read: str, quals: np.ndarray | list[int], hap: str, gap_open: float,
                     gap_extend: float, insert_emission: float = 0.25) -> float:
    """log10 Pr(read | haplotype)."""
    if _KERNEL is None:
        return log10_likelihood_reference(read, quals, hap, gap_open, gap_extend, insert_emission)
    _check(read, quals, hap, gap_open, gap_extend)
    err = 10.0 ** (-np.asarray(quals, dtype=float) / 10.0)
    return _KERNEL(_encode(read), _encode(hap), err, float(gap_open), float(gap_extend), float(insert_emission))


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
