"""Exact closed-form power: sum of signal-fires probabilities over the Fisher
noncentral hypergeometric PMF.

This is the authoritative reference value against which the Wald-analytic and
Monte-Carlo paths are validated. Given fixed margins (n_drug, n_event, N) and
a true odds ratio OR, the only random quantity is the (drug, event) cell `a`,
distributed as nchypergeom_fisher(M=N, n=n_drug, N=n_event, odds_ratio=OR).

Power = sum_{a in support} pmf(a; OR) * 1{signal_fires(a, n_drug, n_event, N)}.

Used internally by tests/test_power_validation.py for G2; not used in the main
MDE search (which uses wald_power for speed; the wald-vs-MC agreement check is
the G2 gate).
"""
from __future__ import annotations

from typing import Callable

import numpy as np
from scipy.stats import nchypergeom_fisher

from .signal import signal_fires


def exact_power(
    n_drug: int,
    n_event: int,
    N: int,
    odds_ratio: float,
    *,
    signal_fn: Callable[[int, int, int, int], bool] = signal_fires,
) -> float:
    """Exact P(signal | margins, true OR) via PMF summation.

    Support: a in [max(0, n_drug+n_event-N), min(n_drug, n_event)].
    """
    a_lo = max(0, n_drug + n_event - N)
    a_hi = min(n_drug, n_event)
    if a_hi < a_lo:
        return 0.0
    ks = np.arange(a_lo, a_hi + 1, dtype=np.int64)
    # scipy nchypergeom_fisher signature: (M, n, N, odds_ratio) positional
    pmf = nchypergeom_fisher.pmf(ks, N, n_drug, n_event, odds_ratio)
    pmf = np.asarray(pmf, dtype=np.float64)
    # Guard against catastrophic underflow at extreme ORs:
    pmf = np.where(np.isfinite(pmf), pmf, 0.0)
    s = pmf.sum()
    if s > 0:
        pmf = pmf / s  # normalise if scipy returned an off-by-eps mass
    fires = np.fromiter(
        (signal_fn(int(k), n_drug, n_event, N) for k in ks),
        dtype=bool,
        count=ks.size,
    )
    return float(pmf[fires].sum())
