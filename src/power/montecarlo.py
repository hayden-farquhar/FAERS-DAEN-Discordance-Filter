"""Monte-Carlo power via Fisher noncentral hypergeometric simulation
(protocol Section 6.1, step 3 'Monte-Carlo' path).

For given margins and a true odds ratio, draw B samples of the (drug, event)
cell `a` from `scipy.stats.nchypergeom_fisher(M=N, n=n_drug, N=n_event,
odds_ratio=OR)` and report the fraction firing the signal rule.

The seed is supplied by the caller (protocol Section 15.3: seed=94 for the
project-default generator).
"""
from __future__ import annotations

import numpy as np
from scipy.stats import nchypergeom_fisher

from .signal import signal_fires

_DEFAULT_B = 5000


def mc_power(
    n_drug: int,
    n_event: int,
    N: int,
    OR_alt: float,
    *,
    B: int = _DEFAULT_B,
    rng: np.random.Generator | None = None,
) -> float:
    """Monte-Carlo estimate of P(signal | margins, true OR_alt).

    Caller-supplied rng (or default_rng(94) if None) governs reproducibility.
    """
    if rng is None:
        rng = np.random.default_rng(94)
    if n_drug <= 0 or n_event <= 0 or N <= 0:
        return 0.0
    # scipy.stats.nchypergeom_fisher.rvs requires shape args as positional
    # (M, n, N, odds_ratio); only loc/size/random_state accept keywords.
    samples = nchypergeom_fisher.rvs(
        N, n_drug, n_event, OR_alt,
        size=B, random_state=rng,
    )
    fires = np.fromiter(
        (signal_fires(int(a), n_drug, n_event, N) for a in samples),
        dtype=bool,
        count=samples.size,
    )
    return float(fires.mean())
