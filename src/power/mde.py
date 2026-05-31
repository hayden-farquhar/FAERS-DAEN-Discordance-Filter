"""Per-pair minimum-detectable-effect (MDE) search (protocol Section 6.1).

For a pair with DAEN marginals (n_drug, n_event, N), the MDE is the smallest
true OR at which DAEN has at least the target power (default 0.80) to fire the
pre-registered signal rule (LB > 1 AND a >= 3).

Search grid: {1.10, 1.15, ..., 5.00} step 0.05 (79 candidates).
"""
from __future__ import annotations

import math

import numpy as np

from .exact import exact_power
from .wald import wald_power  # noqa: F401 — exported for documented sensitivity use


def _or_grid(start: float = 1.10, stop: float = 5.00, step: float = 0.05) -> np.ndarray:
    return np.round(np.arange(start, stop + step / 2, step), 4)


_DEFAULT_GRID = _or_grid()


def daen_mde(
    n_drug: int,
    n_event: int,
    N: int,
    *,
    power_target: float = 0.80,
    power_fn=exact_power,
    grid: np.ndarray = _DEFAULT_GRID,
) -> float:
    """Return the smallest OR in the grid at which power_fn >= power_target.

    Returns math.inf if no OR in the grid reaches the target (DAEN cannot
    detect any effect in this range at this pair's marginals).
    """
    if n_drug <= 0 or n_event <= 0 or N <= 0:
        return math.inf
    # Power is monotone-non-decreasing in OR (for OR >= 1) under nchypergeom,
    # so we can binary-search; but the grid is small, so a linear scan is fine
    # and avoids a monotonicity assumption corner-case.
    for OR in grid:
        p = power_fn(n_drug, n_event, N, float(OR))
        if p >= power_target:
            return float(OR)
    return math.inf


def daen_mde_unconditional(
    n_drug: int,
    n_event: int,
    N: int,
    *,
    power_target: float = 0.80,
    grid: np.ndarray = _DEFAULT_GRID,
    rng: np.random.Generator | None = None,
    B: int = 2000,
) -> float:
    """Poisson-margin sensitivity MDE (protocol Section 6.4).

    Treats the event-side margin as Poisson(lambda = n_event) rather than fixed.
    For each candidate OR:
      1. Sample event-count column from Poisson(n_event), B replicates.
      2. Conditional on each simulated event count, sample `a` from
         Binomial(n_drug, p) where p is set by OR via p_event_marginal.
      3. Apply signal rule to each (a, n_drug, simulated-n_event, N).
      4. Empirical power = proportion firing.
    Return smallest OR achieving power_target.
    """
    if rng is None:
        rng = np.random.default_rng(94)
    if n_drug <= 0 or n_event <= 0 or N <= 0:
        return math.inf
    from .signal import signal_fires
    for OR in grid:
        # p_event = expected event-cases-in-drug / n_drug under OR
        # Approximation: drug-side event rate ~ (n_event/N * OR) / (1 + (OR-1) * n_event/N)
        p_event_marg = n_event / N
        p_event_in_drug = (p_event_marg * OR) / (1.0 + (OR - 1.0) * p_event_marg)
        sim_n_event = rng.poisson(lam=n_event, size=B)
        # For each simulated n_event, simulate a from Binomial(n_drug, p_event_in_drug)
        sim_a = rng.binomial(n_drug, p_event_in_drug, size=B)
        fires_count = 0
        for i in range(B):
            n_e = int(sim_n_event[i])
            a = int(sim_a[i])
            # Clamp a to feasible region given the (random) n_e
            if a > min(n_drug, n_e):
                a = min(n_drug, n_e)
            if a < 0:
                a = 0
            if signal_fires(a, n_drug, n_e, N):
                fires_count += 1
        power = fires_count / B
        if power >= power_target:
            return float(OR)
    return math.inf
