"""Vectorised exact Fisher-NCHG power and MDE search (production fast path).

Numerically identical (to floating-point tolerance) to :func:`power.exact.exact_power`
and :func:`power.mde.daen_mde`, but evaluates the whole OR grid in a single
pure-numpy pass instead of one ``scipy.stats.nchypergeom_fisher`` call per OR.

Why this is exact, not an approximation. Fisher's noncentral hypergeometric
weight factorises as

    w(a; OR) = C(n_drug, a) * C(N - n_drug, n_event - a) * OR**a,

so in log space  log w(a; OR) = g(a) + a*log(OR)  with the OR-free part

    g(a) = -gammaln(a+1) - gammaln(n_drug-a+1)
           - gammaln(n_event-a+1) - gammaln(N-n_drug-n_event+a+1)

(the two a-independent gammaln terms cancel under normalisation). The signal
rejection region {a : signal_fires(a, ...)} depends only on the margins, not on
OR, so it too is computed once. Power at any OR is then

    P(reject; OR) = exp( logsumexp(log w[fires]) - logsumexp(log w[all]) ),

a couple of vectorised reductions over the support. ``exact.py`` remains the
authoritative reference; ``tests/test_power_fast.py`` pins this path to it.
"""
from __future__ import annotations

import math

import numpy as np
from scipy.special import gammaln, logsumexp

from .signal import signal_fires, _Z
from .mde import _DEFAULT_GRID


def _fires_mask(a: np.ndarray, n_drug: int, n_event: int, N: int,
                *, min_cases: int = 3) -> np.ndarray:
    """Vectorised twin of ``signal.signal_fires`` evaluated over a support array.

    Computes the pre-registered signal indicator (lower 95% Wald CI of ROR > 1,
    Haldane-Anscombe cc=0.5 applied only when a cell is zero, AND raw a >=
    min_cases) for every ``a`` in one numpy pass, replacing the per-element Python
    loop. The arithmetic mirrors ``signal_fires`` exactly; ``a`` is assumed to lie
    on the Fisher-NCHG support, where b, c, d are non-negative by construction, so
    the only admissibility test left is the min-cases floor. ``daen_mde_fast``'s
    equality to the scalar path over random margins is asserted in
    ``tests/test_power_fast.py``.
    """
    af = a.astype(np.float64)
    b = n_drug - af
    c = n_event - af
    d = (N - n_drug - n_event) + af
    anyzero = (np.minimum.reduce([af, b, c, d]) == 0.0)
    a2 = np.where(anyzero, af + 0.5, af)
    b2 = np.where(anyzero, b + 0.5, b)
    c2 = np.where(anyzero, c + 0.5, c)
    d2 = np.where(anyzero, d + 0.5, d)
    with np.errstate(divide="ignore", invalid="ignore"):
        log_ror = np.log((a2 * d2) / (b2 * c2))
        se = np.sqrt(1.0 / a2 + 1.0 / b2 + 1.0 / c2 + 1.0 / d2)
        ci_low = np.exp(log_ror - _Z * se)
    return (af >= min_cases) & np.isfinite(ci_low) & (ci_low > 1.0)


def _support_logweights(n_drug: int, n_event: int, N: int):
    """Return (a, g, fires) over the Fisher-NCHG support, or None if empty.

    ``g`` is the OR-free log-weight; ``fires`` is the boolean rejection mask.
    """
    a_lo = max(0, n_drug + n_event - N)
    a_hi = min(n_drug, n_event)
    if a_hi < a_lo:
        return None
    a = np.arange(a_lo, a_hi + 1, dtype=np.int64)
    g = -(
        gammaln(a + 1.0)
        + gammaln(n_drug - a + 1.0)
        + gammaln(n_event - a + 1.0)
        + gammaln(N - n_drug - n_event + a + 1.0)
    )
    fires = _fires_mask(a, n_drug, n_event, N)
    return a, g, fires


def exact_power_fast(n_drug: int, n_event: int, N: int, odds_ratio: float) -> float:
    """Single-OR exact power; same value as ``exact.exact_power`` to fp tolerance."""
    sup = _support_logweights(n_drug, n_event, N)
    if sup is None:
        return 0.0
    a, g, fires = sup
    if not fires.any():
        return 0.0
    logw = g + a * math.log(odds_ratio)
    return float(np.exp(logsumexp(logw[fires]) - logsumexp(logw)))


def daen_mde_fast(
    n_drug: int,
    n_event: int,
    N: int,
    *,
    power_target: float = 0.80,
    grid: np.ndarray = _DEFAULT_GRID,
) -> float:
    """Smallest OR in ``grid`` with exact power >= ``power_target`` (else inf).

    Drop-in replacement for ``power.mde.daen_mde``: computes the per-pair
    log-weights and rejection mask once, then powers the whole grid at once.
    """
    if n_drug <= 0 or n_event <= 0 or N <= 0:
        return math.inf
    sup = _support_logweights(n_drug, n_event, N)
    if sup is None:
        return math.inf
    a, g, fires = sup
    if not fires.any():
        return math.inf
    a_f = a.astype(np.float64)
    g_fires = g[fires]
    a_fires = a_f[fires]

    def power_at(idx: int) -> float:
        lo = math.log(grid[idx])
        den = logsumexp(g + a_f * lo)
        num = logsumexp(g_fires + a_fires * lo)
        return math.exp(num - den)

    # Power is monotone non-decreasing in OR (fixed upper-set rejection region,
    # MLR in the odds parameter), so binary-search the grid for the smallest OR
    # reaching the target — ~log2(G) probes instead of the full G-column matrix.
    hi = len(grid) - 1
    if power_at(hi) < power_target:
        return math.inf
    lo = 0
    while lo < hi:
        mid = (lo + hi) // 2
        if power_at(mid) >= power_target:
            hi = mid
        else:
            lo = mid + 1
    return float(grid[lo])
