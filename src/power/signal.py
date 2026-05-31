"""Deterministic signal-rule evaluator for a single 2x2 cell.

The pre-registered signal rule (protocol Section 9.3):
    Signal iff (a) lower bound of 95% Wald CI of ROR > 1 with Haldane-Anscombe
    correction (cc=0.5) applied only when at least one cell is zero,
    AND (b) raw `a` cell count >= 3.

This is a deterministic function of (a, n_drug, n_event, N). Given the margins,
it tells us at which integer values of `a` the signal rule fires. The power
computations then aggregate this indicator under the relevant distribution.
"""
from __future__ import annotations

import math
from statistics import NormalDist

_Z = NormalDist().inv_cdf(0.975)


def signal_fires(a: int, n_drug: int, n_event: int, N: int, *, min_cases: int = 3) -> bool:
    """Return True iff the 2x2 (a, n_drug-a, n_event-a, N-n_drug-n_event+a) fires
    the pre-registered signal rule.

    Mirrors the contract of `metrics.ror.ror(...)` for the LB > 1 check, with
    Haldane-Anscombe correction (cc=0.5) applied only when at least one cell is zero.
    The min-cases test uses the RAW a (the protocol's anti-Haldane-manufactured-signal
    safeguard).
    """
    if a < 0 or a > min(n_drug, n_event):
        return False
    b = n_drug - a
    c = n_event - a
    d = N - n_drug - n_event + a
    if min(b, c, d) < 0:
        return False
    if a < min_cases:
        return False
    a2, b2, c2, d2 = a, b, c, d
    if min(a, b, c, d) == 0:
        a2, b2, c2, d2 = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    try:
        log_ror = math.log((a2 * d2) / (b2 * c2))
        se = math.sqrt(1.0 / a2 + 1.0 / b2 + 1.0 / c2 + 1.0 / d2)
    except (ValueError, ZeroDivisionError):
        return False
    ci_low = math.exp(log_ror - _Z * se)
    return ci_low > 1.0
