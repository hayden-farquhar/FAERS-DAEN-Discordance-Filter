"""Reporting Odds Ratio (ROR) -- the pre-registered PRIMARY metric.

ROR = (a*d) / (b*c);  SE(ln ROR) = sqrt(1/a + 1/b + 1/c + 1/d);  Wald log-CI.
Signalling rule (pre-registration Section 9.3): lower 95% CI bound > 1 AND a >= 3.
"""
from __future__ import annotations

import math
from statistics import NormalDist
from typing import NamedTuple


class RORResult(NamedTuple):
    ror: float
    ci_low: float
    ci_high: float
    a: float
    signal: bool


def ror(
    a: float,
    b: float,
    c: float,
    d: float,
    alpha: float = 0.05,
    cc: float = 0.5,
    min_cases: int = 3,
) -> RORResult:
    """Reporting odds ratio with a two-sided (1-alpha) Wald confidence interval.

    A Haldane-Anscombe correction (add `cc` to all four cells) is applied only when a
    cell is zero, so non-degenerate tables are computed exactly. `signal` uses the raw
    `a` for the `a >= min_cases` test (the correction must not manufacture a signal).
    """
    a_raw = a
    if min(a, b, c, d) == 0:
        a, b, c, d = a + cc, b + cc, c + cc, d + cc

    point = (a * d) / (b * c)
    se = math.sqrt(1.0 / a + 1.0 / b + 1.0 / c + 1.0 / d)
    z = NormalDist().inv_cdf(1.0 - alpha / 2.0)
    ln_point = math.log(point)
    ci_low = math.exp(ln_point - z * se)
    ci_high = math.exp(ln_point + z * se)

    signal = (ci_low > 1.0) and (a_raw >= min_cases)
    return RORResult(point, ci_low, ci_high, a_raw, signal)
