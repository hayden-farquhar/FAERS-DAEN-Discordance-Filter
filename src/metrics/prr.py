"""Proportional Reporting Ratio (PRR) -- pre-registered sensitivity metric.

PRR = [a/(a+b)] / [c/(c+d)];  SE(ln PRR) = sqrt(1/a - 1/(a+b) + 1/c - 1/(c+d)).
Signalling rule (Evans et al. 2001): PRR >= 2 AND chi-square >= 4 AND a >= 3.
The chi-square is the Pearson 2x2 statistic on the RAW table (Yates' correction optional).
"""
from __future__ import annotations

import math
from statistics import NormalDist
from typing import NamedTuple


class PRRResult(NamedTuple):
    prr: float
    ci_low: float
    ci_high: float
    chi2: float
    a: float
    signal: bool


def _chi2_2x2(a: float, b: float, c: float, d: float, yates: bool = False) -> float:
    n = a + b + c + d
    if n == 0:
        return float("nan")
    num = abs(a * d - b * c)
    if yates:
        num = max(0.0, num - n / 2.0)
    denom = (a + b) * (c + d) * (a + c) * (b + d)
    if denom == 0:
        return 0.0
    return n * num * num / denom


def prr(
    a: float,
    b: float,
    c: float,
    d: float,
    alpha: float = 0.05,
    cc: float = 0.5,
    min_cases: int = 3,
    yates: bool = False,
) -> PRRResult:
    a_raw = a
    chi2 = _chi2_2x2(a, b, c, d, yates=yates)  # computed on the raw table

    if min(a, b, c, d) == 0:
        a, b, c, d = a + cc, b + cc, c + cc, d + cc

    r1 = a + b
    r2 = c + d
    point = (a / r1) / (c / r2)
    se = math.sqrt(1.0 / a - 1.0 / r1 + 1.0 / c - 1.0 / r2)
    z = NormalDist().inv_cdf(1.0 - alpha / 2.0)
    ln_point = math.log(point)
    ci_low = math.exp(ln_point - z * se)
    ci_high = math.exp(ln_point + z * se)

    signal = (point >= 2.0) and (chi2 >= 4.0) and (a_raw >= min_cases)
    return PRRResult(point, ci_low, ci_high, chi2, a_raw, signal)
