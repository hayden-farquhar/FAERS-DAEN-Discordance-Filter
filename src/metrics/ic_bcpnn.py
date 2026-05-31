"""Information Component (IC), BCPNN family -- pre-registered sensitivity metric.

IC = log2( p11_hat / (p1._hat * p.1_hat) ), where each probability is the posterior
mean of a Beta with a Jeffreys (+0.5) prior:  p_hat = (count + 0.5) / (N + 1).

Credibility interval by the delta method: for Beta(al, be), Var[ln p] ~= be / (al*(al+be+1)),
and the three log-probabilities are treated as independent (the standard operational BCPNN
approximation; conservative for a sensitivity metric). Signalling rule: IC025 > 0.

References: Bate et al. (1998) Eur J Clin Pharmacol; Noren et al. (2006) Stat Med.
The exact prior (Jeffreys 0.5) and the delta-method interval are specified here so the
metric is fully reproducible.
"""
from __future__ import annotations

import math
from statistics import NormalDist
from typing import NamedTuple

_LN2 = math.log(2.0)


class ICResult(NamedTuple):
    ic: float
    ic025: float
    ic975: float
    signal: bool


def _var_ln_p(count: float, n_total: float) -> float:
    """Var[ln p] for p ~ Beta(count+0.5, n_total-count+0.5) via Var[ln X] ~= Var[X]/E[X]^2."""
    al = count + 0.5
    be = n_total - count + 0.5
    return be / (al * (al + be + 1.0))


def ic(a: float, n_drug: float, n_event: float, n_total: float, alpha: float = 0.05) -> ICResult:
    N = n_total
    p11 = (a + 0.5) / (N + 1.0)
    p1 = (n_drug + 0.5) / (N + 1.0)
    p2 = (n_event + 0.5) / (N + 1.0)
    point = math.log2(p11 / (p1 * p2))

    var_ic = (_var_ln_p(a, N) + _var_ln_p(n_drug, N) + _var_ln_p(n_event, N)) / (_LN2 ** 2)
    sd = math.sqrt(var_ic)
    z = NormalDist().inv_cdf(1.0 - alpha / 2.0)
    ic025 = point - z * sd
    ic975 = point + z * sd
    return ICResult(point, ic025, ic975, ic025 > 0.0)


def ic_from_table(a: float, b: float, c: float, d: float, alpha: float = 0.05) -> ICResult:
    """Convenience wrapper accepting the four 2x2 cells."""
    return ic(a, a + b, a + c, a + b + c + d, alpha=alpha)
