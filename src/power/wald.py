"""Wald-analytic power (protocol Section 6.1, step 3 'Analytic' path).

Approximates P(signal | margins, true OR) by:
  1. Computing the expected cell counts under the alternative (a_alt such that
     a_alt * d_alt / (b_alt * c_alt) = OR_alt, with margins held fixed).
  2. Treating log(ROR_hat) as approximately N(log(OR_alt), SE_alt^2) where
     SE_alt^2 = 1/a_alt + 1/b_alt + 1/c_alt + 1/d_alt (Haldane-corrected when
     any expected cell is below a small floor).
  3. P(LB of 95% Wald CI > 1) = 1 - Phi(z) where z = (1.96 - log(OR_alt) / SE_alt).
  4. Correcting for the discrete a >= 3 threshold by multiplying by the
     Fisher-noncentral-hypergeometric tail mass P(a >= 3 | OR_alt).

This is "Wald-based power on the log-ROR scale, accounting for the
minimum-case-count threshold via a binomial mass below a = 3" per the protocol.
The "binomial mass" language is a slight imprecision in the locked protocol
text; the correct distribution is the Fisher noncentral hypergeometric (because
DAEN's two margins are conditioned on), and this implementation uses that.
"""
from __future__ import annotations

import math
from statistics import NormalDist

from scipy.stats import nchypergeom_fisher

_Z = NormalDist().inv_cdf(0.975)
_NORMAL = NormalDist()
_CC = 0.5  # Haldane-Anscombe continuity correction


def _expected_a_under_OR(n_drug: int, n_event: int, N: int, OR_alt: float) -> float:
    """Solve a/(n_drug-a) * (N-n_drug-n_event+a)/(n_event-a) = OR for a, with
    a in (max(0, n_drug+n_event-N), min(n_drug, n_event)).

    Use the noncentral-hypergeometric mean as the expected value under the
    alternative (this is exactly what 'expected a under OR_alt' means for fixed
    margins).
    """
    return float(nchypergeom_fisher.mean(N, n_drug, n_event, OR_alt))


def wald_power(
    n_drug: int,
    n_event: int,
    N: int,
    OR_alt: float,
    *,
    min_cases: int = 3,
) -> float:
    """Wald-analytic estimate of P(signal | margins, true OR_alt)."""
    if n_drug <= 0 or n_event <= 0 or N <= 0:
        return 0.0
    if n_drug + n_event > N + min(n_drug, n_event):
        return 0.0

    # 1. Expected cell counts under the alternative.
    a_alt = _expected_a_under_OR(n_drug, n_event, N, OR_alt)
    if not math.isfinite(a_alt) or a_alt <= 0:
        # tiny-effect or pathological margins: collapse to the null tail mass
        a_alt = max(a_alt, 1e-9)
    b_alt = n_drug - a_alt
    c_alt = n_event - a_alt
    d_alt = N - n_drug - n_event + a_alt

    # Haldane on expected cells if any falls below 0.5 (defensive).
    cells = [a_alt, b_alt, c_alt, d_alt]
    if min(cells) < 0.5:
        a_alt, b_alt, c_alt, d_alt = (x + _CC for x in cells)

    # 2. SE of log(ROR) under the alternative.
    se = math.sqrt(1.0 / a_alt + 1.0 / b_alt + 1.0 / c_alt + 1.0 / d_alt)
    if se <= 0 or not math.isfinite(se):
        return 0.0

    # 3. P(LB > 1) = P(log(ROR_hat) - Z*SE > 0) = P(log(ROR_hat)/SE > Z).
    # Under N(log(OR_alt), SE^2), log(ROR_hat)/SE ~ N(log(OR_alt)/SE, 1).
    mu = math.log(OR_alt)
    z_lb = _Z - mu / se
    p_lb = 1.0 - _NORMAL.cdf(z_lb)

    # 4. Correct for the discrete a >= min_cases threshold.
    if min_cases <= 0:
        p_min = 1.0
    else:
        p_min = float(nchypergeom_fisher.sf(min_cases - 1, N, n_drug, n_event, OR_alt))

    # The two events (LB>1 and a>=min_cases) are positively correlated in the
    # FNCHG; the product is an upper-bound-ish proxy that the G2 validation
    # checks against MC truth. The protocol acknowledges this and authorises
    # falling back to MC where the approximation breaks.
    return p_lb * p_min
