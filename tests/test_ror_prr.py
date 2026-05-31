"""ROR and PRR validated against independently hand-computed worked-example values."""
import math
from statistics import NormalDist

import pytest

from metrics import ror, prr

# Worked example 2x2: a strong, non-degenerate signal table.
A, B, C, D = 25, 1000, 50, 10000
Z = NormalDist().inv_cdf(0.975)


def test_ror_point_is_exact():
    # ROR = ad/bc = (25*10000)/(1000*50) = 5.0 exactly
    assert ror(A, B, C, D).ror == pytest.approx(5.0, rel=1e-12)


def test_ror_ci_matches_wald_formula():
    se = math.sqrt(1 / A + 1 / B + 1 / C + 1 / D)
    lo = math.exp(math.log(5.0) - Z * se)
    hi = math.exp(math.log(5.0) + Z * se)
    r = ror(A, B, C, D)
    assert r.ci_low == pytest.approx(lo, rel=1e-12)
    assert r.ci_high == pytest.approx(hi, rel=1e-12)


def test_ror_signal_rule():
    # lower CI ~3.08 > 1 and a=25 >= 3  -> signal
    assert ror(A, B, C, D).signal is True
    # same CI but a below the minimum-case threshold -> no signal
    assert ror(2, 10, 50, 100000, min_cases=3).signal is False


def test_ror_zero_cell_correction_runs_and_uses_raw_a():
    r = ror(0, 100, 5, 100000)
    assert math.isfinite(r.ror) and math.isfinite(r.ci_low)
    assert r.a == 0 and r.signal is False  # raw a=0 fails a>=3


def test_prr_point_and_chi2_match_formula():
    prr_point = (A / (A + B)) / (C / (C + D))
    N = A + B + C + D
    chi2 = N * (A * D - B * C) ** 2 / ((A + B) * (C + D) * (A + C) * (B + D))
    p = prr(A, B, C, D)
    assert p.prr == pytest.approx(prr_point, rel=1e-12)
    assert p.chi2 == pytest.approx(chi2, rel=1e-12)


def test_prr_ci_matches_wald_formula():
    se = math.sqrt(1 / A - 1 / (A + B) + 1 / C - 1 / (C + D))
    point = (A / (A + B)) / (C / (C + D))
    lo = math.exp(math.log(point) - Z * se)
    hi = math.exp(math.log(point) + Z * se)
    p = prr(A, B, C, D)
    assert p.ci_low == pytest.approx(lo, rel=1e-12)
    assert p.ci_high == pytest.approx(hi, rel=1e-12)


def test_prr_evans_signal_rule():
    # PRR ~4.9 >= 2, chi2 large >= 4, a >= 3 -> signal
    assert prr(A, B, C, D).signal is True
    # PRR < 2 -> no signal even with many cases
    assert prr(30, 1000, 1000, 30000).signal is False


def test_ror_prr_converge_for_rare_event():
    # When the event is rare, ROR and PRR coincide.
    a, b, c, d = 8, 200000, 4, 400000
    assert ror(a, b, c, d).ror == pytest.approx(prr(a, b, c, d).prr, rel=1e-2)
