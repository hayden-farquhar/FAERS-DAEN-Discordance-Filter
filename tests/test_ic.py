"""IC (BCPNN) validated against its closed-form definition and known properties."""
import math
from statistics import NormalDist

import pytest

from metrics import ic, ic_from_table


def _ic_reference(a, n_drug, n_event, N, alpha=0.05):
    """Independent re-derivation of the documented IC + delta-method interval."""
    p11 = (a + 0.5) / (N + 1)
    p1 = (n_drug + 0.5) / (N + 1)
    p2 = (n_event + 0.5) / (N + 1)
    point = math.log(p11 / (p1 * p2)) / math.log(2)  # log2 via change of base (independent path)

    def vln(k):
        al = k + 0.5
        be = N - k + 0.5
        return be / (al * (al + be + 1))

    var = (vln(a) + vln(n_drug) + vln(n_event)) / (math.log(2) ** 2)
    sd = math.sqrt(var)
    z = NormalDist().inv_cdf(1 - alpha / 2)
    return point, point - z * sd, point + z * sd


def test_ic_matches_closed_form():
    a, nd, ne, N = 10, 100, 200, 10000
    point, lo, hi = _ic_reference(a, nd, ne, N)
    r = ic(a, nd, ne, N)
    assert r.ic == pytest.approx(point, rel=1e-12)
    assert r.ic025 == pytest.approx(lo, rel=1e-12)
    assert r.ic975 == pytest.approx(hi, rel=1e-12)


def test_ic_credibility_interval_ordering():
    r = ic(10, 100, 200, 10000)
    assert r.ic025 < r.ic < r.ic975


def test_ic_from_table_equivalent():
    # a=10, b=90, c=190, d=9710  -> n_drug=100, n_event=200, N=10000
    assert ic_from_table(10, 90, 190, 9710).ic == pytest.approx(ic(10, 100, 200, 10000).ic, rel=1e-12)


def test_ic_large_sample_approaches_log2_relative_reporting_ratio():
    a, nd, ne, N = 100, 200, 300, 100000
    rrr = a * N / (nd * ne)               # relative reporting ratio
    assert ic(a, nd, ne, N).ic == pytest.approx(math.log2(rrr), rel=0.02)


def test_ic_strong_signal_flagged_and_null_not():
    assert ic(100, 200, 300, 100000).signal is True       # IC025 > 0
    # no co-reports despite exposure -> IC025 < 0, no signal
    assert ic(0, 500, 500, 100000).signal is False
