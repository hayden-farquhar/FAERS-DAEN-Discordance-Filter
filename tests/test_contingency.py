import math

import pytest

from metrics import Contingency


def test_marginals_and_expected():
    t = Contingency(a=25, b=1000, c=50, d=10000)
    assert t.n_drug == 1025
    assert t.n_event == 75
    assert t.n_total == 11075
    # expected = n_drug * n_event / N
    assert t.expected == pytest.approx(1025 * 75 / 11075)


def test_from_marginals_roundtrip():
    t = Contingency.from_marginals(a=25, n_drug=1025, n_event=75, n_total=11075)
    assert t.cells() == (25, 1000, 50, 10000)


def test_negative_cell_rejected():
    with pytest.raises(ValueError):
        Contingency(a=-1, b=1, c=1, d=1)
