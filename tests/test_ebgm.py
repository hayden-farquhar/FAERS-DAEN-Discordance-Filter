"""EBGM/MGPS validated against the gamma posterior closed form, shrinkage properties,
and a seeded hyperparameter-recovery simulation."""
import math

import numpy as np
import pytest
from scipy import special, stats

from metrics import MGPSHyperparams, ebgm_cell, mgps


def _single(alpha, beta):
    """A hyperparameter set whose two components are identical -> a single gamma prior."""
    return MGPSHyperparams(alpha, beta, alpha, beta, 0.5, loglik=float("nan"), converged=True)


def test_ebgm_single_component_matches_gamma_posterior():
    A, B, n, E = 2.0, 1.0, 5.0, 2.0
    r = ebgm_cell(n, E, _single(A, B))
    # identical components -> posterior is Gamma(A+n, rate B+E); geometric mean = exp(psi(a)-ln(b))
    assert r.q1 == pytest.approx(0.5, abs=1e-9)
    assert r.ebgm == pytest.approx(math.exp(special.digamma(A + n) - math.log(B + E)), rel=1e-9)
    assert r.eb05 == pytest.approx(stats.gamma.ppf(0.05, A + n, scale=1 / (B + E)), rel=1e-6)
    assert r.eb95 == pytest.approx(stats.gamma.ppf(0.95, A + n, scale=1 / (B + E)), rel=1e-6)


def test_ebgm_interval_ordering():
    r = ebgm_cell(10.0, 3.0, _single(2.0, 1.0))
    assert r.eb05 < r.ebgm < r.eb95


def test_ebgm_shrinkage_behaviour():
    # comp1: mean 1 (background); comp2: mean 5 (signal)
    h = MGPSHyperparams(2.0, 2.0, 5.0, 1.0, 0.85, loglik=float("nan"), converged=True)
    low = ebgm_cell(0.0, 10.0, h)    # no co-reports despite E=10 -> shrink well below 1
    high = ebgm_cell(100.0, 10.0, h)  # observed RR=10, large n -> EBGM approaches 10
    assert low.ebgm < 1.0
    assert 8.0 < high.ebgm < 10.5
    assert low.ebgm < high.ebgm
    assert high.signal is True        # EB05 > 2
    assert low.signal is False


def test_mgps_recovery_seeded():
    rng = np.random.default_rng(93)
    M = 800
    E = rng.gamma(2.0, 3.0, M) + 0.5
    is_bg = rng.random(M) < 0.85
    lam = np.where(
        is_bg,
        rng.gamma(2.0, 1.0 / 2.0, M),   # background: mean 1
        rng.gamma(5.0, 1.0 / 1.0, M),   # signal: mean 5
    )
    n = rng.poisson(lam * E)

    h, results = mgps(n, E, restarts=2, seed=93)
    assert math.isfinite(h.loglik)
    assert h.alpha1 / h.beta1 <= h.alpha2 / h.beta2     # component ordering enforced

    ebgm = np.array([r.ebgm for r in results])
    rho = stats.spearmanr(ebgm, lam).correlation
    assert rho > 0.5                                    # EBGM tracks the true relative risk
