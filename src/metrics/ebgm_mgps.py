"""EBGM via the Multi-item Gamma-Poisson Shrinker (MGPS) -- sensitivity metric.

DuMouchel (1999, Am Stat). Model for each drug-event cell:
    n | mu ~ Poisson(mu),  mu = lambda * E,  E = expected count under independence,
    lambda ~ P * Gamma(a1, rate b1) + (1-P) * Gamma(a2, rate b2).
Marginally, n | E follows a 2-component negative-binomial mixture; the five
hyperparameters (a1, b1, a2, b2, P) are fit by maximising the marginal log-likelihood
across ALL cells (hence the two-stage API). For each cell the posterior of lambda is a
2-component gamma mixture, from which:
    EBGM = exp(E[ln lambda | n])         (the empirical-Bayes geometric mean)
    EB05, EB95 = 5th / 95th posterior percentiles of lambda.
Signalling rule (pre-registration Section 9.3): EB05 > 2.
"""
from __future__ import annotations

import math
from typing import NamedTuple, Sequence

import numpy as np
from scipy import optimize, special, stats

_DEFAULT_SEED = 93
_DEFAULT_INIT = (0.2, 0.1, 2.0, 4.0, 1.0 / 3.0)  # DuMouchel-style starting point


class MGPSHyperparams(NamedTuple):
    alpha1: float
    beta1: float
    alpha2: float
    beta2: float
    p: float
    loglik: float
    converged: bool


class EBGMResult(NamedTuple):
    ebgm: float
    eb05: float
    eb95: float
    n: float
    expected: float
    q1: float          # posterior weight on component 1 (the lower-mean / background component)
    signal: bool


def _sigmoid(x: float) -> float:
    return 1.0 / (1.0 + math.exp(-x))


def _logit(p: float) -> float:
    return math.log(p / (1.0 - p))


def _log_nb(n, alpha, beta, E):
    """log P(n) for n | E ~ NegBinom(r=alpha, prob=beta/(beta+E)). Vectorised over n, E."""
    n = np.asarray(n, dtype=float)
    E = np.asarray(E, dtype=float)
    log_p = np.log(beta) - np.log(beta + E)          # log of success prob
    log_1mp = np.log(E) - np.log(beta + E)           # log(1 - success prob) = log(E/(beta+E))
    return (
        special.gammaln(n + alpha)
        - special.gammaln(alpha)
        - special.gammaln(n + 1.0)
        + alpha * log_p
        + n * log_1mp
    )


def _neg_loglik(params: np.ndarray, n: np.ndarray, E: np.ndarray) -> float:
    a1, b1, a2, b2 = np.exp(params[:4])
    p = _sigmoid(params[4])
    f1 = _log_nb(n, a1, b1, E)
    f2 = _log_nb(n, a2, b2, E)
    mix = np.logaddexp(math.log(p) + f1, math.log1p(-p) + f2)
    val = -np.sum(mix)
    return val if np.isfinite(val) else 1e18


def fit_hyperparameters(
    n: Sequence[float],
    E: Sequence[float],
    init: tuple[float, float, float, float, float] = _DEFAULT_INIT,
    restarts: int = 4,
    seed: int = _DEFAULT_SEED,
) -> MGPSHyperparams:
    """Maximum-marginal-likelihood fit of the 5 MGPS hyperparameters across all cells.

    Parameters are optimised in an unconstrained space (log for the gamma a/b, logit for P)
    with Nelder-Mead plus a few seeded random restarts to reduce local-minimum risk.
    Components are returned ordered so component 1 has the smaller prior mean (a/b).
    """
    n = np.asarray(n, dtype=float)
    E = np.asarray(E, dtype=float)
    rng = np.random.default_rng(seed)

    x0 = np.array(
        [math.log(init[0]), math.log(init[1]), math.log(init[2]), math.log(init[3]), _logit(init[4])]
    )
    starts = [x0] + [x0 + rng.normal(scale=0.75, size=5) for _ in range(restarts)]

    best = None
    for s in starts:
        res = optimize.minimize(
            _neg_loglik, s, args=(n, E), method="Nelder-Mead",
            options={"maxiter": 20000, "maxfev": 20000, "xatol": 1e-7, "fatol": 1e-7},
        )
        if best is None or res.fun < best.fun:
            best = res

    a1, b1, a2, b2 = np.exp(best.x[:4])
    p = _sigmoid(best.x[4])
    if a1 / b1 > a2 / b2:  # enforce component-1 = lower-mean component
        a1, b1, a2, b2, p = a2, b2, a1, b1, 1.0 - p
    return MGPSHyperparams(float(a1), float(b1), float(a2), float(b2), float(p), float(-best.fun), bool(best.success))


def _posterior_weight_comp1(n: float, E: float, h: MGPSHyperparams) -> float:
    lf1 = float(_log_nb(n, h.alpha1, h.beta1, E))
    lf2 = float(_log_nb(n, h.alpha2, h.beta2, E))
    la = math.log(h.p) + lf1
    lb = math.log1p(-h.p) + lf2
    return math.exp(la - np.logaddexp(la, lb))


def _mixture_cdf(x: float, n: float, E: float, h: MGPSHyperparams, q1: float) -> float:
    c1 = stats.gamma.cdf(x, h.alpha1 + n, scale=1.0 / (h.beta1 + E))
    c2 = stats.gamma.cdf(x, h.alpha2 + n, scale=1.0 / (h.beta2 + E))
    return q1 * c1 + (1.0 - q1) * c2


def _posterior_quantile(q: float, n: float, E: float, h: MGPSHyperparams, q1: float) -> float:
    f = lambda x: _mixture_cdf(x, n, E, h, q1) - q
    hi = 10.0
    while f(hi) < 0.0 and hi < 1e12:
        hi *= 10.0
    return float(optimize.brentq(f, 1e-12, hi, xtol=1e-10, maxiter=500))


def ebgm_cell(n: float, E: float, h: MGPSHyperparams, min_eb05_signal: float = 2.0) -> EBGMResult:
    """Posterior EBGM / EB05 / EB95 for a single cell given fitted hyperparameters `h`."""
    q1 = _posterior_weight_comp1(n, E, h)
    e_ln = (
        q1 * (special.digamma(h.alpha1 + n) - math.log(h.beta1 + E))
        + (1.0 - q1) * (special.digamma(h.alpha2 + n) - math.log(h.beta2 + E))
    )
    ebgm = math.exp(e_ln)
    eb05 = _posterior_quantile(0.05, n, E, h, q1)
    eb95 = _posterior_quantile(0.95, n, E, h, q1)
    return EBGMResult(ebgm, eb05, eb95, float(n), float(E), q1, eb05 > min_eb05_signal)


def mgps(
    n: Sequence[float],
    E: Sequence[float],
    init: tuple[float, float, float, float, float] = _DEFAULT_INIT,
    restarts: int = 4,
    seed: int = _DEFAULT_SEED,
) -> tuple[MGPSHyperparams, list[EBGMResult]]:
    """Fit the prior on all cells, then return per-cell EBGM results."""
    n = np.asarray(n, dtype=float)
    E = np.asarray(E, dtype=float)
    h = fit_hyperparameters(n, E, init=init, restarts=restarts, seed=seed)
    return h, [ebgm_cell(float(ni), float(Ei), h) for ni, Ei in zip(n, E)]
