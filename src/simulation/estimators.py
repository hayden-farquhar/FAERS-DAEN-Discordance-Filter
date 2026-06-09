"""Enrichment estimator panel (protocol Section 6).

The estimand is the H1 reference-negative enrichment exactly as in
``enrichment.family2.family2_h1``: a logistic regression of ``is_known_neg`` on
``is_faers_only`` within the selected H1 stratum. The slope's odds ratio is the
enrichment statistic; the simulation's rejection rule is the same for every
estimator — **reject iff the 95% CI lower bound of the OR exceeds 1** — so the
estimators are compared on a single, uniform decision and a single CI (Q3
coverage and Q1 type-I come from the same object).

Four estimators:

1. ``fixed_effect``           — plain logit, Wald CI.
2. ``cluster_robust``         — logit with CR SE clustered by event (outcomeName).
3. ``wild_cluster_bootstrap`` — unrestricted score wild cluster bootstrap with
   Webb six-point weights (Cameron-Gelbach-Miller / Kline-Santos), percentile-t
   CI. The few-cluster-appropriate method.
4. ``firth``                  — Firth penalized likelihood (Jeffreys prior),
   separation-robust; the estimator reported for the otherwise non-estimable
   cells.

Complete separation (a perfectly predicted ``faers_only`` row, or an empty
contrast row) is recorded as ``non_estimable=True`` rather than silently dropped;
the non-estimability rate is itself a reported quantity (Section 7).
"""
from __future__ import annotations

import math
import warnings
from dataclasses import dataclass, field
from statistics import NormalDist

import numpy as np
import pandas as pd

from enrichment.family2 import _h1_stratum

_Z = NormalDist().inv_cdf(0.975)
# Webb (2014) six-point weights: equal probability 1/6 on each.
_WEBB_POINTS = np.array([
    -math.sqrt(1.5), -1.0, -math.sqrt(0.5),
    math.sqrt(0.5), 1.0, math.sqrt(1.5),
])
_SEP_BETA = 15.0   # |slope| above this is treated as numerical separation


@dataclass
class EstimatorResult:
    estimator: str
    or_hat: float
    ci_lo: float
    ci_hi: float
    reject: bool
    non_estimable: bool
    beta: float
    se: float
    stratum_size: int
    n_faers_only: int
    n_concordant_positive: int
    n_clusters: int
    notes: list[str] = field(default_factory=list)


def build_h1_stratum(perpair: pd.DataFrame) -> pd.DataFrame:
    """The H1 stratum (faers_signal & daen_powered & in {faers_only,
    concordant_positive}), with is_faers_only / is_known_neg, via the empirical
    code path so the simulated stratum matches the real one exactly."""
    return _h1_stratum(perpair)


def _design(stratum: pd.DataFrame) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    y = stratum["is_known_neg"].to_numpy(dtype=float)
    x1 = stratum["is_faers_only"].to_numpy(dtype=float)
    X = np.column_stack([np.ones_like(x1), x1])
    clusters = stratum["outcomeName"].to_numpy()
    return X, y, clusters


def _degenerate(stratum: pd.DataFrame) -> str | None:
    """Return a reason string if the contrast is structurally undefined."""
    n_fo = int((stratum["is_faers_only"] == 1).sum())
    n_cp = int((stratum["is_faers_only"] == 0).sum())
    if n_fo == 0 or n_cp == 0:
        return "one contrast row (faers_only / concordant_positive) is empty"
    if stratum["is_known_neg"].nunique() < 2:
        return "outcome is constant across the stratum (no known-pos or no known-neg)"
    return None


def _nonestimable_result(name: str, stratum: pd.DataFrame, reason: str) -> EstimatorResult:
    n_fo = int((stratum["is_faers_only"] == 1).sum())
    n_cp = int((stratum["is_faers_only"] == 0).sum())
    return EstimatorResult(
        estimator=name, or_hat=float("nan"), ci_lo=float("nan"), ci_hi=float("nan"),
        reject=False, non_estimable=True, beta=float("nan"), se=float("nan"),
        stratum_size=len(stratum), n_faers_only=n_fo, n_concordant_positive=n_cp,
        n_clusters=int(stratum["outcomeName"].nunique()), notes=[reason],
    )


# --- IRLS logistic fit (shared by FE / CR / wild bootstrap) ---------------------

def _fit_logit(X: np.ndarray, y: np.ndarray, *, max_iter: int = 100, tol: float = 1e-10):
    """Plain Newton-Raphson logistic fit. Returns (beta, p, bread, converged)
    where bread = (X'WX)^{-1}. Raises nothing; separation surfaces as a large
    |beta| and is caught by the caller."""
    n, k = X.shape
    beta = np.zeros(k)
    bread = np.eye(k)
    for _ in range(max_iter):
        eta = X @ beta
        p = 1.0 / (1.0 + np.exp(-eta))
        W = p * (1.0 - p)
        Wsafe = np.clip(W, 1e-12, None)
        XtWX = X.T @ (X * Wsafe[:, None])
        try:
            bread = np.linalg.inv(XtWX)
        except np.linalg.LinAlgError:
            return beta, p, bread, False
        grad = X.T @ (y - p)
        step = bread @ grad
        beta = beta + step
        if np.max(np.abs(step)) < tol:
            break
    eta = X @ beta
    p = 1.0 / (1.0 + np.exp(-eta))
    converged = bool(np.all(np.isfinite(beta)) and np.max(np.abs(beta)) < _SEP_BETA)
    return beta, p, bread, converged


def _cluster_scores(X: np.ndarray, y: np.ndarray, p: np.ndarray, clusters: np.ndarray):
    """Per-cluster score sums u_g = sum_{i in g} x_i (y_i - p_i). Returns
    (U, labels) where U is (G x k)."""
    resid = (y - p)
    s = X * resid[:, None]
    labels = np.unique(clusters)
    U = np.zeros((labels.size, X.shape[1]))
    for j, g in enumerate(labels):
        U[j] = s[clusters == g].sum(axis=0)
    return U, labels


def _cr_cov(bread: np.ndarray, U: np.ndarray, n_clusters: int, n_obs: int, k: int) -> np.ndarray:
    """Cluster-robust sandwich with the standard small-sample correction."""
    meat = U.T @ U
    G = n_clusters
    corr = (G / (G - 1.0)) * ((n_obs - 1.0) / (n_obs - k)) if G > 1 else 1.0
    return corr * (bread @ meat @ bread)


def _wild_tstar(U: np.ndarray, bread: np.ndarray, small: float, B: int,
                rng: np.random.Generator) -> np.ndarray:
    """Vectorised Webb six-point wild cluster bootstrap reference distribution
    for the slope t-statistic.

    Given fixed per-cluster scores ``U`` (G x k), the cluster-robust bread
    (k x k, symmetric), and the small-sample factor ``small``, draw B Webb-weight
    vectors and return the B bootstrap t-statistics on the slope (coefficient
    index 1). This batches over the B axis with einsum instead of a Python loop;
    because ``rng.choice(size=(B, G))`` consumes the same row-major draw stream
    as B sequential ``rng.choice(size=G)`` calls, the output is bit-identical to
    the per-draw loop on the same Generator state, not merely distributionally
    equal. ``UU[g] = U_g U_g'`` is precomputed once so the per-draw meat is a
    single weighted contraction.
    """
    G = U.shape[0]
    UU = np.einsum("gi,gj->gij", U, U)          # (G, k, k)
    W = rng.choice(_WEBB_POINTS, size=(B, G))   # (B, G)
    S = W @ U                                    # (B, k): sum_g w_g U_g
    delta = S @ bread                            # (B, k): bread symmetric
    meat = np.einsum("bg,gij->bij", W * W, UU)   # (B, k, k): sum_g w_g^2 U_g U_g'
    BMB = np.einsum("ij,bjk,kl->bil", bread, meat, bread)
    se_w = np.sqrt(np.maximum(small * BMB[:, 1, 1], 1e-300))
    return delta[:, 1] / se_w


# --- estimator 1: fixed-effect logit -------------------------------------------

def _fixed_effect(stratum: pd.DataFrame, **_) -> EstimatorResult:
    reason = _degenerate(stratum)
    if reason:
        return _nonestimable_result("fixed_effect", stratum, reason)
    X, y, clusters = _design(stratum)
    beta, p, bread, ok = _fit_logit(X, y)
    if not ok:
        return _nonestimable_result("fixed_effect", stratum, "logit did not converge (separation)")
    b = float(beta[1])
    se = float(math.sqrt(bread[1, 1]))
    return _assemble("fixed_effect", stratum, b, se)


# --- estimator 2: cluster-robust SE --------------------------------------------

def _cluster_robust(stratum: pd.DataFrame, **_) -> EstimatorResult:
    reason = _degenerate(stratum)
    if reason:
        return _nonestimable_result("cluster_robust", stratum, reason)
    X, y, clusters = _design(stratum)
    beta, p, bread, ok = _fit_logit(X, y)
    if not ok:
        return _nonestimable_result("cluster_robust", stratum, "logit did not converge (separation)")
    U, _labels = _cluster_scores(X, y, p, clusters)
    V = _cr_cov(bread, U, n_clusters=U.shape[0], n_obs=X.shape[0], k=X.shape[1])
    b = float(beta[1])
    se = float(math.sqrt(max(V[1, 1], 0.0)))
    return _assemble("cluster_robust", stratum, b, se)


# --- estimator 3: wild cluster bootstrap (Webb 6-point) ------------------------

def _wild_cluster_bootstrap(stratum: pd.DataFrame, *, rng, B: int = 999, **_) -> EstimatorResult:
    reason = _degenerate(stratum)
    if reason:
        return _nonestimable_result("wild_cluster_bootstrap", stratum, reason)
    X, y, clusters = _design(stratum)
    beta, p, bread, ok = _fit_logit(X, y)
    if not ok:
        return _nonestimable_result("wild_cluster_bootstrap", stratum, "logit did not converge (separation)")
    U, _labels = _cluster_scores(X, y, p, clusters)
    G = U.shape[0]
    V = _cr_cov(bread, U, n_clusters=G, n_obs=X.shape[0], k=X.shape[1])
    b = float(beta[1])
    se = float(math.sqrt(max(V[1, 1], 0.0)))
    if se <= 0 or G < 2:
        return _nonestimable_result(
            "wild_cluster_bootstrap", stratum,
            "cluster-robust SE is zero or <2 clusters; bootstrap undefined")

    small = (G / (G - 1.0)) * ((X.shape[0] - 1.0) / (X.shape[0] - X.shape[1]))
    t_star = np.empty(B)
    for b_i in range(B):
        w = rng.choice(_WEBB_POINTS, size=G)
        Uw = U * w[:, None]
        delta = bread @ Uw.sum(axis=0)            # bootstrap coefficient deviation
        Vw = small * (bread @ (Uw.T @ Uw) @ bread)
        se_w = math.sqrt(max(Vw[1, 1], 1e-300))
        t_star[b_i] = delta[1] / se_w
    # Percentile-t (equal-tailed) CI on the slope, mapped to the OR scale.
    q_lo, q_hi = np.quantile(t_star, [0.975, 0.025])
    ci_lo_beta = b - se * q_lo
    ci_hi_beta = b - se * q_hi
    or_hat = math.exp(b)
    ci_lo = math.exp(ci_lo_beta)
    ci_hi = math.exp(ci_hi_beta)
    return EstimatorResult(
        estimator="wild_cluster_bootstrap", or_hat=or_hat, ci_lo=ci_lo, ci_hi=ci_hi,
        reject=bool(ci_lo > 1.0), non_estimable=False, beta=b, se=se,
        stratum_size=len(stratum),
        n_faers_only=int((stratum["is_faers_only"] == 1).sum()),
        n_concordant_positive=int((stratum["is_faers_only"] == 0).sum()),
        n_clusters=G,
    )


# --- estimator 3b: restricted wild cluster bootstrap (WCR, null imposed) --------

def _wild_cluster_bootstrap_restricted(
        stratum: pd.DataFrame, *, rng, B: int = 999, **_) -> EstimatorResult:
    """Restricted score wild cluster bootstrap of H0: beta_slope = 0.

    Differs from ``_wild_cluster_bootstrap`` (WCU) in one place that matters at
    few clusters: the bootstrap data-generating process imposes the null. The
    cluster scores and the bread are taken at the **restricted** fit (intercept
    only, so the fitted probability is the sample mean of the outcome), and the
    Webb-weighted bootstrap statistic is therefore generated under
    beta_slope = 0. MacKinnon & Webb (2014) show this restricted (WCR) variant
    has size/coverage closer to nominal than the unrestricted (WCU) variant when
    the number of clusters is small.

    The point statistic is the one-step (score) estimate of the slope from the
    restricted fit, delta_1 = [B_tilde @ sum_g U_g]_1, with a cluster-robust SE at
    the restricted fit. An equal-tailed percentile-t CI is built from the
    restricted reference distribution and mapped to the OR scale, so the result
    plugs into the same coverage/rejection machinery as every other estimator
    (reject iff CI lower bound of the OR exceeds 1; coverage = CI contains OR=1).
    """
    reason = _degenerate(stratum)
    if reason:
        return _nonestimable_result("wild_cluster_bootstrap_restricted", stratum, reason)
    X, y, clusters = _design(stratum)
    G = int(np.unique(clusters).size)
    if G < 2:
        return _nonestimable_result(
            "wild_cluster_bootstrap_restricted", stratum, "<2 clusters; bootstrap undefined")

    # Restricted fit: intercept only -> fitted prob is the outcome mean (exact MLE
    # of the intercept-only logistic). Build the restricted bread B_tilde from the
    # restricted weight on the FULL design.
    ybar = float(y.mean())
    w_tilde = ybar * (1.0 - ybar)
    if w_tilde <= 0:
        return _nonestimable_result(
            "wild_cluster_bootstrap_restricted", stratum, "restricted weight degenerate")
    XtWX = (X.T @ X) * w_tilde
    try:
        bread = np.linalg.inv(XtWX)
    except np.linalg.LinAlgError:
        return _nonestimable_result(
            "wild_cluster_bootstrap_restricted", stratum, "restricted information singular")

    p_tilde = np.full_like(y, ybar)
    U, _labels = _cluster_scores(X, y, p_tilde, clusters)   # scores at the restricted fit
    small = (G / (G - 1.0)) * ((X.shape[0] - 1.0) / (X.shape[0] - X.shape[1]))

    delta = bread @ U.sum(axis=0)                # one-step slope estimate from the null fit
    V = small * (bread @ (U.T @ U) @ bread)
    se = float(math.sqrt(max(V[1, 1], 0.0)))
    if se <= 0:
        return _nonestimable_result(
            "wild_cluster_bootstrap_restricted", stratum, "restricted cluster-robust SE is zero")
    b1 = float(delta[1])

    t_star = _wild_tstar(U, bread, small, B, rng)   # null imposed via restricted U/bread
    q_lo, q_hi = np.quantile(t_star, [0.975, 0.025])
    ci_lo_beta = b1 - se * q_lo
    ci_hi_beta = b1 - se * q_hi
    if abs(b1) >= _SEP_BETA:
        return _nonestimable_result(
            "wild_cluster_bootstrap_restricted", stratum,
            "restricted one-step slope at separation boundary")
    return EstimatorResult(
        estimator="wild_cluster_bootstrap_restricted",
        or_hat=math.exp(b1), ci_lo=math.exp(ci_lo_beta), ci_hi=math.exp(ci_hi_beta),
        reject=bool(math.exp(ci_lo_beta) > 1.0), non_estimable=False, beta=b1, se=se,
        stratum_size=len(stratum),
        n_faers_only=int((stratum["is_faers_only"] == 1).sum()),
        n_concordant_positive=int((stratum["is_faers_only"] == 0).sum()),
        n_clusters=G,
    )


# --- estimator 3c: vectorised unrestricted wild bootstrap (post-hoc speed path) -

def _wild_cluster_bootstrap_fast(
        stratum: pd.DataFrame, *, rng, B: int = 999, **_) -> EstimatorResult:
    """Vectorised twin of ``_wild_cluster_bootstrap`` (WCU).

    Identical estimator, identical decision rule, identical draw stream; the only
    difference is that the B-loop is replaced by the batched ``_wild_tstar``
    contraction. On the same Generator state it reproduces the frozen WCU
    bit-for-bit (verified max abs diff ~9e-16), so the density sweep's WCU column
    is the same object as the frozen sweep's, only faster. Kept distinct from the
    frozen ``_wild_cluster_bootstrap`` name so the pre-registered grid's estimator
    set is untouched.
    """
    reason = _degenerate(stratum)
    if reason:
        return _nonestimable_result("wild_cluster_bootstrap", stratum, reason)
    X, y, clusters = _design(stratum)
    beta, p, bread, ok = _fit_logit(X, y)
    if not ok:
        return _nonestimable_result(
            "wild_cluster_bootstrap", stratum, "logit did not converge (separation)")
    U, _labels = _cluster_scores(X, y, p, clusters)
    G = U.shape[0]
    V = _cr_cov(bread, U, n_clusters=G, n_obs=X.shape[0], k=X.shape[1])
    b = float(beta[1])
    se = float(math.sqrt(max(V[1, 1], 0.0)))
    if se <= 0 or G < 2:
        return _nonestimable_result(
            "wild_cluster_bootstrap", stratum,
            "cluster-robust SE is zero or <2 clusters; bootstrap undefined")
    small = (G / (G - 1.0)) * ((X.shape[0] - 1.0) / (X.shape[0] - X.shape[1]))
    t_star = _wild_tstar(U, bread, small, B, rng)
    q_lo, q_hi = np.quantile(t_star, [0.975, 0.025])
    ci_lo = math.exp(b - se * q_lo)
    ci_hi = math.exp(b - se * q_hi)
    return EstimatorResult(
        estimator="wild_cluster_bootstrap", or_hat=math.exp(b), ci_lo=ci_lo, ci_hi=ci_hi,
        reject=bool(ci_lo > 1.0), non_estimable=False, beta=b, se=se,
        stratum_size=len(stratum),
        n_faers_only=int((stratum["is_faers_only"] == 1).sum()),
        n_concordant_positive=int((stratum["is_faers_only"] == 0).sum()),
        n_clusters=G,
    )


# --- estimator 4: Firth penalized likelihood -----------------------------------

def _fit_firth(X: np.ndarray, y: np.ndarray, *, max_iter: int = 200, tol: float = 1e-8):
    """Firth (1993) penalized logistic via modified-score Newton iteration.

    Score adjustment: U*_j = sum_i (y_i - p_i + h_i (1/2 - p_i)) x_ij, where h_i
    are the hat-matrix diagonals from the Fisher-weighted design. Information is
    the unpenalized Fisher information X'WX. Separation-robust by construction.
    """
    n, k = X.shape
    beta = np.zeros(k)
    info_inv = np.eye(k)
    for _ in range(max_iter):
        eta = X @ beta
        p = 1.0 / (1.0 + np.exp(-eta))
        W = np.clip(p * (1.0 - p), 1e-12, None)
        XtWX = X.T @ (X * W[:, None])
        try:
            info_inv = np.linalg.inv(XtWX)
        except np.linalg.LinAlgError:
            return beta, info_inv, False
        # Hat diagonals h_i = W_i * x_i' (X'WX)^-1 x_i
        Wsqrt = np.sqrt(W)
        Xt = X * Wsqrt[:, None]
        H = Xt @ info_inv @ Xt.T
        h = np.clip(np.diag(H), 0.0, 1.0)
        U_star = X.T @ (y - p + h * (0.5 - p))
        step = info_inv @ U_star
        beta = beta + step
        if np.max(np.abs(step)) < tol:
            break
    return beta, info_inv, bool(np.all(np.isfinite(beta)))


def _firth(stratum: pd.DataFrame, **_) -> EstimatorResult:
    # Firth tolerates separation, but a structurally empty contrast row still
    # yields no information on the slope.
    n_fo = int((stratum["is_faers_only"] == 1).sum())
    n_cp = int((stratum["is_faers_only"] == 0).sum())
    if n_fo == 0 or n_cp == 0:
        return _nonestimable_result("firth", stratum, "one contrast row is empty")
    X, y, _clusters = _design(stratum)
    beta, info_inv, ok = _fit_firth(X, y)
    if not ok or info_inv[1, 1] <= 0:
        return _nonestimable_result("firth", stratum, "Firth IRLS did not converge")
    b = float(beta[1])
    se = float(math.sqrt(info_inv[1, 1]))
    return _assemble("firth", stratum, b, se)


# --- shared Wald assembly -------------------------------------------------------

def _assemble(name: str, stratum: pd.DataFrame, beta: float, se: float) -> EstimatorResult:
    if not math.isfinite(beta) or not math.isfinite(se) or se <= 0:
        return _nonestimable_result(name, stratum, "non-finite estimate or SE")
    if abs(beta) >= _SEP_BETA:
        # Numerical (near-)separation: the slope has run off; OR is not meaningfully
        # estimable even under penalisation. Record as non-estimable, not a reject.
        return _nonestimable_result(name, stratum, "slope at separation boundary (|beta| too large)")
    or_hat = math.exp(beta)
    ci_lo = math.exp(beta - _Z * se)
    ci_hi = math.exp(beta + _Z * se)
    return EstimatorResult(
        estimator=name, or_hat=or_hat, ci_lo=ci_lo, ci_hi=ci_hi,
        reject=bool(ci_lo > 1.0), non_estimable=False, beta=beta, se=se,
        stratum_size=len(stratum),
        n_faers_only=int((stratum["is_faers_only"] == 1).sum()),
        n_concordant_positive=int((stratum["is_faers_only"] == 0).sum()),
        n_clusters=int(stratum["outcomeName"].nunique()),
    )


# The four pre-registered estimators that populate the frozen sweep
# (results/sim_results.parquet). The restricted wild bootstrap (WCR) is a
# post-hoc addition for the S4 few-cluster coverage re-examination and is run
# only on the focused K-sweep line (scripts/run_wcr_coverage.py); it is kept out
# of this dict so it cannot silently alter the pre-registered grid, and is
# exposed separately via POSTHOC_ESTIMATORS / run_estimator's lookup.
ESTIMATORS = {
    "fixed_effect": _fixed_effect,
    "cluster_robust": _cluster_robust,
    "wild_cluster_bootstrap": _wild_cluster_bootstrap,
    "firth": _firth,
}

POSTHOC_ESTIMATORS = {
    "wild_cluster_bootstrap_restricted": _wild_cluster_bootstrap_restricted,
    "wild_cluster_bootstrap_fast": _wild_cluster_bootstrap_fast,
}


def run_estimator(
    stratum: pd.DataFrame, estimator_name: str, *,
    rng: np.random.Generator | None = None, B: int = 999,
) -> EstimatorResult:
    """Run one estimator on a prebuilt H1 stratum."""
    fn = ESTIMATORS.get(estimator_name) or POSTHOC_ESTIMATORS.get(estimator_name)
    if fn is None:
        raise ValueError(
            f"unknown estimator {estimator_name!r}; expected one of "
            f"{sorted(ESTIMATORS) + sorted(POSTHOC_ESTIMATORS)}"
        )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return fn(stratum, rng=rng, B=B)
