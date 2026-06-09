"""Data-generating process for the enrichment type-I / power simulation.

Protocol Section 3. A population of drug-event pairs is generated, each carrying
a known true label and known true cross-database behaviour, then passed through
the *same* signal rule (``power.signal.signal_fires``) used in the empirical
pipeline. The output is a per-pair frame with the columns the empirical
enrichment code (``enrichment.family2``) consumes, plus the DAEN/FAERS marginals
the selection rules need.

Two regimes separate the two simulation questions (Section 3.3):

- ``Q1_null``  — discordance is uninformative: FAERS-only status arises *only*
  because DAEN is underpowered at that pair's marginals, a mechanism independent
  of the label. True negative-enrichment among ``faers_only`` is zero by
  construction, so any H1 rejection is a false positive. This is the primary
  (type-I) regime.
- ``Q2_alt``   — discordance is informative: artefactual nulls are given a
  genuine extra propensity to be FAERS-only (their DAEN association is suppressed
  below 1), so true negative-enrichment is positive by construction. This is the
  power regime.

The calibrated marginal-distribution parameters are NOT fixed here; they are read
from ``sim_config.yaml`` and remain TODO placeholders (structural defaults) until
G-SIM check 1 (empirical-anchor sanity) is closed. Nothing in this file hard-codes
a value claimed to be calibrated to the empirical histograms.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy.stats import nchypergeom_fisher

from power.signal import signal_fires

# Columns the downstream code requires (kept explicit so a schema drift is loud).
PERPAIR_COLUMNS = (
    "exposureName", "outcomeName", "groundTruth",
    "faers_signal", "daen_signal", "cross_db_class",
    "n_drug_F", "n_event_F", "N_F", "a_F", "OR_F",
    "n_drug_D", "n_event_D", "N_D", "a_D", "OR_D",
    "ror_F", "ror_F_ci_lo", "true_label", "is_artefact",
)


@dataclass(frozen=True)
class MarginalParams:
    """Structural marginal-draw parameters (uncalibrated until G-SIM check 1).

    Drug- and event-side FAERS marginals are drawn lognormal on the count scale;
    DAEN marginals are the FAERS marginals scaled by the size ratio ``rho`` and
    Poisson-resampled to inject independent small-sample noise. These defaults let
    the simulator run for the G-SIM null-of-the-null and estimator unit tests; the
    calibrated values overwrite them at lock (see ``sim_config.yaml``).
    """
    log_mu_drug: float
    log_sd_drug: float
    log_mu_event: float
    log_sd_event: float
    calibrated: bool = False  # flipped to True only when fitted to empirical data


@dataclass(frozen=True)
class DGPCell:
    """One point in the DGP grid (protocol Section 4)."""
    regime: str                 # "Q1_null" | "Q2_alt"
    pi: float                   # true-null fraction
    lambda_inflate: float       # FAERS inflation factor on the artefactual subset
    phi: float                  # artefact prevalence among nulls
    rho: float                  # DAEN/FAERS size ratio
    K: int                      # number of event clusters
    icc: float                  # within-event ICC on the label logit
    OR_true: float | str        # fixed value, or "lognormal" for the continuous arm
    N_F: int                    # FAERS-scale total report count
    n_pairs: int                # pairs generated per replicate
    marginals: MarginalParams
    q2_daen_suppression: float = 0.5   # DAEN OR multiplier for informative nulls (Q2)
    or_lognormal_mu: float = math.log(2.0)
    or_lognormal_sigma: float = 0.5
    label: str = field(default="")     # free-text tag for the sweep manifest


def _sigma_u_from_icc(icc: float) -> float:
    """Random-intercept SD giving the requested latent-logistic ICC.

    Latent-variable ICC for a logistic model with a Normal(0, sigma_u^2) intercept
    is sigma_u^2 / (sigma_u^2 + pi^2/3). Invert for sigma_u.
    """
    if icc <= 0:
        return 0.0
    if icc >= 1:
        raise ValueError("icc must be < 1")
    resid = math.pi ** 2 / 3.0
    return math.sqrt(icc * resid / (1.0 - icc))


def _draw_or_true(cell: DGPCell, n: int, rng: np.random.Generator) -> np.ndarray:
    if cell.OR_true == "lognormal":
        return rng.lognormal(mean=cell.or_lognormal_mu, sigma=cell.or_lognormal_sigma, size=n)
    return np.full(n, float(cell.OR_true))


def _draw_a_cell(
    n_drug: np.ndarray, n_event: np.ndarray, N: int, OR: np.ndarray,
    rng: np.random.Generator,
) -> np.ndarray:
    """Vectorised Fisher noncentral hypergeometric a-cell draw (one per pair).

    scipy's nchypergeom_fisher takes shape args positionally (M, n, N, odds_ratio);
    they broadcast, so a single rvs call draws the whole population.
    """
    M = np.full(n_drug.shape, int(N), dtype=np.int64)
    return nchypergeom_fisher.rvs(M, n_drug, n_event, OR, random_state=rng)


def generate_population(cell: DGPCell, rng: np.random.Generator) -> pd.DataFrame:
    """Generate one replicate population for ``cell`` and return a per-pair frame."""
    if cell.regime not in ("Q1_null", "Q2_alt"):
        raise ValueError(f"unknown regime {cell.regime!r}")
    n = cell.n_pairs

    # --- event clusters + within-event random intercept on the label logit ---
    cluster = rng.integers(0, cell.K, size=n)
    sigma_u = _sigma_u_from_icc(cell.icc)
    u = rng.normal(0.0, sigma_u, size=cell.K) if sigma_u > 0 else np.zeros(cell.K)
    beta0 = math.log((1.0 - cell.pi) / cell.pi)         # P(real_positive) = 1 - pi at u=0
    logit = beta0 + u[cluster]
    p_pos = 1.0 / (1.0 + np.exp(-logit))
    is_real_positive = rng.random(n) < p_pos
    true_label = np.where(is_real_positive, "real_positive", "true_null")
    ground_truth = is_real_positive.astype(int)         # 1 = positive, 0 = known-negative

    # --- true effect sizes ---------------------------------------------------
    or_pos = _draw_or_true(cell, n, rng)
    OR_true = np.where(is_real_positive, or_pos, 1.0)

    # --- artefact subset: a fraction phi of NULLS receives FAERS inflation ---
    is_null = ~is_real_positive
    is_artefact = is_null & (rng.random(n) < cell.phi)
    OR_F = np.where(is_artefact, OR_true * cell.lambda_inflate, OR_true)
    if cell.regime == "Q1_null":
        # Null regime (protocol Section 3.3): discordance must be uninformative, so
        # DAEN detection is label-blind. Setting OR_D = 1 for EVERY pair (positives
        # included) makes DAEN replication depend only on marginals/power, not on L
        # -> faers_only is independent of the label -> true enrichment is zero by
        # construction. Real positives carry their FAERS effect (OR_F = OR_true) but
        # leave no DAEN footprint here; their realistic DAEN replication belongs to
        # the Q2 regime. Without this, positives replicate into concordant_positive
        # and mechanically enrich faers_only with negatives even at lambda=1, phi=0.
        OR_D = np.ones(n)
    else:
        # Q2_alt — informative discordance: positives carry their real DAEN effect
        # (OR_D = OR_true, so they replicate), while artefactual nulls have a
        # genuinely suppressed DAEN association (below power-alone expectation) ->
        # real negative enrichment among faers_only by construction.
        OR_D = OR_true.copy()
        OR_D = np.where(is_artefact, OR_true * cell.q2_daen_suppression, OR_D)

    # --- marginals (FAERS drawn lognormal; DAEN scaled by rho + resampled) ---
    mp = cell.marginals
    n_drug_F = np.maximum(1, np.round(
        rng.lognormal(mp.log_mu_drug, mp.log_sd_drug, size=n)).astype(np.int64))
    n_event_F = np.maximum(1, np.round(
        rng.lognormal(mp.log_mu_event, mp.log_sd_event, size=n)).astype(np.int64))
    N_F = int(cell.N_F)
    N_D = max(1, int(round(N_F * cell.rho)))
    n_drug_D = np.maximum(1, rng.poisson(np.maximum(n_drug_F * cell.rho, 1e-6))).astype(np.int64)
    n_event_D = np.maximum(1, rng.poisson(np.maximum(n_event_F * cell.rho, 1e-6))).astype(np.int64)
    # Keep marginals feasible against their database total.
    n_drug_F = np.minimum(n_drug_F, N_F - 1)
    n_event_F = np.minimum(n_event_F, N_F - 1)
    n_drug_D = np.minimum(n_drug_D, N_D - 1)
    n_event_D = np.minimum(n_event_D, N_D - 1)

    # --- observed a-cells via Fisher noncentral hypergeometric ---------------
    a_F = _draw_a_cell(n_drug_F, n_event_F, N_F, OR_F, rng)
    a_D = _draw_a_cell(n_drug_D, n_event_D, N_D, OR_D, rng)

    # --- registered signal rule (unchanged) ----------------------------------
    faers_signal = np.fromiter(
        (signal_fires(int(a_F[i]), int(n_drug_F[i]), int(n_event_F[i]), N_F) for i in range(n)),
        dtype=bool, count=n)
    daen_signal = np.fromiter(
        (signal_fires(int(a_D[i]), int(n_drug_D[i]), int(n_event_D[i]), N_D) for i in range(n)),
        dtype=bool, count=n)

    # --- FAERS point + lower-bound ROR (drives the FAERS-derived selection) --
    ror_F, ror_F_lo = _faers_ror(a_F, n_drug_F, n_event_F, N_F)

    cross_db_class = np.select(
        [
            (faers_signal & daen_signal),
            (faers_signal & ~daen_signal),
            (~faers_signal & daen_signal),
        ],
        ["concordant_positive", "faers_only", "daen_only"],
        default="concordant_negative",
    )

    return pd.DataFrame({
        "exposureName": [f"P{i:06d}" for i in range(n)],
        "outcomeName": [f"E{int(k):03d}" for k in cluster],
        "groundTruth": ground_truth,
        "faers_signal": faers_signal,
        "daen_signal": daen_signal,
        "cross_db_class": cross_db_class,
        "n_drug_F": n_drug_F, "n_event_F": n_event_F, "N_F": N_F, "a_F": a_F, "OR_F": OR_F,
        "n_drug_D": n_drug_D, "n_event_D": n_event_D, "N_D": N_D, "a_D": a_D, "OR_D": OR_D,
        "ror_F": ror_F, "ror_F_ci_lo": ror_F_lo,
        "true_label": true_label, "is_artefact": is_artefact,
    })


def _faers_ror(a, n_drug, n_event, N):
    """FAERS ROR point estimate and 95% Wald lower bound per pair (Haldane cc on
    zero cells only), matching the signal-rule convention. Used by the
    FAERS-derived selection rules as the OR they plug into the DAEN power model.
    """
    a = a.astype(float)
    b = n_drug - a
    c = n_event - a
    d = N - n_drug - n_event + a
    z = 1.959963984540054
    a2, b2, c2, d2 = a.copy(), b.copy(), c.copy(), d.copy()
    zero = (np.minimum.reduce([a, b, c, d]) == 0)
    a2[zero] += 0.5; b2[zero] += 0.5; c2[zero] += 0.5; d2[zero] += 0.5
    with np.errstate(divide="ignore", invalid="ignore"):
        log_ror = np.log((a2 * d2) / (b2 * c2))
        se = np.sqrt(1.0 / a2 + 1.0 / b2 + 1.0 / c2 + 1.0 / d2)
        ror = np.exp(log_ror)
        ror_lo = np.exp(log_ror - z * se)
    bad = ~np.isfinite(ror)
    ror[bad] = 1.0
    ror_lo[~np.isfinite(ror_lo)] = 1.0
    return ror, ror_lo
