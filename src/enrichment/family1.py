"""Family 1 — primary descriptive headline (protocol Section 10.1).

H2: power-conditioned non-replication rate with cluster-bootstrap 95% CI
    (B=2,000; seed=94; clusters = events) within the Arm-1 universe restricted
    to FAERS-positive AND daen_powered.

Concordance reporting: Cohen's kappa with explicit marginal-imbalance caveat,
positive-specific agreement (PA+), and power-conditioned PA+ (subset to
mutually-powered pairs).

This module computes substrate-level descriptives only; no formal hypothesis
test. The output dict is the substrate for Section 15.5 Table 1.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

_SEED = 94
_B_BOOT = 2_000


@dataclass
class Family1Result:
    n_total_pairs: int
    n_faers_powered_universe: int
    n_daen_powered_universe: int
    n_faers_positive_daen_powered: int
    h2_non_replication_rate: float
    h2_ci_low: float
    h2_ci_high: float
    h2_underpowered_discordant_count: int
    pa_plus: float
    pa_plus_ci_low: float
    pa_plus_ci_high: float
    pa_plus_power_conditioned: float
    pa_plus_pc_ci_low: float
    pa_plus_pc_ci_high: float
    cohen_kappa: float
    kappa_caveat: str
    cross_db_2x2: dict[str, int] = field(default_factory=dict)
    notes: list[str] = field(default_factory=list)


def _cluster_bootstrap_proportion(
    df: pd.DataFrame, value_col: str, cluster_col: str,
    *, B: int = _B_BOOT, rng: np.random.Generator | None = None,
) -> tuple[float, float]:
    """Cluster-bootstrap 95% CI on the mean of a binary indicator column,
    resampling at the level of the cluster_col then computing the proportion
    over the assembled resample (all rows within each resampled cluster)."""
    if rng is None:
        rng = np.random.default_rng(_SEED)
    clusters = df[cluster_col].unique()
    if len(clusters) == 0 or len(df) == 0:
        return (float("nan"), float("nan"))
    grouped = {c: df.loc[df[cluster_col] == c, value_col].values for c in clusters}
    estimates = np.empty(B, dtype=float)
    for b in range(B):
        sampled_clusters = rng.choice(clusters, size=len(clusters), replace=True)
        vals = np.concatenate([grouped[c] for c in sampled_clusters])
        estimates[b] = float(vals.mean()) if vals.size else float("nan")
    estimates = estimates[~np.isnan(estimates)]
    if estimates.size == 0:
        return (float("nan"), float("nan"))
    return (float(np.percentile(estimates, 2.5)), float(np.percentile(estimates, 97.5)))


def _cohen_kappa(a: int, b: int, c: int, d: int) -> float:
    """2x2 Cohen's kappa where (a,b,c,d) lay out:
       rows = rater 1 (FAERS signal yes/no), cols = rater 2 (DAEN signal yes/no).
       a = both positive; b = FAERS+ only; c = DAEN+ only; d = both negative.
    """
    n = a + b + c + d
    if n == 0:
        return float("nan")
    po = (a + d) / n
    p_yes_1 = (a + b) / n
    p_yes_2 = (a + c) / n
    pe = p_yes_1 * p_yes_2 + (1 - p_yes_1) * (1 - p_yes_2)
    if 1 - pe == 0:
        return float("nan")
    return (po - pe) / (1 - pe)


def _pa_plus(a: int, b: int, c: int) -> float:
    """Positive-specific agreement = 2a / (2a + b + c)."""
    denom = 2 * a + b + c
    return (2 * a) / denom if denom else float("nan")


def family1_headline(perpair: pd.DataFrame) -> Family1Result:
    """Compute Family 1 descriptive headline on the augmented per-pair output."""
    # H2: among FAERS-positive AND daen_powered Arm-1 pairs, fraction faers_only
    h2_universe = perpair[(perpair["faers_signal"]) & (perpair["daen_powered"])].copy()
    h2_universe["is_faers_only"] = (h2_universe["cross_db_class"] == "faers_only").astype(int)
    n_h2 = len(h2_universe)
    if n_h2 == 0:
        h2_point, h2_lo, h2_hi = float("nan"), float("nan"), float("nan")
    else:
        h2_point = float(h2_universe["is_faers_only"].mean())
        rng = np.random.default_rng(_SEED)
        h2_lo, h2_hi = _cluster_bootstrap_proportion(
            h2_universe, "is_faers_only", "outcomeName", B=_B_BOOT, rng=rng,
        )

    # Underpowered-discordant count: FAERS-positive AND NOT daen_powered
    n_underpow_disc = int(
        ((perpair["faers_signal"]) & (~perpair["daen_powered"])).sum()
    )

    # Cross-DB 2x2 on the daen_powered universe (the headline conditioning)
    sub = perpair[perpair["daen_powered"]].copy()
    sub_a = int(((sub["faers_signal"]) & (sub["daen_signal"])).sum())
    sub_b = int(((sub["faers_signal"]) & (~sub["daen_signal"])).sum())
    sub_c = int(((~sub["faers_signal"]) & (sub["daen_signal"])).sum())
    sub_d = int(((~sub["faers_signal"]) & (~sub["daen_signal"])).sum())
    cross_db_2x2 = {
        "faers_pos_daen_pos": sub_a,
        "faers_pos_daen_neg": sub_b,
        "faers_neg_daen_pos": sub_c,
        "faers_neg_daen_neg": sub_d,
    }

    # PA+ on the daen_powered universe
    pa_plus = _pa_plus(sub_a, sub_b, sub_c)
    # Cluster-bootstrap CI on PA+ via resampling the powered universe at event-cluster level
    sub["pa_plus_pair"] = ((sub["faers_signal"]) & (sub["daen_signal"])).astype(int)
    # PA+ is not a simple mean of a binary col; do the bootstrap manually
    def _pa_plus_resample(rng: np.random.Generator) -> float:
        clusters = sub["outcomeName"].unique()
        if len(clusters) == 0:
            return float("nan")
        sampled_clusters = rng.choice(clusters, size=len(clusters), replace=True)
        rows = pd.concat([sub[sub["outcomeName"] == c] for c in sampled_clusters], ignore_index=True)
        a = int(((rows["faers_signal"]) & (rows["daen_signal"])).sum())
        b = int(((rows["faers_signal"]) & (~rows["daen_signal"])).sum())
        c2 = int(((~rows["faers_signal"]) & (rows["daen_signal"])).sum())
        return _pa_plus(a, b, c2)

    rng = np.random.default_rng(_SEED)
    pa_boot = np.array([_pa_plus_resample(rng) for _ in range(_B_BOOT)])
    pa_boot = pa_boot[~np.isnan(pa_boot)]
    pa_lo = float(np.percentile(pa_boot, 2.5)) if pa_boot.size else float("nan")
    pa_hi = float(np.percentile(pa_boot, 97.5)) if pa_boot.size else float("nan")

    # PA+ on the mutually-powered subset (faers_powered AND daen_powered)
    sub_mp = perpair[(perpair["daen_powered"]) & (perpair["faers_powered"])].copy()
    mp_a = int(((sub_mp["faers_signal"]) & (sub_mp["daen_signal"])).sum())
    mp_b = int(((sub_mp["faers_signal"]) & (~sub_mp["daen_signal"])).sum())
    mp_c = int(((~sub_mp["faers_signal"]) & (sub_mp["daen_signal"])).sum())
    pa_plus_pc = _pa_plus(mp_a, mp_b, mp_c)

    def _pa_plus_pc_resample(rng: np.random.Generator) -> float:
        clusters = sub_mp["outcomeName"].unique()
        if len(clusters) == 0:
            return float("nan")
        sampled_clusters = rng.choice(clusters, size=len(clusters), replace=True)
        rows = pd.concat([sub_mp[sub_mp["outcomeName"] == c] for c in sampled_clusters], ignore_index=True)
        a = int(((rows["faers_signal"]) & (rows["daen_signal"])).sum())
        b = int(((rows["faers_signal"]) & (~rows["daen_signal"])).sum())
        c2 = int(((~rows["faers_signal"]) & (rows["daen_signal"])).sum())
        return _pa_plus(a, b, c2)

    rng = np.random.default_rng(_SEED + 1)
    pa_pc_boot = np.array([_pa_plus_pc_resample(rng) for _ in range(_B_BOOT)])
    pa_pc_boot = pa_pc_boot[~np.isnan(pa_pc_boot)]
    pa_pc_lo = float(np.percentile(pa_pc_boot, 2.5)) if pa_pc_boot.size else float("nan")
    pa_pc_hi = float(np.percentile(pa_pc_boot, 97.5)) if pa_pc_boot.size else float("nan")

    # Cohen's kappa on the daen_powered universe, with caveat
    kappa = _cohen_kappa(sub_a, sub_b, sub_c, sub_d)
    n_pos_marg_1 = sub_a + sub_b
    n_pos_marg_2 = sub_a + sub_c
    if min(n_pos_marg_1, n_pos_marg_2) / max(len(sub), 1) < 0.10:
        kappa_caveat = (
            f"Cohen's kappa = {kappa:.3f}. SKEWED MARGINS: "
            f"FAERS-positive marginal = {n_pos_marg_1}/{len(sub)} "
            f"({100*n_pos_marg_1/max(len(sub),1):.1f}%), "
            f"DAEN-positive marginal = {n_pos_marg_2}/{len(sub)} "
            f"({100*n_pos_marg_2/max(len(sub),1):.1f}%). "
            f"kappa is well-known to be misleading under skewed marginals; "
            f"PA+ (= {pa_plus:.3f}) is the preferred concordance measure for this design."
        )
    else:
        kappa_caveat = (
            f"Cohen's kappa = {kappa:.3f}. Marginal balance is reasonable "
            f"({100*n_pos_marg_1/max(len(sub),1):.1f}% / "
            f"{100*n_pos_marg_2/max(len(sub),1):.1f}%); kappa interpretable but "
            f"PA+ (= {pa_plus:.3f}) remains the pre-registered primary concordance statistic."
        )

    return Family1Result(
        n_total_pairs=len(perpair),
        n_faers_powered_universe=int(perpair["faers_powered"].sum()),
        n_daen_powered_universe=int(perpair["daen_powered"].sum()),
        n_faers_positive_daen_powered=n_h2,
        h2_non_replication_rate=h2_point,
        h2_ci_low=h2_lo,
        h2_ci_high=h2_hi,
        h2_underpowered_discordant_count=n_underpow_disc,
        pa_plus=pa_plus,
        pa_plus_ci_low=pa_lo,
        pa_plus_ci_high=pa_hi,
        pa_plus_power_conditioned=pa_plus_pc,
        pa_plus_pc_ci_low=pa_pc_lo,
        pa_plus_pc_ci_high=pa_pc_hi,
        cohen_kappa=kappa,
        kappa_caveat=kappa_caveat,
        cross_db_2x2=cross_db_2x2,
        notes=[
            f"H2 universe (FAERS-positive AND daen_powered): n={n_h2}",
            f"Underpowered-discordant (FAERS-positive AND NOT daen_powered): n={n_underpow_disc}",
            f"daen_powered universe (for concordance + kappa): n={len(sub)}",
            f"Mutually-powered universe (daen_powered AND faers_powered): n={len(sub_mp)}",
        ],
    )
