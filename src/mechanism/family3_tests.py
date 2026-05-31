"""Family 3 confirmatory tests on the Arm-2 daen_powered universe.

For each H5a-d: difference-of-proportions test (faers_only vs concordant_positive)
within the Arm-2 daen_powered universe; Newcombe 95% CI; Holm-Bonferroni at
family-wise α=0.05, k=4.

H5a: proportion with consumer_lift_pp >= 10 pp
H5b: proportion with lawyer_lift_pp >= 5 pp
H5c: proportion with aligned_alert_F == True (among pairs with >= 16 quarters of data)
H5d: proportion with mass_tort_drug == True

Decision rule per Section 4: reject Holm-adjusted p < 0.05 within Family 3 AND
difference-of-proportions >= corresponding threshold (10 pp / 5 pp for H5a/b;
positive lift for H5c/d).
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from statistics import NormalDist

import numpy as np
import pandas as pd

_Z = NormalDist().inv_cdf(0.975)
_FAMILY_ALPHA = 0.05


@dataclass
class H5Result:
    name: str
    description: str
    n_faers_only: int
    n_concordant_positive: int
    p_faers_only: float       # proportion exhibiting mechanism in faers_only
    p_concordant_positive: float
    diff_pp: float             # (p_faers_only - p_concordant_positive) * 100
    ci_low_pp: float           # Newcombe 95% CI lower (percentage points)
    ci_high_pp: float          # Newcombe 95% CI upper
    pvalue: float              # two-sided z-test on difference of proportions
    notes: list[str] = field(default_factory=list)


def _newcombe_diff_ci(p1: float, n1: int, p2: float, n2: int) -> tuple[float, float]:
    """Newcombe-Wilson 95% CI for the difference p1 - p2."""
    if n1 == 0 or n2 == 0:
        return (float("nan"), float("nan"))
    z = _Z
    # Wilson interval on each
    def _wilson(p, n):
        denom = 1 + z * z / n
        centre = (p + z * z / (2 * n)) / denom
        half = (z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n))) / denom
        return centre - half, centre + half
    l1, u1 = _wilson(p1, n1)
    l2, u2 = _wilson(p2, n2)
    # Newcombe method 10
    delta = (p1 - p2)
    lo = delta - math.sqrt((p1 - l1) ** 2 + (u2 - p2) ** 2)
    hi = delta + math.sqrt((u1 - p1) ** 2 + (p2 - l2) ** 2)
    return lo, hi


def _diff_test(p1: float, n1: int, p2: float, n2: int) -> float:
    """Two-sided z-test on difference of proportions (pooled SE)."""
    if n1 == 0 or n2 == 0:
        return float("nan")
    pooled = (p1 * n1 + p2 * n2) / (n1 + n2)
    se = math.sqrt(pooled * (1 - pooled) * (1.0 / n1 + 1.0 / n2))
    if se == 0:
        return float("nan")
    z = (p1 - p2) / se
    # two-sided p
    return 2.0 * (1.0 - NormalDist().cdf(abs(z)))


def _compute_one(stratum: pd.DataFrame, mech_col: str, name: str, desc: str,
                  *, threshold_kind: str = "boolean",
                  threshold_value: float = 0.0) -> H5Result:
    """Compute one H5 test.

    threshold_kind:
      'boolean'      -> proportion with mech_col == True (no value threshold)
      'lift_ge'      -> proportion with mech_col >= threshold_value
    """
    fo = stratum[stratum["cross_db_class"] == "faers_only"]
    cp = stratum[stratum["cross_db_class"] == "concordant_positive"]
    n1, n2 = len(fo), len(cp)
    if threshold_kind == "boolean":
        x1 = int(fo[mech_col].fillna(False).astype(bool).sum())
        x2 = int(cp[mech_col].fillna(False).astype(bool).sum())
    else:
        x1 = int((fo[mech_col].fillna(-1) >= threshold_value).sum())
        x2 = int((cp[mech_col].fillna(-1) >= threshold_value).sum())
    p1 = x1 / n1 if n1 else float("nan")
    p2 = x2 / n2 if n2 else float("nan")
    diff_pp = (p1 - p2) * 100
    ci_lo_pp, ci_hi_pp = _newcombe_diff_ci(p1, n1, p2, n2)
    ci_lo_pp *= 100; ci_hi_pp *= 100
    pval = _diff_test(p1, n1, p2, n2)
    return H5Result(
        name=name, description=desc,
        n_faers_only=n1, n_concordant_positive=n2,
        p_faers_only=p1, p_concordant_positive=p2,
        diff_pp=diff_pp, ci_low_pp=ci_lo_pp, ci_high_pp=ci_hi_pp,
        pvalue=pval,
    )


def family3_holm(arm2_features: pd.DataFrame) -> tuple[list[H5Result], dict[str, str]]:
    """Run H5a-d on the Arm-2 daen_powered universe; apply Holm at α=0.05/k=4."""
    stratum = arm2_features[arm2_features["daen_powered"]].copy()
    # Per Section 4: H5c restricts to pairs with >= 16 quarters of data
    h5c_stratum = stratum[stratum["n_quarters_with_data"] >= 16]

    results = [
        _compute_one(stratum, "consumer_lift_pp",
                       "H5a", "Consumer-share lift ≥ 10 pp vs ingredient baseline",
                       threshold_kind="lift_ge", threshold_value=10.0),
        _compute_one(stratum, "lawyer_lift_pp",
                       "H5b", "Lawyer-share lift ≥ 5 pp vs ingredient baseline",
                       threshold_kind="lift_ge", threshold_value=5.0),
        _compute_one(h5c_stratum, "aligned_alert_F",
                       "H5c", "Bai-Perron break within ±2 quarters of FDA alert (≥16 quarters required)",
                       threshold_kind="boolean"),
        _compute_one(stratum, "mass_tort_drug",
                       "H5d", "Drug on US JPML mass-tort MDL registry",
                       threshold_kind="boolean"),
    ]

    # Holm-Bonferroni at family-wise α=0.05, k=4
    ranked = sorted(enumerate(results), key=lambda x: (x[1].pvalue if not math.isnan(x[1].pvalue) else 999))
    holm_decision: dict[str, str] = {}
    k = 4
    rejected_so_far = True
    for rank, (orig_idx, r) in enumerate(ranked):
        threshold = _FAMILY_ALPHA / (k - rank)
        if rejected_so_far and not math.isnan(r.pvalue) and r.pvalue < threshold:
            holm_decision[r.name] = "REJECT_AT_HOLM"
        else:
            holm_decision[r.name] = "FAIL_TO_REJECT"
            rejected_so_far = False
    return results, holm_decision
