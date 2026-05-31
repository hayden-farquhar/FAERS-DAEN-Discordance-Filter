"""Phase 6 sensitivity arm — alternative signal definition: IC025 > 0 AND a >= 3.

Pre-registered Section 10.7 robustness exhibit (NOT a confirmatory test; not
multiplicity-counted). Re-runs Family 1 (H2 + concordance) and Family 2 (H1, H3,
within-event MH) on the Arm-1 universe under the BCPNN Information-Component
signal rule, holding everything else fixed.

Design decisions (documented for the manuscript):
  * The variation in Section 10.7 is the *signal definition* only. The power
    model (`daen_powered = daen_mde <= 1.5`) is the G2-locked ROR-anchored MDE
    and is held fixed; we vary only which pairs the signal rule flags as
    positive. This isolates the effect of the disproportionality metric on
    cross-database classification from any change in the powered universe.
  * IC025 is computed by the deposited BCPNN module (Jeffreys 0.5 prior, delta-
    method interval) directly from the per-pair 2x2 cells already in the
    substrate parquet. No cell counts are recomputed.
  * Signal rule matches the primary's structure: IC025 > 0 AND raw a >= 3.

Outputs: results/phase6_ic025_sensitivity.{json,md}
"""
from __future__ import annotations

import json
import sys
from dataclasses import asdict
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from metrics.ic_bcpnn import ic_from_table  # noqa: E402
from enrichment.family1 import family1_headline  # noqa: E402
from enrichment.family2 import (  # noqa: E402
    family2_h1,
    family2_h3,
    h1_mantel_haenszel_within_event,
)

RESULTS = ROOT / "results"
MIN_CASES = 3


def _ic025_signal(a: float, b: float, c: float, d: float) -> tuple[float, bool]:
    """Return (ic025, signal) under the rule IC025 > 0 AND raw a >= 3."""
    res = ic_from_table(a, b, c, d)
    return res.ic025, bool((res.ic025 > 0.0) and (a >= MIN_CASES))


def build_ic025_frame(perpair: pd.DataFrame) -> pd.DataFrame:
    """Return a copy of the substrate with faers_signal / daen_signal /
    cross_db_class recomputed under the IC025 rule. All power columns,
    groundTruth, and identifiers are preserved unchanged."""
    df = perpair.copy()

    ic025_F, sig_F, ic025_D, sig_D = [], [], [], []
    for _, r in df.iterrows():
        f_ic, f_sig = _ic025_signal(r["a_F"], r["b_F"], r["c_F"], r["d_F"])
        d_ic, d_sig = _ic025_signal(r["a_D"], r["b_D"], r["c_D"], r["d_D"])
        ic025_F.append(f_ic); sig_F.append(f_sig)
        ic025_D.append(d_ic); sig_D.append(d_sig)

    df["ic025_F"] = ic025_F
    df["ic025_D"] = ic025_D
    # Preserve the ROR-based signals for the flip diagnostic, then overwrite.
    df["faers_signal_ror"] = df["faers_signal"]
    df["daen_signal_ror"] = df["daen_signal"]
    df["faers_signal"] = sig_F
    df["daen_signal"] = sig_D
    df["cross_db_class"] = np.select(
        [
            (df["faers_signal"] & df["daen_signal"]),
            (df["faers_signal"] & ~df["daen_signal"]),
            (~df["faers_signal"] & df["daen_signal"]),
        ],
        ["concordant_positive", "faers_only", "daen_only"],
        default="concordant_negative",
    )
    return df


def main() -> None:
    perpair = pd.read_parquet(RESULTS / "perpair_arm1.parquet")
    ic = build_ic025_frame(perpair)

    # --- flip diagnostic (how much does the metric change classification?) ---
    f_flip = int((ic["faers_signal"] != ic["faers_signal_ror"]).sum())
    d_flip = int((ic["daen_signal"] != ic["daen_signal_ror"]).sum())
    powered = ic[ic["daen_powered"]]
    class_primary = perpair.set_index(["exposureName", "outcomeName", "ref_set"])["cross_db_class"]
    class_ic = ic.set_index(["exposureName", "outcomeName", "ref_set"])["cross_db_class"]
    aligned = class_primary.to_frame("primary").join(class_ic.to_frame("ic025"))
    n_class_flip = int((aligned["primary"] != aligned["ic025"]).sum())

    # --- re-run Family 1 + Family 2 under IC025 ---
    f1 = family1_headline(ic)
    h1 = family2_h1(ic)
    h3 = family2_h3(ic)
    mh = h1_mantel_haenszel_within_event(ic)

    out = {
        "signal_rule": "IC025 > 0 AND a >= 3 (BCPNN, Jeffreys 0.5 prior, delta-method CI)",
        "power_model": "held fixed at primary daen_powered = (daen_mde <= 1.5); ROR-anchored, G2-locked",
        "flip_diagnostic": {
            "faers_signal_flips_vs_ror": f_flip,
            "daen_signal_flips_vs_ror": d_flip,
            "n_pairs_total": int(len(ic)),
            "cross_db_class_flips_vs_primary": n_class_flip,
            "n_daen_powered": int(ic["daen_powered"].sum()),
        },
        "family1": asdict(f1),
        "h1": asdict(h1),
        "h3": asdict(h3),
        "mantel_haenszel": asdict(mh),
    }
    (RESULTS / "phase6_ic025_sensitivity.json").write_text(json.dumps(out, indent=2))

    # --- markdown report ---
    md = []
    md.append("# Phase 6 sensitivity — alternative signal definition (IC025)\n")
    md.append("**Pre-registered Section 10.7 robustness exhibit.** Reported for H2, H1, H3. "
              "Not a confirmatory test; not multiplicity-counted (shares data with the primary).\n")
    md.append("**Signal rule:** `IC025 > 0 AND a >= 3` (BCPNN, Jeffreys 0.5 prior, delta-method "
              "interval) vs the primary `ROR 95% LB > 1 AND a >= 3`.\n")
    md.append("**Power model:** held fixed at the primary `daen_powered = (daen_mde <= 1.5)` "
              "(ROR-anchored, G2-locked). The Section 10.7 variation is the signal definition only.\n")

    md.append("\n## Classification flip diagnostic\n")
    md.append(f"- FAERS signal-status flips vs ROR rule: **{f_flip} / {len(ic)}** "
              f"({100*f_flip/len(ic):.1f}%)")
    md.append(f"- DAEN signal-status flips vs ROR rule: **{d_flip} / {len(ic)}** "
              f"({100*d_flip/len(ic):.1f}%)")
    md.append(f"- `cross_db_class` flips vs primary: **{n_class_flip} / {len(ic)}** "
              f"({100*n_class_flip/len(ic):.1f}%)\n")

    md.append("\n## H2 — power-conditioned non-replication rate (Family 1)\n")
    md.append(f"- H2 (IC025) = **{f1.h2_non_replication_rate:.3f}** "
              f"(95% cluster-bootstrap CI {f1.h2_ci_low:.3f}–{f1.h2_ci_high:.3f}) "
              f"on n={f1.n_faers_positive_daen_powered} FAERS-positive ∧ daen_powered pairs")
    md.append(f"- Underpowered-discordant count reported alongside: {f1.h2_underpowered_discordant_count}")
    md.append(f"- PA+ = **{f1.pa_plus:.3f}** (95% CI {f1.pa_plus_ci_low:.3f}–{f1.pa_plus_ci_high:.3f})")
    md.append(f"- Cohen's κ = {f1.cohen_kappa:.3f}")
    md.append(f"- {f1.kappa_caveat}")
    md.append(f"- Cross-DB 2×2 on daen_powered universe: {f1.cross_db_2x2}\n")

    md.append("\n## H1 — reference-negative enrichment (Family 2)\n")
    md.append(f"- Decision: **{h1.decision}**")
    md.append(f"- Stratum: {h1.stratum_size} pairs ({h1.n_faers_only} faers_only, "
              f"{h1.n_concordant_positive} concordant_positive) across {h1.n_events_in_stratum} events, "
              f"{h1.n_drugs_in_stratum} drugs")
    md.append(f"- Known-negatives: {h1.n_known_neg_faers_only} of {h1.n_faers_only} faers_only; "
              f"{h1.n_known_neg_concordant_pos} of {h1.n_concordant_positive} concordant_positive")
    md.append(f"- Fixed-effect OR = {h1.or_fixed_effect:.3f} "
              f"(95% CI {h1.or_fe_ci_low:.3f}–{h1.or_fe_ci_high:.3f}; p = {h1.or_fe_pvalue:.4f})")
    md.append(f"- Cluster-robust OR = {h1.or_cluster_robust:.3f} "
              f"(95% CI {h1.or_cr_ci_low:.3f}–{h1.or_cr_ci_high:.3f}; p = {h1.or_cr_pvalue:.4f})")
    md.append(f"- Realised MDE = {h1.realised_mde_or:.2f}")
    for n in h1.notes:
        md.append(f"  - {n}")

    md.append("\n## H3 — discordance-filter LR+ (Family 2)\n")
    md.append(f"- Decision: **{h3.decision}**")
    md.append(f"- sensitivity = {h3.sensitivity}, specificity = {h3.specificity}, "
              f"PPV = {h3.ppv}, NPV = {h3.npv}")
    md.append(f"- LR+ = {h3.lr_plus} (95% CI {h3.lr_plus_ci_low}–{h3.lr_plus_ci_high})")
    for n in h3.notes:
        md.append(f"  - {n}")

    md.append("\n## Within-event Mantel–Haenszel (Section 10.5)\n")
    md.append(f"- Strata total: {mh.n_strata_total}; informative: {mh.n_strata_informative}")
    md.append(f"- MH OR = {mh.mh_or} (95% CI {mh.mh_ci_low}–{mh.mh_ci_high})")
    for n in mh.notes:
        md.append(f"  - {n}")

    (RESULTS / "phase6_ic025_sensitivity.md").write_text("\n".join(md) + "\n")
    print("Wrote results/phase6_ic025_sensitivity.{json,md}")
    print(f"H2(IC025) = {f1.h2_non_replication_rate:.3f} "
          f"[{f1.h2_ci_low:.3f}, {f1.h2_ci_high:.3f}]; H1 decision = {h1.decision}; "
          f"H3 decision = {h3.decision}")
    print(f"class flips vs primary: {n_class_flip}/{len(ic)}")


if __name__ == "__main__":
    main()
