"""Phase 4 runner: Family 1 descriptive headline + Family 2 confirmatory.

Produces:
  - results/phase4_results.json   — structured numbers
  - results/phase4_results.md     — human-readable report

Once G2 is closed (Phase 2), G3 lifts and these headline numbers are
authoritative per the protocol. Per Section 16.3, no post-hoc analyses
are included; this script only computes the pre-registered Family 1 + 2
endpoints.
"""
from __future__ import annotations

import json
import sys
import time
from dataclasses import asdict
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from enrichment import (  # noqa: E402
    family1_headline, family2_h1, family2_h3, h1_mantel_haenszel_within_event,
)


RESULTS = ROOT / "results"


def _fmt(v, n=3):
    if v is None:
        return "—"
    try:
        if v != v:  # NaN
            return "—"
        if v == float("inf"):
            return "∞"
        return f"{v:.{n}f}"
    except (TypeError, ValueError):
        return str(v)


def _holm_decision(h1_p: float, h3_p: float, alpha: float = 0.05) -> tuple[str, str]:
    """Apply Holm-Bonferroni at family-wise alpha across {H1, H3} (k=2).

    Returns (h1_status, h3_status) where each is 'reject@Holm' or 'fail'.
    """
    import math
    if math.isnan(h1_p) and math.isnan(h3_p):
        return ("N/A", "N/A")
    ps = []
    if not math.isnan(h1_p):
        ps.append(("H1", h1_p))
    if not math.isnan(h3_p):
        ps.append(("H3", h3_p))
    ps.sort(key=lambda x: x[1])
    decision = {}
    k = len(ps)
    for i, (name, p) in enumerate(ps):
        threshold = alpha / (k - i)
        decision[name] = "reject@Holm" if p < threshold else "fail"
        if decision[name] == "fail":
            # Holm step-down: once a test fails, all subsequent fail
            for j in range(i + 1, k):
                decision[ps[j][0]] = "fail"
            break
    return decision.get("H1", "N/A"), decision.get("H3", "N/A")


def main() -> int:
    t0 = time.time()
    in_path = RESULTS / "perpair_arm1.parquet"
    df = pd.read_parquet(in_path)
    print(f"[Phase 4] loaded {len(df)} rows from {in_path.name}", flush=True)

    print("[Phase 4] computing Family 1 (descriptive headline) ...", flush=True)
    f1 = family1_headline(df)

    print("[Phase 4] computing Family 2 H1 (mixed-effects + cluster-robust) ...", flush=True)
    h1 = family2_h1(df)

    print("[Phase 4] computing Family 2 H3 (LR+) ...", flush=True)
    h3 = family2_h3(df)

    print("[Phase 4] computing Mantel-Haenszel within-event robustness exhibit ...", flush=True)
    mh = h1_mantel_haenszel_within_event(df)

    # Holm decision within Family 2
    h1_holm, h3_holm = _holm_decision(h1.or_fe_pvalue, h3.lr_plus_pvalue_one_sided)

    # Aggregate
    out = {
        "phase": 4,
        "produced_at": time.strftime("%Y-%m-%d %H:%M:%S"),
        "input_file": str(in_path.name),
        "n_pairs": int(len(df)),
        "family_1": asdict(f1),
        "family_2": {
            "h1": asdict(h1),
            "h3": asdict(h3),
            "holm": {"h1": h1_holm, "h3": h3_holm, "family_wise_alpha": 0.05, "k": 2},
            "robustness_within_event_mh": asdict(mh),
        },
    }

    json_path = RESULTS / "phase4_results.json"
    with open(json_path, "w") as f:
        json.dump(out, f, indent=2, default=str)
    print(f"[Phase 4] wrote {json_path.name}", flush=True)

    # Human-readable report
    md = []
    md.append("# Phase 4 Results — Family 1 Headline + Family 2 Confirmatory\n")
    md.append(f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
    md.append(f"Input substrate: `results/{in_path.name}` ({len(df)} pairs)\n")
    md.append("\n---\n\n")

    # Family 1
    md.append("## Family 1 — Primary descriptive headline (no formal test)\n\n")
    md.append(f"**H2 (power-conditioned non-replication rate).** Within the FAERS-positive AND `daen_powered` Arm-1 universe (n={f1.n_faers_positive_daen_powered}), "
              f"the proportion classified `faers_only` is **{_fmt(f1.h2_non_replication_rate)}** "
              f"(95% cluster-bootstrap CI: {_fmt(f1.h2_ci_low)} – {_fmt(f1.h2_ci_high)}; "
              f"B=2,000; seed=94; clusters=events).\n\n")
    md.append(f"**Underpowered-discordant pairs** (FAERS-positive but DAEN underpowered): {f1.h2_underpowered_discordant_count}. "
              f"Reported alongside per Section 4.\n\n")
    md.append("**Concordance (on the `daen_powered` universe):**\n\n")
    md.append("| Statistic | Value | 95% CI |\n|---|---|---|\n")
    md.append(f"| PA+ (positive-specific agreement) | {_fmt(f1.pa_plus)} | {_fmt(f1.pa_plus_ci_low)} – {_fmt(f1.pa_plus_ci_high)} |\n")
    md.append(f"| PA+ (mutually-powered subset) | {_fmt(f1.pa_plus_power_conditioned)} | {_fmt(f1.pa_plus_pc_ci_low)} – {_fmt(f1.pa_plus_pc_ci_high)} |\n")
    md.append(f"| Cohen's κ | {_fmt(f1.cohen_kappa)} | — |\n\n")
    md.append(f"**Cohen's κ caveat:** {f1.kappa_caveat}\n\n")
    md.append("**Cross-DB 2×2 (on `daen_powered` universe):**\n\n")
    md.append("| | DAEN+ | DAEN− |\n|---|---|---|\n")
    md.append(f"| FAERS+ | {f1.cross_db_2x2['faers_pos_daen_pos']} | {f1.cross_db_2x2['faers_pos_daen_neg']} |\n")
    md.append(f"| FAERS− | {f1.cross_db_2x2['faers_neg_daen_pos']} | {f1.cross_db_2x2['faers_neg_daen_neg']} |\n\n")
    md.append("**Substrate notes:**\n")
    for note in f1.notes:
        md.append(f"- {note}\n")
    md.append("\n---\n\n")

    # Family 2
    md.append("## Family 2 — Secondary confirmatory (Holm-Bonferroni, k=2, family-wise α=0.05)\n\n")
    md.append("### H1 — Reference-negative enrichment\n\n")
    md.append(f"**Stratum:** {h1.stratum_size} pairs (FAERS-positive AND `daen_powered` AND in {{faers_only, concordant_positive}}). "
              f"Cluster structure: {h1.n_events_in_stratum} unique events; {h1.n_drugs_in_stratum} unique drugs.\n\n")
    md.append("**Per-stratum cells:**\n\n")
    md.append("| | faers_only | concordant_positive |\n|---|---|---|\n")
    md.append(f"| known-negative | {h1.n_known_neg_faers_only} | {h1.n_known_neg_concordant_pos} |\n")
    md.append(f"| known-positive | {h1.n_faers_only - h1.n_known_neg_faers_only} | "
              f"{h1.n_concordant_positive - h1.n_known_neg_concordant_pos} |\n\n")
    md.append("**Estimates:**\n\n")
    md.append("| Specification | OR | 95% CI | p-value |\n|---|---|---|---|\n")
    md.append(f"| Fixed-effect logistic | {_fmt(h1.or_fixed_effect)} | "
              f"{_fmt(h1.or_fe_ci_low)} – {_fmt(h1.or_fe_ci_high)} | {_fmt(h1.or_fe_pvalue, 4)} |\n")
    md.append(f"| Cluster-robust SE (cluster=event) | {_fmt(h1.or_cluster_robust)} | "
              f"{_fmt(h1.or_cr_ci_low)} – {_fmt(h1.or_cr_ci_high)} | {_fmt(h1.or_cr_pvalue, 4)} |\n\n")
    md.append(f"**Realised MDE (OR at 80% power, α=0.05):** {_fmt(h1.realised_mde_or, 2)}\n\n")
    md.append(f"**Pre-registered decision:** `{h1.decision}` (fragility flag: `{h1.fragility_flag}`)\n\n")
    md.append(f"**Holm-Bonferroni within Family 2:** `{h1_holm}`\n\n")
    for note in h1.notes:
        md.append(f"- {note}\n")
    md.append("\n")

    md.append("### H3 — Discordance-filter positive likelihood ratio\n\n")
    md.append(f"**Stratum:** {h3.stratum_size} pairs.\n\n")
    md.append("**Operating characteristics:**\n\n")
    md.append("| Metric | Value |\n|---|---|\n")
    md.append(f"| Sensitivity | {_fmt(h3.sensitivity)} |\n")
    md.append(f"| Specificity | {_fmt(h3.specificity)} |\n")
    md.append(f"| PPV | {_fmt(h3.ppv)} |\n")
    md.append(f"| NPV | {_fmt(h3.npv)} |\n")
    md.append(f"| LR+ | {_fmt(h3.lr_plus)} (95% CI: {_fmt(h3.lr_plus_ci_low)} – {_fmt(h3.lr_plus_ci_high)}) |\n")
    md.append(f"| One-sided p (LR+ > 1) | {_fmt(h3.lr_plus_pvalue_one_sided, 4)} |\n\n")
    md.append(f"**Realised LR+ MDE (one-sided 80% power):** {_fmt(h3.realised_mde_lrplus, 2)}\n\n")
    md.append(f"**Pre-registered decision:** `{h3.decision}`\n\n")
    md.append(f"**Holm-Bonferroni within Family 2:** `{h3_holm}`\n\n")
    for note in h3.notes:
        md.append(f"- {note}\n")
    md.append("\n---\n\n")

    md.append("## Within-event matched robustness (Mantel-Haenszel, Section 10.5)\n\n")
    md.append(f"- Pairs in analytic universe: {mh.n_pairs}\n")
    md.append(f"- Total event strata: {mh.n_strata_total}\n")
    md.append(f"- **Informative strata** (both rows AND both columns non-zero): **{mh.n_strata_informative}**\n")
    md.append(f"- MH pooled OR: {_fmt(mh.mh_or)} (95% CI: {_fmt(mh.mh_ci_low)} – {_fmt(mh.mh_ci_high)})\n\n")
    for note in mh.notes:
        md.append(f"- {note}\n")
    md.append("\n---\n\n")

    # Pivot status
    md.append("## Pre-registered pivot status (Section 5)\n\n")
    if h1.decision == "EXPLORATORY_DOWNGRADE":
        md.append("**Family 2 pivot for H1 FIRES (MDE > 3.0):** ")
        md.append("the realised MDE exceeds the pre-registered downgrade threshold; H1 is reported as exploratory; ")
        md.append("the headline falls back to Family 1 descriptive (concordance + power-conditioned non-replication rate). ")
        md.append("Pivot framing per protocol Section 5: *Calibrated negative-result paper.*\n\n")
    elif h1.decision == "DEGENERATE":
        md.append("**Family 2 H1 DEGENERATE:** the H1 stratum has too few pairs in at least one cell to estimate an OR. ")
        md.append("Family 1 descriptive remains the headline.\n\n")
    elif h1.decision == "REJECT_H0":
        md.append("**Family 2 H1 REJECTS H0:** reference-negative enrichment confirmed. ")
        md.append("Family 4 (Arm-1 FAERS-only mechanism intersection) determines whether the headline can use the word 'artefact'.\n\n")
    elif h1.decision == "REJECT_H0_FRAGILE":
        md.append("**Family 2 H1 REJECTS H0 (cluster-sensitive):** headline degraded to 'directional enrichment with cluster-sensitive significance' per Section 4.\n\n")
    else:
        md.append(f"**Family 2 H1: `{h1.decision}`.**\n\n")

    md_path = RESULTS / "phase4_results.md"
    with open(md_path, "w") as f:
        f.write("".join(md))
    print(f"[Phase 4] wrote {md_path.name}", flush=True)

    dt = time.time() - t0
    print(f"[Phase 4] complete in {dt:.1f}s", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
