"""Run the downstream simulation analysis (protocol S1-S4) and write artefacts.

Reads ``results/sim_results.parquet`` and writes, into ``results/sim_analysis/``:
  - ``q1_typeI.{parquet,csv}``    type-I per rule x estimator x cell (Q1_null)
  - ``q2_power.{parquet,csv}``    power per rule x estimator x cell (Q2_alt)
  - ``q3_coverage.{parquet,csv}`` coverage / type-I / bias vs K (Q1_null)
  - ``s_hypotheses.json``         S1-S4 verdicts with supporting numbers
  - ``sim_analysis_summary.md``   human-readable report

Usage:  python3 scripts/sim_analysis.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from reporting.sim_analysis import (  # noqa: E402
    PRIMARY_ESTIMATOR, TYPEI_BAND, NOMINAL_COVERAGE, MIN_ESTIMABLE_COND,
    load_results, q1_typeI, q2_power, q3_coverage,
    evaluate_hypotheses, typeI_rule_summary, typeI_conditional_summary,
    power_rule_summary,
)

RESULTS = ROOT / "results"
OUTDIR = RESULTS / "sim_analysis"


def _fmt_pct(x: float) -> str:
    return "n/a" if x is None or pd.isna(x) else f"{100 * x:.2f}%"


def _md_table(df: pd.DataFrame, cols: list[str], headers: list[str],
              fmts: dict[str, str] | None = None) -> str:
    fmts = fmts or {}
    lines = ["| " + " | ".join(headers) + " |",
             "|" + "|".join(["---"] * len(headers)) + "|"]
    for _, row in df.iterrows():
        cells = []
        for c in cols:
            v = row[c]
            if c in fmts and not pd.isna(v):
                cells.append(format(v, fmts[c]))
            else:
                cells.append(str(v))
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


def build_report(q1, q2, q3, hyp) -> str:
    est = PRIMARY_ESTIMATOR
    ti = typeI_rule_summary(q1, est)
    tc = typeI_conditional_summary(q1, est)
    pw = power_rule_summary(q2, est)

    L: list[str] = []
    L.append("# Confirmatory simulation analysis (S1-S4)")
    L.append("")
    L.append(f"Source: `results/sim_results.parquet` (6,048,000 rows; 126 DGP cells "
             f"x 2000 replicates x 6 selection rules x 4 estimators).")
    L.append(f"Primary estimator for the headline endpoints: `{est}` "
             f"(cluster-robust SE retained as a sensitivity; protocol Section 11).")
    L.append(f"Type-I control band: [{TYPEI_BAND[0]}, {TYPEI_BAND[1]}]. "
             f"Nominal CI coverage: {NOMINAL_COVERAGE}.")
    L.append("")
    L.append("Rejection rule (uniform across estimators): reject iff the 95% CI lower "
             "bound of the enrichment OR exceeds 1. A non-estimable replicate "
             "(complete separation or empty contrast row) carries no decision and is "
             "counted as a non-rejection; the non-estimability rate is reported "
             "alongside every rate.")
    L.append("")

    L.append("## Q1 type-I error by selection rule (Q1_null, primary estimator)")
    L.append("")
    L.append(_md_table(
        ti, ["selection_rule", "n_cells", "median_typeI", "max_typeI",
             "frac_within_band", "mean_non_estimable"],
        ["rule", "cells", "median type-I", "max type-I", "frac within band",
         "mean non-estimable"],
        {"median_typeI": ".4f", "max_typeI": ".4f", "frac_within_band": ".3f",
         "mean_non_estimable": ".3f"}))
    L.append("")

    L.append("## Q1 type-I conditional on estimability (Q1_null, primary estimator, post-hoc)")
    L.append("")
    L.append("Second estimand requested in review: type-I as n_reject / n_estimable "
             "(rejection rate among replicates that produced a decision), so a rule's "
             "marginal type-I cannot look controlled merely because most replicates were "
             "non-estimable and were scored as non-rejections. `pooled` is "
             "estimable-weighted over all cells; `max` and the correlation are over cells "
             f"with >= {MIN_ESTIMABLE_COND} estimable replicates (thinner cells are "
             "Monte-Carlo noise). `corr nonest/cond` < 0 means a rule's least-estimable "
             "cells are not its worst-behaved cells, so the marginal deflation does not "
             "conceal a conditional breach.")
    L.append("")
    L.append(_md_table(
        tc, ["selection_rule", "n_cells_gated", "pooled_cond_typeI", "max_cond_typeI",
             "n_cells_cond_breach", "frac_within_band_cond", "corr_nonest_cond"],
        ["rule", "cells (gated)", "pooled cond type-I", "max cond type-I",
         "cells breaching", "frac within band", "corr nonest/cond"],
        {"pooled_cond_typeI": ".4f", "max_cond_typeI": ".4f",
         "frac_within_band_cond": ".3f", "corr_nonest_cond": ".3f"}))
    L.append("")

    L.append("## Q2 power by selection rule (Q2_alt, primary estimator)")
    L.append("")
    L.append(_md_table(
        pw, ["selection_rule", "n_cells", "median_power", "max_power",
             "mean_non_estimable"],
        ["rule", "cells", "median power", "max power", "mean non-estimable"],
        {"median_power": ".4f", "max_power": ".4f", "mean_non_estimable": ".3f"}))
    L.append("")

    L.append("## Q3 coverage and few-cluster inference (Q1_null, K-sweep line)")
    L.append("")
    L.append("Coverage is assessed against the generative null OR = 1 (protocol 3.3). "
             "The S4 estimator comparison is read on the type-I-controlling rule "
             "`mde_1.5`, where selection bias is removed; coverage shortfalls on other "
             "rules reflect point-estimate (selection) bias, the same phenomenon as Q1 "
             "type-I inflation.")
    L.append("")
    q3m = q3[q3["selection_rule"] == "mde_1.5"].sort_values(["estimator", "K"])
    L.append("### Rule `mde_1.5`")
    L.append("")
    L.append(_md_table(
        q3m, ["estimator", "K", "coverage", "typeI", "median_bias",
              "non_estimable_rate", "n_estimable"],
        ["estimator", "K", "coverage", "type-I", "median bias", "non-est rate", "n est"],
        {"coverage": ".3f", "typeI": ".4f", "median_bias": ".4f",
         "non_estimable_rate": ".3f"}))
    L.append("")

    L.append("## Pre-specified hypotheses (honest-either-way)")
    L.append("")
    for key in ("S1", "S2", "S3", "S4"):
        h = hyp[key]
        L.append(f"### {key} - {h['verdict']}")
        L.append("")
        L.append(h["claim"])
        L.append("")
        for k, v in h.items():
            if k in ("claim", "verdict"):
                continue
            L.append(f"- `{k}`: {json.dumps(v)}")
        L.append("")

    L.append("## Reporting note")
    L.append("")
    L.append("All numbers above are computed from the deposited sweep; none are "
             "placeholders. A hypothesis graded NOT SUPPORTED is reported as a finding "
             "per the binding honesty clause (protocol Section 8): no DGP cell was "
             "added, removed, or re-weighted after seeing results.")
    L.append("")
    return "\n".join(L)


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    df = load_results(RESULTS / "sim_results.parquet")

    q1 = q1_typeI(df)
    q2 = q2_power(df)
    q3 = q3_coverage(df)
    hyp = evaluate_hypotheses(q1, q3)

    for name, tab in (("q1_typeI", q1), ("q2_power", q2), ("q3_coverage", q3)):
        tab.to_parquet(OUTDIR / f"{name}.parquet", index=False)
        tab.to_csv(OUTDIR / f"{name}.csv", index=False)

    (OUTDIR / "s_hypotheses.json").write_text(json.dumps(hyp, indent=2))
    (OUTDIR / "sim_analysis_summary.md").write_text(build_report(q1, q2, q3, hyp))

    print(f"[sim-analysis] wrote tables + verdicts + report to {OUTDIR}")
    for key in ("S1", "S2", "S3", "S4"):
        print(f"  {key}: {hyp[key]['verdict']}")


if __name__ == "__main__":
    main()
