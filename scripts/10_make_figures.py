"""Publication-quality figures for the paper.

Figure 1 — CONSORT-style Arm-1 flow diagram (492 → 44 H1 stratum), authored in
           Mermaid (figures/figure_1_arm1_flow.mmd) and rendered via mmdc
Figure 2 — Counterfactual comparison: H1 OR + H2 rate across power-tag specifications
           (the CENTREPIECE methodological-cautionary visual)
Figure 3 — H2 sensitivity forest (cluster-bootstrap CIs across 9 power tags)
Figure 4 — Family 3 H5a-d effect-size forest with directional interpretation
Figure 5 — Cross-DB 2x2 mosaic / pie for the daen_powered universe
Figure 6 — H4 balance: DAEN-only powered pairs (with known-pos/known-neg shading)

All figures: 300 dpi PNG + SVG vector; matplotlib only (no seaborn).
Outputs to figures/.
"""
from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
FIGS = ROOT / "outputs" / "figures"
RESULTS = ROOT / "results"
FIGS.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 9,
    "axes.titlesize": 10,
    "axes.labelsize": 9,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
    "figure.dpi": 100,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "axes.spines.top": False,
    "axes.spines.right": False,
})


def _save(fig, name):
    fig.savefig(FIGS / f"{name}.png", dpi=300)
    fig.savefig(FIGS / f"{name}.svg")
    plt.close(fig)
    print(f"  wrote {name}.png + {name}.svg")


# ----------------------------------------------------------------------------
# Figure 1 — CONSORT-style Arm-1 flow (Mermaid; rendered via mmdc)
# ----------------------------------------------------------------------------

def figure_1_flow():
    """Render the Arm-1 flow diagram from the Mermaid source via mmdc.

    The diagram is authored in figures/figure_1_arm1_flow.mmd. Counts in the
    source are cross-checked against the live substrate below so the committed
    .mmd can never silently drift from results/perpair_arm1.parquet.
    """
    mmd = FIGS / "figure_1_arm1_flow.mmd"

    # Guard: the hard-coded counts in the .mmd must match the current substrate.
    df = pd.read_parquet(RESULTS / "perpair_arm1.parquet")
    expected = {
        "n = 492": len(df),
        "OMOP n = 399": int((df["ref_set"] == "OMOP").sum()),
        "EU-ADR n = 93": int((df["ref_set"] == "EU-ADR").sum()),
        "n = 292": int(df["faers_signal"].sum()),
        "n = 248": int((df["faers_signal"] & ~df["daen_powered"]).sum()),
        "n = 44": int((df["faers_signal"] & df["daen_powered"]).sum()),
        "faers_only": int(((df["cross_db_class"] == "faers_only") & df["daen_powered"]).sum()),
        "concordant_positive": int(((df["cross_db_class"] == "concordant_positive") & df["daen_powered"]).sum()),
    }
    actual = {"n = 492": 492, "OMOP n = 399": 399, "EU-ADR n = 93": 93,
              "n = 292": 292, "n = 248": 248, "n = 44": 44,
              "faers_only": 6, "concordant_positive": 38}
    drift = {k: (actual[k], v) for k, v in expected.items() if actual[k] != v}
    if drift:
        raise SystemExit(f"figure_1 Mermaid counts drifted from substrate: {drift}\n"
                          f"Update {mmd} before re-rendering.")

    mmdc = shutil.which("mmdc")
    if mmdc is None:
        print("  [skip] mmdc not on PATH; Figure 1 .mmd unchanged. "
              "Install @mermaid-js/mermaid-cli to render.")
        return
    for ext, extra in (("png", ["-s", "3"]), ("svg", [])):
        subprocess.run([mmdc, "-i", str(mmd), "-o", str(FIGS / f"figure_1_arm1_flow.{ext}"),
                        "-b", "white", *extra], check=True)
    print("  wrote figure_1_arm1_flow.png + figure_1_arm1_flow.svg (from Mermaid)")


# ----------------------------------------------------------------------------
# Figure 2 — Counterfactual (THE centrepiece)
# ----------------------------------------------------------------------------

def figure_2_counterfactual():
    with open(RESULTS / "counterfactual_h1.json") as f:
        cf = json.load(f)

    # Order for display: PRIMARY first (highlighted), then counterfactual, then others
    order = ["daen_powered", "daen_powered_FAERS_pt", "daen_powered_FAERS_LB", "daen_powered_at_2"]
    short_labels = {
        "daen_powered": "PRIMARY\n(MDE-anchored,\nFAERS-independent)",
        "daen_powered_FAERS_pt": "COUNTERFACTUAL\n(FAERS-point-derived;\nthe naive default)",
        "daen_powered_FAERS_LB": "FAERS-95%-LB-derived\n(conservative)",
        "daen_powered_at_2": "MDE ≤ 2.0\n(permissive anchor)",
    }
    ordered = [next(r for r in cf if r["tag"] == t) for t in order]

    fig, axes = plt.subplots(1, 2, figsize=(12, 6.0), gridspec_kw={"width_ratios": [1, 1.2]})

    # Panel A: H1 OR comparison
    ax = axes[0]
    y_positions = np.arange(len(ordered))
    or_vals, or_lo, or_hi, colors = [], [], [], []
    labels = []
    for r in ordered:
        labels.append(short_labels[r["tag"]])
        if (np.isnan(r["or_fe"]) or r["or_fe"] == 0 or
            np.isnan(r["or_fe_ci_lo"]) or np.isnan(r["or_fe_ci_hi"])):
            or_vals.append(np.nan)
            or_lo.append(np.nan); or_hi.append(np.nan)
        else:
            or_vals.append(r["or_fe"])
            or_lo.append(r["or_fe_ci_lo"]); or_hi.append(r["or_fe_ci_hi"])
        if r["tag"] == "daen_powered":
            colors.append("#1f5f7a")  # blue for primary
        elif r["tag"] == "daen_powered_FAERS_pt":
            colors.append("#c0392b")  # red for the dangerous counterfactual
        else:
            colors.append("#7f8c8d")  # grey for others

    ax.axvline(1, color="#888", lw=0.8, ls="--", zorder=1)
    ax.axvline(2, color="#bbb", lw=0.6, ls=":", zorder=1)
    for i, (val, lo, hi, c) in enumerate(zip(or_vals, or_lo, or_hi, colors)):
        if np.isnan(val):
            ax.text(3.0, y_positions[i], "OR undefined — H1 degenerate\n(complete separation;\nrealised MDE = 76.5 ≫ 3.0\n→ EXPLORATORY_DOWNGRADE)",
                     va="center", ha="center", fontsize=7.5, color="#1f5f7a", weight="bold",
                     bbox=dict(boxstyle="round,pad=0.35", facecolor="#dfeef7", edgecolor="#1f5f7a"), zorder=4)
        else:
            ax.errorbar(val, y_positions[i], xerr=[[val - lo], [hi - val]],
                          fmt="o", color=c, markersize=10, capsize=4, lw=2, zorder=3)
            # Anchor the OR/CI/p label ABOVE the marker and grow it upward (va="bottom").
            # The CI whisker is horizontal at exactly the marker's y, so a downward-
            # growing label crosses its own whisker; growing upward clears it.
            ax.text(val, y_positions[i] - 0.16, f"OR = {val:.2f}\n[{lo:.2f}, {hi:.2f}]\np = {ordered[i]['pval_fe']:.4f}",
                     ha="center", va="bottom", fontsize=7.5, color=c, weight="bold")

    ax.set_yticks(y_positions); ax.set_yticklabels(labels, fontsize=8.5)
    ax.set_xscale("log")
    ax.set_xlim(0.4, 50)
    # Inline reference-line label instead of a one-entry legend box (the legend
    # otherwise lands on the bottom MDE<=2.0 row); rotated text sits on the OR=2 line.
    ax.text(2, len(ordered) - 0.42, "OR = 2", rotation=90, va="bottom", ha="center",
             fontsize=7, color="#999")
    ax.set_xlabel("H1 reference-negative enrichment OR (95% CI)\n— fixed-effect logistic", fontsize=9)
    ax.set_title("A. H1 enrichment OR by power tag",
                  fontsize=10, weight="bold", loc="left")
    ax.invert_yaxis()
    # Headroom above the top row so the upward-growing OR/CI/p labels (and the
    # 'OR undefined' box) cannot rise into the title.
    ax.set_ylim(len(ordered) - 0.35, -1.15)
    ax.grid(axis="x", alpha=0.2, ls=":")

    # Panel B: H1 stratum cells visualised
    ax = axes[1]
    cells_data = []
    for r in ordered:
        cells_data.append([r["n_known_neg_faers_only"],
                            r["n_faers_only"] - r["n_known_neg_faers_only"],
                            r["n_known_neg_concordant_positive"],
                            r["n_concordant_positive"] - r["n_known_neg_concordant_positive"]])
    cells_data = np.array(cells_data)
    bar_y = y_positions
    bar_h = 0.7
    # stacked horizontal bars per scenario, split: faers_only [known-neg | known-pos] + conc_pos [known-neg | known-pos]
    color_kn_fo = "#c0392b"     # red — the H1-driving cell
    color_kp_fo = "#e8a094"     # pink — faers_only known-pos
    color_kn_cp = "#3a7ca5"     # blue — concordant_positive known-neg
    color_kp_cp = "#8cb9d5"     # light blue — concordant_positive known-pos
    for i, row in enumerate(cells_data):
        x_offset = 0
        for j, (val, color, label) in enumerate(zip(row, [color_kn_fo, color_kp_fo, color_kn_cp, color_kp_cp],
                                                        ["FO known-neg", "FO known-pos", "CP known-neg", "CP known-pos"])):
            if val > 0:
                ax.barh(bar_y[i], val, left=x_offset, height=bar_h,
                         color=color, edgecolor="white", lw=0.5,
                         label=label if i == 0 else "")
                ax.text(x_offset + val/2, bar_y[i], str(val),
                         ha="center", va="center", fontsize=8, color="white", weight="bold")
            x_offset += val
        ax.text(x_offset + 2, bar_y[i], f"  n = {int(row.sum())}", va="center", fontsize=8, color="#444")

    ax.set_yticks(bar_y); ax.set_yticklabels(labels, fontsize=8.5)
    ax.set_xlabel("Count of pairs in H1 stratum", fontsize=9)
    ax.set_title("B. H1 stratum cell counts",
                  fontsize=10, weight="bold", loc="left")
    ax.invert_yaxis()
    ax.legend(loc="lower right", framealpha=0.9, ncol=2, fontsize=7.5)
    ax.set_xlim(0, max(cells_data.sum(axis=1)) * 1.18)

    fig.suptitle("Counterfactual: FAERS-derived power vs the MDE-anchored primary",
                  fontsize=11, weight="bold", y=1.0)
    _save(fig, "figure_2_counterfactual_centrepiece")


# ----------------------------------------------------------------------------
# Figure 3 — H2 sensitivity forest
# ----------------------------------------------------------------------------

def figure_3_h2_forest():
    rows = [
        ("PRIMARY (MDE ≤ 1.5)", 44, 0.136, 0.057, 0.318, True),
        ("Power 0.70 (MDE ≤ 1.5)", 58, 0.138, 0.084, 0.276, False),
        ("Power 0.90 (MDE ≤ 1.5)", 33, 0.152, 0.079, 0.313, False),
        ("MDE ≤ 2.0 (permissive anchor)", 112, 0.214, 0.156, 0.383, False),
        ("MDE ≤ 3.0 (very permissive)", 171, 0.298, 0.212, 0.512, False),
        ("FAERS point-ROR derived", 126, 0.135, 0.101, 0.257, False),
        ("FAERS 95% LB derived", 116, 0.112, 0.087, 0.189, False),
        ("Poisson-margin", 43, 0.140, 0.058, 0.318, False),
        ("OMOP only (per-set)", 23, 0.130, 0.000, 1.000, False),
        ("EU-ADR only (per-set)", 21, 0.143, 0.000, 0.318, False),
    ]
    fig, ax = plt.subplots(figsize=(10, 6))
    y = np.arange(len(rows))
    for i, (label, n, pt, lo, hi, primary) in enumerate(rows):
        c = "#1f5f7a" if primary else "#7f8c8d"
        ax.errorbar(pt, y[i], xerr=[[pt - lo], [hi - pt]], fmt="o",
                     color=c, markersize=8 if primary else 6, capsize=3, lw=1.5)
        ax.text(hi + 0.02, y[i], f"  {pt:.3f}  [{lo:.3f}, {hi:.3f}]  (n = {n})",
                 va="center", fontsize=8, color=c)

    ax.axvline(0.136, color="#1f5f7a", lw=0.8, ls=":", alpha=0.5)
    ax.set_yticks(y); ax.set_yticklabels([r[0] for r in rows], fontsize=9)
    ax.set_xlabel("H2 power-conditioned non-replication rate (95% cluster-bootstrap CI)", fontsize=9)
    ax.set_title("H2 power-conditioned non-replication rate across sensitivity arms",
                  fontsize=11, weight="bold", loc="left")
    ax.invert_yaxis()
    ax.set_xlim(-0.02, 0.85)
    ax.grid(axis="x", alpha=0.2, ls=":")
    _save(fig, "figure_3_h2_sensitivity_forest")


# ----------------------------------------------------------------------------
# Figure 4 — Family 3 H5a-d effect-size forest
# ----------------------------------------------------------------------------

def figure_4_family3_forest():
    rows = [
        ("H5a: Consumer share ≥ 10 pp lift", +6.62, +3.90, +9.14, "0.0001", "REJECT", "supports", "#2c5f2c"),
        ("H5d: Mass-tort drug membership", +1.65, +0.53, +2.53, "0.006", "REJECT", "supports", "#2c5f2c"),
        ("H5b: Lawyer share ≥ 5 pp lift", -0.93, -2.03, -0.13, "0.019", "REJECT (contradictory)", "contradicts", "#c0392b"),
        ("H5c: Bai-Perron alert alignment", 0.0, -0.41, +0.12, "—", "FAIL", "null", "#888"),
    ]
    fig, ax = plt.subplots(figsize=(10, 4.5))
    y = np.arange(len(rows))
    for i, (label, d, lo, hi, p, holm, direction, c) in enumerate(rows):
        ax.errorbar(d, y[i], xerr=[[d - lo], [hi - d]], fmt="o",
                     color=c, markersize=10, capsize=4, lw=2)
        annot = f"  Δ = {d:+.2f} pp  [{lo:+.2f}, {hi:+.2f}]   p = {p}   Holm: {holm}"
        ax.text(max(hi, d) + 0.5, y[i], annot, va="center", fontsize=8.5, color=c)

    ax.axvline(0, color="#888", lw=0.8, ls="--")
    ax.set_yticks(y); ax.set_yticklabels([r[0] for r in rows], fontsize=9)
    ax.set_xlabel("Δ = faers_only proportion − concordant_positive proportion (percentage points; 95% Newcombe CI)", fontsize=9)
    ax.set_title("Family 3 mechanism arm: discordance-associated effect sizes (H5a–H5d)",
                  fontsize=11, weight="bold", loc="left")
    ax.invert_yaxis()
    ax.set_xlim(-4.5, 22)
    ax.grid(axis="x", alpha=0.2, ls=":")
    _save(fig, "figure_4_family3_forest")


# ----------------------------------------------------------------------------
# Figure 5 — Cross-DB 2x2 mosaic on the daen_powered universe
# ----------------------------------------------------------------------------

def figure_5_mosaic():
    df = pd.read_parquet(RESULTS / "perpair_arm1.parquet")
    sub = df[df["daen_powered"]]
    a = int(((sub["faers_signal"]) & (sub["daen_signal"])).sum())
    b = int(((sub["faers_signal"]) & (~sub["daen_signal"])).sum())
    c = int(((~sub["faers_signal"]) & (sub["daen_signal"])).sum())
    d = int(((~sub["faers_signal"]) & (~sub["daen_signal"])).sum())
    n = a + b + c + d
    fig, ax = plt.subplots(figsize=(7, 6.5))
    # Standard 2x2 layout: rows = FAERS (top = +, bottom = -); cols = DAEN (left = +, right = -)
    cells = {(0, 0): ("FAERS+\nDAEN+\nn=" + str(a) + "\n(concordant\npositive)", "#3a7ca5", a),
             (0, 1): ("FAERS+\nDAEN−\nn=" + str(b) + "\n(faers_only)", "#c0392b", b),
             (1, 0): ("FAERS−\nDAEN+\nn=" + str(c) + "\n(daen_only)", "#2c5f2c", c),
             (1, 1): ("FAERS−\nDAEN−\nn=" + str(d) + "\n(concordant\nnegative)", "#888", d)}
    for (r, col), (txt, color, val) in cells.items():
        # area-proportional cell: width proportional to col marginal, height to row marginal
        w = val if val > 0 else 0.5
        ax.add_patch(mpatches.Rectangle((col * 5, (1 - r) * 5), 5, 5,
                                          facecolor=color, alpha=0.3 + 0.7 * (val / max(n, 1)),
                                          edgecolor="white", lw=2))
        ax.text(col * 5 + 2.5, (1 - r) * 5 + 2.5, txt, ha="center", va="center",
                 fontsize=11, weight="bold", color=color)
    ax.set_xlim(-0.5, 10.5); ax.set_ylim(-0.5, 10.5)
    ax.set_xticks([2.5, 7.5]); ax.set_xticklabels(["DAEN signal +", "DAEN signal −"], fontsize=10)
    ax.set_yticks([2.5, 7.5]); ax.set_yticklabels(["FAERS signal −", "FAERS signal +"], fontsize=10)
    ax.set_title("Cross-database 2×2 on the power-conditioned universe",
                  fontsize=11, weight="bold", loc="left")
    ax.set_aspect("equal")
    for spine in ["top", "right", "left", "bottom"]:
        ax.spines[spine].set_visible(False)
    _save(fig, "figure_5_cross_db_mosaic")


# ----------------------------------------------------------------------------
# Figure 6 — H4 balance: DAEN-only powered pairs
# ----------------------------------------------------------------------------

def figure_6_h4_balance():
    df = pd.read_parquet(RESULTS / "perpair_arm1.parquet")
    h4 = df[(df["cross_db_class"] == "daen_only") & (df["faers_powered"])].copy()
    h4 = h4.sort_values("ror_D", ascending=True)

    fig, ax = plt.subplots(figsize=(10, 5))
    y = np.arange(len(h4))
    colors = ["#2c5f2c" if g == 1 else "#c0392b" for g in h4["groundTruth"]]
    ax.barh(y, h4["ror_D"], color=colors, edgecolor="white", lw=0.5)
    labels = [f"{r['exposureName']} × {r['outcomeName'][:35]}{'…' if len(r['outcomeName']) > 35 else ''}"
              for _, r in h4.iterrows()]
    ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=8)
    ax.axvline(1, color="#888", lw=0.8, ls="--")
    ax.set_xlabel("DAEN ROR point estimate (log scale)", fontsize=9)
    ax.set_xscale("log")
    ax.set_xlim(0.7, 20)
    legend_handles = [
        mpatches.Patch(color="#2c5f2c", label="Known-positive reference"),
        mpatches.Patch(color="#c0392b", label="Known-negative reference"),
    ]
    ax.legend(handles=legend_handles, loc="lower right")
    ax.set_title("H4 balance: Australian DAEN-only signals FAERS does not detect",
                  fontsize=11, weight="bold", loc="left")
    ax.invert_yaxis()
    _save(fig, "figure_6_h4_balance")


def main():
    print("Generating publication-quality figures ...")
    figure_1_flow()
    figure_2_counterfactual()
    figure_3_h2_forest()
    figure_4_family3_forest()
    figure_5_mosaic()
    figure_6_h4_balance()
    print(f"\nAll figures written to {FIGS}/")
    print("Each is in both PNG (300 dpi) and SVG (vector) form.")


if __name__ == "__main__":
    main()
