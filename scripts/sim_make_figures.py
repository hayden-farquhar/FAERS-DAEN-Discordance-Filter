"""Figures for the confirmatory simulation (protocol S1-S4).

Reads the analysis tables written by ``scripts/sim_analysis.py`` into
``results/sim_analysis/`` (plus the post-hoc tables in ``results/sim_wcr/`` and
``results/sim_density/``) and writes four figures to ``outputs/figures/``:

  sim_1_typeI_surface   type-I across lambda x phi, faceted by selection rule
                        (Q1_null, primary estimator), with the [0,0.075] band
  sim_2_power_curve     power vs true enrichment OR per rule on the OR sweep
                        line (Q2_alt, primary estimator)
  sim_3_coverage_K      CI coverage vs cluster count K per estimator on the
                        type-I-controlling rule mde_1.5 (S4 panel); the post-hoc
                        restricted wild bootstrap (WCR) is overlaid as a dashed
                        curve read from results/sim_wcr/
  sim_4_density_sweep   post-hoc reference-negative-density sweep: type-I (both
                        estimands) and non-estimability vs n_pairs for mde_1.5
                        and none, read from results/sim_density/

All figures: 300 dpi PNG + SVG; matplotlib only (no seaborn). The numbers are
read from the deposited tables, never re-derived here, so the figures cannot
drift from the analysis.

Usage:  python3 scripts/sim_make_figures.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from reporting.sim_analysis import (  # noqa: E402
    PRIMARY_ESTIMATOR, TYPEI_BAND, NOMINAL_COVERAGE, SWEEP_BASE,
)

FIGS = ROOT / "outputs" / "figures"
ANALYSIS = ROOT / "results" / "sim_analysis"
WCR = ROOT / "results" / "sim_wcr"
DENSITY = ROOT / "results" / "sim_density"
FIGS.mkdir(parents=True, exist_ok=True)

# Tier-1 portfolio figure standard: SciencePlots (nature/science, no-latex) if
# available, else a curated fallback that matches its look on this toolchain.
# Either way the colour identities come from the shared Okabe-Ito palette and the
# fonts are embedded (TrueType type-42) so the vector exports are press-ready.
try:
    import scienceplots  # noqa: F401
    plt.style.use(["science", "nature", "no-latex"])
except Exception:
    pass

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
    "pdf.fonttype": 42,   # embed TrueType (no Type-3) for journal vector files
    "ps.fonttype": 42,
})

RULE_ORDER = ["none", "faers_point", "faers_lb", "mde_3.0", "mde_2.0", "mde_1.5"]

# Estimator display names + plot order for the coverage panel (S4). The
# restricted wild bootstrap (WCR) is a post-hoc addition, drawn dashed.
EST_NICE = {
    "cluster_robust": "cluster-robust SE",
    "wild_cluster_bootstrap": "wild bootstrap (WCU)",
    "wild_cluster_bootstrap_restricted": "wild bootstrap, restricted (WCR)",
    "fixed_effect": "fixed-effect",
    "firth": "Firth",
}
EST_ORDER = [
    "cluster_robust", "wild_cluster_bootstrap",
    "wild_cluster_bootstrap_restricted", "fixed_effect", "firth",
]
# Okabe-Ito colourblind-safe palette for the density panel (three selection
# rules: the registered MDE anchor, the contaminated FAERS-derived rule, and the
# unconditioned collider baseline).
DENS_COLOR = {"mde_1.5": "#0072B2", "faers_point": "#009E73", "none": "#D55E00"}


def _save(fig, name: str) -> None:
    fig.savefig(FIGS / f"{name}.png", dpi=300)
    fig.savefig(FIGS / f"{name}.svg")
    plt.close(fig)
    print(f"  wrote {name}.png + {name}.svg")


def figure_typeI_surface(q1: pd.DataFrame) -> None:
    """Type-I heatmaps over lambda x phi, one panel per rule, primary estimator.

    Averaged over pi (the third primary-slice factor). The [0,0.075] band edge
    is drawn as a contour so band exceedance is visible at a glance."""
    d = q1[(q1["estimator"] == PRIMARY_ESTIMATOR) & (q1["label"] == "primary")]
    lambdas = sorted(d["lambda_inflate"].unique())
    phis = sorted(d["phi"].unique())

    rules = [r for r in RULE_ORDER if r in set(d["selection_rule"])]
    ncol = 3
    nrow = int(np.ceil(len(rules) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.1 * ncol, 2.7 * nrow),
                             constrained_layout=True)
    axes = np.atleast_1d(axes).ravel()

    vmax = max(float(d["reject_rate"].max()), TYPEI_BAND[1])
    grids = {}
    for rule in rules:
        dr = d[d["selection_rule"] == rule]
        g = (dr.groupby(["phi", "lambda_inflate"], observed=True)["reject_rate"]
               .mean().reset_index())
        grid = (g.pivot(index="phi", columns="lambda_inflate", values="reject_rate")
                  .reindex(index=phis, columns=lambdas))
        grids[rule] = grid.to_numpy()

    im = None
    for ax, rule in zip(axes, rules):
        arr = grids[rule]
        im = ax.imshow(arr, origin="lower", aspect="auto", cmap="magma",
                       vmin=0.0, vmax=vmax)
        ax.set_xticks(range(len(lambdas)), [str(x) for x in lambdas])
        ax.set_yticks(range(len(phis)), [str(x) for x in phis])
        ax.set_title(rule)
        ax.set_xlabel(r"$\lambda$ (notoriety inflation)")
        ax.set_ylabel(r"$\varphi$ (artefact share)")
        for i in range(len(phis)):
            for j in range(len(lambdas)):
                v = arr[i, j]
                if np.isnan(v):
                    continue
                ax.text(j, i, f"{v:.3f}", ha="center", va="center",
                        fontsize=6,
                        color="white" if v < 0.6 * vmax else "black")
        # band-edge contour (type-I = 0.075)
        if np.nanmax(arr) >= TYPEI_BAND[1] >= np.nanmin(arr):
            ax.contour(arr, levels=[TYPEI_BAND[1]], colors="cyan",
                       linewidths=1.2)

    for ax in axes[len(rules):]:
        ax.set_visible(False)

    if im is not None:
        cb = fig.colorbar(im, ax=axes.tolist(), shrink=0.6, location="right")
        cb.set_label("type-I error rate")
        cb.ax.axhline(TYPEI_BAND[1], color="cyan", lw=1.2)
    fig.suptitle(
        f"Q1 type-I across $\\lambda\\times\\varphi$ (Q1_null, {PRIMARY_ESTIMATOR}; "
        f"cyan = {TYPEI_BAND[1]:.3f} band edge)", fontsize=10)
    _save(fig, "sim_1_typeI_surface")


def _numeric_or(s: pd.Series) -> pd.Series:
    return pd.to_numeric(s, errors="coerce")


def figure_power_curve(q2: pd.DataFrame) -> None:
    """Power vs true enrichment OR per rule, primary estimator.

    Read on the OR_true sweep line (base pi/lambda/phi), where OR_true is varied
    and parses to a number; the 'lognormal' mixture level is excluded from the
    line because it has no single x-coordinate."""
    d = q2[q2["estimator"] == PRIMARY_ESTIMATOR].copy()
    base = ((d["pi"] == SWEEP_BASE["pi"])
            & (d["lambda_inflate"] == SWEEP_BASE["lambda_inflate"])
            & (d["phi"] == SWEEP_BASE["phi"]))
    line = d[(d["label"] == "sweep_OR_true") | ((d["label"] == "primary") & base)]
    line = line.assign(or_num=_numeric_or(line["OR_true"]))
    line = line[line["or_num"].notna()].sort_values("or_num")

    fig, ax = plt.subplots(figsize=(5.4, 3.8), constrained_layout=True)
    rules = [r for r in RULE_ORDER if r in set(line["selection_rule"])]
    cmap = plt.get_cmap("viridis", len(rules))
    for i, rule in enumerate(rules):
        dr = (line[line["selection_rule"] == rule]
              .groupby("or_num", observed=True)["reject_rate"].mean())
        ax.plot(dr.index, dr.values, marker="o", ms=4, lw=1.4,
                color=cmap(i), label=rule)
    ax.set_xlabel("true enrichment OR (artefact concentration in FAERS-only)")
    ax.set_ylabel("power (correct-rejection rate)")
    ax.set_ylim(-0.02, 1.0)
    ax.set_title(f"Q2 power vs effect size (Q2_alt, {PRIMARY_ESTIMATOR})")
    ax.legend(title="selection rule", ncol=2, frameon=False)
    _save(fig, "sim_2_power_curve")


def figure_coverage_K(q3: pd.DataFrame, wcr: pd.DataFrame | None = None) -> None:
    """CI coverage vs cluster count K per estimator on mde_1.5 (S4 panel).

    The four pre-registered estimators are read from the frozen q3 table; the
    post-hoc restricted wild bootstrap (WCR), if supplied, is appended from
    results/sim_wcr/ and drawn dashed to mark it as a post-hoc addition. The
    WCR file's cluster_robust / wild_cluster_bootstrap rows reproduce the frozen
    q3 exactly, so only the restricted estimator is taken from it to avoid
    double-drawing the shared curves."""
    d = q3[q3["selection_rule"] == "mde_1.5"][["estimator", "K", "coverage"]]
    frames = [d]
    if wcr is not None:
        wr = wcr[wcr["estimator"] == "wild_cluster_bootstrap_restricted"]
        frames.append(wr[["estimator", "K", "coverage"]])
    d = pd.concat(frames, ignore_index=True).sort_values(["estimator", "K"])

    present = set(d["estimator"])
    ests = [e for e in EST_ORDER if e in present] + \
           [e for e in sorted(present) if e not in EST_ORDER]
    cmap = plt.get_cmap("tab10")

    fig, ax = plt.subplots(figsize=(5.8, 3.9), constrained_layout=True)
    for i, est in enumerate(ests):
        de = d[d["estimator"] == est].sort_values("K")
        is_wcr = est == "wild_cluster_bootstrap_restricted"
        ax.plot(de["K"], de["coverage"],
                marker="D" if is_wcr else "o", ms=5, lw=1.6,
                ls="--" if is_wcr else "-",
                color=cmap(i), label=EST_NICE.get(est, est),
                zorder=4 if is_wcr else 3)
    ax.axhline(NOMINAL_COVERAGE, color="black", ls=":", lw=1.0,
               label=f"nominal {NOMINAL_COVERAGE:.2f}")
    ax.set_xscale("log")
    ax.set_xticks(sorted(d["K"].unique()))
    ax.get_xaxis().set_major_formatter(plt.ScalarFormatter())
    ax.set_xlabel("clusters K (log scale)")
    ax.set_ylabel("empirical 95% CI coverage (vs generative null OR=1)")
    ax.set_title("S4 few-cluster coverage on mde_1.5 (Q1_null)")
    ax.legend(frameon=False, ncol=1, fontsize=7)
    _save(fig, "sim_3_coverage_K")


def figure_density_sweep(dens: pd.DataFrame) -> None:
    """Post-hoc reference-negative-density sweep (base null cell, primary est).

    Reads the Q1_null rows of the two-regime density summary. Panel A: type-I
    (marginal solid, conditional-on-estimable dashed) vs n_pairs for the
    registered MDE anchor, the contaminated FAERS-derived rule, and the
    unconditioned baseline, with the [0,0.075] acceptance band shaded. Panel B:
    non-estimable replicate rate vs n_pairs, showing the powered strata filling
    in as density grows. The faers_point curve answers the reviewer's
    estimability objection directly: if the contaminated rule stays in band even
    at 8x density, its in-band type-I is not an artefact of thin-stratum
    non-estimability."""
    if "regime" in dens.columns:
        dens = dens[dens["regime"] == "Q1_null"]
    rules = [r for r in ["mde_1.5", "faers_point", "none"] if r in set(dens["rule"])]
    xs = sorted(dens["n_pairs"].unique())

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(8.8, 3.9),
                                   constrained_layout=True)

    axA.axhspan(TYPEI_BAND[0], TYPEI_BAND[1], color="0.90", zorder=0)
    axA.axhline(TYPEI_BAND[1], color="0.55", lw=0.9, ls="-", zorder=1,
                label=f"band edge {TYPEI_BAND[1]:.3f}")
    for rule in rules:
        dr = dens[dens["rule"] == rule].sort_values("n_pairs")
        c = DENS_COLOR.get(rule, "black")
        axA.plot(dr["n_pairs"], dr["marginal_reject"], marker="o", ms=5, lw=1.6,
                 color=c, label=f"{rule} (marginal)", zorder=3)
        axA.plot(dr["n_pairs"], dr["conditional_reject"], marker="s", ms=4,
                 lw=1.3, ls="--", color=c, label=f"{rule} (conditional)",
                 zorder=3)
    axA.set_xscale("log")
    axA.set_xticks(xs)
    axA.get_xaxis().set_major_formatter(plt.ScalarFormatter())
    axA.minorticks_off()
    axA.set_xlabel("reference-negative density, n_pairs (log scale)")
    axA.set_ylabel("type-I error rate")
    axA.set_title("Type-I vs density")
    axA.legend(frameon=False, fontsize=7)

    for rule in rules:
        dr = dens[dens["rule"] == rule].sort_values("n_pairs")
        axB.plot(dr["n_pairs"], dr["non_estimable_rate"], marker="o", ms=5,
                 lw=1.6, color=DENS_COLOR.get(rule, "black"), label=rule)
    axB.set_xscale("log")
    axB.set_xticks(xs)
    axB.get_xaxis().set_major_formatter(plt.ScalarFormatter())
    axB.minorticks_off()
    axB.set_ylim(-0.03, 1.0)
    axB.set_xlabel("reference-negative density, n_pairs (log scale)")
    axB.set_ylabel("non-estimable replicate rate")
    axB.set_title("Estimability vs density")
    axB.legend(frameon=False, title="selection rule")

    fig.suptitle(
        f"Post-hoc reference-negative-density sweep "
        f"(base null cell; {PRIMARY_ESTIMATOR})", fontsize=10)
    _save(fig, "sim_4_density_sweep")


def main() -> None:
    q1 = pd.read_parquet(ANALYSIS / "q1_typeI.parquet")
    q2 = pd.read_parquet(ANALYSIS / "q2_power.parquet")
    q3 = pd.read_parquet(ANALYSIS / "q3_coverage.parquet")
    wcr_path = WCR / "wcr_coverage.parquet"
    dens_path = DENSITY / "density_summary.parquet"
    wcr = pd.read_parquet(wcr_path) if wcr_path.exists() else None
    print("[sim-figures] writing to", FIGS)
    figure_typeI_surface(q1)
    figure_power_curve(q2)
    figure_coverage_K(q3, wcr)
    if dens_path.exists():
        figure_density_sweep(pd.read_parquet(dens_path))
    else:
        print("  [skip] density_summary.parquet not found; sim_4 not written")


if __name__ == "__main__":
    main()
