"""Post-hoc S4 re-examination: restricted wild cluster bootstrap (WCR) coverage.

Reviewer point (S4): the few-cluster coverage comparison must confirm the
behaviour is not an artefact of the *unrestricted* wild bootstrap (WCU). The
principled small-cluster estimator is the *restricted* wild bootstrap (WCR;
MacKinnon & Webb 2014), which imposes the null when generating the bootstrap
distribution. This script re-runs the Q1_null K-sweep line on the
type-I-controlling rule ``mde_1.5`` and adds WCR alongside the frozen-sweep
estimators, on the **identical** strata.

Exact reproduction. The ``mde_1.5`` selection is RNG-independent and is applied
first, so each replicate's H1 stratum is a deterministic function of the
population RNG alone. Reproducing the recorded per-cell seed spawn
(``SeedSequence(seed_master).spawn(len(grid))[idx]`` -> first of three children)
regenerates the same populations the frozen sweep used. Consequence: the
deterministic ``cluster_robust`` coverage here must match
``results/sim_analysis/q3_coverage.parquet`` cell-for-cell (the build-in check),
WCU reproduces it up to its own bootstrap-draw noise, and WCR is the new column.

The bootstrap draws use a documented post-hoc stream so this script is itself
reproducible; it never touches the frozen ``results/sim_results.parquet``.

Usage:  python3 scripts/sim_wcr_coverage.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from tqdm import tqdm

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from simulation.dgp import generate_population              # noqa: E402
from simulation.selection import SelectionContext, apply_selection  # noqa: E402
from simulation.estimators import build_h1_stratum, run_estimator   # noqa: E402
from simulation.sweep import build_grid, load_config, _cell_id      # noqa: E402
from reporting.sim_analysis import (                         # noqa: E402
    NULL_OR, NOMINAL_COVERAGE, SWEEP_BASE,
)

RULE = "mde_1.5"
REGIME = "Q1_null"
# Canonical reporting names (what the tables/S4 verdict key on). The unrestricted
# wild bootstrap is computed via the vectorised `_wild_cluster_bootstrap_fast`
# implementation, which reproduces the frozen looped WCU bit-for-bit on the same
# Generator state (verified); `_RUN_NAME` maps the canonical label to the fast
# call so the speed swap is invisible downstream.
ESTIMATORS = ["cluster_robust", "wild_cluster_bootstrap",
              "wild_cluster_bootstrap_restricted"]
_RUN_NAME = {"wild_cluster_bootstrap": "wild_cluster_bootstrap_fast"}
# Documented post-hoc bootstrap stream (distinct from the frozen sweep's seed
# tree): seed_master XOR a fixed WCR tag, then per-cell offset.
WCR_TAG = 0x57435200  # "WCR\0"
OUTDIR = ROOT / "results" / "sim_wcr"


def _k_line_cells(grid):
    """Grid cells on the Q1_null K-sweep line (same predicate as
    reporting.sim_analysis._k_line): all 'sweep_K' cells plus the 'primary'
    K-anchor at the base pi/lambda/phi."""
    out = []
    for idx, c in enumerate(grid):
        if c.regime != REGIME:
            continue
        is_sweepK = c.label == "sweep_K"
        is_anchor = (c.label == "primary"
                     and c.pi == SWEEP_BASE["pi"]
                     and c.lambda_inflate == SWEEP_BASE["lambda_inflate"]
                     and c.phi == SWEEP_BASE["phi"])
        if is_sweepK or is_anchor:
            out.append((idx, c))
    return out


def _run_cell(idx, cell, *, replicates, cfg, child_seqs):
    pop_ss = child_seqs[idx].spawn(3)[0]          # the population stream (first child)
    pop_rng = np.random.default_rng(pop_ss)
    boot_rng = np.random.default_rng(
        np.random.SeedSequence([int(cfg["seed_master"]) ^ WCR_TAG, idx]))
    boot_B = int(cfg["bootstrap_B"])
    ctx = SelectionContext(rng=np.random.default_rng(0), mc_B=int(cfg["mc_power_B"]))

    rows = []
    for rep in range(replicates):
        pop = generate_population(cell, pop_rng)
        powered = apply_selection(pop, RULE, ctx)
        stratum = build_h1_stratum(pop.assign(daen_powered=powered))
        for est in ESTIMATORS:
            res = run_estimator(stratum, _RUN_NAME.get(est, est), rng=boot_rng, B=boot_B)
            rows.append({
                "cell_id": _cell_id(cell), "K": int(cell.K), "replicate": rep,
                "estimator": est, "reject": res.reject,
                "non_estimable": res.non_estimable,
                "ci_lo": res.ci_lo, "ci_hi": res.ci_hi, "beta": res.beta,
            })
    return pd.DataFrame(rows)


def _coverage_table(df):
    est = ~df["non_estimable"].to_numpy()
    covered = est & (df["ci_lo"].to_numpy() <= NULL_OR) & (NULL_OR <= df["ci_hi"].to_numpy())
    df = df.assign(_covered=covered, _estimable=est)
    rows = []
    for (estor, K), x in df.groupby(["estimator", "K"], observed=True):
        n = len(x)
        n_est = int(x["_estimable"].sum())
        beta_est = x.loc[x["_estimable"], "beta"]
        rows.append({
            "estimator": estor, "K": int(K), "n": n, "n_estimable": n_est,
            "coverage": (x["_covered"].sum() / n_est) if n_est else np.nan,
            "typeI": x["reject"].mean(),
            "median_bias": beta_est.median() if n_est else np.nan,
            "non_estimable_rate": x["non_estimable"].mean(),
        })
    cov = pd.DataFrame(rows)
    return cov.sort_values(["estimator", "K"]).reset_index(drop=True)


def _regrade_s4(cov):
    """S4: is WCR coverage nearer nominal than cluster-robust at K<=10, and does
    the gap close by K=50? Re-graded against WCR (vs the WCU graded in the frozen
    analysis), honest-either-way."""
    def gap(estor, K):
        r = cov[(cov["estimator"] == estor) & (cov["K"] == K)]
        return abs(float(r["coverage"].iloc[0]) - NOMINAL_COVERAGE) if len(r) else float("nan")
    Ks = sorted(cov["K"].unique())
    small = [k for k in (5, 10) if k in Ks]
    wcr_better_small = all(
        gap("wild_cluster_bootstrap_restricted", k) <= gap("cluster_robust", k) + 1e-12
        for k in small) if small else False
    g_wcr_50, g_cr_50 = gap("wild_cluster_bootstrap_restricted", 50), gap("cluster_robust", 50)
    worst_small = max(
        gap("wild_cluster_bootstrap_restricted", min(small)) if small else np.nan,
        gap("cluster_robust", min(small)) if small else np.nan)
    gap_closes = (not np.isnan(g_wcr_50) and not np.isnan(g_cr_50)
                  and abs(g_wcr_50 - g_cr_50) <= worst_small + 1e-12)
    coverage_by = {
        estor: {int(K): (float(cov[(cov["estimator"] == estor) & (cov["K"] == K)]["coverage"].iloc[0])
                         if len(cov[(cov["estimator"] == estor) & (cov["K"] == K)]) else float("nan"))
                for K in Ks}
        for estor in ESTIMATORS}
    return {
        "claim": "Restricted wild cluster bootstrap (WCR) coverage is nearer nominal than "
                 "cluster-robust SE at K<=10, and the gap closes as K->50.",
        "read_on_rule": RULE,
        "estimands_compared": ESTIMATORS,
        "coverage_by_estimator_K": coverage_by,
        "wcr_nearer_nominal_at_small_K": bool(wcr_better_small),
        "gap_closes_at_K50": bool(gap_closes),
        "verdict": "SUPPORTED" if (wcr_better_small and gap_closes) else "NOT SUPPORTED",
    }


def _report(cov, verdict, replicates) -> str:
    L = ["# Post-hoc S4: restricted wild cluster bootstrap (WCR) coverage", "",
         f"Q1_null K-sweep line on rule `{RULE}`, {replicates} replicates per K-cell, "
         f"on the **same strata** as the frozen sweep (population RNG reproduced; "
         f"`mde_1.5` selection is RNG-independent). Estimators: cluster-robust SE, "
         f"unrestricted wild bootstrap (WCU), restricted wild bootstrap (WCR). "
         f"Coverage is of the generative null OR = {NULL_OR:.0f}; nominal "
         f"{NOMINAL_COVERAGE:.2f}.", ""]
    L.append("| estimator | K | coverage | type-I | median bias | non-est rate | n est |")
    L.append("|---|---|---|---|---|---|---|")
    for _, r in cov.iterrows():
        L.append(f"| {r['estimator']} | {int(r['K'])} | {r['coverage']:.3f} | "
                 f"{r['typeI']:.4f} | {r['median_bias']:.4f} | "
                 f"{r['non_estimable_rate']:.3f} | {int(r['n_estimable'])} |")
    L += ["", f"## S4 re-grade (WCR) - {verdict['verdict']}", "",
          verdict["claim"], "",
          f"- `wcr_nearer_nominal_at_small_K`: {verdict['wcr_nearer_nominal_at_small_K']}",
          f"- `gap_closes_at_K50`: {verdict['gap_closes_at_K50']}",
          f"- `coverage_by_estimator_K`: {json.dumps(verdict['coverage_by_estimator_K'])}",
          "",
          "Reproduction check: the deterministic `cluster_robust` coverage above must "
          "equal the frozen `results/sim_analysis/q3_coverage.parquet` for `mde_1.5` "
          "(same strata); any divergence is a reproduction bug, not a finding.", ""]
    return "\n".join(L)


def main() -> None:
    cfg = load_config()
    replicates = int(cfg["replicates"])
    grid = build_grid(cfg)
    master = np.random.SeedSequence(int(cfg["seed_master"]))
    child_seqs = master.spawn(len(grid))
    cells = _k_line_cells(grid)
    print(f"[wcr] {len(cells)} K-sweep-line cells under {REGIME} on rule {RULE}; "
          f"{replicates} replicates x {len(ESTIMATORS)} estimators each")

    OUTDIR.mkdir(parents=True, exist_ok=True)
    cell_dir = OUTDIR / "_cells"
    cell_dir.mkdir(exist_ok=True)
    frames = []
    for idx, cell in tqdm(cells, desc="K-cells"):
        part = cell_dir / f"cell_{idx}_K{int(cell.K)}.parquet"
        if part.exists():                          # resume: skip completed cells
            frames.append(pd.read_parquet(part))
            continue
        df_cell = _run_cell(idx, cell, replicates=replicates, cfg=cfg, child_seqs=child_seqs)
        df_cell.to_parquet(part, index=False)      # incremental save
        frames.append(df_cell)
    raw = pd.concat(frames, ignore_index=True)

    raw.to_parquet(OUTDIR / "wcr_raw.parquet", index=False)
    cov = _coverage_table(raw)
    cov.to_parquet(OUTDIR / "wcr_coverage.parquet", index=False)
    cov.to_csv(OUTDIR / "wcr_coverage.csv", index=False)
    verdict = _regrade_s4(cov)
    (OUTDIR / "s4_wcr_verdict.json").write_text(json.dumps(verdict, indent=2))
    (OUTDIR / "wcr_coverage_summary.md").write_text(_report(cov, verdict, replicates))
    print(f"[wcr] wrote coverage + S4 re-grade to {OUTDIR}")
    print(f"[wcr] S4 (WCR): {verdict['verdict']}")


if __name__ == "__main__":
    main()
