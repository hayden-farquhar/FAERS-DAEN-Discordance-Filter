"""Post-hoc: reference-negative-density / estimability sweep over n_pairs.

Reviewer point (#4): the confirmatory grid pins the population size at the
calibrated ``n_pairs = 492``, where the H1 stratum is thin and non-estimability
is high (~0.88 on ``mde_1.5``). The worry is that the in-band type-I is an
artefact of that "degenerate corner": almost nothing is estimable, and
non-estimable replicates are scored as non-rejections. This sweep moves the test
out of the corner by growing ``n_pairs`` (1x, 2x, 4x, 8x the calibrated value) at
the base cell, holding every other factor fixed, and reports BOTH type-I
estimands (marginal and conditional-on-estimable) plus the non-estimability rate
as the stratum fills in.

Two regimes are swept:

- **Q1_null** (type-I): the check is whether each rule stays inside the
  [0, 0.075] band on BOTH estimands as estimability rises. ``faers_point`` is now
  included alongside ``none`` and ``mde_1.5``: if the contaminated FAERS-derived
  rule's type-I stays in band even at 8x density, the S1 (inflation) hypothesis
  is dead beyond resuscitation, not merely masked by thin-stratum
  non-estimability.
- **Q2_alt** (power): the check is whether the registered ``mde_1.5`` rule ever
  becomes usefully powered as the reference set grows, or whether it is
  permanently dominated on power by the contaminated ``faers_point`` rule. Here
  the rejection rate is power (the enrichment alternative holds), so higher is
  better and there is no band.

Reproducibility: this is one coherent, fully deterministic run (fixed
``SeedSequence`` per regime x n_pairs x stream). Populations are drawn once per
replicate before the rule loop, so every rule is compared on identical
populations; the selection and bootstrap generators are threaded across
replicates and shared by the rules within a level. Because the bootstrap stream
is threaded, adding ``faers_point`` to the rule loop perturbs the later-replicate
bootstrap draws that ``none``/``mde_1.5`` also consume, so the Q1_null
``none``/``mde_1.5`` rates here are Monte-Carlo-consistent with (not bit-identical
to) the earlier two-rule sweep: the 1x/2x rows match exactly and the 4x/8x rows
differ by <= 0.002, inside the stated +/- ~0.006 precision. The manuscript's
density table, figure, and text are regenerated from this run. Q2_alt uses an
independent regime-salted seed stream. Never touches the frozen confirmatory
results in ``results/sim_results.parquet``.

Usage:  python3 scripts/sim_density_sweep.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from tqdm import tqdm

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from simulation.dgp import DGPCell, generate_population         # noqa: E402
from simulation.sweep import load_config, _marginals           # noqa: E402
from simulation.selection import SelectionContext, apply_selection  # noqa: E402
from simulation.estimators import build_h1_stratum, run_estimator   # noqa: E402
from reporting.sim_analysis import TYPEI_BAND, SWEEP_BASE        # noqa: E402

REGIMES = ["Q1_null", "Q2_alt"]
# Collider baseline (`none`), the registered type-I-controlling rule (`mde_1.5`),
# and the winner's-curse-contaminated FAERS-derived rule (`faers_point`).
# `faers_point` is evaluated last (it is the only rule that consumes selection
# RNG, via mc_power); the threaded bootstrap stream is shared across rules, so
# the run is reproduced as a whole rather than rule-by-rule (see module docstring).
RULES = ["none", "mde_1.5", "faers_point"]
# Vectorised WCU (bit-identical to the frozen looped wild bootstrap on the same
# Generator state; reported under its canonical name via res.estimator).
ESTIMATOR = "wild_cluster_bootstrap_fast"
DENSITY_TAG = 0x44454e53   # "DENS"
# Regime salt keeps Q2_alt's seed stream independent of Q1_null. Q1_null uses the
# 3-element seed (no salt) so its already-reported numbers reproduce exactly;
# appending a 4th element would change the SeedSequence, so the salt is applied
# only to the alternative regime, which has no prior published numbers.
REGIME_SALT = {"Q1_null": None, "Q2_alt": 0x514132}  # "QA2"
N_PAIRS_MULTIPLIERS = [1, 2, 4, 8]
REPLICATES = 1000          # post-hoc; type-I precision +/- ~0.006 at p~0.03
OUTDIR = ROOT / "results" / "sim_density"


def _seed(cfg: dict, regime: str, n_pairs: int, stream: int) -> np.random.SeedSequence:
    base = [int(cfg["seed_master"]) ^ DENSITY_TAG, n_pairs, stream]
    salt = REGIME_SALT[regime]
    if salt is not None:
        base.append(salt)
    return np.random.SeedSequence(base)


def _base_cell(cfg: dict, regime: str, n_pairs: int) -> DGPCell:
    s = cfg["primary_slice"]
    return DGPCell(
        regime=regime, pi=SWEEP_BASE["pi"], lambda_inflate=SWEEP_BASE["lambda_inflate"],
        phi=SWEEP_BASE["phi"], rho=s["rho"], K=s["K"], icc=s["icc"],
        OR_true=s["OR_true"], N_F=int(cfg["population"]["N_F"]), n_pairs=n_pairs,
        marginals=_marginals(cfg), q2_daen_suppression=float(cfg["q2_daen_suppression"]),
        label=f"density_{regime}_npairs{n_pairs}",
    )


def _run_level(regime: str, n_pairs: int, *, cfg: dict, replicates: int) -> pd.DataFrame:
    cell = _base_cell(cfg, regime, n_pairs)
    pop_rng = np.random.default_rng(_seed(cfg, regime, n_pairs, 1))
    sel_rng = np.random.default_rng(_seed(cfg, regime, n_pairs, 2))
    boot_rng = np.random.default_rng(_seed(cfg, regime, n_pairs, 3))
    boot_B = int(cfg["bootstrap_B"])
    ctx = SelectionContext(rng=sel_rng, mc_B=int(cfg["mc_power_B"]))

    rows = []
    for rep in range(replicates):
        pop = generate_population(cell, pop_rng)
        for rule in RULES:
            powered = apply_selection(pop, rule, ctx)
            stratum = build_h1_stratum(pop.assign(daen_powered=powered))
            res = run_estimator(stratum, ESTIMATOR, rng=boot_rng, B=boot_B)
            rows.append({
                "regime": regime, "n_pairs": n_pairs, "rule": rule, "replicate": rep,
                "reject": res.reject, "non_estimable": res.non_estimable,
                "stratum_size": res.stratum_size, "n_faers_only": res.n_faers_only,
                "n_clusters": res.n_clusters,
            })
    return pd.DataFrame(rows)


def _summarise(raw: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (regime, n_pairs, rule), x in raw.groupby(
            ["regime", "n_pairs", "rule"], observed=True):
        n = len(x)
        n_reject = int(x["reject"].sum())
        n_nonest = int(x["non_estimable"].sum())
        n_est = n - n_nonest
        rows.append({
            "regime": regime, "rule": rule, "n_pairs": int(n_pairs), "n": n,
            "marginal_reject": n_reject / n,
            "conditional_reject": (n_reject / n_est) if n_est > 0 else np.nan,
            "non_estimable_rate": n_nonest / n,
            "n_estimable": n_est,
            "median_stratum_size": float(x["stratum_size"].median()),
            "median_faers_only": float(x["n_faers_only"].median()),
            "marginal_within_band": (n_reject / n) <= TYPEI_BAND[1],
            "conditional_within_band": (n_est > 0 and (n_reject / n_est) <= TYPEI_BAND[1]),
        })
    return (pd.DataFrame(rows)
            .sort_values(["regime", "rule", "n_pairs"]).reset_index(drop=True))


def _table(L: list, summ: pd.DataFrame, regime: str, value_label: str) -> None:
    sub = summ[summ["regime"] == regime]
    L.append(f"| rule | n_pairs | non-est rate | n est | marginal {value_label} | "
             f"conditional {value_label} | median stratum | median faers-only |")
    L.append("|---|---|---|---|---|---|---|---|")
    for _, r in sub.iterrows():
        L.append(f"| {r['rule']} | {int(r['n_pairs'])} | {r['non_estimable_rate']:.3f} | "
                 f"{int(r['n_estimable'])} | {r['marginal_reject']:.4f} | "
                 f"{r['conditional_reject']:.4f} | {r['median_stratum_size']:.0f} | "
                 f"{r['median_faers_only']:.0f} |")
    L.append("")


def _report(summ: pd.DataFrame, replicates: int) -> str:
    L = ["# Post-hoc: reference-negative-density (n_pairs) sweep", "",
         f"Base cell (pi={SWEEP_BASE['pi']}, lambda={SWEEP_BASE['lambda_inflate']}, "
         f"phi={SWEEP_BASE['phi']}), primary estimator `{ESTIMATOR}`, {replicates} "
         f"replicates per (regime x rule x n_pairs). n_pairs grows from the "
         f"calibrated 492 (1x) to 8x; everything else is held at the slice values. "
         f"Type-I band [{TYPEI_BAND[0]}, {TYPEI_BAND[1]}].", "",
         "Both estimands are shown: `marginal` scores non-estimable replicates as "
         "non-rejections (the pre-registered endpoint); `conditional` divides "
         "rejections by estimable replicates only. As n_pairs grows the stratum "
         "fills in, non-estimability falls, and the two estimands converge.", ""]

    L += ["## Q1_null (type-I)", "",
          "The rejection rate is the type-I error. The question is whether each "
          "rule stays in band as the stratum fills. `faers_point` is the "
          "winner's-curse-contaminated rule: if its type-I stays in band even at "
          "8x density, the inflation hypothesis (S1) is dead, not merely masked by "
          "thin-stratum non-estimability.", ""]
    _table(L, summ, "Q1_null", "type-I")

    L += ["## Q2_alt (power)", "",
          "The enrichment alternative holds, so the rejection rate is power "
          "(higher is better; no band). The question is whether the registered "
          "`mde_1.5` rule ever becomes usefully powered as the reference set "
          "grows, or is permanently dominated on power by the contaminated "
          "`faers_point` rule.", ""]
    _table(L, summ, "Q2_alt", "power")

    # ---- honest-either-way roll-ups (cited verbatim in the manuscript) ----
    q1 = summ[summ["regime"] == "Q1_null"]
    mde1 = q1[q1["rule"] == "mde_1.5"]
    fp1 = q1[q1["rule"] == "faers_point"]
    none1 = q1[q1["rule"] == "none"]
    mde_ok = bool(mde1["conditional_within_band"].all() and mde1["marginal_within_band"].all())
    fp_ok = bool(fp1["conditional_within_band"].all() and fp1["marginal_within_band"].all())

    q2 = summ[summ["regime"] == "Q2_alt"]
    mde2 = q2[q2["rule"] == "mde_1.5"].set_index("n_pairs")
    fp2 = q2[q2["rule"] == "faers_point"].set_index("n_pairs")
    common = sorted(set(mde2.index) & set(fp2.index))
    fp_dominates_mde_power = bool(
        all(fp2.loc[n, "marginal_reject"] >= mde2.loc[n, "marginal_reject"] for n in common))

    L += ["## Roll-up", "",
          f"- `Q1_mde_1.5_in_band_all_densities_both_estimands`: {mde_ok}",
          f"- `Q1_faers_point_in_band_all_densities_both_estimands`: {fp_ok}",
          f"- `Q1_none_conditional_typeI_range`: "
          f"[{none1['conditional_reject'].min():.4f}, {none1['conditional_reject'].max():.4f}]",
          f"- `Q2_mde_1.5_marginal_power_range`: "
          f"[{mde2['marginal_reject'].min():.4f}, {mde2['marginal_reject'].max():.4f}]",
          f"- `Q2_faers_point_marginal_power_range`: "
          f"[{fp2['marginal_reject'].min():.4f}, {fp2['marginal_reject'].max():.4f}]",
          f"- `Q2_faers_point_dominates_mde_1.5_on_power_all_densities`: {fp_dominates_mde_power}",
          ""]
    return "\n".join(L)


def main() -> None:
    cfg = load_config()
    n_pairs_base = int(cfg["population"]["n_pairs"])
    levels = [n_pairs_base * m for m in N_PAIRS_MULTIPLIERS]
    print(f"[density] regimes {REGIMES}; n_pairs levels {levels}; rules {RULES}; "
          f"{REPLICATES} replicates each on {ESTIMATOR}")

    OUTDIR.mkdir(parents=True, exist_ok=True)
    level_dir = OUTDIR / "_levels"
    level_dir.mkdir(exist_ok=True)
    frames = []
    for regime in REGIMES:
        for n in tqdm(levels, desc=f"{regime} n_pairs"):
            part = level_dir / f"{regime}_npairs_{n}.parquet"
            if part.exists():                      # resume: skip completed levels
                frames.append(pd.read_parquet(part))
                continue
            df_lvl = _run_level(regime, n, cfg=cfg, replicates=REPLICATES)
            df_lvl.to_parquet(part, index=False)   # incremental save
            frames.append(df_lvl)
    raw = pd.concat(frames, ignore_index=True)

    raw.to_parquet(OUTDIR / "density_raw.parquet", index=False)
    summ = _summarise(raw)
    summ.to_parquet(OUTDIR / "density_summary.parquet", index=False)
    summ.to_csv(OUTDIR / "density_summary.csv", index=False)
    (OUTDIR / "density_summary.md").write_text(_report(summ, REPLICATES))
    print(f"[density] wrote summary to {OUTDIR}")


if __name__ == "__main__":
    main()
