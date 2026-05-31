"""Phase 2.3: augment results/perpair_arm1.parquet with the daen_powered_* family.

For each pair, compute (using the pinned metrics module + the G2-validated power
module):
  - daen_mde            : smallest OR with exact_power(...) >= 0.80
  - daen_powered        : daen_mde <= 1.5  (primary; protocol Section 6.1)
  - daen_powered_at_2   : daen_mde <= 2.0  (permissive anchor)
  - daen_powered_at_3   : daen_mde <= 3.0  (very permissive)
  - daen_powered_07     : MDE @ 0.70 power <= 1.5  (permissive threshold)
  - daen_powered_09     : MDE @ 0.90 power <= 1.5  (stringent threshold)
  - daen_powered_FAERS_pt : power at OR = FAERS point estimate >= 0.80
  - daen_powered_FAERS_LB : power at OR = FAERS 95% LB >= 0.80
  - daen_powered_unconditional : Poisson-margin daen_mde <= 1.5 (Section 6.4)
  - daen_powered_min    : daen_powered AND daen_powered_unconditional
  - faers_powered       : symmetric on FAERS margins (single OR=1.5 check; used by H4)

Architecture: caches by (n_drug, n_event, N). DAEN side uses early-terminating MDE
search (exact PMF). FAERS side computes only the single OR=1.5 power needed for
faers_powered. Large-margin combos (min(n_drug, n_event) > 10,000) fall back to
Monte-Carlo (B=5000, seed=94) since exact-PMF over a 10k+ element support is slow.
"""
from __future__ import annotations

import math
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from power import exact_power, mc_power, daen_mde_unconditional  # noqa: E402


RESULTS = ROOT / "results"
LARGE_MARGIN_THRESHOLD = 10_000  # use MC above this; exact below
MC_B = 5000
OR_GRID = np.round(np.arange(1.10, 5.00 + 0.025, 0.05), 4)


def _log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def _power_at(n_drug: int, n_event: int, N: int, OR: float) -> float:
    """Route to exact or MC based on PMF-support size."""
    if n_drug <= 0 or n_event <= 0 or N <= 0 or not math.isfinite(OR) or OR <= 0:
        return 0.0
    support_size = min(n_drug, n_event)
    if support_size <= LARGE_MARGIN_THRESHOLD:
        return exact_power(n_drug, n_event, N, OR)
    rng = np.random.default_rng(94)
    return mc_power(n_drug, n_event, N, OR, B=MC_B, rng=rng)


def _mde_search(n_drug: int, n_event: int, N: int, *, power_target: float = 0.80) -> float:
    """Linear scan over OR_GRID with early termination at first OR with power >= target."""
    if n_drug <= 0 or n_event <= 0 or N <= 0:
        return math.inf
    for OR in OR_GRID:
        if _power_at(n_drug, n_event, N, float(OR)) >= power_target:
            return float(OR)
    return math.inf


def augment(in_path: Path, out_path: Path) -> pd.DataFrame:
    _log(f"Loading {in_path.name} ...")
    df = pd.read_parquet(in_path).copy()
    _log(f"  {len(df)} pairs; columns: {len(df.columns)}")

    # Margin helper columns
    df["n_drug_D"] = (df["a_D"] + df["b_D"]).astype(int)
    df["n_event_D"] = (df["a_D"] + df["c_D"]).astype(int)
    df["N_D"] = (df["a_D"] + df["b_D"] + df["c_D"] + df["d_D"]).astype(int)
    df["n_drug_F"] = (df["a_F"] + df["b_F"]).astype(int)
    df["n_event_F"] = (df["a_F"] + df["c_F"]).astype(int)
    df["N_F"] = (df["a_F"] + df["b_F"] + df["c_F"] + df["d_F"]).astype(int)

    # ---------- DAEN MDE at three power targets (early-terminating) ----------
    daen_keys = df[["n_drug_D", "n_event_D", "N_D"]].drop_duplicates()
    _log(f"DAEN unique margin-combos: {len(daen_keys)}")
    _log("Computing DAEN MDE @ powers 0.70 / 0.80 / 0.90 (exact, early-terminated) ...")
    mde_08: dict = {}
    mde_07: dict = {}
    mde_09: dict = {}
    t0 = time.time()
    for i, (_, r) in enumerate(daen_keys.iterrows(), 1):
        nd, ne, N = int(r["n_drug_D"]), int(r["n_event_D"]), int(r["N_D"])
        key = (nd, ne, N)
        mde_08[key] = _mde_search(nd, ne, N, power_target=0.80)
        mde_07[key] = _mde_search(nd, ne, N, power_target=0.70)
        mde_09[key] = _mde_search(nd, ne, N, power_target=0.90)
        if i % 50 == 0:
            _log(f"  {i}/{len(daen_keys)}  ({time.time()-t0:.0f}s)")
    _log(f"  done DAEN MDE in {time.time()-t0:.0f}s")

    # ---------- Poisson-margin sensitivity (MC) ----------
    _log("Computing Poisson-margin MDE (Section 6.4) ...")
    poisson_mde: dict = {}
    t0 = time.time()
    for i, (_, r) in enumerate(daen_keys.iterrows(), 1):
        nd, ne, N = int(r["n_drug_D"]), int(r["n_event_D"]), int(r["N_D"])
        key = (nd, ne, N)
        if nd > 0 and ne > 0 and N > 0:
            rng = np.random.default_rng(94)
            poisson_mde[key] = daen_mde_unconditional(nd, ne, N, power_target=0.80, rng=rng, B=2000)
        else:
            poisson_mde[key] = math.inf
    _log(f"  done Poisson MDE in {time.time()-t0:.0f}s")

    # ---------- FAERS-side: single OR=1.5 power check ----------
    faers_keys = df[["n_drug_F", "n_event_F", "N_F"]].drop_duplicates()
    _log(f"FAERS unique margin-combos: {len(faers_keys)}")
    _log("Computing FAERS power at OR=1.5 (single eval for faers_powered tag) ...")
    faers_pow15: dict = {}
    t0 = time.time()
    for i, (_, r) in enumerate(faers_keys.iterrows(), 1):
        nd, ne, N = int(r["n_drug_F"]), int(r["n_event_F"]), int(r["N_F"])
        key = (nd, ne, N)
        faers_pow15[key] = _power_at(nd, ne, N, 1.5)
        if i % 50 == 0:
            _log(f"  {i}/{len(faers_keys)}  ({time.time()-t0:.0f}s)")
    _log(f"  done FAERS power in {time.time()-t0:.0f}s")

    # ---------- Apply to each row ----------
    _log("Applying tags ...")
    d_keys = list(zip(df["n_drug_D"], df["n_event_D"], df["N_D"]))
    f_keys = list(zip(df["n_drug_F"], df["n_event_F"], df["N_F"]))

    df["daen_mde"] = [mde_08.get(k, math.inf) for k in d_keys]
    df["daen_mde_07"] = [mde_07.get(k, math.inf) for k in d_keys]
    df["daen_mde_09"] = [mde_09.get(k, math.inf) for k in d_keys]
    df["daen_mde_unconditional"] = [poisson_mde.get(k, math.inf) for k in d_keys]
    df["faers_power_at_15"] = [faers_pow15.get(k, 0.0) for k in f_keys]

    # FAERS-derived DAEN power tags
    _log("Computing FAERS-derived DAEN power tags (per-pair) ...")
    daen_pow_FAERS_pt = []
    daen_pow_FAERS_LB = []
    for _, r in df.iterrows():
        nd, ne, N = int(r["n_drug_D"]), int(r["n_event_D"]), int(r["N_D"])
        daen_pow_FAERS_pt.append(_power_at(nd, ne, N, float(r["ror_F"])) >= 0.80)
        daen_pow_FAERS_LB.append(_power_at(nd, ne, N, float(r["ror_F_ci_lo"])) >= 0.80)

    df["daen_powered_FAERS_pt"] = daen_pow_FAERS_pt
    df["daen_powered_FAERS_LB"] = daen_pow_FAERS_LB

    # Threshold tags
    df["daen_powered"] = df["daen_mde"] <= 1.5
    df["daen_powered_at_2"] = df["daen_mde"] <= 2.0
    df["daen_powered_at_3"] = df["daen_mde"] <= 3.0
    df["daen_powered_07"] = df["daen_mde_07"] <= 1.5
    df["daen_powered_09"] = df["daen_mde_09"] <= 1.5
    df["daen_powered_unconditional"] = df["daen_mde_unconditional"] <= 1.5
    df["daen_powered_min"] = df["daen_powered"] & df["daen_powered_unconditional"]
    df["faers_powered"] = df["faers_power_at_15"] >= 0.80

    out_path.parent.mkdir(exist_ok=True)
    df.to_parquet(out_path, index=False)
    _log(f"Wrote {len(df)} pairs with {len(df.columns)} columns -> {out_path}")
    return df


def main() -> int:
    in_path = RESULTS / "perpair_arm1.parquet"
    backup = RESULTS / "perpair_arm1.phase1.parquet"
    if not backup.exists():
        import shutil
        shutil.copy(in_path, backup)
        _log(f"Backed up Phase-1 substrate to {backup.name}")

    t0 = time.time()
    df = augment(in_path, in_path)
    dt = time.time() - t0

    print(f"\nWall time: {dt:.1f}s")
    print()
    print("=== daen_powered_* distributions ===")
    for col in ["daen_powered", "daen_powered_at_2", "daen_powered_at_3",
                "daen_powered_07", "daen_powered_09",
                "daen_powered_FAERS_pt", "daen_powered_FAERS_LB",
                "daen_powered_unconditional", "daen_powered_min", "faers_powered"]:
        n_true = int(df[col].sum())
        print(f"  {col:32s}  TRUE = {n_true:4d}  ({n_true/len(df)*100:5.1f}%)")
    print()
    print("=== daen_mde distribution (primary; exact path) ===")
    finite = df["daen_mde"].replace([math.inf], np.nan).dropna()
    print(f"  finite: {len(finite)}/{len(df)} pairs ({len(finite)/len(df)*100:.1f}%)")
    if len(finite) > 0:
        print(f"  median={finite.median():.2f}  p25={finite.quantile(0.25):.2f}  p75={finite.quantile(0.75):.2f}")
        print(f"  min={finite.min():.2f}  max={finite.max():.2f}")
    print()
    print("=== HEADLINE-FORBIDDEN cross-tab (G3 in effect; this is substrate not result) ===")
    print(df.groupby("cross_db_class", observed=True)["daen_powered"].agg(["sum", "count"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
