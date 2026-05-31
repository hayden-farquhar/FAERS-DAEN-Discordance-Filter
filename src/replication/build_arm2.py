"""Phase 5 Arm-2 discovery substrate.

Universe per protocol Section 8.2: all FAERS drug-event pairs meeting the signal
threshold (LB of 95% CI of ROR > 1 AND raw a >= 3) — across the full FAERS
ingredient × MedDRA-PT space, not restricted to the OMOP+EU-ADR ground-truth
universe.

Algorithm:
  1. Build per-ingredient case-set inverted index (key = prod_ai if non-empty,
     else lowercase-stripped drugname; values = sorted unique primaryid arrays).
  2. Build per-PT case-set inverted index.
  3. Construct sparse matrices: drug_matrix is (n_drugs × n_cases) with 1s where
     a drug appears in a case; pt_matrix is (n_pts × n_cases) similarly.
  4. Co-occurrence: a_matrix = drug_matrix @ pt_matrix.T is the (drug, PT) cell
     `a`. Keep cells with a >= 3 — this is the candidate set.
  5. For each candidate pair, compute the full 2x2 from drug-side and pt-side
     marginals, apply the pre-registered signal rule via the pinned metrics
     module, and tag faers_signal.
  6. For pairs with faers_signal=True, compute the DAEN 2x2 from the DAEN
     substrate (active_ingredient + PT exact match after bullet-stripping +
     lowercasing) and the DAEN signal flag.
  7. For all FAERS-positive Arm-2 pairs, compute daen_mde (exact-PMF MDE search,
     same as Phase 2) and the daen_powered tag.
  8. Persist to results/perpair_arm2.parquet.

This is the substrate for Phase 5 Family 3 (Arm-2 mechanism arm H5a-d).
"""
from __future__ import annotations

import math
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from metrics import ror as _ror  # noqa: E402
from power import exact_power  # noqa: E402

RESULTS = ROOT / "results"
DATA_DAEN = ROOT / "data" / "daen"
FAERS_SUBSTRATE = Path(
    os.environ.get("FAERS_SUBSTRATE", str(ROOT / "data" / "raw" / "faers_substrate"))
)
MIN_CO_REPORT = 3       # raw a >= 3 (protocol signal rule)
OR_GRID = np.round(np.arange(1.10, 5.00 + 0.025, 0.05), 4)


def _log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def _ingredient_key(prod_ai: str, drugname: str) -> str:
    pa = (prod_ai or "").strip().lower()
    if pa:
        return pa
    return (drugname or "").strip().lower()


def load_faers_long() -> tuple[pd.DataFrame, pd.DataFrame, int]:
    """Long-form (primaryid, ingredient) and (primaryid, pt) tables + N_total.

    Filters both drug and reac to the global surviving-primaryid set, so that
    n_drug, n_pt, a, and N are all computed against the same case base.
    """
    _log("FAERS: loading surviving_cases.parquet ...")
    surv = pd.read_parquet(FAERS_SUBSTRATE / "surviving_cases.parquet", columns=["primaryid"])
    surviving = np.unique(surv["primaryid"].astype(np.int64).values)
    N = int(surviving.size)
    _log(f"  N_FAERS (surviving) = {N:,}")

    _log("FAERS: scanning 88 drug-quarter parquets (with surviving-id filter) ...")
    parts = []
    for p in sorted((FAERS_SUBSTRATE / "drug").glob("*.parquet")):
        df = pd.read_parquet(p, columns=["primaryid", "drugname", "prod_ai"])
        df["primaryid"] = pd.to_numeric(df["primaryid"], errors="coerce")
        df = df.dropna(subset=["primaryid"])
        df["primaryid"] = df["primaryid"].astype(np.int64)
        # Filter to surviving cases (must be done BEFORE building indexes so n_drug/n_pt/N match)
        df = df[np.isin(df["primaryid"].values, surviving, assume_unique=False)]
        pa = df["prod_ai"].astype(str).str.lower().str.strip()
        dn = df["drugname"].astype(str).str.lower().str.strip()
        ingr = pa.where(pa != "", dn)
        parts.append(pd.DataFrame({
            "primaryid": df["primaryid"].values,
            "ingredient": ingr.values,
        }))
    drug = pd.concat(parts, ignore_index=True)
    drug = drug[drug["ingredient"] != ""]
    _log(f"  drug rows (post surviving-filter): {len(drug):,}")

    _log("FAERS: scanning 88 reac-quarter parquets (with surviving-id filter) ...")
    parts = []
    for p in sorted((FAERS_SUBSTRATE / "reac").glob("*.parquet")):
        df = pd.read_parquet(p, columns=["primaryid", "pt"])
        df["primaryid"] = pd.to_numeric(df["primaryid"], errors="coerce")
        df = df.dropna(subset=["primaryid"])
        df["primaryid"] = df["primaryid"].astype(np.int64)
        df = df[np.isin(df["primaryid"].values, surviving, assume_unique=False)]
        parts.append(pd.DataFrame({
            "primaryid": df["primaryid"].values,
            "pt": df["pt"].astype(str).str.lower().str.strip().values,
        }))
    reac = pd.concat(parts, ignore_index=True)
    reac = reac[reac["pt"] != ""]
    _log(f"  reac rows (post surviving-filter): {len(reac):,}")

    return drug, reac, N


def build_cooccurrence(drug: pd.DataFrame, reac: pd.DataFrame, N: int,
                       *, min_a: int = MIN_CO_REPORT,
                       min_drug_cases: int = 3,
                       min_pt_cases: int = 3,
                       ) -> pd.DataFrame:
    """Build the candidate (ingredient, PT) pair set with a >= min_a via sparse matmul.

    Returns a DataFrame with columns: ingredient, pt, a, n_drug, n_pt.
    """
    _log("Building unique primaryid → case-index map ...")
    all_ids = np.unique(np.concatenate([drug["primaryid"].values, reac["primaryid"].values]))
    id_to_idx = {pid: i for i, pid in enumerate(all_ids)}
    n_cases_effective = len(all_ids)
    _log(f"  unique primaryids (drug ∪ reac): {n_cases_effective:,}")
    _log(f"  N_total (protocol denominator from surviving_cases.parquet): {N:,}")

    # Map drug.primaryid -> case-index
    _log("Mapping drug rows to case-indexes ...")
    drug["case_idx"] = drug["primaryid"].map(id_to_idx)
    reac["case_idx"] = reac["primaryid"].map(id_to_idx)

    _log("Filtering tiny ingredients and tiny PTs ...")
    drug_counts = drug.groupby("ingredient")["primaryid"].nunique()
    keep_drugs = drug_counts[drug_counts >= min_drug_cases].index
    drug = drug[drug["ingredient"].isin(keep_drugs)]
    pt_counts = reac.groupby("pt")["primaryid"].nunique()
    keep_pts = pt_counts[pt_counts >= min_pt_cases].index
    reac = reac[reac["pt"].isin(keep_pts)]
    _log(f"  drug rows after ingredient-prune: {len(drug):,}  ({drug['ingredient'].nunique():,} ingredients)")
    _log(f"  reac rows after PT-prune:        {len(reac):,}  ({reac['pt'].nunique():,} PTs)")

    # Build sparse boolean matrices
    _log("Building sparse drug-case matrix ...")
    ingr_unique = sorted(drug["ingredient"].unique())
    ingr_to_idx = {g: i for i, g in enumerate(ingr_unique)}
    drug_rows = drug["ingredient"].map(ingr_to_idx).values
    drug_cols = drug["case_idx"].values
    drug_data = np.ones(len(drug), dtype=np.int8)
    drug_mat = sp.csr_matrix(
        (drug_data, (drug_rows, drug_cols)),
        shape=(len(ingr_unique), n_cases_effective),
    )
    # Deduplicate (multiple drug-rows for same case): use sum_duplicates then clip
    drug_mat.sum_duplicates()
    drug_mat.data = np.minimum(drug_mat.data, 1).astype(np.int8)
    _log(f"  drug matrix: {drug_mat.shape}  nnz={drug_mat.nnz:,}")

    _log("Building sparse PT-case matrix ...")
    pt_unique = sorted(reac["pt"].unique())
    pt_to_idx = {p: i for i, p in enumerate(pt_unique)}
    pt_rows = reac["pt"].map(pt_to_idx).values
    pt_cols = reac["case_idx"].values
    pt_data = np.ones(len(reac), dtype=np.int8)
    pt_mat = sp.csr_matrix(
        (pt_data, (pt_rows, pt_cols)),
        shape=(len(pt_unique), n_cases_effective),
    )
    pt_mat.sum_duplicates()
    pt_mat.data = np.minimum(pt_mat.data, 1).astype(np.int8)
    _log(f"  PT matrix: {pt_mat.shape}  nnz={pt_mat.nnz:,}")

    _log("Computing co-occurrence (sparse matmul) ...")
    # a_matrix = drug_mat @ pt_mat.T -- (n_drugs × n_pts) sparse with a-counts
    drug_mat_int = drug_mat.astype(np.int32)
    pt_mat_int = pt_mat.astype(np.int32)
    a_mat = drug_mat_int @ pt_mat_int.T
    _log(f"  co-occurrence matrix: {a_mat.shape}  nnz={a_mat.nnz:,}")

    # Filter to a >= min_a
    _log(f"Filtering co-occurrences to a >= {min_a} ...")
    a_coo = a_mat.tocoo()
    mask = a_coo.data >= min_a
    rows = a_coo.row[mask]
    cols = a_coo.col[mask]
    vals = a_coo.data[mask]
    _log(f"  candidate pairs (a >= {min_a}): {len(vals):,}")

    # n_drug and n_pt margins
    n_drug_per = np.asarray(drug_mat.sum(axis=1)).flatten()  # cases per ingredient
    n_pt_per = np.asarray(pt_mat.sum(axis=1)).flatten()      # cases per PT

    candidates = pd.DataFrame({
        "ingredient": [ingr_unique[r] for r in rows],
        "pt": [pt_unique[c] for c in cols],
        "a": vals.astype(np.int64),
        "n_drug": n_drug_per[rows].astype(np.int64),
        "n_pt": n_pt_per[cols].astype(np.int64),
    })
    _log(f"  candidate frame built: {len(candidates):,} rows")
    return candidates


def apply_faers_signal_rule(candidates: pd.DataFrame, N: int) -> pd.DataFrame:
    """Compute b, c, d, ROR + CI, faers_signal for each candidate pair."""
    _log("Computing FAERS 2x2 + ROR + signal rule ...")
    df = candidates.copy()
    df["b_F"] = df["n_drug"] - df["a"]
    df["c_F"] = df["n_pt"] - df["a"]
    df["d_F"] = N - df["n_drug"] - df["n_pt"] + df["a"]

    # Defensive: drop any rows where the 2x2 is structurally infeasible
    # (this should be impossible after the surviving-primaryid filter, but guard anyway)
    bad = (df["b_F"] < 0) | (df["c_F"] < 0) | (df["d_F"] < 0)
    if bad.any():
        _log(f"  WARNING: dropping {int(bad.sum())} rows with negative 2x2 cells (post-filter sanity)")
        df = df[~bad].reset_index(drop=True)

    # Apply ROR via the pinned metrics module
    rors, lo, hi, sig = [], [], [], []
    for a, b, c, d in zip(df["a"], df["b_F"], df["c_F"], df["d_F"]):
        try:
            r = _ror(int(a), int(b), int(c), int(d))
            rors.append(r.ror); lo.append(r.ci_low); hi.append(r.ci_high); sig.append(bool(r.signal))
        except (ValueError, ZeroDivisionError):
            rors.append(float("nan")); lo.append(float("nan")); hi.append(float("nan")); sig.append(False)
    df["ror_F"] = rors
    df["ror_F_ci_lo"] = lo
    df["ror_F_ci_hi"] = hi
    df["faers_signal"] = sig
    df = df.rename(columns={"a": "a_F", "n_drug": "n_drug_F", "n_pt": "n_event_F"})
    df["N_F"] = N
    _log(f"  FAERS-positive Arm-2 pairs: {df['faers_signal'].sum():,} / {len(df):,}")
    return df


def add_daen_2x2(faers_pos: pd.DataFrame) -> pd.DataFrame:
    """For each FAERS-positive Arm-2 pair, compute the DAEN 2x2 + signal."""
    _log("Loading DAEN substrate ...")
    cases = pd.read_csv(DATA_DAEN / "daen_cases.csv", usecols=["case_number"], low_memory=False)
    cases["case_number"] = pd.to_numeric(cases["case_number"], errors="coerce").dropna()
    N_D = int(cases["case_number"].nunique())
    _log(f"  N_DAEN = {N_D:,}")

    dd = pd.read_csv(DATA_DAEN / "daen_case_drugs.csv",
                      usecols=["case_number", "active_ingredient"], low_memory=False)
    dd["case_number"] = pd.to_numeric(dd["case_number"], errors="coerce")
    dd = dd.dropna(subset=["case_number"])
    dd["ai"] = dd["active_ingredient"].astype(str).str.lower().str.strip()
    dd = dd[dd["ai"] != ""][["case_number", "ai"]]
    drug_idx = {k: np.unique(g["case_number"].astype(np.int64).values)
                for k, g in dd.groupby("ai")}
    _log(f"  DAEN drug-index keys: {len(drug_idx):,}")

    dr = pd.read_csv(DATA_DAEN / "daen_case_reactions.csv",
                      usecols=["case_number", "reaction"], low_memory=False)
    dr["case_number"] = pd.to_numeric(dr["case_number"], errors="coerce")
    dr = dr.dropna(subset=["case_number"])
    dr["pt"] = dr["reaction"].astype(str).str.strip().str.removeprefix("•").str.strip().str.lower()
    dr = dr[dr["pt"] != ""][["case_number", "pt"]]
    pt_idx = {k: np.unique(g["case_number"].astype(np.int64).values)
              for k, g in dr.groupby("pt")}
    _log(f"  DAEN PT-index keys: {len(pt_idx):,}")

    _log("Computing DAEN 2x2 per FAERS-positive Arm-2 pair ...")
    a_D, n_drug_D, n_event_D = [], [], []
    for ingr, pt in zip(faers_pos["ingredient"], faers_pos["pt"]):
        ds = drug_idx.get(ingr, np.empty(0, dtype=np.int64))
        ps = pt_idx.get(pt, np.empty(0, dtype=np.int64))
        if ds.size == 0 or ps.size == 0:
            a_D.append(0)
        else:
            a_D.append(int(np.intersect1d(ds, ps, assume_unique=True).size))
        n_drug_D.append(int(ds.size))
        n_event_D.append(int(ps.size))

    out = faers_pos.copy()
    out["a_D"] = a_D
    out["n_drug_D"] = n_drug_D
    out["n_event_D"] = n_event_D
    out["b_D"] = out["n_drug_D"] - out["a_D"]
    out["c_D"] = out["n_event_D"] - out["a_D"]
    out["d_D"] = N_D - out["n_drug_D"] - out["n_event_D"] + out["a_D"]
    out["N_D"] = N_D

    rors_D, lo_D, hi_D, sig_D = [], [], [], []
    for a, b, c, d in zip(out["a_D"], out["b_D"], out["c_D"], out["d_D"]):
        if min(a, b, c, d) < 0:
            rors_D.append(float("nan")); lo_D.append(float("nan"))
            hi_D.append(float("nan")); sig_D.append(False); continue
        r = _ror(int(a), int(b), int(c), int(d))
        rors_D.append(r.ror); lo_D.append(r.ci_low); hi_D.append(r.ci_high); sig_D.append(bool(r.signal))
    out["ror_D"] = rors_D
    out["ror_D_ci_lo"] = lo_D
    out["ror_D_ci_hi"] = hi_D
    out["daen_signal"] = sig_D

    out["cross_db_class"] = np.select(
        [(out["faers_signal"] & out["daen_signal"]),
         (out["faers_signal"] & ~out["daen_signal"]),
         (~out["faers_signal"] & out["daen_signal"])],
        ["concordant_positive", "faers_only", "daen_only"],
        default="concordant_negative",
    )
    return out


def add_daen_mde(df: pd.DataFrame) -> pd.DataFrame:
    """Compute daen_powered tag per pair via a single power check at OR=1.5
    (the protocol's primary MDE threshold).

    For Family 3 we only need the boolean daen_powered flag, not the full MDE
    search. One PMF eval per unique DAEN margin tuple. With ~260k unique
    tuples this completes in ~5-10 minutes, vs >6 hours for the full 79-grid
    search.

    Sensitivity tags (daen_powered_at_2, daen_powered_at_3) are computed in
    the same pass with one extra eval each.
    """
    _log("Computing daen_powered (single-OR check at 1.5/2.0/3.0) per unique margin tuple ...")
    keys = list(zip(df["n_drug_D"].astype(int), df["n_event_D"].astype(int), df["N_D"].astype(int)))
    unique_keys = list(set(keys))
    _log(f"  unique DAEN margin tuples: {len(unique_keys):,}")

    pow_at_15: dict[tuple[int, int, int], bool] = {}
    pow_at_20: dict[tuple[int, int, int], bool] = {}
    pow_at_30: dict[tuple[int, int, int], bool] = {}

    t0 = time.time()
    for i, key in enumerate(unique_keys, 1):
        nd, ne, N = key
        if nd <= 0 or ne <= 0 or N <= 0:
            pow_at_15[key] = False
            pow_at_20[key] = False
            pow_at_30[key] = False
            continue
        # Skip clearly underpowered combos quickly: if even OR=3 doesn't fire, OR=1.5/2 won't either
        p3 = exact_power(nd, ne, N, 3.0)
        pow_at_30[key] = p3 >= 0.80
        if pow_at_30[key]:
            p2 = exact_power(nd, ne, N, 2.0)
            pow_at_20[key] = p2 >= 0.80
            if pow_at_20[key]:
                p15 = exact_power(nd, ne, N, 1.5)
                pow_at_15[key] = p15 >= 0.80
            else:
                pow_at_15[key] = False
        else:
            pow_at_20[key] = False
            pow_at_15[key] = False
        if i % 5000 == 0:
            _log(f"    {i:,}/{len(unique_keys):,} margins powered (elapsed {time.time()-t0:.0f}s)")
    _log(f"  power-check done in {time.time()-t0:.0f}s")

    df["daen_powered"] = [pow_at_15[k] for k in keys]
    df["daen_powered_at_2"] = [pow_at_20[k] for k in keys]
    df["daen_powered_at_3"] = [pow_at_30[k] for k in keys]
    _log(f"  daen_powered:       {df['daen_powered'].sum():,} / {len(df):,}  ({df['daen_powered'].mean()*100:.1f}%)")
    _log(f"  daen_powered_at_2:  {df['daen_powered_at_2'].sum():,} / {len(df):,}  ({df['daen_powered_at_2'].mean()*100:.1f}%)")
    _log(f"  daen_powered_at_3:  {df['daen_powered_at_3'].sum():,} / {len(df):,}  ({df['daen_powered_at_3'].mean()*100:.1f}%)")
    return df


def build(out_path: Path | str | None = None) -> pd.DataFrame:
    if out_path is None:
        RESULTS.mkdir(exist_ok=True)
        out_path = RESULTS / "perpair_arm2.parquet"
    else:
        out_path = Path(out_path)

    drug, reac, N = load_faers_long()
    candidates = build_cooccurrence(drug, reac, N)
    with_signal = apply_faers_signal_rule(candidates, N)

    # Restrict to FAERS-positive before adding DAEN 2x2 (saves a lot of work)
    faers_pos = with_signal[with_signal["faers_signal"]].copy()
    _log(f"FAERS-positive Arm-2 pairs: {len(faers_pos):,}")

    with_daen = add_daen_2x2(faers_pos)
    with_mde = add_daen_mde(with_daen)

    with_mde.to_parquet(out_path, index=False)
    _log(f"Wrote {len(with_mde):,} FAERS-positive Arm-2 rows to {out_path}")
    return with_mde


if __name__ == "__main__":
    build()
