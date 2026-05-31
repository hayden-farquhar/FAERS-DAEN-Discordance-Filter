"""Per-pair H5 mechanism feature computation (Arm-2 substrate).

For each FAERS-positive Arm-2 pair (and only those that are also daen_powered for
H5a-d cross-pair comparison), compute:

  - H5a: consumer_share_F = proportion of cases with occp_cod in {CN, CN-LWY}
         (Consumer / Consumer-Lawyer; FDA DEMO field codes)
  - H5b: lawyer_share_F   = proportion of cases with occp_cod in {LW}
  - H5c: time_breakpoint_F (BIC-selected single break on quarterly case counts;
         only if >= 16 quarters per protocol Section 10.3); aligned_alert
         (boolean: break within +/- 2 quarters of any alert for this drug in the
         alert registry)
  - H5d: mass_tort_drug (boolean)

The drug-class baseline for the H5a/b enrichment threshold is computed at the
drug level (proportion of consumer/lawyer reports across ALL cases of the drug,
not pair-specific).

FAERS occp_cod values per FDA DEMO docs:
  MD  = Physician (medical doctor)
  PH  = Pharmacist
  OT  = Other health professional
  LW  = Lawyer
  CN  = Consumer
  (blank / null / variants)
"""
from __future__ import annotations

import math
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import gammaln

ROOT = Path(__file__).resolve().parents[2]
DATA_REF = ROOT / "data" / "reference"
FAERS_SUBSTRATE = Path(
    os.environ.get("FAERS_SUBSTRATE", str(ROOT / "data" / "raw" / "faers_substrate"))
)

_CONSUMER_CODES = {"cn", "consumer"}
_LAWYER_CODES = {"lw", "lawyer"}
_MIN_QUARTERS_H5C = 16
_BREAK_ALIGNMENT_QUARTERS = 2


def _log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] [features] {msg}", flush=True)


def _ingredient_key(prod_ai: str, drugname: str) -> str:
    pa = (prod_ai or "").strip().lower()
    return pa if pa else (drugname or "").strip().lower()


def _quarter_to_index(q: str) -> int:
    """Convert 'YYYYQn' string to integer index (years from 2004Q1 baseline * 4 + quarter offset)."""
    if not isinstance(q, str) or len(q) < 6:
        return -1
    try:
        y = int(q[:4]); qn = int(q[5])
        return (y - 2004) * 4 + (qn - 1)
    except (ValueError, IndexError):
        return -1


# --- BIC-selected single-break test on Poisson quarterly counts ---

def _loglik_poisson(counts: np.ndarray, lam: float) -> float:
    if lam <= 0:
        return -np.inf
    counts = counts.astype(np.float64)
    return float((counts * np.log(lam) - lam - gammaln(counts + 1)).sum())


def _fit_one_break(counts: np.ndarray, k: int) -> float:
    """counts[:k] ~ Pois(lam1); counts[k:] ~ Pois(lam2)."""
    if k < 1 or k >= len(counts):
        return -np.inf
    left = counts[:k]; right = counts[k:]
    if left.sum() == 0 or right.sum() == 0:
        return -np.inf
    return _loglik_poisson(left, left.mean()) + _loglik_poisson(right, right.mean())


def find_break_bic(counts: np.ndarray, *, min_segment: int = 4) -> tuple[int | None, bool]:
    """Returns (best_k, has_break) per BIC selection.

    has_break = True iff BIC of one-break model strictly less than BIC of no-break.
    """
    n = len(counts)
    if n < 2 * min_segment:
        return (None, False)
    lam_null = max(counts.mean(), 1e-9)
    ll_null = _loglik_poisson(counts, lam_null)
    best_k = None; best_ll = -np.inf
    for k in range(min_segment, n - min_segment + 1):
        ll = _fit_one_break(counts, k)
        if ll > best_ll:
            best_ll = ll; best_k = k
    if best_k is None:
        return (None, False)
    bic_null = -2 * ll_null
    bic_alt = -2 * best_ll + np.log(n)  # one extra parameter
    has_break = bic_alt < bic_null
    return (best_k, has_break)


# --- Per-case feature lookup from FAERS demo (occp_cod, quarter) ---

def load_case_features() -> pd.DataFrame:
    """Returns per-case (primaryid, occp_cod_lc, quarter)."""
    _log("Loading per-case occp_cod + quarter from FAERS demo parquets ...")
    parts = []
    for p in sorted((FAERS_SUBSTRATE / "demo").glob("*.parquet")):
        df = pd.read_parquet(p, columns=["primaryid", "occp_cod", "quarter"])
        df["primaryid"] = pd.to_numeric(df["primaryid"], errors="coerce")
        df = df.dropna(subset=["primaryid"])
        df["primaryid"] = df["primaryid"].astype(np.int64)
        df["occp_lc"] = df["occp_cod"].astype(str).str.lower().str.strip()
        parts.append(df[["primaryid", "occp_lc", "quarter"]])
    out = pd.concat(parts, ignore_index=True)
    # primaryid -> last seen value (dedup; surviving_cases ensures one per case)
    out = out.drop_duplicates(subset="primaryid", keep="last").set_index("primaryid")
    _log(f"  per-case feature table: {len(out):,} cases")
    return out


# --- Drug-class baseline for H5a/b (drug-level consumer/lawyer share across all cases) ---

def compute_drug_baselines(drug_long: pd.DataFrame,
                            case_features: pd.DataFrame) -> dict[str, dict[str, float]]:
    """Per-ingredient baselines: {ingredient: {'consumer': p, 'lawyer': p, 'n': N}}."""
    _log("Computing per-ingredient consumer + lawyer baselines ...")
    # drug_long: (primaryid, ingredient). Join with case_features for occp.
    merged = drug_long.drop_duplicates(["primaryid", "ingredient"])
    merged = merged.join(case_features[["occp_lc"]], on="primaryid", how="left")
    merged["is_consumer"] = merged["occp_lc"].isin(_CONSUMER_CODES).astype(int)
    merged["is_lawyer"] = merged["occp_lc"].isin(_LAWYER_CODES).astype(int)
    g = merged.groupby("ingredient").agg(
        n=("primaryid", "nunique"),
        consumer=("is_consumer", "mean"),
        lawyer=("is_lawyer", "mean"),
    )
    out = {row.Index: {"n": int(row.n), "consumer": float(row.consumer), "lawyer": float(row.lawyer)}
            for row in g.itertuples()}
    _log(f"  baselines for {len(out):,} ingredients")
    return out


# --- Per-pair feature loop ---

def build_arm2_features(perpair_arm2: pd.DataFrame,
                          drug_long: pd.DataFrame,
                          reac_long: pd.DataFrame,
                          alert_registry: pd.DataFrame,
                          mass_tort_drugs: pd.DataFrame,
                          ) -> pd.DataFrame:
    """For each Arm-2 pair, compute the four mechanism features.

    Inputs:
      perpair_arm2  : pairs to enrich (must have 'ingredient' and 'pt' columns)
      drug_long     : (primaryid, ingredient) FAERS long table
      reac_long     : (primaryid, pt) FAERS long table
      alert_registry: {alert_date, drug_ingredient, event_category, ...}
      mass_tort_drugs: {drug_ingredient, mdl_number, ...}
    """
    df = perpair_arm2.copy()
    # Build case-features lookup
    case_feat = load_case_features()
    # Drug baselines
    baselines = compute_drug_baselines(drug_long, case_feat)
    # Mass-tort drug set (lowercased)
    mt_drugs = set(mass_tort_drugs["drug_ingredient"].astype(str).str.lower())
    # Alert registry per drug (list of (quarter_index, event_category))
    alert_registry["alert_date"] = pd.to_datetime(alert_registry["alert_date"])
    alert_registry["alert_q"] = (
        (alert_registry["alert_date"].dt.year - 2004) * 4
        + (alert_registry["alert_date"].dt.quarter - 1)
    )
    alerts_by_drug = (
        alert_registry.groupby(alert_registry["drug_ingredient"].str.lower())["alert_q"]
        .apply(list).to_dict()
    )
    _log(f"  alert registry: {len(alert_registry)} entries; {len(alerts_by_drug)} drugs covered")
    _log(f"  mass-tort drugs: {len(mt_drugs)}")

    # Build drug-pt case-set lookup (slow but per-pair we need the case set)
    _log("Building per-(ingredient,pt) case-set lookup for cohort-restricted features ...")
    drug_index: dict[str, np.ndarray] = {}
    for k, g in drug_long.groupby("ingredient"):
        drug_index[k] = np.unique(g["primaryid"].astype(np.int64).values)
    pt_index: dict[str, np.ndarray] = {}
    for k, g in reac_long.groupby("pt"):
        pt_index[k] = np.unique(g["primaryid"].astype(np.int64).values)
    _log(f"  drug_index: {len(drug_index):,} keys; pt_index: {len(pt_index):,} keys")

    # Quarterly time-series per ingredient — for H5c we need the drug-pt cohort
    # specifically (not the whole drug); build per-pair quarterly counts inline.
    _log("Computing per-pair features (RPSR + breakpoint + mass-tort) ...")

    n_total = len(df)
    consumer_share = np.full(n_total, np.nan)
    lawyer_share = np.full(n_total, np.nan)
    baseline_cons = np.full(n_total, np.nan)
    baseline_lawy = np.full(n_total, np.nan)
    n_cases_pair = np.zeros(n_total, dtype=np.int64)
    n_quarters_with_data = np.zeros(n_total, dtype=np.int64)
    break_k = np.full(n_total, -1, dtype=np.int64)
    has_break = np.zeros(n_total, dtype=bool)
    aligned_alert = np.zeros(n_total, dtype=bool)
    mass_tort_flag = np.zeros(n_total, dtype=bool)

    t0 = time.time()
    for i, (ingr, pt) in enumerate(zip(df["ingredient"], df["pt"])):
        ds = drug_index.get(ingr, np.empty(0, dtype=np.int64))
        ps = pt_index.get(pt, np.empty(0, dtype=np.int64))
        if ds.size == 0 or ps.size == 0:
            continue
        cohort = np.intersect1d(ds, ps, assume_unique=True)
        n = cohort.size
        n_cases_pair[i] = n
        if n == 0:
            continue

        # H5a/b: cohort RPSR
        cohort_idx = case_feat.index.intersection(cohort)
        cohort_feat = case_feat.loc[cohort_idx]
        if len(cohort_feat):
            cs = cohort_feat["occp_lc"].isin(_CONSUMER_CODES).mean()
            ls = cohort_feat["occp_lc"].isin(_LAWYER_CODES).mean()
            consumer_share[i] = cs
            lawyer_share[i] = ls
        base = baselines.get(ingr, {"consumer": np.nan, "lawyer": np.nan})
        baseline_cons[i] = base["consumer"]
        baseline_lawy[i] = base["lawyer"]

        # H5c: quarterly counts
        if len(cohort_feat) >= _MIN_QUARTERS_H5C:
            q_indices = cohort_feat["quarter"].map(_quarter_to_index)
            q_indices = q_indices[q_indices >= 0]
            if len(q_indices):
                # Quarterly count vector across 2004Q1 -> 2025Q4 (88 quarters)
                counts = np.bincount(q_indices.values, minlength=88)[:88]
                n_q_with_data = int((counts > 0).sum())
                n_quarters_with_data[i] = n_q_with_data
                if n_q_with_data >= _MIN_QUARTERS_H5C:
                    k, hb = find_break_bic(counts)
                    if k is not None:
                        break_k[i] = k
                        has_break[i] = hb
                        if hb:
                            drug_alerts = alerts_by_drug.get(ingr, [])
                            for aq in drug_alerts:
                                if abs(k - aq) <= _BREAK_ALIGNMENT_QUARTERS:
                                    aligned_alert[i] = True
                                    break

        # H5d: mass-tort drug membership
        mass_tort_flag[i] = ingr in mt_drugs

        if (i + 1) % 500 == 0:
            _log(f"  {i+1}/{n_total}  ({time.time()-t0:.0f}s)")

    df["consumer_share_F"] = consumer_share
    df["lawyer_share_F"] = lawyer_share
    df["baseline_consumer_share_F"] = baseline_cons
    df["baseline_lawyer_share_F"] = baseline_lawy
    df["consumer_lift_pp"] = (df["consumer_share_F"] - df["baseline_consumer_share_F"]) * 100
    df["lawyer_lift_pp"] = (df["lawyer_share_F"] - df["baseline_lawyer_share_F"]) * 100
    df["n_cases_pair_F"] = n_cases_pair
    df["n_quarters_with_data"] = n_quarters_with_data
    df["break_k"] = break_k
    df["has_break"] = has_break
    df["aligned_alert_F"] = aligned_alert
    df["mass_tort_drug"] = mass_tort_flag

    return df
