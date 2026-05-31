"""Phase 1 Arm-1 per-pair 2x2 build (vectorised).

For each drug-event pair in the pooled OMOP+EU-ADR universe, compute the 2x2
contingency table (a, b, c, d) in FAERS and in DAEN, apply the pre-registered
signal rule (LB of 95% CI of ROR > 1 AND raw a >= 3), and persist the per-pair
output to results/perpair_arm1.parquet.

Drug resolution uses the OMOP/EU-ADR ingredient-aliases map plus the US-AU drug
name crosswalk; event resolution uses the pre-deposited OMOP-outcome -> MedDRA-PT
mapping (OSF Amendment 1, deposited 2026-05-30).

Performance: the per-quarter substrate parquets are already deduplicated to
surviving cases (one row per (case, drug-string) and one row per (case, PT)).
We build per-drug-string and per-PT case-set indexes in a single O(N log N)
groupby, then resolve each pair via dict lookups. End-to-end runtime is dominated
by parquet I/O.
"""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

from metrics import ror as _ror

ROOT = Path(__file__).resolve().parents[2]
DATA_REF = ROOT / "data" / "reference"
DATA_DAEN = ROOT / "data" / "daen"
FAERS_SUBSTRATE = Path(
    os.environ.get("FAERS_SUBSTRATE", str(ROOT / "data" / "raw" / "faers_substrate"))
)
RESULTS = ROOT / "results"


def _log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


# --- universe + mappings ---------------------------------------------------------

def load_reference_universe() -> pd.DataFrame:
    omop = pd.read_csv(DATA_REF / "omop_reference_set.csv")
    euadr = pd.read_csv(DATA_REF / "euadr_reference_set.csv")
    omop["ref_set"] = "OMOP"
    euadr["ref_set"] = "EU-ADR"
    df = pd.concat([omop, euadr], ignore_index=True)
    df["exposureName_lc"] = df["exposureName"].astype(str).str.lower().str.strip()
    return df


def load_pt_mapping() -> dict[str, frozenset[str]]:
    df = pd.read_csv(DATA_REF / "omop_outcome_to_meddra_pt_map.csv")
    df["pt_lc"] = df["meddra_pt"].astype(str).str.lower().str.strip()
    return {
        outcome: frozenset(g["pt_lc"].tolist())
        for outcome, g in df.groupby("outcomeName")
    }


def load_drug_search_strings(exposures_lc: list[str]) -> dict[str, frozenset[str]]:
    with open(DATA_REF / "omop_euadr_ingredient_map.json") as f:
        aliases = json.load(f)
    aliases_lc = {k.strip().lower(): [v.strip().lower() for v in vs] for k, vs in aliases.items()}
    cw = pd.read_csv(DATA_REF / "drug_name_crosswalk.csv")
    us_au = dict(zip(cw["us_name"].str.lower(), cw["au_name"].str.lower()))

    out: dict[str, frozenset[str]] = {}
    for e in exposures_lc:
        s = {e}
        for v in aliases_lc.get(e, []):
            if v:
                s.add(v)
        if e in us_au and us_au[e]:
            s.add(us_au[e])
        out[e] = frozenset(s)
    return out


# --- substrate -> per-key case-set index ----------------------------------------

def _build_index(df: pd.DataFrame, key_col: str, id_col: str) -> dict[str, np.ndarray]:
    """Return {key_lc -> sorted unique int64 array of ids}, dropping empty keys."""
    df = df[df[key_col] != ""]
    grouped = df.groupby(key_col, sort=False)[id_col]
    out: dict[str, np.ndarray] = {}
    for k, ids in grouped:
        u = np.unique(ids.values.astype(np.int64))
        out[k] = u
    return out


def _union_caseset(index: dict[str, np.ndarray], keys: frozenset[str]) -> np.ndarray:
    parts = [index[k] for k in keys if k in index]
    if not parts:
        return np.empty(0, dtype=np.int64)
    if len(parts) == 1:
        return parts[0]
    return np.unique(np.concatenate(parts))


# --- FAERS --------------------------------------------------------------------

def load_faers() -> tuple[dict[str, np.ndarray], dict[str, np.ndarray], int]:
    """Returns (drug_index, pt_index, N_total)."""
    _log("FAERS: loading surviving_cases.parquet for N ...")
    surv = pd.read_parquet(FAERS_SUBSTRATE / "surviving_cases.parquet", columns=["primaryid"])
    N = int(surv["primaryid"].nunique())
    _log(f"  N_FAERS = {N:,}")

    _log("FAERS: scanning 88 drug quarter parquets ...")
    drug_parts = []
    for p in sorted((FAERS_SUBSTRATE / "drug").glob("*.parquet")):
        df = pd.read_parquet(p, columns=["primaryid", "drugname", "prod_ai"])
        df["primaryid"] = pd.to_numeric(df["primaryid"], errors="coerce")
        df = df.dropna(subset=["primaryid"])
        # Long-form: emit two normalized strings per row (drugname + prod_ai),
        # collected into a single index keyed by either field.
        dn = pd.DataFrame({
            "primaryid": df["primaryid"].astype(np.int64).values,
            "key": df["drugname"].astype(str).str.lower().str.strip().values,
        })
        pa = pd.DataFrame({
            "primaryid": df["primaryid"].astype(np.int64).values,
            "key": df["prod_ai"].astype(str).str.lower().str.strip().values,
        })
        drug_parts.append(pd.concat([dn, pa], ignore_index=True))
    drug_long = pd.concat(drug_parts, ignore_index=True)
    _log(f"  drug_long rows: {len(drug_long):,}")
    _log("FAERS: building drug-string -> case-set index ...")
    drug_index = _build_index(drug_long, key_col="key", id_col="primaryid")
    _log(f"  drug_index keys: {len(drug_index):,}")
    del drug_parts, drug_long

    _log("FAERS: scanning 88 reac quarter parquets ...")
    reac_parts = []
    for p in sorted((FAERS_SUBSTRATE / "reac").glob("*.parquet")):
        df = pd.read_parquet(p, columns=["primaryid", "pt"])
        df["primaryid"] = pd.to_numeric(df["primaryid"], errors="coerce")
        df = df.dropna(subset=["primaryid"])
        reac_parts.append(pd.DataFrame({
            "primaryid": df["primaryid"].astype(np.int64).values,
            "key": df["pt"].astype(str).str.lower().str.strip().values,
        }))
    reac_long = pd.concat(reac_parts, ignore_index=True)
    _log(f"  reac_long rows: {len(reac_long):,}")
    _log("FAERS: building PT -> case-set index ...")
    pt_index = _build_index(reac_long, key_col="key", id_col="primaryid")
    _log(f"  pt_index keys: {len(pt_index):,}")
    del reac_parts, reac_long

    return drug_index, pt_index, N


# --- DAEN ---------------------------------------------------------------------

def load_daen() -> tuple[dict[str, np.ndarray], dict[str, np.ndarray], int]:
    _log("DAEN: loading cases ...")
    cases = pd.read_csv(DATA_DAEN / "daen_cases.csv", usecols=["case_number"], low_memory=False)
    cases["case_number"] = pd.to_numeric(cases["case_number"], errors="coerce").dropna()
    N = int(cases["case_number"].nunique())
    _log(f"  N_DAEN = {N:,}")

    _log("DAEN: loading case_drugs ...")
    drug = pd.read_csv(DATA_DAEN / "daen_case_drugs.csv",
                        usecols=["case_number", "active_ingredient"], low_memory=False)
    drug["case_number"] = pd.to_numeric(drug["case_number"], errors="coerce")
    drug = drug.dropna(subset=["case_number"])
    drug_long = pd.DataFrame({
        "case_number": drug["case_number"].astype(np.int64).values,
        "key": drug["active_ingredient"].astype(str).str.lower().str.strip().values,
    })
    _log(f"  drug rows: {len(drug_long):,}")
    _log("DAEN: building active-ingredient -> case-set index ...")
    drug_index = _build_index(drug_long, key_col="key", id_col="case_number")
    _log(f"  drug_index keys: {len(drug_index):,}")
    del drug, drug_long

    _log("DAEN: loading case_reactions ...")
    reac = pd.read_csv(DATA_DAEN / "daen_case_reactions.csv",
                        usecols=["case_number", "reaction"], low_memory=False)
    reac["case_number"] = pd.to_numeric(reac["case_number"], errors="coerce")
    reac = reac.dropna(subset=["case_number"])
    rxn_norm = reac["reaction"].astype(str).str.strip()
    # Strip a leading bullet ("•" + optional spaces) before lowercasing.
    rxn_norm = rxn_norm.str.removeprefix("•").str.strip().str.lower()
    reac_long = pd.DataFrame({
        "case_number": reac["case_number"].astype(np.int64).values,
        "key": rxn_norm.values,
    })
    _log(f"  reac rows: {len(reac_long):,}")
    _log("DAEN: building PT -> case-set index ...")
    pt_index = _build_index(reac_long, key_col="key", id_col="case_number")
    _log(f"  pt_index keys: {len(pt_index):,}")

    return drug_index, pt_index, N


# --- 2x2 + signal -------------------------------------------------------------

def _intersect_size(a: np.ndarray, b: np.ndarray) -> int:
    if a.size == 0 or b.size == 0:
        return 0
    # both are already sorted unique int64 arrays from _build_index / _union_caseset
    return int(np.intersect1d(a, b, assume_unique=True).size)


def two_by_two(drug_set: np.ndarray, event_set: np.ndarray, N: int) -> tuple[int, int, int, int]:
    n_drug = int(drug_set.size)
    n_event = int(event_set.size)
    a = _intersect_size(drug_set, event_set)
    b = n_drug - a
    c = n_event - a
    d = N - a - b - c
    return a, b, c, d


def signal_row(a: int, b: int, c: int, d: int) -> dict:
    r = _ror(a, b, c, d)
    return {
        "ror": r.ror,
        "ror_ci_low": r.ci_low,
        "ror_ci_high": r.ci_high,
        "signal": bool(r.signal),
    }


# --- main ---------------------------------------------------------------------

def build(out_path: Path | str | None = None) -> pd.DataFrame:
    if out_path is None:
        RESULTS.mkdir(exist_ok=True)
        out_path = RESULTS / "perpair_arm1.parquet"
    else:
        out_path = Path(out_path)

    _log("Loading reference universe + mappings ...")
    ref = load_reference_universe()
    pt_map = load_pt_mapping()
    _log(f"  reference pairs: {len(ref)}  "
         f"({(ref['ref_set']=='OMOP').sum()} OMOP + {(ref['ref_set']=='EU-ADR').sum()} EU-ADR)")
    _log(f"  unique exposures: {ref['exposureName_lc'].nunique()}; outcomes: {ref['outcomeName'].nunique()}")

    search_for = load_drug_search_strings(ref["exposureName_lc"].unique().tolist())

    faers_drug_idx, faers_pt_idx, N_F = load_faers()
    daen_drug_idx, daen_pt_idx, N_D = load_daen()

    _log("Resolving per-exposure case-sets ...")
    faers_drug_sets: dict[str, np.ndarray] = {}
    daen_drug_sets: dict[str, np.ndarray] = {}
    for exp in ref["exposureName_lc"].unique():
        keys = search_for[exp]
        faers_drug_sets[exp] = _union_caseset(faers_drug_idx, keys)
        daen_drug_sets[exp] = _union_caseset(daen_drug_idx, keys)

    _log("Resolving per-outcome case-sets ...")
    faers_event_sets: dict[str, np.ndarray] = {}
    daen_event_sets: dict[str, np.ndarray] = {}
    for outcome in ref["outcomeName"].unique():
        pts = pt_map[outcome]
        faers_event_sets[outcome] = _union_caseset(faers_pt_idx, pts)
        daen_event_sets[outcome] = _union_caseset(daen_pt_idx, pts)

    _log("Building per-pair 2x2 + applying signal rule ...")
    rows = []
    for _, r in ref.iterrows():
        exp = r["exposureName_lc"]
        out = r["outcomeName"]
        af, bf, cf, df_ = two_by_two(faers_drug_sets[exp], faers_event_sets[out], N_F)
        fsig = signal_row(af, bf, cf, df_)
        ad, bd, cd, dd_ = two_by_two(daen_drug_sets[exp], daen_event_sets[out], N_D)
        dsig = signal_row(ad, bd, cd, dd_)
        rows.append({
            "ref_set": r["ref_set"],
            "exposureName": r["exposureName"],
            "outcomeName": r["outcomeName"],
            "groundTruth": int(r["groundTruth"]),
            "indicationName": r["indicationName"],
            "comparatorName": r["comparatorName"],
            "comparatorType": r["comparatorType"],
            "a_F": af, "b_F": bf, "c_F": cf, "d_F": df_,
            "ror_F": fsig["ror"], "ror_F_ci_lo": fsig["ror_ci_low"], "ror_F_ci_hi": fsig["ror_ci_high"],
            "faers_signal": fsig["signal"],
            "a_D": ad, "b_D": bd, "c_D": cd, "d_D": dd_,
            "ror_D": dsig["ror"], "ror_D_ci_lo": dsig["ror_ci_low"], "ror_D_ci_hi": dsig["ror_ci_high"],
            "daen_signal": dsig["signal"],
        })

    out = pd.DataFrame(rows)
    out["cross_db_class"] = np.select(
        [
            (out["faers_signal"] & out["daen_signal"]),
            (out["faers_signal"] & ~out["daen_signal"]),
            (~out["faers_signal"] & out["daen_signal"]),
        ],
        ["concordant_positive", "faers_only", "daen_only"],
        default="concordant_negative",
    )

    out.to_parquet(out_path, index=False)
    _log(f"Wrote {len(out)} per-pair rows -> {out_path}")
    return out


if __name__ == "__main__":
    build()
