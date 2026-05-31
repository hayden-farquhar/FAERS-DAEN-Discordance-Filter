"""Phase 6 sensitivity — event-granularity substitute for the HLT arm.

WHY A SUBSTITUTE (documented deviation from Section 10.7).
The pre-registered HLT-level aggregation sensitivity requires the licensed
MedDRA v27.1 PT->HLT hierarchy, which is not deposited in the repository (and
fabricating it would breach the protocol's ontology-pin commitments). Moreover,
the Arm-1 events are OMOP/EU-ADR *baskets* of 3-41 PTs that already span their
HLTs, so a literal HLT roll-up in Arm 1 is close to a no-op by construction.

This exhibit instead probes the SAME underlying concern -- "does the cross-
database classification depend on event-coding granularity?" -- in the only
direction executable without licensed data: it makes the event definition
*finer*, not coarser. Each Arm-1 event is re-defined as its single INDEX PT
(the most-reported PT in FAERS within the basket; deterministic, no subjective
concept-matching), the 2x2 tables are rebuilt in both databases, the G2-locked
power model is re-run on the finer marginals (so daen_powered is recomputed,
not borrowed), and Family 1 (H2) and Family 2 (H1) are re-estimated. Reported
for H2 and H1 per Section 10.7.

Outputs: results/perpair_arm1_indexpt.parquet, results/phase6_granularity_substitute.{json,md}
"""
from __future__ import annotations

import json
import sys
import time
from dataclasses import asdict
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from replication.build_arm1 import (  # noqa: E402
    load_reference_universe,
    load_pt_mapping,
    load_drug_search_strings,
    load_faers,
    load_daen,
    two_by_two,
    signal_row,
    _union_caseset,
)
from enrichment.family1 import family1_headline  # noqa: E402
from enrichment.family2 import family2_h1  # noqa: E402

RESULTS = ROOT / "results"


def _log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def pick_index_pts(
    pt_map: dict[str, frozenset[str]], faers_pt_idx: dict[str, np.ndarray]
) -> dict[str, str]:
    """For each outcome basket, the index PT = the basket PT with the largest
    FAERS case-set. Deterministic; ties broken by lexical order of the PT."""
    chosen: dict[str, str] = {}
    for outcome, pts in pt_map.items():
        best_pt, best_n = None, -1
        for pt in sorted(pts):  # sorted -> deterministic tie-break
            n = int(faers_pt_idx[pt].size) if pt in faers_pt_idx else 0
            if n > best_n:
                best_pt, best_n = pt, n
        chosen[outcome] = best_pt
    return chosen


def build_index_pt_substrate() -> pd.DataFrame:
    ref = load_reference_universe()
    pt_map = load_pt_mapping()
    search_for = load_drug_search_strings(ref["exposureName_lc"].unique().tolist())

    faers_drug_idx, faers_pt_idx, N_F = load_faers()
    daen_drug_idx, daen_pt_idx, N_D = load_daen()

    _log("Selecting index PT per outcome (max FAERS case-set within basket) ...")
    index_pt = pick_index_pts(pt_map, faers_pt_idx)
    for oc, pt in index_pt.items():
        nF = int(faers_pt_idx[pt].size) if pt in faers_pt_idx else 0
        nD = int(daen_pt_idx[pt].size) if pt in daen_pt_idx else 0
        _log(f"  {oc!r:55s} -> index PT {pt!r}  (FAERS n={nF}, DAEN n={nD})")

    _log("Resolving per-exposure case-sets ...")
    faers_drug_sets, daen_drug_sets = {}, {}
    for exp in ref["exposureName_lc"].unique():
        keys = search_for[exp]
        faers_drug_sets[exp] = _union_caseset(faers_drug_idx, keys)
        daen_drug_sets[exp] = _union_caseset(daen_drug_idx, keys)

    _log("Resolving per-outcome SINGLE-INDEX-PT case-sets ...")
    faers_event_sets, daen_event_sets = {}, {}
    for outcome in ref["outcomeName"].unique():
        pt = index_pt[outcome]
        faers_event_sets[outcome] = faers_pt_idx.get(pt, np.empty(0, dtype=np.int64))
        daen_event_sets[outcome] = daen_pt_idx.get(pt, np.empty(0, dtype=np.int64))

    _log("Building per-pair 2x2 (index-PT events) + signal rule ...")
    rows = []
    for _, r in ref.iterrows():
        exp, out = r["exposureName_lc"], r["outcomeName"]
        af, bf, cf, df_ = two_by_two(faers_drug_sets[exp], faers_event_sets[out], N_F)
        fsig = signal_row(af, bf, cf, df_)
        ad, bd, cd, dd_ = two_by_two(daen_drug_sets[exp], daen_event_sets[out], N_D)
        dsig = signal_row(ad, bd, cd, dd_)
        rows.append({
            "ref_set": r["ref_set"], "exposureName": r["exposureName"],
            "outcomeName": r["outcomeName"], "groundTruth": int(r["groundTruth"]),
            "indicationName": r["indicationName"], "comparatorName": r["comparatorName"],
            "comparatorType": r["comparatorType"], "index_pt": index_pt[out],
            "a_F": af, "b_F": bf, "c_F": cf, "d_F": df_,
            "ror_F": fsig["ror"], "ror_F_ci_lo": fsig["ror_ci_low"],
            "ror_F_ci_hi": fsig["ror_ci_high"], "faers_signal": fsig["signal"],
            "a_D": ad, "b_D": bd, "c_D": cd, "d_D": dd_,
            "ror_D": dsig["ror"], "ror_D_ci_lo": dsig["ror_ci_low"],
            "ror_D_ci_hi": dsig["ror_ci_high"], "daen_signal": dsig["signal"],
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
    return out


def main() -> None:
    t0 = time.time()
    substrate = build_index_pt_substrate()
    raw_path = RESULTS / "perpair_arm1_indexpt.raw.parquet"
    substrate.to_parquet(raw_path, index=False)
    _log(f"Built {len(substrate)} index-PT pairs ({time.time()-t0:.0f}s)")

    # Recompute the power model on the finer marginals (reuse the G2-locked augment).
    from run_power_augment import augment  # noqa: E402
    aug_path = RESULTS / "perpair_arm1_indexpt.parquet"
    _log("Recomputing daen_powered on finer (index-PT) marginals ...")
    aug = augment(raw_path, aug_path)

    # Re-run H2 (Family 1) and H1 (Family 2) on the finer-event substrate.
    f1 = family1_headline(aug)
    h1 = family2_h1(aug)

    # Compare against the primary (basket) substrate.
    primary = pd.read_parquet(RESULTS / "perpair_arm1.parquet")
    out = {
        "deviation_note": (
            "Substitute for the Section 10.7 HLT arm: licensed MedDRA v27.1 PT->HLT "
            "hierarchy not deposited; Arm-1 events are supra-PT baskets making literal "
            "HLT roll-up a near-no-op. This exhibit instead makes events FINER (single "
            "index PT = max-FAERS-count PT per basket) and re-runs H2 + H1 with the "
            "power model recomputed on the finer marginals."
        ),
        "index_pt_per_outcome": dict(zip(substrate["outcomeName"], substrate["index_pt"]))
        if "index_pt" in substrate else {},
        "n_daen_powered_indexpt": int(aug["daen_powered"].sum()),
        "n_daen_powered_primary": int(primary["daen_powered"].sum()),
        "family1_indexpt": asdict(f1),
        "h1_indexpt": asdict(h1),
    }
    (RESULTS / "phase6_granularity_substitute.json").write_text(json.dumps(out, indent=2, default=str))

    md = []
    md.append("# Phase 6 sensitivity — event-granularity substitute (index-PT)\n")
    md.append("**Documented deviation from the Section 10.7 HLT arm.** The literal HLT-level "
              "aggregation requires the licensed MedDRA v27.1 PT→HLT hierarchy, which is not "
              "deposited (fabricating it would breach the protocol ontology pins); and the Arm-1 "
              "events are OMOP/EU-ADR baskets of 3–41 PTs that already span their HLTs, so a "
              "literal roll-up is a near-no-op in Arm 1. This exhibit probes the same "
              "granularity-robustness concern in the executable direction: it makes events "
              "**finer** — each event re-defined as its single index PT (max-FAERS-count PT in "
              "the basket) — and re-runs H2 + H1 with the G2-locked power model **recomputed** on "
              "the finer marginals. Reported for H2, H1 per Section 10.7.\n")
    md.append("\n## Index PT chosen per outcome\n")
    md.append("| Outcome (basket) | Index PT (finer event) |")
    md.append("|---|---|")
    for oc, pt in zip(substrate["outcomeName"].drop_duplicates(),
                      substrate.drop_duplicates("outcomeName")["index_pt"]):
        md.append(f"| {oc} | {pt} |")
    md.append("")
    md.append(f"\n**daen_powered universe:** index-PT n={out['n_daen_powered_indexpt']} "
              f"vs primary basket n={out['n_daen_powered_primary']} "
              f"(finer events reduce DAEN event-side counts, so fewer pairs are powered).\n")
    md.append("\n## H2 — power-conditioned non-replication rate (finer events)\n")
    md.append(f"- H2 (index-PT) = **{f1.h2_non_replication_rate:.3f}** "
              f"(95% cluster-bootstrap CI {f1.h2_ci_low:.3f}–{f1.h2_ci_high:.3f}) "
              f"on n={f1.n_faers_positive_daen_powered} FAERS-positive ∧ daen_powered pairs")
    md.append(f"- (primary basket H2 = 0.136 [0.057–0.318] on n=44)")
    md.append(f"- Underpowered-discordant: {f1.h2_underpowered_discordant_count}")
    md.append(f"- PA+ = {f1.pa_plus:.3f} (95% CI {f1.pa_plus_ci_low:.3f}–{f1.pa_plus_ci_high:.3f})")
    md.append(f"- Cohen's κ = {f1.cohen_kappa:.3f}")
    md.append(f"- Cross-DB 2×2 (daen_powered): {f1.cross_db_2x2}\n")
    md.append("\n## H1 — reference-negative enrichment (finer events)\n")
    md.append(f"- Decision: **{h1.decision}**")
    md.append(f"- Stratum: {h1.stratum_size} pairs ({h1.n_faers_only} faers_only, "
              f"{h1.n_concordant_positive} concordant_positive); known-neg "
              f"{h1.n_known_neg_faers_only}/{h1.n_faers_only} vs "
              f"{h1.n_known_neg_concordant_pos}/{h1.n_concordant_positive}")
    md.append(f"- Realised MDE = {h1.realised_mde_or:.2f}")
    for n in h1.notes:
        md.append(f"  - {n}")
    (RESULTS / "phase6_granularity_substitute.md").write_text("\n".join(md) + "\n")

    print("\n=== SUMMARY ===")
    print(f"H2(index-PT) = {f1.h2_non_replication_rate:.3f} "
          f"[{f1.h2_ci_low:.3f}, {f1.h2_ci_high:.3f}] on n={f1.n_faers_positive_daen_powered}")
    print(f"daen_powered: index-PT={out['n_daen_powered_indexpt']} vs primary={out['n_daen_powered_primary']}")
    print(f"H1 decision = {h1.decision} (MDE={h1.realised_mde_or:.2f})")
    print(f"total wall: {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
