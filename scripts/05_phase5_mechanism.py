"""Phase 5 runner: build Family 3 mechanism features on the Arm-2 substrate +
run the H5a-d Holm-Bonferroni tests.

Prerequisite: results/perpair_arm2.parquet exists (run scripts/run_arm2_build.py first).
"""
from __future__ import annotations

import json
import math
import sys
import time
from dataclasses import asdict
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from mechanism import build_arm2_features, family3_holm  # noqa: E402
from replication.build_arm2 import load_faers_long  # noqa: E402


RESULTS = ROOT / "results"
DATA_REF = ROOT / "data" / "reference"


def _fmt(v, n=3):
    if v is None: return "—"
    try:
        if v != v: return "—"
        if v == float("inf"): return "∞"
        return f"{v:.{n}f}"
    except (TypeError, ValueError):
        return str(v)


def main() -> int:
    t0 = time.time()
    arm2_full = pd.read_parquet(RESULTS / "perpair_arm2.parquet")
    print(f"[Phase 5] loaded {len(arm2_full):,} Arm-2 pairs total", flush=True)

    # Filter to daen_powered for Family 3 test universe; mechanism features
    # only matter for pairs that enter the H5 contrast.
    arm2 = arm2_full[arm2_full["daen_powered"]].copy()
    print(f"[Phase 5] Family 3 universe (daen_powered): {len(arm2):,} pairs", flush=True)
    print(f"  cross_db_class distribution:")
    print(arm2["cross_db_class"].value_counts().to_string())

    # Need FAERS long tables for cohort feature lookups; reload (the build script discarded them)
    print("[Phase 5] reloading FAERS long tables for feature extraction ...", flush=True)
    drug_long, reac_long, _ = load_faers_long()

    alert_registry = pd.read_csv(DATA_REF / "alert_registry.csv")
    mass_tort = pd.read_csv(DATA_REF / "mass_tort_drugs.csv")

    print("[Phase 5] computing per-pair H5 mechanism features (on daen_powered subset) ...", flush=True)
    feats = build_arm2_features(arm2, drug_long, reac_long, alert_registry, mass_tort)
    feat_path = RESULTS / "perpair_arm2_with_features.parquet"
    feats.to_parquet(feat_path, index=False)
    print(f"[Phase 5] wrote {feat_path.name}  ({len(feats):,} rows)", flush=True)

    print("[Phase 5] running Family 3 H5a-d tests + Holm-Bonferroni ...", flush=True)
    results, holm = family3_holm(feats)

    out = {
        "phase": 5,
        "produced_at": time.strftime("%Y-%m-%d %H:%M:%S"),
        "arm2_n_pairs": int(len(feats)),
        "arm2_daen_powered_n": int(feats["daen_powered"].sum()),
        "family_3": {
            "tests": [asdict(r) for r in results],
            "holm": holm,
            "family_wise_alpha": 0.05,
            "k": 4,
        },
    }
    json_path = RESULTS / "phase5_results.json"
    with open(json_path, "w") as f:
        json.dump(out, f, indent=2, default=str)
    print(f"[Phase 5] wrote {json_path.name}", flush=True)

    # Markdown
    md = []
    md.append("# Phase 5 Results — Family 3 mechanism arm (Arm-2 discovery universe)\n\n")
    md.append(f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
    md.append(f"Input: `results/perpair_arm2_with_features.parquet` ({len(feats):,} pairs)\n")
    md.append(f"Arm-2 daen_powered universe: {int(feats['daen_powered'].sum()):,} pairs\n")
    md.append(f"  — cross_db_class distribution within daen_powered:\n\n")
    md.append("```\n")
    md.append(feats[feats["daen_powered"]]["cross_db_class"].value_counts().to_string())
    md.append("\n```\n\n")

    md.append("## Family 3 — Holm-Bonferroni at family-wise α = 0.05 (k = 4)\n\n")
    md.append("| Test | Description | n (faers_only) | n (conc_pos) | p (faers_only) | p (conc_pos) | Δ (pp) | 95% CI (pp) | p-value | Holm decision |\n")
    md.append("|---|---|---:|---:|---:|---:|---:|---|---:|---|\n")
    for r in results:
        md.append(f"| **{r.name}** | {r.description} | {r.n_faers_only} | {r.n_concordant_positive} | "
                  f"{_fmt(r.p_faers_only)} | {_fmt(r.p_concordant_positive)} | "
                  f"{_fmt(r.diff_pp, 2)} | {_fmt(r.ci_low_pp, 2)} – {_fmt(r.ci_high_pp, 2)} | "
                  f"{_fmt(r.pvalue, 4)} | `{holm[r.name]}` |\n")
    md.append("\n")

    md.append("## Pivot status (Section 5 Family 3 pivot)\n\n")
    all_fail = all(d == "FAIL_TO_REJECT" for d in holm.values())
    if all_fail:
        md.append("**Family 3 cluster NULL** — none of H5a-d rejects at Holm-Bonferroni. ")
        md.append("Per Section 5 Family 3 pivot framing: *Cross-database discordance without a canonical artefact signature: the unexplained heterogeneity of FAERS-only signals.* ")
        md.append("Combined with the H5c one-sidedness limitation (US-only alert registry), this is informative — ")
        md.append("the cross-database discordance pattern on the broader Arm-2 universe is not explained by the ")
        md.append("canonical US-side reporting-artefact catalogue, motivating mechanism-specific follow-up.\n\n")
    else:
        rejected = [name for name, d in holm.items() if d == "REJECT_AT_HOLM"]
        md.append(f"**Family 3 rejects at Holm:** {', '.join(rejected)}. ")
        md.append("Per Section 4 H5 + Section 4 'reservation rule', the headline use of the word 'artefact' ")
        md.append("ALSO requires Family 4 positivity on the Arm-1 intersection. Since Family 4 was insufficient-n, ")
        md.append("the reservation rule still binds the headline.\n\n")

    md_path = RESULTS / "phase5_results.md"
    with open(md_path, "w") as f:
        f.write("".join(md))
    print(f"[Phase 5] wrote {md_path.name}", flush=True)
    print(f"[Phase 5] complete in {time.time()-t0:.1f}s", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
