"""H2 sensitivity arms (Section 10.7) — quick re-run of the power-conditioned
non-replication rate under each pre-registered sensitivity arm.

For each daen_powered_* tag in the augmented substrate, recompute:
  - n in H2 universe (FAERS-positive AND <tag>)
  - H2 = proportion with cross_db_class == 'faers_only'
  - cluster-bootstrap 95% CI (B=2,000; seed=94; clusters=events)

Produces results/phase6_h2_sensitivity.md.
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

RESULTS = ROOT / "results"

SENSITIVITY_TAGS = [
    ("daen_powered",                "PRIMARY (MDE ≤ 1.5)",                                ),
    ("daen_powered_07",             "MDE ≤ 1.5 at 70% power (permissive power threshold)" ),
    ("daen_powered_09",             "MDE ≤ 1.5 at 90% power (stringent power threshold)"  ),
    ("daen_powered_at_2",           "MDE ≤ 2.0 (permissive anchor)"                       ),
    ("daen_powered_at_3",           "MDE ≤ 3.0 (very permissive anchor)"                  ),
    ("daen_powered_FAERS_pt",       "Power ≥ 0.80 vs FAERS point ROR (the contested fix)" ),
    ("daen_powered_FAERS_LB",       "Power ≥ 0.80 vs FAERS 95% LB"                        ),
    ("daen_powered_unconditional",  "Poisson-margin MDE ≤ 1.5"                            ),
    ("daen_powered_min",            "primary AND Poisson (conservative composite)"        ),
]


def cluster_boot_rate(df: pd.DataFrame, indicator_col: str, cluster_col: str,
                       B: int = 2000, seed: int = 94) -> tuple[float, float, float]:
    rng = np.random.default_rng(seed)
    clusters = df[cluster_col].unique()
    if len(clusters) == 0 or len(df) == 0:
        return (float("nan"), float("nan"), float("nan"))
    grouped = {c: df.loc[df[cluster_col] == c, indicator_col].values for c in clusters}
    point = float(df[indicator_col].mean())
    boot = np.empty(B, dtype=float)
    for b in range(B):
        sampled = rng.choice(clusters, size=len(clusters), replace=True)
        vals = np.concatenate([grouped[c] for c in sampled])
        boot[b] = float(vals.mean()) if vals.size else float("nan")
    boot = boot[~np.isnan(boot)]
    if boot.size == 0:
        return (point, float("nan"), float("nan"))
    return (point, float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5)))


def main() -> int:
    t0 = time.time()
    df = pd.read_parquet(RESULTS / "perpair_arm1.parquet")
    print(f"Loaded {len(df)} pairs", flush=True)

    rows = []
    for tag, label in SENSITIVITY_TAGS:
        sub = df[(df["faers_signal"]) & (df[tag])].copy()
        if len(sub) == 0:
            rows.append({"tag": tag, "label": label, "n": 0,
                          "h2_point": float("nan"), "h2_lo": float("nan"), "h2_hi": float("nan"),
                          "n_underpow": int((df["faers_signal"] & ~df[tag]).sum())})
            continue
        sub["is_faers_only"] = (sub["cross_db_class"] == "faers_only").astype(int)
        point, lo, hi = cluster_boot_rate(sub, "is_faers_only", "outcomeName")
        rows.append({
            "tag": tag, "label": label,
            "n": len(sub),
            "h2_point": point, "h2_lo": lo, "h2_hi": hi,
            "n_underpow": int((df["faers_signal"] & ~df[tag]).sum()),
        })

    table = pd.DataFrame(rows)
    print()
    print(table.to_string(index=False))

    md = ["# Phase 6 (partial) — H2 power-conditioned non-replication rate, sensitivity arms\n\n",
          f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}\n",
          "Source: `results/perpair_arm1.parquet` (492 pairs)\n",
          "Endpoint: H2 = proportion of FAERS-positive AND <power-tag> Arm-1 pairs classified `faers_only`.\n",
          "Cluster-bootstrap CIs (B=2000, seed=94, clusters=events).\n\n",
          "## Sensitivity table\n\n",
          "| Power tag | Description | n (H2 universe) | H2 point | 95% CI | underpowered-discordant alongside |\n",
          "|---|---|---:|---:|---|---:|\n"]
    for r in rows:
        ci = f"{r['h2_lo']:.3f} – {r['h2_hi']:.3f}" if not (r["h2_lo"] != r["h2_lo"]) else "—"
        h2 = f"{r['h2_point']:.3f}" if not (r["h2_point"] != r["h2_point"]) else "—"
        md.append(f"| `{r['tag']}` | {r['label']} | {r['n']} | {h2} | {ci} | {r['n_underpow']} |\n")
    md.append("\n## Interpretation\n\n")
    md.append("**Primary (PRIMARY MDE ≤ 1.5):** the headline number.\n\n")
    md.append("**Power-threshold sensitivities (0.7 / 0.9):** test whether the headline is sensitive to the power cutoff.\n\n")
    md.append("**MDE-anchor sensitivities (at_2 / at_3):** test whether the headline moves as the implicit ")
    md.append("'real-effect' threshold rises. **A meaningful gap here is the H2 construct caveat (Section 4) made concrete** ")
    md.append("— pairs with true effects between the primary 1.5 and the sensitivity anchor inflate apparent non-replication.\n\n")
    md.append("**FAERS-derived sensitivities (FAERS_pt / FAERS_LB):** the contested earlier-design approach. ")
    md.append("**The gap from the primary quantifies the FAERS-derived-power circularity** that the Round-2 must-fix addressed. ")
    md.append("If the FAERS-derived H2 is markedly higher than the primary, the circularity would have inflated the headline.\n\n")
    md.append("**Poisson-margin sensitivities (unconditional / min):** the fixed-margin caveat (Section 6.3).\n\n")
    with open(RESULTS / "phase6_h2_sensitivity.md", "w") as f:
        f.write("".join(md))
    print(f"\nWrote results/phase6_h2_sensitivity.md")
    print(f"Wall time: {time.time()-t0:.1f}s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
