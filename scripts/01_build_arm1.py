"""Phase 1 Arm-1 substrate build runner.

Builds per-pair 2x2 contingency tables for the OMOP+EU-ADR Arm-1 universe
in both FAERS and DAEN, applies the pre-registered signal rule (LB > 1 AND
a >= 3), and persists to results/perpair_arm1.parquet.

Run from project root:
    python scripts/run_arm1_build.py
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

# Allow `from replication ...` when run as a script
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from replication.build_arm1 import build  # noqa: E402


def main() -> int:
    t0 = time.time()
    out = build()
    dt = time.time() - t0
    print(f"\nWall time: {dt:.1f}s")
    print(f"Output rows: {len(out)}")
    print()
    print("=== Cross-DB class distribution (Arm-1, signal rule applied) ===")
    print(out["cross_db_class"].value_counts().to_string())
    print()
    print("=== Class x groundTruth ===")
    print(out.groupby(["cross_db_class", "groundTruth"]).size().unstack(fill_value=0))
    print()
    print("=== Sample rows ===")
    cols = ["ref_set", "exposureName", "outcomeName", "groundTruth",
            "a_F", "ror_F", "faers_signal", "a_D", "ror_D", "daen_signal", "cross_db_class"]
    print(out[cols].head(10).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
