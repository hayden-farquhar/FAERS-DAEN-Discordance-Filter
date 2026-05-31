"""Phase 5.2 runner: build Arm-2 discovery substrate."""
from __future__ import annotations

import sys, time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from replication.build_arm2 import build  # noqa: E402

t0 = time.time()
df = build()
dt = time.time() - t0
print(f"\nWall time: {dt:.1f}s ({dt/60:.1f} min)")
print(f"Output rows (FAERS-positive Arm-2 pairs): {len(df):,}")
print()
print("=== cross_db_class distribution ===")
print(df["cross_db_class"].value_counts().to_string())
print()
print("=== daen_powered distribution ===")
print(df["daen_powered"].value_counts().to_string())
print()
print(f"=== H5 Family 3 universe (FAERS+ AND daen_powered): {((df['faers_signal']) & (df['daen_powered'])).sum():,} pairs ===")
print(df[df["daen_powered"]]["cross_db_class"].value_counts().to_string())
