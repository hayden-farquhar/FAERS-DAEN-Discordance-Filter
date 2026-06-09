"""Compute specific DGP cells by grid index, byte-identically to the sweep.

Used to rebalance a lopsided tail: when one shard-worker is left grinding
several cells serially, stop it and launch one of these per remaining cell so
they run in parallel. Reproduces ``run_sweep``'s seed architecture exactly --
``SeedSequence(seed_master).spawn(len(grid))`` indexed by grid position -- so a
shard written here is identical to one written by the worker that owned it.
Idempotent: skips any cell whose shard already exists.

Usage:  python3 run_cells.py 107 113 119 125
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))

import numpy as np
from simulation.sweep import load_config, build_grid, run_cell, _cell_id, SHARD_DIR


def main(indices: list[int]) -> None:
    cfg = load_config()
    grid = build_grid(cfg)
    reps = int(cfg["replicates"])
    child_seqs = np.random.SeedSequence(int(cfg["seed_master"])).spawn(len(grid))
    SHARD_DIR.mkdir(parents=True, exist_ok=True)
    for idx in indices:
        cell = grid[idx]
        shard = SHARD_DIR / f"{_cell_id(cell)}.parquet"
        if shard.exists():
            print(f"[run_cells] skip idx={idx} (shard exists)", flush=True)
            continue
        print(f"[run_cells] computing idx={idx} cid={_cell_id(cell)} "
              f"({reps} reps)", flush=True)
        df = run_cell(cell, replicates=reps, cfg=cfg, seed_seq=child_seqs[idx])
        df.to_parquet(shard, index=False)
        print(f"[run_cells] done idx={idx} rows={len(df)}", flush=True)


if __name__ == "__main__":
    main([int(x) for x in sys.argv[1:]])
