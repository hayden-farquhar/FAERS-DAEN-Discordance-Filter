"""Sweep driver (protocol Sections 4, 7, 10).

Builds the DGP grid from ``sim_config.yaml`` (primary slice + one-at-a-time
marginal sweeps), runs ``replicates`` populations per cell, applies every
selection rule and every estimator, and writes one tidy row per
(replicate x selection rule x estimator x DGP cell) to
``results/sim_results.parquet``.

Robustness (per the portfolio notebook/long-running-code rule):
- ``numpy.random.SeedSequence(master).spawn`` gives an independent, recorded
  child stream per cell; the spawn map is written alongside the results.
- Each cell's rows are checkpointed to ``results/sim_shards/<cell>.parquet`` as
  it finishes, so an interrupted sweep resumes by skipping existing shards.
- ``tqdm`` reports cell-level and replicate-level progress.

The full grid is large; ``--quick`` runs a tiny smoke configuration for wiring
checks, and ``--replicates`` / ``--cells`` override the config for partial runs.
G-SIM (``tests/test_simulation_gsim.py``) must pass before any reported number is
drawn from a full run.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import time
from dataclasses import asdict, replace
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from tqdm import tqdm

from .dgp import DGPCell, MarginalParams, generate_population
from .selection import SELECTION_RULES, SelectionContext, apply_selection
from .estimators import ESTIMATORS, build_h1_stratum, run_estimator

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG = ROOT / "sim_config.yaml"
RESULTS = ROOT / "results"
SHARD_DIR = RESULTS / "sim_shards"


def _log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def load_config(path: Path | str = DEFAULT_CONFIG) -> dict:
    with open(path) as f:
        return yaml.safe_load(f)


def _marginals(cfg: dict) -> MarginalParams:
    m = cfg["marginals"]
    return MarginalParams(
        log_mu_drug=float(m["log_mu_drug"]), log_sd_drug=float(m["log_sd_drug"]),
        log_mu_event=float(m["log_mu_event"]), log_sd_event=float(m["log_sd_event"]),
        calibrated=bool(m.get("calibrated", False)),
    )


def _cell_id(cell: DGPCell) -> str:
    """Stable short id for a cell (used as the shard filename + seed-spawn key)."""
    d = {k: v for k, v in asdict(cell).items() if k != "marginals"}
    blob = json.dumps(d, sort_keys=True, default=str)
    return hashlib.sha1(blob.encode()).hexdigest()[:16]


def build_grid(cfg: dict) -> list[DGPCell]:
    """Primary slice (sweep lambda x phi x pi) + marginal sweeps (one factor off
    the slice at a time), each crossed with both regimes. De-duplicated by id."""
    f = cfg["factors"]
    slice_ = cfg["primary_slice"]
    mp = _marginals(cfg)
    n_pairs = int(cfg["population"]["n_pairs"])
    N_F = int(cfg["population"]["N_F"])
    q2 = float(cfg["q2_daen_suppression"])
    regimes = cfg["regimes"]

    def make(regime, pi, lam, phi, rho, K, icc, OR_true, tag) -> DGPCell:
        return DGPCell(
            regime=regime, pi=pi, lambda_inflate=lam, phi=phi, rho=rho, K=K, icc=icc,
            OR_true=OR_true, N_F=N_F, n_pairs=n_pairs, marginals=mp,
            q2_daen_suppression=q2, label=tag,
        )

    cells: dict[str, DGPCell] = {}

    def add(cell: DGPCell) -> None:
        cells.setdefault(_cell_id(cell), cell)

    # (a) primary slice: lambda x phi x pi
    for regime in regimes:
        for pi in f["pi"]:
            for lam in f["lambda_inflate"]:
                for phi in f["phi"]:
                    add(make(regime, pi, lam, phi, slice_["rho"], slice_["K"],
                             slice_["icc"], slice_["OR_true"], "primary"))

    # (b) marginal sweeps: vary one factor off the slice at its central pi/lambda/phi
    base_pi, base_lam, base_phi = 0.5, 2, 0.5
    for factor in cfg["marginal_sweeps"]:
        for level in f[factor]:
            kw = dict(rho=slice_["rho"], K=slice_["K"], icc=slice_["icc"],
                      OR_true=slice_["OR_true"])
            kw[factor] = level
            for regime in regimes:
                add(make(regime, base_pi, base_lam, base_phi,
                         kw["rho"], kw["K"], kw["icc"], kw["OR_true"], f"sweep_{factor}"))

    return list(cells.values())


def run_cell(cell: DGPCell, *, replicates: int, cfg: dict,
             seed_seq: np.random.SeedSequence) -> pd.DataFrame:
    """Run all replicates for one cell; return the tidy row frame."""
    # Three independent streams off this cell's seed: population, selection, bootstrap.
    pop_ss, sel_ss, boot_ss = seed_seq.spawn(3)
    pop_rng = np.random.default_rng(pop_ss)
    sel_rng = np.random.default_rng(sel_ss)
    boot_rng = np.random.default_rng(boot_ss)
    mc_B = int(cfg["mc_power_B"])
    boot_B = int(cfg["bootstrap_B"])

    rows: list[dict] = []
    base = {k: v for k, v in asdict(cell).items() if k != "marginals"}
    for rep in range(replicates):
        pop = generate_population(cell, pop_rng)
        sel_ctx = SelectionContext(rng=sel_rng, mc_B=mc_B)
        for rule in SELECTION_RULES:
            pop_r = pop.assign(daen_powered=apply_selection(pop, rule, sel_ctx))
            stratum = build_h1_stratum(pop_r)
            for est in ESTIMATORS:
                res = run_estimator(stratum, est, rng=boot_rng, B=boot_B)
                rows.append({
                    **base, "cell_id": _cell_id(cell), "replicate": rep,
                    "selection_rule": rule, "estimator": est,
                    "reject": res.reject, "non_estimable": res.non_estimable,
                    "or_hat": res.or_hat, "ci_lo": res.ci_lo, "ci_hi": res.ci_hi,
                    "beta": res.beta, "se": res.se,
                    "stratum_size": res.stratum_size,
                    "n_faers_only": res.n_faers_only,
                    "n_concordant_positive": res.n_concordant_positive,
                    "n_clusters": res.n_clusters,
                })
    return pd.DataFrame(rows)


def _write_spawn_map(grid: list[DGPCell], child_seqs, seed_master: int) -> None:
    """Record the full per-cell seed spawn map (partition-invariant)."""
    spawn_map = {}
    for cell, ss in zip(grid, child_seqs):
        cid = _cell_id(cell)
        spawn_map[cid] = list(ss.entropy) if hasattr(ss.entropy, "__iter__") else ss.entropy
    (RESULTS / "sim_seed_spawn_map.json").write_text(
        json.dumps({"seed_master": seed_master, "cells": spawn_map}, indent=2,
                   default=str))


def _combine_shards(grid: list[DGPCell], out: Path) -> Path:
    """Concatenate all per-cell shards into the tidy results parquet (only when
    every cell's shard is present)."""
    present = [c for c in grid if (SHARD_DIR / f"{_cell_id(c)}.parquet").exists()]
    if len(present) < len(grid):
        _log(f"Combine skipped: {len(present)}/{len(grid)} shards present.")
        return out
    combined = pd.concat(
        [pd.read_parquet(SHARD_DIR / f"{_cell_id(c)}.parquet") for c in grid],
        ignore_index=True)
    # OR_true mixes numeric levels (1.5/2/3/5) with the categorical 'lognormal'
    # mixture level, so the concatenated column is object-typed and unwritable as
    # a single Arrow type. Normalise to a string label; downstream derives a
    # numeric value where the label parses.
    combined = combined.astype({"OR_true": str})
    combined.to_parquet(out, index=False)
    _log(f"Wrote {len(combined):,} rows -> {out}")
    return out


def run_sweep(*, config: Path | str = DEFAULT_CONFIG, replicates: int | None = None,
              max_cells: int | None = None, out: Path | None = None,
              quick: bool = False, shard_index: int | None = None,
              num_shards: int | None = None, combine_only: bool = False) -> Path:
    cfg = load_config(config)
    if quick:
        cfg = _quickify(cfg)
    reps = replicates if replicates is not None else int(cfg["replicates"])
    grid = build_grid(cfg)
    if max_cells is not None:
        grid = grid[:max_cells]

    SHARD_DIR.mkdir(parents=True, exist_ok=True)
    master = np.random.SeedSequence(int(cfg["seed_master"]))
    # One child SeedSequence per cell, keyed by *position in the full grid*. This
    # is what makes partitioned parallel runs reproducible: every worker spawns
    # the same full set and picks its assigned indices, so a cell's seed never
    # depends on how the work was split.
    child_seqs = master.spawn(len(grid))

    out = Path(out) if out else (RESULTS / "sim_results.parquet")
    is_worker = shard_index is not None and num_shards is not None

    if not combine_only:
        assigned = [i for i in range(len(grid))
                    if (not is_worker) or (i % num_shards == shard_index)]
        tag = f"[shard {shard_index}/{num_shards}] " if is_worker else ""
        _log(f"{tag}Sweep: {len(assigned)}/{len(grid)} cells x {reps} replicates "
             f"x {len(SELECTION_RULES)} rules x {len(ESTIMATORS)} estimators")
        for idx in tqdm(assigned, desc=f"{tag}cells"):
            cell = grid[idx]
            cid = _cell_id(cell)
            shard = SHARD_DIR / f"{cid}.parquet"
            if shard.exists():
                continue   # resumable: skip already-computed cells
            df = run_cell(cell, replicates=reps, cfg=cfg, seed_seq=child_seqs[idx])
            df.to_parquet(shard, index=False)

    # A worker never finalises (other workers may still be running). The full run
    # (no partition) and the explicit --combine-only pass write the spawn map and
    # concatenate the shards.
    if combine_only or not is_worker:
        _write_spawn_map(grid, child_seqs, int(cfg["seed_master"]))
        return _combine_shards(grid, out)
    return SHARD_DIR


def _quickify(cfg: dict) -> dict:
    """Shrink the config for a wiring smoke-test (not for any reported number)."""
    cfg = json.loads(json.dumps(cfg))
    cfg["replicates"] = 5
    cfg["mc_power_B"] = 100
    cfg["bootstrap_B"] = 99
    cfg["population"]["n_pairs"] = 400
    cfg["factors"]["lambda_inflate"] = [1, 5]
    cfg["factors"]["phi"] = [0.0, 1.0]
    cfg["factors"]["pi"] = [0.5]
    cfg["marginal_sweeps"] = []
    return cfg


def _parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Run the enrichment type-I/power sweep.")
    p.add_argument("--config", default=str(DEFAULT_CONFIG))
    p.add_argument("--replicates", type=int, default=None)
    p.add_argument("--cells", type=int, default=None, help="cap number of DGP cells")
    p.add_argument("--out", default=None)
    p.add_argument("--quick", action="store_true", help="tiny smoke configuration")
    p.add_argument("--shard-index", type=int, default=None,
                   help="this worker's index in [0, num-shards) for parallel runs")
    p.add_argument("--num-shards", type=int, default=None,
                   help="total number of parallel workers (cells split by index %% n)")
    p.add_argument("--combine-only", action="store_true",
                   help="skip computation; just write spawn map + combine shards")
    return p.parse_args(argv)


if __name__ == "__main__":
    args = _parse_args()
    run_sweep(config=args.config, replicates=args.replicates,
              max_cells=args.cells, out=args.out, quick=args.quick,
              shard_index=args.shard_index, num_shards=args.num_shards,
              combine_only=args.combine_only)
