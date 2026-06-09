"""Pin the vectorised fast power/MDE path to the scipy-exact reference.

``power.fast`` is the production path used by the simulation sweep; ``power.exact``
+ ``power.mde`` are the authoritative (slow) reference. These tests guarantee the
fast path returns numerically identical power (to fp tolerance) and *exactly* the
same MDE grid point as the reference, so swapping it into the sweep cannot move
the registered ``daen_powered`` classification.
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest

from power.exact import exact_power
from power.mde import daen_mde
from power.fast import exact_power_fast, daen_mde_fast

ROOT = Path(__file__).resolve().parents[1]
ANCHOR_FRAME = ROOT / "results" / "perpair_arm1.parquet"


def _random_margins(rng):
    N = int(rng.integers(200, 700_000))
    n_drug = int(rng.integers(1, max(2, N // 50)))
    n_event = int(rng.integers(1, max(2, N // 5)))
    return n_drug, n_event, N


def test_fast_power_matches_reference():
    rng = np.random.default_rng(20260608)
    max_abs = 0.0
    checked = 0
    for _ in range(400):
        n_drug, n_event, N = _random_margins(rng)
        if n_drug + n_event >= N:
            continue
        OR = float(rng.uniform(1.0, 5.0))
        ref = exact_power(n_drug, n_event, N, OR)
        fast = exact_power_fast(n_drug, n_event, N, OR)
        max_abs = max(max_abs, abs(ref - fast))
        checked += 1
    assert checked > 100
    # 1e-7 is pure logsumexp-vs-scipy normalisation roundoff over large supports;
    # far tighter than anything that can move an MDE decision (target 0.80 on a
    # 0.05-spaced grid). The grid-point equivalence test below is the real lock.
    assert max_abs < 1e-7, f"fast power diverged from reference by {max_abs:.2e}"


def test_fast_mde_matches_reference_grid_point():
    rng = np.random.default_rng(94)
    mism = 0
    checked = 0
    for _ in range(400):
        n_drug, n_event, N = _random_margins(rng)
        if n_drug + n_event >= N:
            continue
        ref = daen_mde(n_drug, n_event, N)
        fast = daen_mde_fast(n_drug, n_event, N)
        same = (ref == fast) or (math.isinf(ref) and math.isinf(fast))
        mism += int(not same)
        checked += 1
    assert checked > 100
    assert mism == 0, f"{mism}/{checked} MDE grid-point mismatches fast-vs-reference"


@pytest.mark.slow
def test_fast_mde_reproduces_empirical_anchor():
    if not ANCHOR_FRAME.exists():
        pytest.skip("empirical anchor frame not present")
    import pandas as pd

    df = pd.read_parquet(ANCHOR_FRAME)
    nd = df["n_drug_D"].to_numpy()
    ne = df["n_event_D"].to_numpy()
    Nd = df["N_D"].to_numpy()
    mde = [daen_mde_fast(int(nd[i]), int(ne[i]), int(Nd[i])) for i in range(len(df))]
    powered = sum(1 for v in mde if v <= 1.5)
    assert powered == 52, f"fast MDE powered count {powered} != registered anchor 52"
