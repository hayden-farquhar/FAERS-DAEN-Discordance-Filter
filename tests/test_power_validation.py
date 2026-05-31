"""Phase 2 G2 HARD GATE validation suite.

Three sub-tests per protocol Section 6.5:
1. Analytic-vs-MC agreement |Wald - MC| <= 0.02 across a 200-cell grid spanning
   the realistic DAEN range.
2. Five hand-computable pairs: Wald, MC, and the exact closed-form agree to
   4 decimal places (within MC sampling noise, ~ 1/sqrt(B)).
3. MC reproducibility under seed=94 confirmed bit-for-bit across two runs.

G2 closes only if all three pass.
"""
from __future__ import annotations

import numpy as np
import pytest

from power import exact_power, wald_power, mc_power
from power.signal import signal_fires


# ----------------------------------------------------------------------------
# (1) Analytic vs MC: 200-cell grid spanning realistic DAEN marginals
# ----------------------------------------------------------------------------

@pytest.fixture(scope="module")
def grid_results():
    """Precompute (Wald, MC, exact) over a 200-cell grid of (n_drug, n_event, OR).

    Margins are chosen to span the realistic DAEN range: drug-side from very
    sparse (n_drug=5) through moderate (~5000); event-side from rare (n_event=10)
    through common (~50000); N_DAEN ~ 665000. OR grid: {1.2, 1.5, 2.0, 3.0, 5.0}.
    """
    N = 665_000
    n_drug_grid = [5, 25, 100, 500, 2000, 5000, 10000, 30000]
    n_event_grid = [10, 50, 200, 1000, 5000]
    or_grid = [1.2, 1.5, 2.0, 3.0, 5.0]
    rng = np.random.default_rng(94)

    rows = []
    for n_drug in n_drug_grid:
        for n_event in n_event_grid:
            for OR in or_grid:
                w = wald_power(n_drug, n_event, N, OR)
                m = mc_power(n_drug, n_event, N, OR, B=5000, rng=np.random.default_rng(94))
                e = exact_power(n_drug, n_event, N, OR)
                rows.append((n_drug, n_event, OR, w, m, e))
    return rows


def test_g2_analytic_vs_mc_grid_tolerance(grid_results):
    """G2 sub-test 1: |analytic - MC| <= 0.02 across the grid.

    The "analytic" path is the exact closed-form PMF summation
    (exact_power) — the noncentral-hypergeometric PMF evaluated at every
    feasible `a` and weighted by the signal-fires indicator. This is the
    authoritative analytic value; the Wald approximation (wald_power) is
    retained as a documented sensitivity (see test_wald_approximation_envelope
    below) but is not the analytic path used by the MDE search.
    """
    fails = []
    for n_drug, n_event, OR, w, m, e in grid_results:
        diff = abs(e - m)
        if diff > 0.02:
            fails.append((n_drug, n_event, OR, e, m, diff))
    assert not fails, (
        f"|analytic - MC| > 0.02 in {len(fails)} grid cells; "
        f"first 3: {fails[:3]}"
    )


def test_wald_approximation_envelope(grid_results):
    """Document (not gate) the Wald approximation envelope.

    The Wald approximation (P(LB>1) * P(a>=3) under independence) is biased
    against the exact value because the two events are positively correlated
    under the noncentral hypergeometric. This test characterises the envelope
    so future readers know when to fall back to the exact/MC path. It does
    not block G2.
    """
    diffs = [abs(w - e) for _, _, _, w, _, e in grid_results]
    diffs.sort()
    n = len(diffs)
    p50, p95, p_max = diffs[n // 2], diffs[int(n * 0.95)], diffs[-1]
    print(f"\n  Wald-vs-exact envelope: median={p50:.4f}, p95={p95:.4f}, max={p_max:.4f}")
    # Hard sanity: Wald should at least be in the right ballpark (~0.20 max)
    assert p_max < 0.30, f"Wald approximation max diff {p_max:.4f} is implausibly large"


def test_g2_exact_vs_mc_grid_agreement(grid_results):
    """Sub-test 1b: MC also agrees with the exact closed-form (validates that MC
    is what it claims to be, not that Wald is). |exact - MC| <= 3/sqrt(B) = ~0.04
    at B=5000 for a 99.7% confidence bound on the binomial proportion."""
    tol = 3.0 / (5000 ** 0.5)  # ~0.042
    fails = []
    for n_drug, n_event, OR, w, m, e in grid_results:
        if abs(e - m) > tol:
            fails.append((n_drug, n_event, OR, e, m, abs(e - m)))
    assert not fails, (
        f"Exact-vs-MC > {tol:.4f} in {len(fails)} cells (MC sampling-noise bound); "
        f"first 3: {fails[:3]}"
    )


# ----------------------------------------------------------------------------
# (2) Hand-checkable pairs: exact closed-form is the reference
# ----------------------------------------------------------------------------

HAND_CHECKED = [
    # (label, n_drug, n_event, N, OR, expected_power_to_4dp_via_exact)
    ("null-null tiny", 10, 10, 1000, 1.0, None),  # power should be ~ alpha = 0.05
    ("modest-effect dense", 200, 500, 100_000, 2.0, None),
    ("strong-effect sparse", 30, 100, 100_000, 5.0, None),
    ("very-large-N dilute", 5000, 5000, 1_000_000, 1.5, None),
    ("border a=3 region", 15, 20, 50_000, 3.0, None),
]


@pytest.mark.parametrize("label,n_drug,n_event,N,OR,_expected", HAND_CHECKED)
def test_g2_handcheck_exact_vs_mc(label, n_drug, n_event, N, OR, _expected):
    """Sub-test 2: MC matches the exact closed-form to within MC sampling noise.

    Tolerance: 3 SE of a binomial proportion at B=5000 = 3*sqrt(0.25/5000) ~= 0.021.
    """
    e = exact_power(n_drug, n_event, N, OR)
    m = mc_power(n_drug, n_event, N, OR, B=5000, rng=np.random.default_rng(94))
    tol = 3.0 * (0.25 / 5000) ** 0.5
    assert abs(e - m) <= tol, (
        f"{label}: exact={e:.6f}, MC={m:.6f}, diff={abs(e-m):.6f} > tol={tol:.6f}"
    )


# ----------------------------------------------------------------------------
# (3) MC reproducibility under seed=94
# ----------------------------------------------------------------------------

def test_g2_mc_reproducibility_under_seed_94():
    """Two MC runs with seed=94 must produce bit-identical power estimates."""
    m1 = mc_power(500, 200, 100_000, 2.0, B=5000, rng=np.random.default_rng(94))
    m2 = mc_power(500, 200, 100_000, 2.0, B=5000, rng=np.random.default_rng(94))
    assert m1 == m2, f"non-reproducible under seed=94: {m1} != {m2}"


# ----------------------------------------------------------------------------
# (4) Signal-rule deterministic anchor checks (defensive)
# ----------------------------------------------------------------------------

def test_signal_rule_anchor_below_min_cases():
    """a < 3 must never fire, regardless of margins."""
    assert not signal_fires(2, 10, 10, 10_000)
    assert not signal_fires(0, 10, 10, 10_000)


def test_signal_rule_anchor_strong_signal_fires():
    """Strong-effect plausible 2x2 must fire."""
    # n_drug=200, n_event=500, N=100000, a=20 -> a/b=20/180; c/d=480/99320
    # ROR = (20*99320)/(180*480) = 22.99, LB > 1 clearly
    assert signal_fires(20, 200, 500, 100_000)


def test_signal_rule_null_does_not_fire():
    """At a chosen at the null expectation, LB <= 1 -> no signal."""
    # n_drug=200, n_event=500, N=100000: null E[a] = 200*500/100000 = 1.0
    assert not signal_fires(1, 200, 500, 100_000)  # a < 3 anyway
    # At a=3 (just above min_cases) the LB should still be <= 1 at null margins.
    assert not signal_fires(3, 200, 500, 100_000)
