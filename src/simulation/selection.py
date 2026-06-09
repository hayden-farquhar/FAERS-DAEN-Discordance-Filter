"""Power-conditioning selection rules (protocol Section 5).

Each rule maps a per-pair frame to a boolean ``daen_powered`` flag. The two
families differ only in *what odds ratio the DAEN power model is asked about*:

- **MDE-anchored** rules use the smaller database's own marginals to find the
  minimum detectable effect (``power.daen_mde``) and admit the pair iff that MDE
  is at or below a FAERS-independent anchor. This is the registered fix.
- **FAERS-derived** rules ask whether DAEN is powered for the *FAERS* effect
  estimate (point or lower bound) via ``power.mc_power``. This is the
  winner's-curse-contaminated default the simulation is designed to indict.
- **none** admits every FAERS-positive pair (the unconditioned baseline).

Only FAERS-positive pairs can enter the H1 stratum, so the (expensive) power
calls are evaluated only on those rows; all other rows are ``False``. MDE values
are memoised per (DAEN marginals) key — they are a deterministic function of the
marginals, so caching does not perturb reproducibility.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable

import numpy as np
import pandas as pd

from power import mc_power
from power.fast import daen_mde_fast

# Default Monte-Carlo budget for the FAERS-derived power check inside the sweep.
# Smaller than power.mc_power's standalone default (5000) because it is called
# once per FAERS-positive pair per replicate; the G-SIM null-of-the-null check
# confirms this budget still controls type-I where it should.
_DEFAULT_MC_B = 400
_POWER_TARGET = 0.80


@dataclass
class SelectionContext:
    """Per-invocation state passed to every rule (carries rng + caches)."""
    rng: np.random.Generator
    mc_B: int = _DEFAULT_MC_B
    power_target: float = _POWER_TARGET
    _mde_cache: dict = field(default_factory=dict)

    def mde(self, n_drug: int, n_event: int, N: int) -> float:
        key = (int(n_drug), int(n_event), int(N))
        v = self._mde_cache.get(key)
        if v is None:
            # Vectorised exact Fisher-NCHG MDE; pinned bit-for-grid-point to the
            # scipy reference (power.mde.daen_mde) by tests/test_power_fast.py.
            v = daen_mde_fast(n_drug, n_event, N, power_target=self.power_target)
            self._mde_cache[key] = v
        return v


def _faers_positive_mask(perpair: pd.DataFrame) -> np.ndarray:
    return perpair["faers_signal"].to_numpy(dtype=bool)


def _mde_rule(anchor: float) -> Callable[[pd.DataFrame, SelectionContext], np.ndarray]:
    def rule(perpair: pd.DataFrame, ctx: SelectionContext) -> np.ndarray:
        out = np.zeros(len(perpair), dtype=bool)
        fp = _faers_positive_mask(perpair)
        nd = perpair["n_drug_D"].to_numpy()
        ne = perpair["n_event_D"].to_numpy()
        Nd = perpair["N_D"].to_numpy()
        for i in np.nonzero(fp)[0]:
            out[i] = ctx.mde(nd[i], ne[i], Nd[i]) <= anchor
        return out
    return rule


def _faers_derived_rule(which: str) -> Callable[[pd.DataFrame, SelectionContext], np.ndarray]:
    col = "ror_F" if which == "point" else "ror_F_ci_lo"

    def rule(perpair: pd.DataFrame, ctx: SelectionContext) -> np.ndarray:
        out = np.zeros(len(perpair), dtype=bool)
        fp = _faers_positive_mask(perpair)
        nd = perpair["n_drug_D"].to_numpy()
        ne = perpair["n_event_D"].to_numpy()
        Nd = perpair["N_D"].to_numpy()
        orf = perpair[col].to_numpy()
        for i in np.nonzero(fp)[0]:
            or_alt = float(orf[i])
            if not np.isfinite(or_alt) or or_alt <= 1.0:
                # No detectable positive effect claimed by FAERS -> not powered.
                out[i] = False
                continue
            p = mc_power(int(nd[i]), int(ne[i]), int(Nd[i]), or_alt,
                         B=ctx.mc_B, rng=ctx.rng)
            out[i] = p >= ctx.power_target
        return out
    return rule


def _none_rule(perpair: pd.DataFrame, ctx: SelectionContext) -> np.ndarray:
    # Baseline: every FAERS-positive pair is admitted (no power conditioning).
    return _faers_positive_mask(perpair)


SELECTION_RULES: dict[str, Callable[[pd.DataFrame, SelectionContext], np.ndarray]] = {
    "mde_1.5": _mde_rule(1.5),
    "mde_2.0": _mde_rule(2.0),
    "mde_3.0": _mde_rule(3.0),
    "faers_point": _faers_derived_rule("point"),
    "faers_lb": _faers_derived_rule("lb"),
    "none": _none_rule,
}


def apply_selection(
    perpair: pd.DataFrame, rule_name: str, ctx: SelectionContext,
) -> np.ndarray:
    """Return the ``daen_powered`` boolean array for ``rule_name``."""
    try:
        rule = SELECTION_RULES[rule_name]
    except KeyError:
        raise ValueError(
            f"unknown selection rule {rule_name!r}; "
            f"expected one of {sorted(SELECTION_RULES)}"
        ) from None
    return rule(perpair, ctx)
