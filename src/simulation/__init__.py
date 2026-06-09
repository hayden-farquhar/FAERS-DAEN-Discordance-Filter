"""Simulation module for the cross-database enrichment type-I/power study.

Implements the locked simulation pre-analysis plan (OSF Amendment 3,
``osf/simulation_protocol.md``, locked 2026-06-08). Four components:

- :mod:`simulation.dgp`        — data-generating process (protocol Section 3)
- :mod:`simulation.selection`  — power-conditioning selection rules (Section 5)
- :mod:`simulation.estimators` — enrichment estimator panel (Section 6)
- :mod:`simulation.sweep`      — sweep driver -> results/sim_results.parquet (Sections 4/7/10)

The DGP applies the registered signal rule (``power.signal.signal_fires``) and
power model (``power.daen_mde`` / ``power.mc_power``) UNCHANGED, so the simulated
estimand is identical to the empirical one. The enrichment statistic mirrors the
fixed-effect / cluster-robust spec of ``enrichment.family2.family2_h1``.
"""
from __future__ import annotations

from .dgp import DGPCell, generate_population
from .selection import SELECTION_RULES, apply_selection
from .estimators import ESTIMATORS, EstimatorResult, run_estimator

__all__ = [
    "DGPCell",
    "generate_population",
    "SELECTION_RULES",
    "apply_selection",
    "ESTIMATORS",
    "EstimatorResult",
    "run_estimator",
]
