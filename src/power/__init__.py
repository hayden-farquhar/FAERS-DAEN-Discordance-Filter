"""Per-pair power model (protocol Section 6).

The primary `daen_powered` tag is MDE-based: daen_powered = (daen_mde <= 1.5),
where daen_mde is the smallest OR at which DAEN has >= 80% power to fire the
signal rule (LB of 95% CI of ROR > 1 AND a >= 3) given the pair's marginals.
"""
from __future__ import annotations

from .signal import signal_fires
from .exact import exact_power
from .wald import wald_power
from .montecarlo import mc_power
from .mde import daen_mde, daen_mde_unconditional

__all__ = [
    "signal_fires",
    "exact_power",
    "wald_power",
    "mc_power",
    "daen_mde",
    "daen_mde_unconditional",
]
