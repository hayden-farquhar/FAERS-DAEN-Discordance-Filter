"""Shared disproportionality-metrics module (single source of truth).

Implements the four metrics pre-registered for the FAERS signal re-analysis audit,
each returning a point estimate plus an interval and a pre-registered signalling flag:

- ROR  (reporting odds ratio)            -- primary metric; signal: lower 95% CI > 1, a >= 3
- PRR  (proportional reporting ratio)     -- Evans 2001: PRR >= 2, chi2 >= 4, a >= 3
- IC   (BCPNN information component)       -- signal: IC025 > 0
- EBGM (MGPS Gamma-Poisson shrinker)      -- signal: EB05 > 2

ROR/PRR/IC are per-cell (one 2x2 + total N). EBGM is database-wide: fit the prior on
all cells with `fit_hyperparameters`, then evaluate cells with `ebgm_cell`.
"""

from .contingency import Contingency
from .ror import ror, RORResult
from .prr import prr, PRRResult
from .ic_bcpnn import ic, ic_from_table, ICResult
from .ebgm_mgps import (
    fit_hyperparameters,
    ebgm_cell,
    mgps,
    MGPSHyperparams,
    EBGMResult,
)

__all__ = [
    "Contingency",
    "ror", "RORResult",
    "prr", "PRRResult",
    "ic", "ic_from_table", "ICResult",
    "fit_hyperparameters", "ebgm_cell", "mgps", "MGPSHyperparams", "EBGMResult",
]
