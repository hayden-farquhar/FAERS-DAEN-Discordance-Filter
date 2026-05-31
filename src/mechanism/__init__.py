"""Family 3 mechanism arm (protocol Section 4 Family 3 + Section 10.3).

H5a — Consumer report source share enrichment
H5b — Lawyer report source share enrichment
H5c — Post-alert temporal-breakpoint alignment
H5d — Mass-tort drug membership

Each H5a-d is tested via difference-of-proportions (faers_only vs concordant_positive)
within the Arm-2 daen_powered universe; Holm-Bonferroni at family-wise α=0.05, k=4.
"""
from .h5_features import build_arm2_features
from .family3_tests import family3_holm

__all__ = ["build_arm2_features", "family3_holm"]
