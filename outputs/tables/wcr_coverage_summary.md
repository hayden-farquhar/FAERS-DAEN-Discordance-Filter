# Post-hoc S4: restricted wild cluster bootstrap (WCR) coverage

Q1_null K-sweep line on rule `mde_1.5`, 2000 replicates per K-cell, on the **same strata** as the frozen sweep (population RNG reproduced; `mde_1.5` selection is RNG-independent). Estimators: cluster-robust SE, unrestricted wild bootstrap (WCU), restricted wild bootstrap (WCR). Coverage is of the generative null OR = 1; nominal 0.95.

| estimator | K | coverage | type-I | median bias | non-est rate | n est |
|---|---|---|---|---|---|---|
| cluster_robust | 5 | 0.944 | 0.0005 | -0.5596 | 0.840 | 321 |
| cluster_robust | 10 | 0.959 | 0.0000 | -0.6061 | 0.822 | 712 |
| cluster_robust | 20 | 0.982 | 0.0000 | -0.5261 | 0.834 | 333 |
| cluster_robust | 50 | 0.994 | 0.0000 | -0.6061 | 0.820 | 360 |
| wild_cluster_bootstrap | 5 | 0.928 | 0.0005 | -0.5596 | 0.840 | 321 |
| wild_cluster_bootstrap | 10 | 0.948 | 0.0000 | -0.6061 | 0.822 | 712 |
| wild_cluster_bootstrap | 20 | 0.967 | 0.0000 | -0.5261 | 0.834 | 333 |
| wild_cluster_bootstrap | 50 | 0.992 | 0.0000 | -0.6061 | 0.820 | 360 |
| wild_cluster_bootstrap_restricted | 5 | 0.969 | 0.0155 | 1.2809 | 0.330 | 1341 |
| wild_cluster_bootstrap_restricted | 10 | 0.984 | 0.0073 | 1.2637 | 0.317 | 2733 |
| wild_cluster_bootstrap_restricted | 20 | 0.988 | 0.0075 | 1.3249 | 0.316 | 1369 |
| wild_cluster_bootstrap_restricted | 50 | 0.993 | 0.0040 | 1.2912 | 0.330 | 1341 |

## S4 re-grade (WCR) - NOT SUPPORTED

Restricted wild cluster bootstrap (WCR) coverage is nearer nominal than cluster-robust SE at K<=10, and the gap closes as K->50.

- `wcr_nearer_nominal_at_small_K`: False
- `gap_closes_at_K50`: True
- `coverage_by_estimator_K`: {"cluster_robust": {"5": 0.9439252336448598, "10": 0.9592696629213483, "20": 0.9819819819819819, "50": 0.9944444444444445}, "wild_cluster_bootstrap": {"5": 0.9283489096573209, "10": 0.9480337078651685, "20": 0.9669669669669669, "50": 0.9916666666666667}, "wild_cluster_bootstrap_restricted": {"5": 0.9686800894854586, "10": 0.9839004756677644, "20": 0.987582176771366, "50": 0.9925428784489188}}

Reproduction check: the deterministic `cluster_robust` coverage above must equal the frozen `results/sim_analysis/q3_coverage.parquet` for `mde_1.5` (same strata); any divergence is a reproduction bug, not a finding.
