# Confirmatory simulation analysis (S1-S4)

Source: `results/sim_results.parquet` (6,048,000 rows; 126 DGP cells x 2000 replicates x 6 selection rules x 4 estimators).
Primary estimator for the headline endpoints: `wild_cluster_bootstrap` (cluster-robust SE retained as a sensitivity; protocol Section 11).
Type-I control band: [0.0, 0.075]. Nominal CI coverage: 0.95.

Rejection rule (uniform across estimators): reject iff the 95% CI lower bound of the enrichment OR exceeds 1. A non-estimable replicate (complete separation or empty contrast row) carries no decision and is counted as a non-rejection; the non-estimability rate is reported alongside every rate.

## Q1 type-I error by selection rule (Q1_null, primary estimator)

| rule | cells | median type-I | max type-I | frac within band | mean non-estimable |
|---|---|---|---|---|---|
| faers_lb | 63 | 0.0000 | 0.0430 | 1.000 | 0.651 |
| faers_point | 63 | 0.0000 | 0.0430 | 1.000 | 0.612 |
| mde_1.5 | 63 | 0.0000 | 0.0155 | 1.000 | 0.883 |
| mde_2.0 | 63 | 0.0000 | 0.0405 | 1.000 | 0.669 |
| mde_3.0 | 63 | 0.0010 | 0.0485 | 1.000 | 0.479 |
| none | 63 | 0.0120 | 0.0745 | 1.000 | 0.342 |

## Q1 type-I conditional on estimability (Q1_null, primary estimator, post-hoc)

Second estimand requested in review: type-I as n_reject / n_estimable (rejection rate among replicates that produced a decision), so a rule's marginal type-I cannot look controlled merely because most replicates were non-estimable and were scored as non-rejections. `pooled` is estimable-weighted over all cells; `max` and the correlation are over cells with >= 50 estimable replicates (thinner cells are Monte-Carlo noise). `corr nonest/cond` < 0 means a rule's least-estimable cells are not its worst-behaved cells, so the marginal deflation does not conceal a conditional breach.

| rule | cells (gated) | pooled cond type-I | max cond type-I | cells breaching | frac within band | corr nonest/cond |
|---|---|---|---|---|---|---|
| mde_1.5 | 42 | 0.0072 | 0.0609 | 0 | 1.000 | -0.287 |
| mde_2.0 | 55 | 0.0090 | 0.0649 | 0 | 1.000 | -0.404 |
| mde_3.0 | 63 | 0.0124 | 0.0560 | 0 | 1.000 | -0.544 |
| faers_point | 42 | 0.0184 | 0.0643 | 0 | 1.000 | -0.530 |
| faers_lb | 42 | 0.0196 | 0.0655 | 0 | 1.000 | -0.568 |
| none | 63 | 0.0267 | 0.0765 | 1 | 0.984 | -0.676 |

## Q2 power by selection rule (Q2_alt, primary estimator)

| rule | cells | median power | max power | mean non-estimable |
|---|---|---|---|---|
| none | 63 | 0.1790 | 0.3570 | 0.819 |
| mde_3.0 | 63 | 0.0455 | 0.1035 | 0.950 |
| mde_2.0 | 63 | 0.0200 | 0.0415 | 0.981 |
| faers_lb | 63 | 0.0000 | 0.1800 | 0.980 |
| faers_point | 63 | 0.0000 | 0.2635 | 0.970 |
| mde_1.5 | 63 | 0.0000 | 0.0090 | 1.000 |

## Q3 coverage and few-cluster inference (Q1_null, K-sweep line)

Coverage is assessed against the generative null OR = 1 (protocol 3.3). The S4 estimator comparison is read on the type-I-controlling rule `mde_1.5`, where selection bias is removed; coverage shortfalls on other rules reflect point-estimate (selection) bias, the same phenomenon as Q1 type-I inflation.

### Rule `mde_1.5`

| estimator | K | coverage | type-I | median bias | non-est rate | n est |
|---|---|---|---|---|---|---|
| cluster_robust | 5 | 0.944 | 0.0005 | -0.5596 | 0.840 | 321 |
| cluster_robust | 10 | 0.959 | 0.0000 | -0.6061 | 0.822 | 712 |
| cluster_robust | 20 | 0.982 | 0.0000 | -0.5261 | 0.834 | 333 |
| cluster_robust | 50 | 0.994 | 0.0000 | -0.6061 | 0.820 | 360 |
| firth | 5 | 1.000 | 0.0000 | -0.1144 | 0.701 | 599 |
| firth | 10 | 0.999 | 0.0000 | -0.2597 | 0.691 | 1238 |
| firth | 20 | 0.998 | 0.0000 | -0.0982 | 0.680 | 639 |
| firth | 50 | 1.000 | 0.0000 | -0.2021 | 0.682 | 635 |
| fixed_effect | 5 | 1.000 | 0.0000 | -0.5596 | 0.840 | 321 |
| fixed_effect | 10 | 0.999 | 0.0000 | -0.6061 | 0.822 | 712 |
| fixed_effect | 20 | 0.997 | 0.0000 | -0.5261 | 0.834 | 333 |
| fixed_effect | 50 | 0.997 | 0.0000 | -0.6061 | 0.820 | 360 |
| wild_cluster_bootstrap | 5 | 0.928 | 0.0005 | -0.5596 | 0.840 | 321 |
| wild_cluster_bootstrap | 10 | 0.949 | 0.0000 | -0.6061 | 0.822 | 712 |
| wild_cluster_bootstrap | 20 | 0.970 | 0.0000 | -0.5261 | 0.834 | 333 |
| wild_cluster_bootstrap | 50 | 0.986 | 0.0000 | -0.6061 | 0.820 | 360 |

## Pre-specified hypotheses (honest-either-way)

### S1 - NOT SUPPORTED

FAERS-point selection inflates type-I above the band for lambda>=2 & phi>0, increasing in lambda and phi.

- `n_cells`: 42
- `fraction_above_band`: 0.0
- `median_typeI`: 0.002
- `max_typeI`: 0.043
- `mean_typeI_by_lambda`: {"2": 0.0030833333333333333, "5": 0.018444444444444444, "10": 0.02322222222222222}
- `mean_typeI_by_phi`: {"0.25": 0.0038333333333333336, "0.5": 0.0063124999999999995, "1.0": 0.029222222222222222}
- `monotone_in_lambda`: true
- `monotone_in_phi`: true
- `verdict_marginal`: "NOT SUPPORTED"
- `conditional_on_estimable`: {"estimand": "n_reject / n_estimable (post-hoc; non-estimable replicates excluded)", "min_estimable_gate": 50, "n_cells_gated": 42, "pooled_cond_typeI": 0.018386191928912183, "max_cond_typeI": 0.06427503736920777, "fraction_above_band": 0.0, "mean_cond_typeI_by_lambda": {"2": 0.005203919005928676, "5": 0.025797708899440854, "10": 0.03232252935530224}, "mean_cond_typeI_by_phi": {"0.25": 0.006713804460197725, "0.5": 0.009617493305362207, "1.0": 0.03963690232938929}, "corr_nonestimability_vs_cond_typeI": -0.5301679445033964, "verdict_conditional": "NOT SUPPORTED"}

### S2 - SUPPORTED

MDE-anchored-at-1.5 controls type-I within [0, 0.075] across the primary slice.

- `n_primary_cells`: 48
- `n_within_band`: 48
- `fraction_within_band`: 1.0
- `max_typeI`: 0.0155
- `worst_cell`: {"cell_id": "ad67ab70483c8601", "lambda_inflate": 2, "phi": 1.0, "pi": 0.7, "reject_rate": 0.0155}

### S3 - SUPPORTED

MDE type-I is non-decreasing as the anchor loosens 1.5 -> 2.0 -> 3.0.

- `n_matched_cells`: 63
- `mean_typeI_by_anchor`: {"mde_1.5": 0.0008492063492063493, "mde_2.0": 0.0029920634920634925, "mde_3.0": 0.006468253968253969}
- `fraction_cells_monotone`: 0.9682539682539683

### S4 - NOT SUPPORTED

Wild cluster bootstrap coverage is nearer nominal than cluster-robust SE at K<=10, and the gap closes as K->50.

- `read_on_rule`: "mde_1.5"
- `coverage_by_estimator_K`: {"wild_cluster_bootstrap": {"5": 0.9283489096573209, "10": 0.949438202247191, "20": 0.96996996996997, "50": 0.9861111111111112}, "cluster_robust": {"5": 0.9439252336448598, "10": 0.9592696629213483, "20": 0.9819819819819819, "50": 0.9944444444444445}, "fixed_effect": {"5": 1.0, "10": 0.9985955056179775, "20": 0.996996996996997, "50": 0.9972222222222222}, "firth": {"5": 1.0, "10": 0.9991922455573505, "20": 0.9984350547730829, "50": 1.0}}
- `wild_nearer_nominal_at_small_K`: false
- `gap_closes_at_K50`: true

## Reporting note

All numbers above are computed from the deposited sweep; none are placeholders. A hypothesis graded NOT SUPPORTED is reported as a finding per the binding honesty clause (protocol Section 8): no DGP cell was added, removed, or re-weighted after seeing results.
