# Post-hoc: reference-negative-density (n_pairs) sweep

Base cell (pi=0.5, lambda=2, phi=0.5), primary estimator `wild_cluster_bootstrap_fast`, 1000 replicates per (regime x rule x n_pairs). n_pairs grows from the calibrated 492 (1x) to 8x; everything else is held at the slice values. Type-I band [0.0, 0.075].

Both estimands are shown: `marginal` scores non-estimable replicates as non-rejections (the pre-registered endpoint); `conditional` divides rejections by estimable replicates only. As n_pairs grows the stratum fills in, non-estimability falls, and the two estimands converge.

## Q1_null (type-I)

The rejection rate is the type-I error. The question is whether each rule stays in band as the stratum fills. `faers_point` is the winner's-curse-contaminated rule: if its type-I stays in band even at 8x density, the inflation hypothesis (S1) is dead, not merely masked by thin-stratum non-estimability.

| rule | n_pairs | non-est rate | n est | marginal type-I | conditional type-I | median stratum | median faers-only |
|---|---|---|---|---|---|---|---|
| faers_point | 492 | 0.500 | 500 | 0.0010 | 0.0020 | 92 | 89 |
| faers_point | 984 | 0.186 | 814 | 0.0040 | 0.0049 | 184 | 179 |
| faers_point | 1968 | 0.034 | 966 | 0.0200 | 0.0207 | 369 | 359 |
| faers_point | 3936 | 0.002 | 998 | 0.0430 | 0.0431 | 740 | 718 |
| mde_1.5 | 492 | 0.829 | 171 | 0.0000 | 0.0000 | 41 | 40 |
| mde_1.5 | 984 | 0.577 | 423 | 0.0010 | 0.0024 | 82 | 79 |
| mde_1.5 | 1968 | 0.244 | 756 | 0.0030 | 0.0040 | 165 | 161 |
| mde_1.5 | 3936 | 0.057 | 943 | 0.0140 | 0.0148 | 330 | 321 |
| none | 492 | 0.044 | 956 | 0.0110 | 0.0115 | 332 | 322 |
| none | 984 | 0.004 | 996 | 0.0370 | 0.0371 | 665 | 646 |
| none | 1968 | 0.000 | 1000 | 0.0570 | 0.0570 | 1329 | 1292 |
| none | 3936 | 0.000 | 1000 | 0.0520 | 0.0520 | 2658 | 2585 |

## Q2_alt (power)

The enrichment alternative holds, so the rejection rate is power (higher is better; no band). The question is whether the registered `mde_1.5` rule ever becomes usefully powered as the reference set grows, or is permanently dominated on power by the contaminated `faers_point` rule.

| rule | n_pairs | non-est rate | n est | marginal power | conditional power | median stratum | median faers-only |
|---|---|---|---|---|---|---|---|
| faers_point | 492 | 1.000 | 0 | 0.0000 | nan | 92 | 33 |
| faers_point | 984 | 0.999 | 1 | 0.0010 | 1.0000 | 184 | 67 |
| faers_point | 1968 | 0.999 | 1 | 0.0010 | 1.0000 | 370 | 132 |
| faers_point | 3936 | 0.998 | 2 | 0.0020 | 1.0000 | 740 | 266 |
| mde_1.5 | 492 | 1.000 | 0 | 0.0000 | nan | 41 | 14 |
| mde_1.5 | 984 | 1.000 | 0 | 0.0000 | nan | 82 | 28 |
| mde_1.5 | 1968 | 0.999 | 1 | 0.0010 | 1.0000 | 166 | 55 |
| mde_1.5 | 3936 | 1.000 | 0 | 0.0000 | nan | 329 | 111 |
| none | 492 | 0.815 | 185 | 0.1850 | 1.0000 | 331 | 227 |
| none | 984 | 0.638 | 362 | 0.3620 | 1.0000 | 662 | 454 |
| none | 1968 | 0.380 | 620 | 0.6200 | 1.0000 | 1326 | 908 |
| none | 3936 | 0.153 | 847 | 0.8470 | 1.0000 | 2654 | 1818 |

## Roll-up

- `Q1_mde_1.5_in_band_all_densities_both_estimands`: True
- `Q1_faers_point_in_band_all_densities_both_estimands`: True
- `Q1_none_conditional_typeI_range`: [0.0115, 0.0570]
- `Q2_mde_1.5_marginal_power_range`: [0.0000, 0.0010]
- `Q2_faers_point_marginal_power_range`: [0.0000, 0.0020]
- `Q2_faers_point_dominates_mde_1.5_on_power_all_densities`: True
