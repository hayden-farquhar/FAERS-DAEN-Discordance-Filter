# Phase 5 Results — Family 3 mechanism arm (Arm-2 discovery universe)

Generated: 2026-05-30 23:21:30
Substrate: `results/perpair_arm2.parquet` (4,300 Arm-2 pairs after filtering)
Family 3 universe (FAERS-positive AND daen_powered): 4,300 pairs

*(Note: input row count is the daen_powered subset of the full Arm-2 substrate. The full Arm-2 substrate is 1,867,523 FAERS-positive pairs in `results/perpair_arm2.parquet`.)*

## Family 3 — Holm-Bonferroni at family-wise α = 0.05 (k = 4)

**Test framing**: for each mechanism, the contrast is the proportion of pairs exhibiting the mechanism in `cross_db_class == 'faers_only'` vs `concordant_positive` within the Arm-2 daen_powered universe. Two-sided z-test on difference of proportions; Newcombe 95% CI on the difference; Holm-Bonferroni step-down at family-wise α = 0.05 across k = 4 tests.

| Test | Mechanism | n (faers_only) | n (conc_pos) | p (faers_only) | p (conc_pos) | Δ (pp) | 95% CI (pp) | p-value | Holm decision | Directional interpretation |
|---|---|---:|---:|---:|---:|---:|---|---:|---|---|
| **H5a** | Consumer-share lift ≥ 10 pp vs ingredient baseline | 3314 | 986 | 0.214 | 0.148 | +6.62 | 3.90 – 9.14 | 0.0000 | `REJECT_AT_HOLM` | **SUPPORTS H5** (FAERS-only enriched) |
| **H5b** | Lawyer-share lift ≥ 5 pp vs ingredient baseline | 3314 | 986 | 0.010 | 0.019 | -0.93 | -2.03 – -0.13 | 0.0188 | `REJECT_AT_HOLM` | **REJECTS H5** (CONTRADICTORY direction — concordant_positive is enriched) |
| **H5c** | Bai-Perron break within ±2 quarters of FDA alert (≥16 quarters required) | 3246 | 944 | 0.000 | 0.000 | +0.00 | -0.41 – 0.12 | — | `FAIL_TO_REJECT` | no enrichment detected |
| **H5d** | Drug on US JPML mass-tort MDL registry | 3314 | 986 | 0.032 | 0.015 | +1.65 | 0.53 – 2.53 | 0.0058 | `REJECT_AT_HOLM` | **SUPPORTS H5** (FAERS-only enriched) |

## Mechanism-by-mechanism summary

**H5a — Consumer-share enrichment.** **SUPPORTS H5a.** Pairs with consumer-report-share lifted ≥ 10 pp above the drug's baseline are more frequent in `faers_only` (21.4%) than in `concordant_positive` (14.8%) by 6.6 pp (95% CI 3.9 – 9.1; p < 0.0001). REJECT_AT_HOLM in the supportive direction. The discordance pattern on the Arm-2 discovery universe is enriched for the consumer-reporting mechanism canonical to FAERS artefacts.

**H5b — Lawyer-share enrichment. REJECTS in the CONTRADICTORY direction.** Pairs with lawyer-share lifted ≥ 5 pp above baseline are *less* frequent in `faers_only` (1.0%) than in `concordant_positive` (1.9%) by −0.9 pp (95% CI −2.0 – −0.1; p = 0.019). The two-sided test technically rejects H0 of equality at the Holm-adjusted threshold, but the direction is the *opposite* of H5b's prediction. Mechanistically: the discordance pattern on this universe is NOT enriched for lawyer-stimulated reporting. Reported as `REJECTS_DIRECTIONALLY_CONTRADICTORY` in the manuscript Methods.

**H5c — Post-alert temporal-breakpoint alignment. FAIL_TO_REJECT, with detection-method caveat.** The test ran on 4,190 / 4,300 pairs (those with ≥ 16 quarters of FAERS data). The BIC-selected single-break test detected a structural break in 100% (4,190 / 4,190) of pairs, reflecting the substantial growth in FAERS reporting volume over 2004 – 2025 (~40× from start to end of window). Of 225 pairs whose drug appears in the alert registry, **zero** had a break aligned within ±2 quarters of any FDA alert for that drug. The H5c null is therefore not surprising given the BIC's leniency on a 88-quarter heteroskedastic count series — the test as pre-specified detects breaks promiscuously but those breaks reflect FAERS-wide volume growth rather than drug-specific post-alert spikes. A more discriminating H5c test (e.g., detrended against the FAERS-wide volume baseline, or focused on residuals after a Poisson-regression secular trend) would be a Phase-6 sensitivity arm.

**H5d — Mass-tort drug membership. SUPPORTS H5d.** Drugs on the US JPML MDL registry appear more often in `faers_only` (3.2%) than in `concordant_positive` (1.5%) by 1.65 pp (95% CI 0.5 – 2.5; p = 0.006). REJECT_AT_HOLM in the supportive direction. The discordance pattern is enriched for drugs that have been mass-tort-litigated, consistent with the H5d-mechanism reading that litigation-stimulated reporting drives FAERS-only signals.

## Family 3 verdict

**Two of four pre-registered mechanisms (H5a consumer reporting, H5d mass-tort membership) confirm in the supportive direction on the Arm-2 discovery universe.** H5b contradicts its predicted direction (concordant-positive is more lawyer-enriched than faers_only). H5c is null with a documented detection-method limitation. Combined, this is meaningful evidence that the canonical artefact mechanisms ARE associated with cross-database discordance at the discovery scale — specifically the consumer-reporting and mass-tort axes.

## Headline reservation rule remains binding (Section 4)

Per the protocol Section 4 operational definition of 'artefact', the headline use of that word requires both (i) H1 confirmation of reference-negative enrichment (Phase 4 result: H1 EXPLORATORY_DOWNGRADE) AND (ii) at least one H5 mechanism positive on the **Arm-1 intersection** (Phase 4 Family 4 result: 6 pairs, insufficient n). Family 3's discovery-arm findings, while supportive of H5a and H5d in the broader universe, **do not satisfy the reservation rule** because the rule requires positivity on the *same* Arm-1 pairs that drive enrichment. The headline therefore reports the Family 1 descriptive (H2 = 13.6%, CI 5.7 – 31.8%) plus the discovery-arm mechanism enrichment for H5a/H5d, but does NOT claim 'artefact' in the headline.

## Pivot framing

Combining Family 2 (H1 calibrated null) + Family 3 (H5a/H5d supportive, H5b contradictory, H5c null) + Family 4 (insufficient n): the paper's contribution is a **calibrated, partially-supportive mechanism story**. The Section 5 Family 2 pivot framing remains the active narrative for the H1 outcome, but Family 3's supportive H5a/H5d findings warrant a Discussion paragraph on the consumer/mass-tort mechanism in the broader discovery universe — with the caveat that the same mechanism could not be validated on the Arm-1 ground-truth subset (Family 4 insufficient n).
