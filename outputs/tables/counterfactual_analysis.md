# Counterfactual analysis — what would the naive FAERS-derived-power approach have shown?

**The single most impactful finding in this paper.** This counterfactual quantifies, on the same data, how much the H1 result depends on the choice between (i) the pre-registered MDE-anchored power model (FAERS-independent) and (ii) the natural-default FAERS-derived power model that any other research team would likely have chosen.

Generated: 2026-05-31. Source: `scripts/03_phase4_confirmatory.py` machinery applied to `results/perpair_arm1.parquet` with the `daen_powered` definition varied. Reported under the Section 16.3 post-hoc-disclosure rule (this counterfactual exhibit was not pre-specified; it is presented here as a methodological-cautionary demonstration, not as a confirmatory result, and is flagged as such in the manuscript Methods).

## Headline comparison

| Specification | H2 rate | H2 n | H1 stratum | known-neg / faers_only | H1 OR (fixed-effect logistic) | H1 OR (cluster-robust by event) | H1 p | LR+ |
|---|---|---|---|---|---|---|---|---|
| **PRIMARY: MDE-anchored (daen_mde ≤ 1.5)** — pre-registered, FAERS-independent | 0.136 (0.057–0.318) | 44 | 6 faers_only + 38 conc_pos | **0 / 6** | **0** (degenerate) | undefined | — | undefined |
| **COUNTERFACTUAL: FAERS-point-ROR-derived** — natural default | 0.135 (0.101–0.257) | 126 | 17 faers_only + 109 conc_pos | 7 / 17 | **5.17** (1.68–15.94) | **5.17** (1.23–21.72) | **0.004 / 0.025** | **3.71** (1.60–8.59) |
| **FAERS-95%-LB-derived** — more conservative FAERS-derived | 0.112 (0.087–0.189) | 116 | 13 faers_only + 103 conc_pos | 4 / 13 | 3.37 (0.90–12.66) | 3.37 (0.81–14.12) | 0.072 / 0.096 | 2.78 (0.97–7.96) |
| **MDE ≤ 2.0** (the earlier-revision anchor) | 0.214 (0.156–0.383) | 112 | 24 faers_only + 88 conc_pos | 10 / 24 | 4.12 (1.51–11.23) | 4.12 (2.13–7.98) | 0.006 / <0.001 | 2.76 (1.42–5.40) |

## The methodological cautionary

**The naive FAERS-point-ROR-derived analysis returns a publication-grade significant confirmation of H1 (OR = 5.17, p = 0.004) supporting the artefact-filter thesis.** Both the fixed-effect logistic and the cluster-robust specification reject H0; the LR+ is 3.71. A reader of a paper using this approach would conclude that cross-database discordance is a reasonably-strong predictor of reference-negative status — a sufficient finding to warrant adoption of discordance-as-filter in practice.

**The pre-registered MDE-anchored primary correctly identifies the same data as H1-degenerate.** The truly-powered FAERS-only stratum (n = 6) contains zero known-negatives; the realised MDE is 76.51, well above the Section 6.6 downgrade threshold of 3.0; H1 is formally downgraded to exploratory and the Family 2 pivot fires.

**The discrepancy is entirely the FAERS-derived-power circularity.** The naive approach asks: "is DAEN powered to detect the FAERS effect?" When FAERS is artefactually inflated for some pairs (the very pairs the test is designed to detect), the FAERS-derived power answer is *yes* — those pairs are classified as `daen_powered` on the basis of an inflated effect that may not be real. They then enter the H1 stratum, where DAEN of course fails to replicate the (non-existent) inflated effect, and they show up as `faers_only` known-negatives. The naive H1 confirms exactly the artefact it was designed to detect — but the confirmation is bookkeeping, not biology.

## The cell-count diagnosis

The H1 stratum cells, comparing MDE-anchored vs FAERS-point-derived:

|  | MDE-anchored (PRIMARY) | FAERS-point-derived (COUNTERFACTUAL) |
|---|---|---|
| faers_only — known-negative | **0** | 7 |
| faers_only — known-positive | 6 | 10 |
| concordant_positive — known-negative | 5 | 13 |
| concordant_positive — known-positive | 33 | 96 |

The 17 FAERS-derived-powered `faers_only` pairs include **7 known-negatives that are not in the MDE-anchored stratum** — these are pairs where FAERS-derived power said "DAEN can detect this" because FAERS's effect estimate was large, but in fact DAEN cannot reliably detect a true effect of 1.5 at these marginals (MDE > 1.5). The 7 known-negatives drive the naive OR = 5.17.

## Methodological lesson for the field

Any future cross-database disproportionality study that asks "did this signal replicate?" must ground the power calculation in **what the second database can actually detect**, not in **what the first database appears to show**. The FAERS-derived-power approach is contaminated by exactly the inflation it purports to detect. The MDE-anchored approach (Fisher noncentral hypergeometric power computed against a fixed clinically-meaningful effect using only the second database's marginals) breaks this circularity by construction.

Concretely, the Section 6.1 specification recommended:

1. For each pair, compute the per-pair minimum-detectable-effect `daen_mde` = smallest OR at which the second database has ≥ 80% power to fire the signal rule, using `scipy.stats.nchypergeom_fisher` PMF summation over the support (validated against Monte-Carlo with `seed = 94` to within sampling noise).
2. The `daen_powered` tag is `daen_mde ≤ threshold` where threshold is chosen to sit just above the signal-rule boundary (here: 1.5 against a LB > 1 rule).
3. Retain the FAERS-derived-power tag as a *sensitivity* arm only, specifically to expose the gap from the primary as the inflation diagnostic.

## Citation

When citing this counterfactual exhibit, the recommended phrasing is:

> "Had the natural-default FAERS-derived-power approach been used in place of the pre-registered MDE-anchored primary, the H1 reference-negative enrichment test on the same data would have returned a significant confirmation (OR = 5.17, 95% CI 1.68–15.94, p = 0.004; cluster-robust OR = 5.17, 95% CI 1.23–21.72, p = 0.025). The MDE-anchored analysis correctly identifies the same data as H1-degenerate. The seven known-negative `faers_only` pairs that drive the naive 'positive' result are pairs where FAERS-derived power classifies the pair as adequately-powered on the basis of an inflated FAERS effect — the FAERS-derived-power circularity that the MDE-anchored approach was designed to break. (Farquhar 2026, post-hoc counterfactual exhibit per Section 16.3, OSF DOI 10.17605/OSF.IO/ZQ97A)."
