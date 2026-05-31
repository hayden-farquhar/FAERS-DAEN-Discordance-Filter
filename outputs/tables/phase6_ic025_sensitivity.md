# Phase 6 sensitivity — alternative signal definition (IC025)

**Pre-registered Section 10.7 robustness exhibit.** Reported for H2, H1, H3. Not a confirmatory test; not multiplicity-counted (shares data with the primary).

**Signal rule:** `IC025 > 0 AND a >= 3` (BCPNN, Jeffreys 0.5 prior, delta-method interval) vs the primary `ROR 95% LB > 1 AND a >= 3`.

**Power model:** held fixed at the primary `daen_powered = (daen_mde <= 1.5)` (ROR-anchored, G2-locked). The Section 10.7 variation is the signal definition only.


## Classification flip diagnostic

- FAERS signal-status flips vs ROR rule: **2 / 492** (0.4%)
- DAEN signal-status flips vs ROR rule: **7 / 492** (1.4%)
- `cross_db_class` flips vs primary: **9 / 492** (1.8%)


## H2 — power-conditioned non-replication rate (Family 1)

- H2 (IC025) = **0.136** (95% cluster-bootstrap CI 0.057–0.318) on n=44 FAERS-positive ∧ daen_powered pairs
- Underpowered-discordant count reported alongside: 250
- PA+ = **0.916** (95% CI 0.811–0.964)
- Cohen's κ = 0.588
- Cohen's kappa = 0.588. Marginal balance is reasonable (84.6% / 75.0%); kappa interpretable but PA+ (= 0.916) remains the pre-registered primary concordance statistic.
- Cross-DB 2×2 on daen_powered universe: {'faers_pos_daen_pos': 38, 'faers_pos_daen_neg': 6, 'faers_neg_daen_pos': 1, 'faers_neg_daen_neg': 7}


## H1 — reference-negative enrichment (Family 2)

- Decision: **EXPLORATORY_DOWNGRADE**
- Stratum: 44 pairs (6 faers_only, 38 concordant_positive) across 7 events, 28 drugs
- Known-negatives: 0 of 6 faers_only; 5 of 38 concordant_positive
- Fixed-effect OR = 0.000 (95% CI 0.000–nan; p = nan)
- Cluster-robust OR = nan (95% CI nan–nan; p = nan)
- Realised MDE = 76.51
  - Separation: known-neg counts in (faers_only, concordant_positive) = (0, 5) of (6, 38). Standard mixed-effects logistic OR is undefined or infinite.
  - Logit fit failed: math range error
  - Realised MDE = 76.51 > 3.0 ⇒ H1 downgraded to exploratory per Section 6.6; Family 2 pivot (Section 5) fires.

## H3 — discordance-filter LR+ (Family 2)

- Decision: **DEGENERATE**
- sensitivity = 0.0, specificity = 0.8461538461538461, PPV = 0.0, NPV = 0.868421052631579
- LR+ = nan (95% CI nan–nan)
  - LR+ degenerate: sensitivity=0.0, specificity=0.8461538461538461, tp/fp/fn/tn=0/6/5/33. One of the operating-characteristic boundaries was hit.

## Within-event Mantel–Haenszel (Section 10.5)

- Strata total: 7; informative: 3
- MH OR = nan (95% CI nan–nan)
  - MH OR undefined: no informative stratum had both rows AND both columns nonzero. This is the expected pattern under severe degeneracy in the H1 stratum.
