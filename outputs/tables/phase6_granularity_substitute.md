# Phase 6 sensitivity — event-granularity substitute (index-PT)

**Documented deviation from the Section 10.7 HLT arm.** The literal HLT-level aggregation requires the licensed MedDRA v27.1 PT→HLT hierarchy, which is not deposited (fabricating it would breach the protocol ontology pins); and the Arm-1 events are OMOP/EU-ADR baskets of 3–41 PTs that already span their HLTs, so a literal roll-up is a near-no-op in Arm 1. This exhibit probes the same granularity-robustness concern in the executable direction: it makes events **finer** — each event re-defined as its single index PT (max-FAERS-count PT in the basket) — and re-runs H2 + H1 with the G2-locked power model **recomputed** on the finer marginals. Reported for H2, H1 per Section 10.7.


## Index PT chosen per outcome

| Outcome (basket) | Index PT (finer event) |
|---|---|
| OMOP Acute Liver Failure 1 | hepatic enzyme increased |
| OMOP Acute Renal Failure 1 | acute kidney injury |
| OMOP Acute myocardial Infarction  1 | myocardial infarction |
| HOI Upper GI #3 | gastrointestinal haemorrhage |
| OMOP Aplastic Anemia 1 | pancytopenia |
| Anaphylaxis #1 | anaphylactic reaction |
| Stevens-Johnson Syndrome #1 | stevens-johnson syndrome |
| Leukopenia Including Neutropenia and Agranulocytosis | neutropenia |
| Rhabdomyolysis #1 | rhabdomyolysis |
| Cardiac Valve Fibrosis #1 | mitral valve incompetence |


**daen_powered universe:** index-PT n=13 vs primary basket n=52 (finer events reduce DAEN event-side counts, so fewer pairs are powered).


## H2 — power-conditioned non-replication rate (finer events)

- H2 (index-PT) = **0.333** (95% cluster-bootstrap CI 0.000–0.500) on n=6 FAERS-positive ∧ daen_powered pairs
- (primary basket H2 = 0.136 [0.057–0.318] on n=44)
- Underpowered-discordant: 240
- PA+ = 0.727 (95% CI 0.667–1.000)
- Cohen's κ = 0.530
- Cross-DB 2×2 (daen_powered): {'faers_pos_daen_pos': 4, 'faers_pos_daen_neg': 2, 'faers_neg_daen_pos': 1, 'faers_neg_daen_neg': 6}


## H1 — reference-negative enrichment (finer events)

- Decision: **EXPLORATORY_DOWNGRADE**
- Stratum: 6 pairs (2 faers_only, 4 concordant_positive); known-neg 0/2 vs 0/4
- Realised MDE = 448.50
  - Separation: known-neg counts in (faers_only, concordant_positive) = (0, 0) of (2, 4). Standard mixed-effects logistic OR is undefined or infinite.
  - Logit fit failed: math range error
  - Realised MDE = 448.50 > 3.0 ⇒ H1 downgraded to exploratory per Section 6.6; Family 2 pivot (Section 5) fires.
