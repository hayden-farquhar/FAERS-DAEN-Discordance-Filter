# Phase 1 Arm-1 Substrate — Build Summary

Pipeline: `src/replication/build_arm1.py`
Output:   `results/perpair_arm1.parquet` (492 rows)
Generated: 2026-05-30 14:39

## Universe (pooled OMOP + EU-ADR per protocol Section 9.1)
- Total per-pair rows:    492
- OMOP subset:            399
- EU-ADR subset:          93
- Unique exposures:       232
- Unique outcomes:        10

## Drug resolution (Phase 1 coverage check)
- Exposures with ZERO FAERS cases (drug never resolved): 4 / 232 (1.7%)
- Exposures with ZERO DAEN cases (drug never resolved):  30 / 232 (12.9%)

## Signal counts (PRE-power-conditioning; G3 in effect)
- FAERS signals:  292  (59.3% of pairs)
- DAEN  signals:  155   (31.5% of pairs)

## Cross-DB class distribution

```
cross_db_class
concordant_negative    185
faers_only             152
concordant_positive    140
daen_only               15
```

## Class x groundTruth (descriptive only; NOT a confirmatory result)

```
groundTruth            0    1
cross_db_class               
concordant_negative  150   35
concordant_positive   28  112
daen_only              9    6
faers_only            97   55
```

## Exposures with ZERO DAEN cases (30)

These are likely drug-name resolution failures (DAEN's active_ingredient column
uses different naming conventions for some entries). Each is flagged for review
at Phase 2 before downstream analysis proceeds. The companion US-AU name crosswalk
in `data/reference/drug_name_crosswalk.csv` covers only 40 mappings; expanding it
may recover some of these.

  - Amoxapine
  - Amylases
  - Capreomycin
  - Carteolol
  - Chlorazepate
  - Cosyntropin
  - Dicyclomine
  - Enalaprilat
  - Endopeptidases
  - Estrogens, Conjugated (USP)
  - Etodolac
  - Factor VIIa
  - Flavoxate
  - Pemoline
  - Propafenone
  - Regular Insulin, Human
  - Salsalate
  - Sodium Phosphate, Monobasic
  - Tolmetin
  - almotriptan
  - benzonatate
  - bromfenac
  - frovatriptan
  - gemifloxacin
  - lithium citrate
  - moexipril
  - oxaprozin
  - ramelteon
  - trovafloxacin
  - valdecoxib

## Hard-gate status

- **G1 (metrics module clean):** PASSED — 20/20 metric unit tests pass against hand-worked examples (`tests/test_*.py`); metrics module SHA-256 hashes verified against the protocol Section 15.2 pin.
- **G2 (power model locked):** Phase 2, not yet attempted.
- **G3 (no headline without power-conditioning):** In effect — no replication-rate number, no H1 enrichment OR, no H3 likelihood ratio is quoted from this substrate. The cross-DB class distribution above is the Phase 1 substrate, not the headline.

## Next step

Phase 2: implement the per-pair MDE-based power model (protocol Section 6.1).
Once G2 closes, the `daen_powered` flag is added to the per-pair output and the
power-conditioned classification (Family 1, Family 2) becomes computable.
