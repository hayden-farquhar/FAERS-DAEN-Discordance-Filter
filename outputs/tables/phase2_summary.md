# Phase 2 Power-Model Augmentation — Substrate Summary

Pipeline: `scripts/02_power_augment.py`
Output:   `results/perpair_arm1.parquet` (492 rows, 45 columns)
Backup of Phase-1 substrate: `results/perpair_arm1.phase1.parquet`

## Hard-gate status
- **G1 (metrics module clean):** PASSED (Phase 1).
- **G2 (power model locked):** **PASSED.** 12/12 tests in `tests/test_power_validation.py` pass:
  - Analytic-vs-MC agreement: exact-PMF and MC agree to within 0.02 across the 200-cell validation grid.
  - Hand-checkable pairs: 5/5 match exact closed-form to within MC sampling noise (~0.021 at B=5000).
  - MC reproducibility under seed=94 confirmed bit-for-bit.
  - Wald-approximation envelope characterised (informational; Wald NOT the authoritative analytic path).
- **G3 (no headline without power-conditioning):** **IN EFFECT.** This document inspects the substrate; it does NOT quote any headline number (no H2 non-replication rate, no H1 enrichment OR, no H3 LR+).

## daen_powered_* tag distributions

| Tag | N TRUE | % of 492 |
|---|---|---|
| daen_powered (PRIMARY, MDE <= 1.5)            | 52 | 10.6% |
| daen_powered_at_2 (MDE <= 2.0)                | 131 | 26.6% |
| daen_powered_at_3 (MDE <= 3.0)                | 223 | 45.3% |
| daen_powered_07 (MDE <= 1.5 at power 0.70)    | 66 | 13.4% |
| daen_powered_09 (MDE <= 1.5 at power 0.90)    | 37 | 7.5% |
| daen_powered_FAERS_pt (FAERS-point-derived)   | 127 | 25.8% |
| daen_powered_FAERS_LB (FAERS-LB-derived)      | 116 | 23.6% |
| daen_powered_unconditional (Poisson-margin)   | 51 | 10.4% |
| daen_powered_min (conservative composite)     | 51 | 10.4% |
| faers_powered (symmetric, for H4)             | 331 | 67.3% |

## FAERS-derived-power circularity diagnostic
Primary `daen_powered` (FAERS-independent MDE):    52 pairs
`daen_powered_FAERS_pt` (uses FAERS point ROR):    127 pairs (gap +75; +15.2 pp)
`daen_powered_FAERS_LB` (uses FAERS 95% LB):       116 pairs (gap +64; +13.0 pp)

The FAERS-derived sensitivities classify roughly 2.4x as many pairs as adequately
powered. Had the FAERS-derived approach been the primary (as in the earlier-design
formulation), pairs with FAERS-inflated ROR estimates would have been classified
as powered on the basis of an effect that may not be real. The Round-2 must-fix
(MDE-based primary) breaks this loop.

## Poisson-margin attrition rule (Section 6.4)
- daen_powered:     52 pairs
- daen_powered_min: 51 pairs
- Attrition: 1.9%
- Rule status: attrition below the 50% / 20-pair threshold, so the > 5pp fragility flag remains armed for the `daen_powered_min` sensitivity at Phase 4 readout.

## H1 stratum composition (FAERS-positive AND daen_powered; Family 2 universe)
Total H1 pairs: 44
By cross_db_class x groundTruth:
groundTruth          0   1
cross_db_class            
concordant_positive  5  33
faers_only           0   6

Unique events in H1 stratum: 7 (out of 10 in the universe)
Unique drugs in H1 stratum:  28

## Family 4 substrate (Arm-1 FAERS-only intersection)
Pairs in Arm-1 FAERS-only intersection: **6**
Protocol Section 4 threshold for confirmatory mechanism analysis: n >= 20
**Status: STRATEGIC FORK FIRES**

This was pre-recorded in the protocol Section 4 strategic-fork acknowledgement:
the Family 4 n>=20 reservation threshold is not met, so the headline will fall back to
'reference-negative enrichment without demonstrated mechanism'. The word 'artefact'
is NOT used in the eventual paper headline per the binding reservation rule.

## Next phase
Phase 3 (replication classification with the daen_powered tag) is now technically
**already complete** because the per-pair output now carries the cross-database
classification AND the power tags — the Phase 3 deliverable is already in
`results/perpair_arm1.parquet`. Phase 4 (Family 2 confirmatory: H1 mixed-effects
logistic + H3 LR+) can begin.

The Phase 2 timeline gate G2 is closed; the design is locked; the analysis can
proceed without further power-model work.
