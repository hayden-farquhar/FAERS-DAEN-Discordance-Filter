# Phase 6.3 — Per-reference-set H2 + H1 sensitivity

Pre-registered sensitivity arm per Section 10.7 (ground-truth set row). Reports H2 (power-conditioned non-replication rate) and H1-stratum composition under each reference-set restriction.

**TGA-action set is correctly EXCLUDED** from this exhibit per Section 9.1 circularity caveat (TGA acted because the signal appears in Australian data, so TGA-action pairs are by construction DAEN-visible).

## H2 power-conditioned non-replication rate by reference set

| Reference set | Total pairs | H2 universe (FAERS+ ∧ daen_powered) | H2 rate | 95% cluster-boot CI |
|---|---|---|---|---|
| Primary (pooled OMOP + EU-ADR) | 492 | 44 | 0.136 | 0.057 – 0.318 |
| OMOP only | 399 | 23 | 0.130 | 0.000 – 1.000 |
| EU-ADR only | 93 | 21 | 0.143 | 0.000 – 0.318 |

## H1 stratum composition by reference set

Cells = pairs in {faers_only, concordant_positive} × {known-neg, known-pos}, restricted to FAERS-positive AND daen_powered.

### Primary (pooled OMOP + EU-ADR)

Stratum n: 44; events represented: 7

```
groundTruth          0   1
cross_db_class            
concordant_positive  5  33
faers_only           0   6
```

### OMOP only

Stratum n: 23; events represented: 3

```
groundTruth          0   1
cross_db_class            
concordant_positive  1  19
faers_only           0   3
```

### EU-ADR only

Stratum n: 21; events represented: 7

```
groundTruth          0   1
cross_db_class            
concordant_positive  4  14
faers_only           0   3
```

## Interpretation

- The headline H2 = 0.136 (CI 0.057–0.318) on the pooled set is reported alongside per-set rates for transparency.
- Differences between OMOP-only and EU-ADR-only H2 estimates characterise whether the pooled headline is driven by one source.
- The pre-registered decision rule per Section 10.7 declares the headline robust if H1 retains direction and significance across the sensitivity cells. Since H1 itself was EXPLORATORY_DOWNGRADE (Phase 4), this exhibit reports H1 stratum composition only — there is no per-set H1 OR to test.
