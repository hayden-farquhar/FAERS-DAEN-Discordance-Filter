# Phase 4 — H4 balance + Family 4 formal report

## H4 — DAEN-only direction (Section 4 H4; Section 10.6)

**Universe:** `cross_db_class == 'daen_only'` AND `faers_powered` (symmetric power check at OR=1.5 on FAERS margins).
**N pairs in H4 stratum:** 11

| Drug | Outcome | groundTruth | ROR_F | a_F | a_D | ROR_D | DAEN MDE |
|---|---|---|---|---|---|---|---|
| Hyoscyamine | OMOP Acute Liver Failure 1 | 0 | 1.10 | 107 | 5 | 2.99 | 4.05 |
| oxybutynin | OMOP Acute Liver Failure 1 | 0 | 0.76 | 305 | 28 | 1.95 | 1.90 |
| Primidone | OMOP Acute Liver Failure 1 | 0 | 1.08 | 203 | 11 | 2.36 | 2.70 |
| Nevirapine | OMOP Acute myocardial Infarction  1 | 0 | 0.96 | 112 | 7 | 6.94 | inf |
| fluticasone | HOI Upper GI #3 | 0 | 1.04 | 1304 | 15 | 1.96 | 2.25 |
| Clozapine | OMOP Acute Liver Failure 1 | 1 | 0.78 | 1823 | 594 | 1.36 | 1.15 |
| Interferon beta-1a | OMOP Acute Liver Failure 1 | 1 | 0.96 | 1656 | 33 | 2.24 | 1.85 |
| pioglitazone | OMOP Acute Liver Failure 1 | 1 | 0.89 | 470 | 19 | 2.90 | 2.25 |
| darbepoetin alfa | OMOP Acute myocardial Infarction  1 | 1 | 0.90 | 405 | 6 | 2.30 | 3.05 |
| Levodopa | OMOP Acute Renal Failure 1 | 0 | 0.77 | 194 | 4 | 3.22 | 4.45 |
| dorzolamide | HOI Upper GI #3 | 0 | 1.18 | 96 | 3 | 5.57 | inf |

## Family 4 — Arm-1 mechanism intersection (Section 4 + Section 10.4)

**Universe:** Arm-1 FAERS-positive AND `daen_powered` AND `cross_db_class == 'faers_only'`.
**N pairs:** 6
**Pre-registered reservation threshold:** n ≥ 20.
**Status: INSUFFICIENT N** ⇒ Family 4 reports 'insufficient n for mechanism characterisation' per Section 4.

### Binding reservation rule fires

Per the pre-registered Section 4 operational definition of 'artefact', the headline can use that word ONLY when (i) H1 confirms reference-negative enrichment AND (ii) at least one H5 mechanism is positive on the Arm-1 intersection. Family 4 is insufficient-n; the rule binds.

The 6 pairs in the intersection (each a known-positive published association):

| Drug | Outcome | groundTruth | ROR_F | DAEN MDE |
|---|---|---|---|---|
| Captopril | OMOP Acute Liver Failure 1 | 1 | 2.09 | 1.30 |
| infliximab | OMOP Acute Liver Failure 1 | 1 | 4.57 | 1.35 |
| Sertraline | HOI Upper GI #3 | 1 | 1.61 | 1.40 |
| Aspirin | Anaphylaxis #1 | 1 | 1.11 | 1.20 |
| Ciprofloxacin | Anaphylaxis #1 | 1 | 1.94 | 1.50 |
| Captopril | Leukopenia Including Neutropenia and Agranulocytosis | 1 | 1.17 | 1.35 |

These are well-established published drug-event associations that DAEN does not signal despite being adequately powered. Directionally informative for the Discussion: the powered FAERS-only stratum is enriched for **known-positive** pairs, not the **known-negative** pairs the H1 thesis predicted. This substrate observation reinforces the Family 2 H1 pivot.
