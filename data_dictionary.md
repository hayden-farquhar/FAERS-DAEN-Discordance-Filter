# Data dictionary

Two classes of data feed this study: small **public reference datasets** (included
in `data/reference/`) and large **raw spontaneous-report substrates** (not
redistributed; see `data/raw/README.md`). This file documents the reference
datasets and the key columns of the per-pair output tables the pipeline produces.

---

## Reference datasets (`data/reference/`)

### `omop_reference_set.csv` and `euadr_reference_set.csv`
Drug–outcome ground-truth reference sets in OMOP format. Each row is a
drug–outcome pair labelled as a known positive or negative control.

| Column | Type | Description |
|---|---|---|
| `exposureId` / `exposureName` | int / str | Drug (RxNorm/OMOP concept id and name) |
| `outcomeId` / `outcomeName` | int / str | Outcome (OMOP outcome cohort id and name) |
| `groundTruth` | int | 1 = known positive association; 0 = negative control |
| `indicationId` / `indicationName` | int / str | Drug indication (for confounding control) |
| `comparatorId` / `comparatorName` | int / str | Comparator drug |
| `comparatorType` | str | Comparator selection rule |

Provenance: OMOP/OHDSI drug-safety reference set (Ryan et al.) and the EU-ADR
reference set (Coloma et al.; Avillach et al.), both published, openly available
control sets widely used to benchmark disproportionality methods.

### `harpaz_2014_reference_set.csv`
Independent drug–event ground-truth set.

| Column | Type | Description |
|---|---|---|
| `EVENT_CONCEPT_NAME` | str | Adverse-event concept |
| `DRUG_CONCEPT_NAME` | str | Drug concept |
| `GROUND_TRUTH` | int | 1 = positive; 0 = negative |

### `harpaz_2014_event_definitions.csv`
MedDRA Preferred-Term definitions for the Harpaz event concepts.

| Column | Type | Description |
|---|---|---|
| `EVENT_CONCEPT_NAME` | str | Event concept |
| `MEDDRA_PT` | str | MedDRA Preferred Term string |
| `DEFINITION_LEVEL` | str | `narrow` or `broad` definition |
| `MDR_CODE` | int | MedDRA PT code |
| `UMLS_CUI` | str | UMLS concept unique identifier |

Provenance (both Harpaz files): Harpaz et al. (2014), *Drug Safety*,
"Performance of pharmacovigilance signal-detection algorithms for the FDA
Adverse Event Reporting System", published reference standard. Only the PT
strings/codes used as outcome definitions are reproduced here; the licensed
MedDRA dictionary and PT→HLT hierarchy are **not** distributed.

### `omop_outcome_to_meddra_pt_map.csv`
Maps each OMOP/EU-ADR outcome cohort to its constituent MedDRA Preferred Terms.

| Column | Type | Description |
|---|---|---|
| `outcomeName` / `outcomeId` | str / int | OMOP outcome cohort |
| `meddra_pt` | str | MedDRA Preferred Term included in the outcome basket |
| `source` | str | Mapping source (e.g. `harpaz_2014_narrow`) |
| `rationale` | str | Free-text justification for the mapping |

### `drug_name_crosswalk.csv`
US↔Australian drug-name reconciliation.

| Column | Type | Description |
|---|---|---|
| `us_name` | str | US ingredient name (FAERS) |
| `au_name` | str | Australian ingredient name (DAEN) |

### `atc4_resolved.json`
Resolved ATC level-4 class memberships. Keys: `classes` (ATC4 → member drugs),
`drugs` (drug → ATC4 classes), `unresolved` (drugs without an ATC4 mapping).

### `omop_euadr_ingredient_map.json`
Ingredient-alias map: canonical ingredient → list of alias strings used to
resolve drug names across databases.

### `alert_registry.csv`
Regulatory safety-alert timeline (for the H5c time-break alignment test).

| Column | Type | Description |
|---|---|---|
| `alert_date` | date | Date of the regulatory action/communication |
| `drug_ingredient` | str | Drug ingredient |
| `event_category` | str | Safety domain (e.g. cardiovascular) |
| `source_type` | str | Source (e.g. `fda_dsc`) |
| `citation` | str | Citation / public reference |
| `notes` | str | Free-text note |

### `mass_tort_drugs.csv`
US mass-tort / multidistrict-litigation drug registry (for the H5d litigation test).

| Column | Type | Description |
|---|---|---|
| `drug_ingredient` | str | Drug ingredient |
| `mdl_number` | str | Multidistrict-litigation docket number |
| `mdl_filing_date` | date | MDL filing date |
| `jurisdiction` | str | Court jurisdiction |
| `event_category` | str | Safety domain |
| `citation` | str | Public case reference |
| `notes` | str | Free-text note |

---

## Per-pair output tables (produced into `results/`)

`perpair_arm1.parquet` (and the Arm-2 analogue) carry one row per drug–event pair.
Key columns:

| Column | Description |
|---|---|
| `ref_set` | `OMOP` or `EU-ADR` |
| `exposureName` / `outcomeName` | drug / outcome |
| `groundTruth` | 1 = known positive; 0 = negative control (the artefact-enrichment label) |
| `a_F,b_F,c_F,d_F` / `a_D,b_D,c_D,d_D` | 2×2 contingency cells in FAERS (`_F`) and DAEN (`_D`) |
| `ror_F`, `ror_F_ci_lo/hi` / `ror_D`, … | ROR point estimate + 95% CI per database |
| `faers_signal` / `daen_signal` | signal rule met (lower 95% CI of ROR > 1 AND a ≥ 3) |
| `cross_db_class` | concordance class: `concordant_positive`, `faers_only`, `daen_only`, `concordant_negative` |
| `n_drug_*`, `n_event_*`, `N_*` | drug, event, and total report margins per database |
| `daen_mde` | minimum detectable OR in DAEN at 80% power (exact noncentral-hypergeometric search) — the pre-registered power anchor |
| `daen_mde_07` / `daen_mde_09` | same at 70% / 90% power (red-team sensitivity) |
| `daen_mde_unconditional` | MDE without conditioning on the FAERS-observed margins |
| `daen_powered` | **primary** power flag: `daen_mde ≤ 1.5` |
| `daen_powered_at_2` / `_at_3` | permissive anchors (`≤ 2.0` / `≤ 3.0`) |
| `daen_powered_FAERS_pt` / `_FAERS_LB` | the **FAERS-derived-power** definitions whose circularity the paper demonstrates |
| `faers_power_at_15`, `faers_powered` | analogous FAERS-side power diagnostics |

---

## Simulation output table (`results/sim_results.parquet`)

Produced by `src/simulation/sweep.py` (driven by `run_local_fill.sh` / `run_cells.py`).
One row per fitted estimator within each replicate of each grid cell: 126 cells
× 2,000 replicates × 6 selection rules × 4 estimators = 6,048,000 rows, 29 columns.
The file reads only `sim_config.yaml`; no external data is required.

| Column | Type | Description |
|---|---|---|
| `regime` | str | data-generating regime: `Q1_null` (type-I) or `Q2_alt` (power) |
| `pi` | float | prevalence of reference-negative status among generated pairs |
| `lambda_inflate` | int | FAERS-only inflation factor applied under the alternative |
| `phi` | float | dispersion / over-reporting parameter of the report-count model |
| `rho` | float | within-cluster correlation target on the latent scale |
| `K` | int | number of clusters (drug families) in the stratum |
| `icc` | float | realised intracluster correlation |
| `OR_true` | str | true enrichment odds ratio for the cell (`null` under Q1) |
| `N_F` | int | FAERS background total report count (fixed at 20,005,192) |
| `n_pairs` | int | target number of drug–event pairs generated per replicate |
| `q2_daen_suppression` | float | DAEN suppression multiplier under the alternative |
| `or_lognormal_mu` / `or_lognormal_sigma` | float | log-normal effect-size distribution parameters |
| `label` | str | human-readable cell label |
| `cell_id` | str | grid-cell identifier (keys the per-cell seed stream) |
| `replicate` | int | replicate index within the cell |
| `selection_rule` | str | one of `mde_1.5`, `mde_2.0`, `mde_3.0`, `faers_point`, `faers_lb`, `none` |
| `estimator` | str | one of `fixed_effect`, `cluster_robust`, `wild_cluster_bootstrap`, `firth` |
| `reject` | bool | 1 if 95% CI lower bound of the enrichment OR > 1 (non-estimable scored 0) |
| `non_estimable` | bool | 1 if the fit failed (e.g. complete separation) |
| `or_hat` | float | estimated enrichment odds ratio |
| `ci_lo` / `ci_hi` | float | 95% confidence-interval bounds for the enrichment OR |
| `beta` / `se` | float | log-OR coefficient on `is_faers_only` and its standard error |
| `stratum_size` | int | number of pairs surviving the selection rule in this replicate |
| `n_faers_only` | int | count of FAERS-only (discordant) pairs in the stratum |
| `n_concordant_positive` | int | count of concordant-positive pairs in the stratum |
| `n_clusters` | int | number of clusters represented after selection |
