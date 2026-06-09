# Power-conditioning selection bias in cross-database pharmacovigilance replication

Code repository for:

> **Power-conditioning selection bias in cross-database pharmacovigilance
> replication: a simulation study of the reference-negative enrichment test, with
> a pre-registered FAERS–DAEN case study.**

Hayden Farquhar, MBBS MPHTM. Independent Researcher, Finley, NSW, Australia.
ORCID: [0009-0002-6226-440X](https://orcid.org/0009-0002-6226-440X)

Pre-registration (OSF): <https://doi.org/10.17605/OSF.IO/ZQ97A> (simulation
hypotheses S1–S4 locked as Amendment 3, 2026-06-08, before any cell was run).

## Overview

A signal present in one spontaneous-reporting database but absent in a second,
independent one is sometimes read as more likely to be a reporting artefact.
Testing that reading requires deciding whether the smaller database *could* have
replicated the larger one's signal, which is a power calculation. The analytic
default of anchoring that power to the larger database's own effect estimate is a
special case of two established phenomena: selection (collider) bias, and the
observed-power fallacy or winner's curse. Neither has been characterised in
cross-database pharmacovigilance.

This repository is simulation-led. The primary contribution is a simulation of
the reference-negative enrichment test across a grid of 126 data-generating cells
(2,000 replicates each, six power-conditioning selection rules, four cluster-aware
estimators, 6,048,000 fitted rows), reporting type-I error under a null where
discordance is uninformative (Q1), power under an alternative where it is
informative (Q2), and confidence-interval coverage as a function of the
event-cluster count K (Q3). Four hypotheses were pre-registered and graded
honest-either-way:

- S2 supported: a minimum-detectable-effect (MDE) rule anchored to the smaller
  database's marginals held type-I within the pre-specified [0, 0.075] band across
  all 48 primary-slice cells (maximum 0.0155).
- S3 supported: MDE type-I rose monotonically as the anchor loosened from 1.5 to
  2.0 to 3.0 (mean 0.0008, 0.0030, 0.0065).
- S1 not supported: the FAERS-point-derived rule did not inflate type-I above the
  band (maximum 0.043). On thin power-conditioned strata the test is conservative
  and frequently non-estimable, not anti-conservative.
- S4 not supported: the wild cluster bootstrap was not nearer nominal coverage
  than cluster-robust standard errors at few clusters (0.928 versus 0.944 at K = 5).

A reference-negative-density sweep then separated the rules: across an eightfold
density increase the MDE anchor held in band (type-I 0.014) while the FAERS-derived
rule rose to 0.043 (a roughly threefold gap) and tracked toward the unconditioned
collider baseline (0.052), even as its non-estimability all but vanished (0.002).
The transferable result is structural: on thin cross-database strata no rule is
both valid and powered, because the pairs that carry detectable enrichment signal
are the same pairs that carry the collider bias. The MDE anchor, computed from the
smaller database's marginals alone, is the principled corner of that trade-off.

The empirical FAERS–DAEN study is demoted to an illustrative case study. On the
genuine 492-pair substrate the FAERS-point default returned a publication-grade
reference-negative enrichment (OR = 5.17, 95% CI 1.68–15.94, p = 0.004) where the
pre-registered MDE-anchored primary found the same test non-estimable (complete
separation; realised MDE = 76.5). The simulation locates this as a finite-sample
fragility of thin strata, not a generic rate inflation.

## Repository layout

Two pipelines. The simulation pipeline is the lead contribution and needs no
external data: its fitted marginals are baked into `sim_config.yaml`, so it runs
from a clean clone. The empirical case-study pipeline reproduces the demoted
FAERS–DAEN illustration and requires the (free, public) source databases.

```
sim_config.yaml          locked simulation grid, seed, and fitted marginals (SSOT)
run_local_fill.sh        compute the 126-cell grid with 6 parallel shard-workers
run_cells.py             optional: recompute specific grid cells by index
src/simulation/          DGP, selection rules, estimators, sweep driver
src/power/               MDE-anchored and FAERS-derived power models
src/metrics/             shared ROR/PRR/IC/EBGM + CIs (vendored; standalone)
src/{replication,enrichment,mechanism,reporting,harmonise}/   case-study modules
scripts/sim_*.py         simulation analysis, sensitivity sweeps, figures
scripts/01_..10_*.py     empirical case-study pipeline (illustration)
outputs/                 committed reference results (figures + tables)
tests/                   worked-example, power-model, and G-SIM validation tests
```

## Data sources

The simulation needs no external data. The empirical case study uses:

| Source | URL | Access | Redistributed here? |
|---|---|---|---|
| FDA FAERS Quarterly Data Extracts | <https://fis.fda.gov/extensions/FPD-QDE-FAERS/FPD-QDE-FAERS.html> | Free, public domain | No, see `data/raw/README.md` |
| TGA DAEN (Medicines) | <https://www.tga.gov.au/safety/safety/safety-monitoring-daen-database-adverse-event-notifications/about-database-adverse-event-notifications-daen> | Free, public | No, see `data/raw/README.md` |
| OMOP & EU-ADR reference sets | published control sets | Free | Yes, `data/reference/` |
| Harpaz et al. (2014) reference set | published (*Drug Safety*) | Free | Yes, `data/reference/` |

Raw line-level reports are not redistributed (the licensed MedDRA hierarchy is also
excluded). Both databases are free to obtain; `data/raw/README.md` gives the exact
expected schema. The small public reference/control datasets are included under
`data/reference/` (see `data_dictionary.md` for provenance).

## Requirements

Python >= 3.11 (results were produced on Python 3.14.3). Install dependencies into
a clean virtual environment:

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

Figure 5 (the case-study flow diagram) is rendered from Mermaid source via the
Mermaid CLI (`mmdc`); if `mmdc` is not on the PATH, the other figures still render.
Install via `npm install -g @mermaid-js/mermaid-cli` (optional). SciencePlots is an
optional styling dependency; the figures fall back to a curated style if it is not
installed.

## Reproduction

### Part A: the simulation (lead contribution, no external data)

Run from the repository root. The sweep checkpoints each cell to
`results/sim_shards/` and is resumable; workers skip any cell already on disk.

```bash
# 1. Compute the 126-cell grid (6 parallel workers, single-threaded BLAS).
#    On a clean clone this computes every cell; re-runs only fill gaps.
bash run_local_fill.sh

# 2. Combine the shards into the analysis frame + write the seed spawn map.
cd src && python -m simulation.sweep --combine-only && cd ..

# 3. Downstream analysis (Q1 type-I, Q2 power, Q3 coverage, S1-S4 verdicts).
python scripts/sim_analysis.py            # -> results/sim_analysis/
python scripts/sim_density_sweep.py       # post-hoc density sweep -> results/sim_density/
python scripts/sim_wcr_coverage.py        # post-hoc restricted bootstrap -> results/sim_wcr/
python scripts/sim_make_figures.py        # Figures 1-4 -> outputs/figures/sim_*
```

The full sweep is the heavy step (6,048,000 fitted rows); it runs in well under an
hour across six laptop cores (CPU only, no GPU). The seed is fixed via
`numpy.random.SeedSequence(94).spawn(126)` keyed by grid position, so a cell's seed
is independent of how the work is partitioned. The fitted marginals in
`sim_config.yaml` were produced once by `src/simulation/calibrate.py` from the
case-study Arm-1 frame; re-fitting is optional and is the only simulation step that
touches the empirical data.

### Part B: the empirical case study (illustration)

First obtain the raw data as described in `data/raw/README.md`, then run the
numbered scripts in order from the repository root. Each reads/writes the
intermediate `results/` directory and is resumable:

```bash
# If your FAERS substrate is not at data/raw/faers_substrate, point to it:
# export FAERS_SUBSTRATE=/path/to/faers_substrate

python scripts/01_build_arm1.py            # Arm-1 per-pair 2x2 (FAERS + DAEN) + signal rule
python scripts/02_power_augment.py         # add daen_mde + daen_powered power family
python scripts/03_phase4_confirmatory.py   # Family 1 (H2 non-replication) + Family 2 (H1)
python scripts/04_build_arm2.py            # Arm-2 discovery substrate (full FAERS space)
python scripts/05_phase5_mechanism.py      # Family 3 mechanism features + H5a-d tests
python scripts/06_h2_sensitivity.py        # H2 across power-tag specifications
python scripts/07_ic025_sensitivity.py     # alternative signal definition (IC025 > 0)
python scripts/08_granularity_substitute.py# event-granularity sensitivity
python scripts/09_h5c_detrended.py         # detrended H5c time-break sensitivity
python scripts/10_make_figures.py          # case-study figures -> outputs/figures/
```

Verify the worked-example, power-model, and simulation-validation tests pass:

```bash
pytest            # uses pyproject pythonpath=src, testpaths=tests
```

The committed reference results live in `outputs/`; compare your regenerated
`results/` against them.

## Script descriptions

### Simulation pipeline

| Script | Description | Key outputs |
|---|---|---|
| `sim_config.yaml` | Locked grid, replicate count, master seed, and fitted marginal parameters (single source of truth the sweep reads) | (config) |
| `run_local_fill.sh` | Runs the sweep as six parallel shard-workers over the 126-cell grid | `results/sim_shards/*.parquet` |
| `python -m simulation.sweep --combine-only` | Combines the shards and writes the seed spawn map | `results/sim_results.parquet`, `results/sim_seed_spawn_map.json` |
| `run_cells.py` | Recomputes specific grid cells by index (rebalances a lopsided tail; byte-identical to the worker) | `results/sim_shards/*.parquet` |
| `scripts/sim_analysis.py` | Q1 type-I, Q2 power, Q3 coverage, and the S1-S4 verdicts | `results/sim_analysis/` (incl. `sim_analysis_summary.md`, `s_hypotheses.json`) |
| `scripts/sim_density_sweep.py` | Post-hoc reference-negative-density sweep (1x to 8x), both regimes and three rules | `results/sim_density/density_summary.md` |
| `scripts/sim_wcr_coverage.py` | Post-hoc restricted (null-imposed) wild bootstrap on the Q3 coverage line | `results/sim_wcr/wcr_coverage_summary.md`, `s4_wcr_verdict.json` |
| `scripts/sim_make_figures.py` | Renders simulation Figures 1-4 | `outputs/figures/sim_*.{png,svg}` |

### Empirical case-study pipeline

| Script | Description | Key inputs | Key outputs |
|---|---|---|---|
| `01_build_arm1.py` | Per-pair 2×2 in FAERS & DAEN for the OMOP+EU-ADR universe; applies signal rule | FAERS substrate, DAEN, reference sets | `results/perpair_arm1.parquet` |
| `02_power_augment.py` | Adds `daen_mde` (exact-PMF MDE search) and the `daen_powered` flag family | `perpair_arm1.parquet` | augmented `perpair_arm1.parquet`, `phase2_summary.md` |
| `03_phase4_confirmatory.py` | Family 1 descriptive headline (H2) + Family 2 confirmatory (H1) | augmented `perpair_arm1.parquet` | `phase4_results.{json,md}`, `phase4_h4_family4.md` |
| `04_build_arm2.py` | Arm-2 discovery substrate across the full FAERS ingredient×PT space | FAERS substrate, DAEN | `results/perpair_arm2.parquet` |
| `05_phase5_mechanism.py` | Family 3 mechanism features (consumer/lawyer share, time-break, mass-tort) + Holm H5a-d | `perpair_arm2.parquet`, alert & mass-tort registries | `perpair_arm2_with_features.parquet`, `phase5_results.{json,md}` |
| `06_h2_sensitivity.py` | Power-conditioned non-replication rate across power-tag specifications | augmented `perpair_arm1.parquet` | `phase6_h2_sensitivity.md` |
| `07_ic025_sensitivity.py` | Re-runs Families 1-2 under the BCPNN IC025 signal definition | `perpair_arm1.parquet` | `phase6_ic025_sensitivity.{json,md}` |
| `08_granularity_substitute.py` | Event-granularity substitute sensitivity arm | `perpair_arm1.parquet` | `phase6_granularity_substitute.{json,md}` |
| `09_h5c_detrended.py` | Volume-detrended H5c time-break alignment sensitivity | FAERS substrate, `perpair_arm2_with_features.parquet`, alert registry | `phase6_h5c_detrended.md` |
| `10_make_figures.py` | Renders the case-study figures | results parquets + `.mmd` source | `outputs/figures/*.{png,svg}` |

Shared statistical-metrics code (ROR/PRR/IC/EBGM with confidence intervals and
signalling rules) lives in `src/metrics/` and is the single source of truth used
across the author's pharmacovigilance pipelines; it is vendored here so this
repository runs standalone. The MDE-anchored and FAERS-derived power models are in
`src/power/`; the simulation DGP, selection rules, and estimators are in
`src/simulation/`.

## Outputs

| File | Paper reference |
|---|---|
| `outputs/figures/sim_1_typeI_surface.{png,svg}` | Figure 1 (simulation type-I surface) |
| `outputs/figures/sim_2_power_curve.{png,svg}` | Figure 2 (simulation power curve) |
| `outputs/figures/sim_3_coverage_K.{png,svg}` | Figure 3 (few-clusters coverage) |
| `outputs/figures/sim_4_density_sweep.{png,svg}` | Figure 4 (reference-negative-density sweep) |
| `outputs/figures/figure_2_counterfactual_centrepiece.{png,svg}` | Figure 5 (case study: FAERS-derived vs MDE-anchored) |
| `outputs/figures/figure_5_cross_db_mosaic.{png,svg}` | Figure 6 (case study: cross-database 2×2 mosaic) |
| `outputs/figures/figure_6_h4_balance.{png,svg}` | Figure 7 (case study: DAEN-only balance) |
| `outputs/figures/figure_1_arm1_flow.{png,svg}` | Case-study Arm-1 flow diagram |
| `outputs/figures/figure_3_h2_sensitivity_forest.{png,svg}` | Case-study sensitivity detail (supplement) |
| `outputs/figures/figure_4_family3_forest.{png,svg}` | Case-study mechanism detail (supplement) |
| `outputs/tables/sim_analysis_summary.md` | Q1-Q3 endpoints and S1-S4 verdicts |
| `outputs/tables/s_hypotheses.json` | Machine-readable S1-S4 verdicts |
| `outputs/tables/density_summary.md` | Reference-negative-density sweep |
| `outputs/tables/wcr_coverage_summary.md`, `s4_wcr_verdict.json` | Restricted wild bootstrap coverage |
| `outputs/tables/phase4_results.md` | Case-study primary + Family 2 confirmatory |
| `outputs/tables/counterfactual_analysis.md` | Case-study FAERS-derived vs MDE-anchored power |
| `outputs/tables/phase5_results.md` | Case-study Family 3 mechanism (H5a-d) |
| `outputs/tables/phase6_*.md` | Case-study sensitivity arms |

## Citation

If you use this code, please cite the accompanying manuscript and the
pre-registration:

```
Farquhar H. Power-conditioning selection bias in cross-database pharmacovigilance
replication: a simulation study of the reference-negative enrichment test, with a
pre-registered FAERS-DAEN case study. Pre-registration:
https://doi.org/10.17605/OSF.IO/ZQ97A
```

## License

Code (`src/`, `scripts/`, `tests/`): MIT. Data and documentation (`data/`,
`outputs/`, `*.md`): CC-BY 4.0. See `LICENSE`.
