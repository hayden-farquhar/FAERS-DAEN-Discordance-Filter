# Cross-database discordance as a drug-safety artefact filter (FAERS–DAEN)

Code repository for:

> **Cross-database discordance as a drug-safety artefact filter: a pre-registered
> FAERS–DAEN test, and the power-conditioning circularity that would have
> manufactured its confirmation.**

Hayden Farquhar, MBBS MPHTM — Independent Researcher, Finley, NSW, Australia
ORCID: [0009-0002-6226-440X](https://orcid.org/0009-0002-6226-440X)

Pre-registration (OSF): <https://doi.org/10.17605/OSF.IO/ZQ97A>

## Overview

This study asks whether **non-replication of a disproportionality signal across
two independent national spontaneous-reporting databases** is informative — i.e.
whether a drug–event signal that is disproportionate in the US FAERS but absent
in the Australian DAEN is enriched for known reporting artefacts rather than real
pharmacology. The entire argument hinges on **power-conditioning**: a signal can
fail to replicate in DAEN simply because DAEN is ~40× smaller, so the
non-replication rate is only ever quoted within drug–event pairs where DAEN has
enough reports to detect the FAERS effect if it were real.

The headline methodological result is a **cautionary counterfactual**: the
natural-default choice of a *FAERS-derived* power model induces a circularity that
manufactures an apparent confirmation (OR ≈ 5) of the artefact-filter hypothesis,
whereas the pre-registered, FAERS-independent **minimum-detectable-effect (MDE)**
power anchor shows the confirmatory enrichment test is non-estimable on this
substrate. This repository reproduces both routes and every sensitivity arm.

## Data sources

| Source | URL | Access | Redistributed here? |
|---|---|---|---|
| FDA FAERS Quarterly Data Extracts | <https://fis.fda.gov/extensions/FPD-QDE-FAERS/FPD-QDE-FAERS.html> | Free, public domain | No — see `data/raw/README.md` |
| TGA DAEN (Medicines) | <https://www.tga.gov.au/safety/safety/safety-monitoring-daen-database-adverse-event-notifications/about-database-adverse-event-notifications-daen> | Free, public | No — see `data/raw/README.md` |
| OMOP & EU-ADR reference sets | published control sets | Free | Yes — `data/reference/` |
| Harpaz et al. (2014) reference set | published (*Drug Safety*) | Free | Yes — `data/reference/` |

Raw line-level reports are **not** redistributed (the licensed MedDRA hierarchy is
also excluded). Both databases are free to obtain; `data/raw/README.md` gives the
exact expected schema. The small public reference/control datasets **are** included
under `data/reference/` (see `data_dictionary.md` for provenance).

## Requirements

Python ≥ 3.11 (results were produced on Python 3.14.3). Install dependencies into
a clean virtual environment:

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

Figure 1 (the flow diagram) is rendered from Mermaid source via the Mermaid CLI
(`mmdc`); if `mmdc` is not on the PATH, the other figures still render and Figure 1
is skipped with a warning. Install via `npm install -g @mermaid-js/mermaid-cli`
(optional).

## Reproduction

First obtain the raw data as described in `data/raw/README.md`, then run the
scripts in numbered order from the repository root. Each script reads/writes the
intermediate `results/` directory (created on first run) and is resumable:

```bash
mkdir -p results

# If your FAERS substrate is not at data/raw/faers_substrate, point to it:
# export FAERS_SUBSTRATE=/path/to/faers_substrate

python scripts/01_build_arm1.py            # Arm-1 per-pair 2x2 (FAERS + DAEN) + signal rule
python scripts/02_power_augment.py         # add daen_mde + daen_powered power family
python scripts/03_phase4_confirmatory.py   # Family 1 (H2 non-replication) + Family 2 (H1)
python scripts/04_build_arm2.py            # Arm-2 discovery substrate (full FAERS space)
python scripts/05_phase5_mechanism.py      # Family 3 mechanism features + H5a–d tests
python scripts/06_h2_sensitivity.py        # H2 across power-tag specifications
python scripts/07_ic025_sensitivity.py     # alternative signal definition (IC025 > 0)
python scripts/08_granularity_substitute.py# event-granularity sensitivity
python scripts/09_h5c_detrended.py         # detrended H5c time-break sensitivity
python scripts/10_make_figures.py          # Figures 1–6 -> outputs/figures/
```

Verify the worked-example and power-model tests pass:

```bash
pytest            # uses pyproject pythonpath=src, testpaths=tests
```

Runtime is dominated by parquet I/O over the FAERS substrate; the full pipeline
runs in well under an hour on a laptop (CPU only — no GPU). Random seed is fixed
to 94 for the Monte-Carlo power validation and cluster bootstraps.

The committed reference results live in `outputs/`; compare your regenerated
`results/` against them.

## Script descriptions

| Script | Description | Key inputs | Key outputs |
|---|---|---|---|
| `01_build_arm1.py` | Per-pair 2×2 in FAERS & DAEN for the OMOP+EU-ADR universe; applies signal rule | FAERS substrate, DAEN, reference sets | `results/perpair_arm1.parquet` |
| `02_power_augment.py` | Adds `daen_mde` (exact-PMF MDE search) and the `daen_powered` flag family | `perpair_arm1.parquet` | augmented `perpair_arm1.parquet`, `phase2_summary.md` |
| `03_phase4_confirmatory.py` | Family 1 descriptive headline (H2) + Family 2 confirmatory (H1) | augmented `perpair_arm1.parquet` | `phase4_results.{json,md}`, `phase4_h4_family4.md` |
| `04_build_arm2.py` | Arm-2 discovery substrate across the full FAERS ingredient×PT space | FAERS substrate, DAEN | `results/perpair_arm2.parquet` |
| `05_phase5_mechanism.py` | Family 3 mechanism features (consumer/lawyer share, time-break, mass-tort) + Holm H5a–d | `perpair_arm2.parquet`, alert & mass-tort registries | `perpair_arm2_with_features.parquet`, `phase5_results.{json,md}` |
| `06_h2_sensitivity.py` | Power-conditioned non-replication rate across power-tag specifications | augmented `perpair_arm1.parquet` | `phase6_h2_sensitivity.md` |
| `07_ic025_sensitivity.py` | Re-runs Families 1–2 under the BCPNN IC025 signal definition | `perpair_arm1.parquet` | `phase6_ic025_sensitivity.{json,md}` |
| `08_granularity_substitute.py` | Event-granularity substitute sensitivity arm | `perpair_arm1.parquet` | `phase6_granularity_substitute.{json,md}` |
| `09_h5c_detrended.py` | Volume-detrended H5c time-break alignment sensitivity | FAERS substrate, `perpair_arm2_with_features.parquet`, alert registry | `phase6_h5c_detrended.md` |
| `10_make_figures.py` | Renders Figures 1–6 | results parquets + `.mmd` source | `outputs/figures/*.{png,svg}` |

Shared statistical-metrics code (ROR/PRR/IC/EBGM with confidence intervals and
signalling rules) lives in `src/metrics/` and is the single source of truth used
across the author's pharmacovigilance pipelines; it is vendored here so this
repository runs standalone. The MDE-anchored power model is in `src/power/`.

## Outputs

| File | Paper reference |
|---|---|
| `outputs/figures/figure_1_arm1_flow.{png,svg}` | Figure 1 (Arm-1 flow) |
| `outputs/figures/figure_2_counterfactual_centrepiece.{png,svg}` | Figure 2 (counterfactual centrepiece) |
| `outputs/figures/figure_3_h2_sensitivity_forest.{png,svg}` | Figure 3 (H2 sensitivity forest) |
| `outputs/figures/figure_4_family3_forest.{png,svg}` | Figure 4 (Family 3 mechanism forest) |
| `outputs/figures/figure_5_cross_db_mosaic.{png,svg}` | Figure 5 (cross-database mosaic) |
| `outputs/figures/figure_6_h4_balance.{png,svg}` | Figure 6 (H4 balance) |
| `outputs/tables/phase4_results.md` | Primary endpoint + Family 2 confirmatory results |
| `outputs/tables/counterfactual_analysis.md` | Counterfactual (FAERS-derived vs MDE-anchored power) |
| `outputs/tables/phase5_results.md` | Family 3 mechanism (H5a–d) results |
| `outputs/tables/phase6_*.md` | Sensitivity arms |

## Citation

If you use this code, please cite the accompanying manuscript and the
pre-registration:

```
Farquhar H. Cross-database discordance as a drug-safety artefact filter:
a pre-registered FAERS–DAEN test, and the power-conditioning circularity that
would have manufactured its confirmation. Pre-registration:
https://doi.org/10.17605/OSF.IO/ZQ97A
```

## License

Code (`src/`, `scripts/`, `tests/`): MIT. Data and documentation (`data/`,
`outputs/`, `*.md`): CC-BY 4.0. See `LICENSE`.
