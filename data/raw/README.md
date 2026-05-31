# Raw data — acquisition instructions

The raw spontaneous-report data underlying this study are **not redistributed**
in this repository. Both source databases are publicly available at no cost;
follow the steps below to assemble the two inputs the pipeline expects. Place
them as described, or point the pipeline at an existing FAERS substrate via the
`FAERS_SUBSTRATE` environment variable.

The small public **reference** datasets (ground-truth sets, crosswalks, alert and
mass-tort registries) **are** included, under `data/reference/`. See
`../../data_dictionary.md` for their provenance and schema.

---

## 1. FAERS substrate

**Source:** FDA Adverse Event Reporting System (FAERS) Quarterly Data Extract
files (ASCII), public domain.
Download: <https://fis.fda.gov/extensions/FPD-QDE-FAERS/FPD-QDE-FAERS.html>

The pipeline does **not** read the raw quarterly ASCII files directly. It reads a
**deduplicated parquet substrate** with the layout below (one quarter per parquet
in each subdirectory). Build this once from the public quarterly extracts:
deduplicate by retaining the latest `CASEVERSION` per `CASEID`, then write the
surviving case set and the long drug/reaction/demographic tables.

Expected layout (default location: `data/raw/faers_substrate/`, or set
`FAERS_SUBSTRATE=/path/to/substrate`):

```
faers_substrate/
├── surviving_cases.parquet      # column: primaryid (int64) — deduplicated case set
├── drug/<quarter>.parquet       # columns: primaryid, drugname, prod_ai
├── reac/<quarter>.parquet       # columns: primaryid, pt          (MedDRA Preferred Term)
└── demo/<quarter>.parquet       # columns: primaryid, quarter, occp_cod  (reporter type)
```

Notes:
- `prod_ai` is the FDA active-ingredient string; where empty the pipeline falls
  back to the lowercased `drugname`.
- `pt` is the MedDRA Preferred Term verbatim string. The licensed MedDRA
  hierarchy (PT→HLT) is **not** required and is **not** distributed here.
- `occp_cod` reporter-type codes follow the FDA DEMO dictionary
  (CN = consumer, LW = lawyer, MD = physician, PH = pharmacist, OT = other).

## 2. DAEN

**Source:** Therapeutic Goods Administration (TGA) Database of Adverse Event
Notifications (DAEN) — Medicines, public.
Portal: <https://www.tga.gov.au/safety/safety/safety-monitoring-daen-database-adverse-event-notifications/about-database-adverse-event-notifications-daen>

Export the line-level case, drug, and reaction tables and place them at
`data/daen/` (relative to the repository root) as:

```
data/daen/
├── daen_cases.csv            # column: case_number
├── daen_case_drugs.csv       # case_number + active-ingredient column
└── daen_case_reactions.csv   # case_number + MedDRA PT reaction column
```

Apply the same case-level deduplication used for FAERS (one row per case per
drug-string, one row per case per PT). Active-ingredient strings are matched to
FAERS after bullet-stripping and lowercasing; the US↔AU spelling crosswalk in
`data/reference/drug_name_crosswalk.csv` reconciles spelling differences
(e.g. acetaminophen ↔ paracetamol).
