# Phase 6.2 — Detrended-H5c sensitivity

Generated: 2026-05-31 12:20:31

## Setup

FAERS-wide quarterly volume baseline computed from 20,005,192 surviving cases across 88 quarters. Volume growth from min-quarter to max-quarter: ~10×.

Per-pair Poisson regression of pair quarterly counts on log-volume baseline; Pearson residuals fed to BIC-selected single-break test.

H5c stratum (pairs with >= 16 quarters of data): 4,190 (3246 faers_only + 944 concordant_positive)

## Comparison: naive H5c vs detrended H5c

| Quantity | Naive H5c (raw counts) | Detrended H5c (residuals) |
|---|---|---|
| Pairs with detected break | 4,190 / 4,190 (100.0%) | 4,162 / 4,190 (99.3%) |
| Proportion faers_only with alert-aligned break | 0.0000 | 0.0000 |
| Proportion concordant_positive with alert-aligned break | 0.0000 | 0.0000 |
| Δ (faers_only − concordant_positive) | +0.00 pp | +0.00 pp |

| Two-sided p-value | nan | nan |

## Interpretation

**Detrending removed 28 of 4,190 naive-detected breaks** (0.7%). These were breaks driven by FAERS-wide volume growth, not by pair-specific reporting spikes.

**Detrended H5c is still null** (0/4190 pairs aligned an alert). The H5c null is robust to the volume-confound concern — neither the naive nor the detrended test finds the post-alert temporal-cluster signature on this 4,300-pair universe. Possible reasons: (i) the alert registry (45 alerts × 37 drugs) covers a small fraction of the FAERS-positive daen_powered universe (225 / 4,300 ≈ 5%); (ii) post-alert spikes may be transient (a few quarters) and not detected by a single-break model; (iii) the daen_powered restriction selects pairs with substantial DAEN volume, which preferentially includes well-established drug-event associations rather than mass-tort-litigated drugs.

## Decision

The detrended-H5c sensitivity is reported alongside the naive H5c result in the Methods, with the difference (or lack thereof) characterised. Per Section 16.3, this is a post-hoc sensitivity arm and flagged as such in the manuscript.
