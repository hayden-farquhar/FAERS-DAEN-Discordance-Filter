# Phase 4 Results — Family 1 Headline + Family 2 Confirmatory
Generated: 2026-05-30 17:11:39
Input substrate: `results/perpair_arm1.parquet` (492 pairs)

---

## Family 1 — Primary descriptive headline (no formal test)

**H2 (power-conditioned non-replication rate).** Within the FAERS-positive AND `daen_powered` Arm-1 universe (n=44), the proportion classified `faers_only` is **0.136** (95% cluster-bootstrap CI: 0.057 – 0.318; B=2,000; seed=94; clusters=events).

**Underpowered-discordant pairs** (FAERS-positive but DAEN underpowered): 248. Reported alongside per Section 4.

**Concordance (on the `daen_powered` universe):**

| Statistic | Value | 95% CI |
|---|---|---|
| PA+ (positive-specific agreement) | 0.916 | 0.811 – 0.964 |
| PA+ (mutually-powered subset) | 0.916 | 0.811 – 0.966 |
| Cohen's κ | 0.588 | — |

**Cohen's κ caveat:** Cohen's kappa = 0.588. Marginal balance is reasonable (84.6% / 75.0%); kappa interpretable but PA+ (= 0.916) remains the pre-registered primary concordance statistic.

**Cross-DB 2×2 (on `daen_powered` universe):**

| | DAEN+ | DAEN− |
|---|---|---|
| FAERS+ | 38 | 6 |
| FAERS− | 1 | 7 |

**Substrate notes:**
- H2 universe (FAERS-positive AND daen_powered): n=44
- Underpowered-discordant (FAERS-positive AND NOT daen_powered): n=248
- daen_powered universe (for concordance + kappa): n=52
- Mutually-powered universe (daen_powered AND faers_powered): n=52

---

## Family 2 — Secondary confirmatory (Holm-Bonferroni, k=2, family-wise α=0.05)

### H1 — Reference-negative enrichment

**Stratum:** 44 pairs (FAERS-positive AND `daen_powered` AND in {faers_only, concordant_positive}). Cluster structure: 7 unique events; 28 unique drugs.

**Per-stratum cells:**

| | faers_only | concordant_positive |
|---|---|---|
| known-negative | 0 | 5 |
| known-positive | 6 | 33 |

**Estimates:**

| Specification | OR | 95% CI | p-value |
|---|---|---|---|
| Fixed-effect logistic | 0.000 | 0.000 – — | — |
| Cluster-robust SE (cluster=event) | — | — – — | — |

**Realised MDE (OR at 80% power, α=0.05):** 76.51

**Pre-registered decision:** `EXPLORATORY_DOWNGRADE` (fragility flag: `N/A`)

**Holm-Bonferroni within Family 2:** `N/A`

- Separation: known-neg counts in (faers_only, concordant_positive) = (0, 5) of (6, 38). Standard mixed-effects logistic OR is undefined or infinite.
- Logit fit failed: math range error
- Realised MDE = 76.51 > 3.0 ⇒ H1 downgraded to exploratory per Section 6.6; Family 2 pivot (Section 5) fires.

### H3 — Discordance-filter positive likelihood ratio

**Stratum:** 44 pairs.

**Operating characteristics:**

| Metric | Value |
|---|---|
| Sensitivity | 0.000 |
| Specificity | 0.846 |
| PPV | 0.000 |
| NPV | 0.868 |
| LR+ | — (95% CI: — – —) |
| One-sided p (LR+ > 1) | — |

**Realised LR+ MDE (one-sided 80% power):** ∞

**Pre-registered decision:** `DEGENERATE`

**Holm-Bonferroni within Family 2:** `N/A`

- LR+ degenerate: sensitivity=0.0, specificity=0.8461538461538461, tp/fp/fn/tn=0/6/5/33. One of the operating-characteristic boundaries was hit.

---

## Within-event matched robustness (Mantel-Haenszel, Section 10.5)

- Pairs in analytic universe: 44
- Total event strata: 7
- **Informative strata** (both rows AND both columns non-zero): **3**
- MH pooled OR: — (95% CI: — – —)

- MH OR undefined: no informative stratum had both rows AND both columns nonzero. This is the expected pattern under severe degeneracy in the H1 stratum.

---

## Pre-registered pivot status (Section 5)

**Family 2 pivot for H1 FIRES (MDE > 3.0):** the realised MDE exceeds the pre-registered downgrade threshold; H1 is reported as exploratory; the headline falls back to Family 1 descriptive (concordance + power-conditioned non-replication rate). Pivot framing per protocol Section 5: *Calibrated negative-result paper.*

