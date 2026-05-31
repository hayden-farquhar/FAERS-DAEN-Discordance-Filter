# Phase 6 (partial) — H2 power-conditioned non-replication rate, sensitivity arms

Generated: 2026-05-30 17:14:12
Source: `results/perpair_arm1.parquet` (492 pairs)
Endpoint: H2 = proportion of FAERS-positive AND <power-tag> Arm-1 pairs classified `faers_only`.
Cluster-bootstrap CIs (B=2000, seed=94, clusters=events).

## Sensitivity table

| Power tag | Description | n (H2 universe) | H2 point | 95% CI | underpowered-discordant alongside |
|---|---|---:|---:|---|---:|
| `daen_powered` | PRIMARY (MDE ≤ 1.5) | 44 | 0.136 | 0.057 – 0.318 | 248 |
| `daen_powered_07` | MDE ≤ 1.5 at 70% power (permissive power threshold) | 58 | 0.138 | 0.084 – 0.276 | 234 |
| `daen_powered_09` | MDE ≤ 1.5 at 90% power (stringent power threshold) | 33 | 0.152 | 0.079 – 0.312 | 259 |
| `daen_powered_at_2` | MDE ≤ 2.0 (permissive anchor) | 112 | 0.214 | 0.156 – 0.383 | 180 |
| `daen_powered_at_3` | MDE ≤ 3.0 (very permissive anchor) | 171 | 0.298 | 0.212 – 0.512 | 121 |
| `daen_powered_FAERS_pt` | Power ≥ 0.80 vs FAERS point ROR (the contested fix) | 126 | 0.135 | 0.101 – 0.257 | 166 |
| `daen_powered_FAERS_LB` | Power ≥ 0.80 vs FAERS 95% LB | 116 | 0.112 | 0.087 – 0.189 | 176 |
| `daen_powered_unconditional` | Poisson-margin MDE ≤ 1.5 | 43 | 0.140 | 0.058 – 0.318 | 249 |
| `daen_powered_min` | primary AND Poisson (conservative composite) | 43 | 0.140 | 0.058 – 0.318 | 249 |

## Interpretation

**Primary (PRIMARY MDE ≤ 1.5):** the headline number.

**Power-threshold sensitivities (0.7 / 0.9):** test whether the headline is sensitive to the power cutoff.

**MDE-anchor sensitivities (at_2 / at_3):** test whether the headline moves as the implicit 'real-effect' threshold rises. **A meaningful gap here is the H2 construct caveat (Section 4) made concrete** — pairs with true effects between the primary 1.5 and the sensitivity anchor inflate apparent non-replication.

**FAERS-derived sensitivities (FAERS_pt / FAERS_LB):** the contested earlier-design approach. **The gap from the primary quantifies the FAERS-derived-power circularity** that the Round-2 must-fix addressed. If the FAERS-derived H2 is markedly higher than the primary, the circularity would have inflated the headline.

**Poisson-margin sensitivities (unconditional / min):** the fixed-margin caveat (Section 6.3).

