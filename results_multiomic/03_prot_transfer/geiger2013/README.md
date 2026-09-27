# Phase 3b — protein transfer to geiger2013 (one sample per tissue)

Built by `scripts/multiomic/03b_prot_transfer_atlases.py` on 2026-09-27 07:29 UTC; numbers from the CSVs here.

- Matching: gene symbol, upper-cased (rat ↔ mouse); 2151 genes after the ≥ 50% completeness filter (`gene_overlap.csv`).
- Target: 29 tissue samples; 9 map to CORTEX;HEART;KIDNEY;LIVER;LUNG;SKM-GN;WAT-SC; without target: none; 20 OOD samples.
- Caveats: pooled mice per tissue; diaphragm and limb muscle both → SKM-GN, kidney cortex and medulla both → KIDNEY; matched by gene symbol across species.

## Accuracy over mapped samples — `accuracy_overall.csv`, `predictions_by_sample.csv`

| model | accuracy | n_mapped | n_ood | chance_1_over_7 |
|---|---|---|---|---|
| full | 1.0 | 9 | 20 | 0.143 |
| k20 | 0.6666666666666666 | 9 | 20 | 0.143 |
| k50 | 0.8888888888888888 | 9 | 20 | 0.143 |

Per sample at k20:

| atlas_tissue | rat_class | prediction | correct | p_max |
|---|---|---|---|---|
| Unnamed: 1 | OOD | WAT-SC | nan | 0.26 |
| Adrenal gland | OOD | SKM-GN | nan | 0.32 |
| Brain cortex | CORTEX | LIVER | False | 0.47 |
| Brain medulla | OOD | CORTEX | nan | 0.89 |
| Brown fat | OOD | KIDNEY | nan | 0.54 |
| Cerebellum | OOD | CORTEX | nan | 0.91 |
| Colon | OOD | LIVER | nan | 0.34 |
| Diaphragm | SKM-GN | SKM-GN | True | 0.82 |
| Duodenum | OOD | WAT-SC | nan | 0.37 |
| Embryonic tissue | OOD | CORTEX | nan | 0.53 |
| Eye | OOD | WAT-SC | nan | 0.31 |
| Heart | HEART | SKM-GN | False | 0.40 |
| Ileum | OOD | HEART | nan | 0.38 |
| Jejunum | OOD | LIVER | nan | 0.34 |
| Kidney cortex | KIDNEY | KIDNEY | True | 0.93 |
| Kidney medulla | KIDNEY | KIDNEY | True | 0.93 |
| Liver | LIVER | LIVER | True | 0.53 |
| Lung | LUNG | LUNG | True | 0.98 |
| Midbrain | OOD | CORTEX | nan | 0.62 |
| Muscle | SKM-GN | KIDNEY | False | 0.36 |
| Olfactory bulb | OOD | CORTEX | nan | 0.88 |
| Ovary | OOD | LUNG | nan | 0.38 |
| Pancreas | OOD | WAT-SC | nan | 0.33 |
| Salivary gland | OOD | WAT-SC | nan | 0.34 |
| Spleeen | OOD | LUNG | nan | 0.84 |
| Stomach | OOD | WAT-SC | nan | 0.46 |
| Thymus | OOD | LIVER | nan | 0.39 |
| Uterus | OOD | LUNG | nan | 0.40 |
| White fat | WAT-SC | WAT-SC | True | 0.50 |

## Conformal sets with MoTrPAC calibration (α = 0.1) — `conformal_transfer.csv`, `recalibration.csv`, `ood_sets.csv`

| model | conformal | n_mapped | coverage_mapped | frac_empty_mapped | avg_set_size_mapped | ood_frac_empty | ood_avg_set_size |
|---|---|---|---|---|---|---|---|
| k20 | marginal | 9 | 0.111 | 0.889 | 0.111 | 1.000 | 0.000 |
| k20 | mondrian | 9 | 0.111 | 0.889 | 0.111 | 1.000 | 0.000 |
| k20 | floored | 9 | 0.111 | 0.889 | 0.111 | 1.000 | 0.000 |
| k50 | marginal | 9 | 0.111 | 0.889 | 0.111 | 1.000 | 0.000 |
| k50 | mondrian | 9 | 0.111 | 0.889 | 0.111 | 1.000 | 0.000 |
| k50 | floored | 9 | 0.111 | 0.889 | 0.111 | 1.000 | 0.000 |
| full | marginal | 9 | 0.000 | 1.000 | 0.000 | 1.000 | 0.000 |
| full | mondrian | 9 | 0.000 | 1.000 | 0.000 | 1.000 | 0.000 |
| full | floored | 9 | 0.000 | 1.000 | 0.000 | 1.000 | 0.000 |

| model | n_recal | n_test_individuals | draws | coverage_recalibrated | coverage_source_cal_same_test | frac_empty_recalibrated | frac_empty_source_cal | set_size_recalibrated |
|---|---|---|---|---|---|---|---|---|
| k20 | 3 | 26 | 14 | 1.000 | 0.119 | 0.000 | 0.881 | 7.000 |
| k20 | 5 | 24 | 19 | 1.000 | 0.129 | 0.000 | 0.871 | 7.000 |
| k50 | 3 | 26 | 12 | 1.000 | 0.112 | 0.000 | 0.888 | 7.000 |
| k50 | 5 | 24 | 16 | 1.000 | 0.128 | 0.000 | 0.872 | 7.000 |
| full | 3 | 26 | 16 | 1.000 | 0.000 | 0.000 | 1.000 | 7.000 |
| full | 5 | 24 | 19 | 1.000 | 0.000 | 0.000 | 1.000 | 7.000 |

## What this does not show

- One sample per tissue: no within-tissue spread, no donor-level statement; recalibration 'individuals' are single tissue samples.

_Run time 0.0 min._