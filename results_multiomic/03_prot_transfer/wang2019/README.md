# Phase 3b — protein transfer to wang2019 (one sample per tissue)

Built by `scripts/multiomic/03b_prot_transfer_atlases.py` on 2026-09-27 07:23 UTC; numbers from the CSVs here.

- Matching: 1:1 rat–human orthologs by Ensembl id; 3037 genes after the ≥ 80% completeness filter (`gene_overlap.csv`).
- Target: 29 tissue samples; 6 map to CORTEX;HEART;KIDNEY;LIVER;LUNG;WAT-SC; without target: SKM-GN; 23 OOD samples.
- Caveats: whole brain → CORTEX; subcutaneous 'Fat' → WAT-SC; no skeletal muscle in the atlas (SKM-GN has no target); one adult donor per tissue, label-free intensities.

## Accuracy over mapped samples — `accuracy_overall.csv`, `predictions_by_sample.csv`

| model | accuracy | n_mapped | n_ood | chance_1_over_7 |
|---|---|---|---|---|
| full | 0.6666666666666666 | 6 | 23 | 0.143 |
| k20 | 0.5 | 6 | 23 | 0.143 |
| k50 | 0.6666666666666666 | 6 | 23 | 0.143 |

Per sample at k20:

| atlas_tissue | rat_class | prediction | correct | p_max |
|---|---|---|---|---|
| Adrenal gland | OOD | CORTEX | nan | 0.35 |
| Appendix | OOD | LUNG | nan | 0.44 |
| Brain | CORTEX | CORTEX | True | 0.99 |
| Colon | OOD | WAT-SC | nan | 0.59 |
| Duodenum | OOD | CORTEX | nan | 0.29 |
| Endometrium | OOD | WAT-SC | nan | 0.24 |
| Esophagus | OOD | LUNG | nan | 0.27 |
| Fallopian tube | OOD | SKM-GN | nan | 0.38 |
| Fat | WAT-SC | LUNG | False | 0.39 |
| Gallbladder | OOD | WAT-SC | nan | 0.65 |
| Heart | HEART | KIDNEY | False | 0.68 |
| Kidney | KIDNEY | KIDNEY | True | 0.93 |
| Liver | LIVER | LIVER | True | 0.55 |
| Lung | LUNG | SKM-GN | False | 0.48 |
| Lymph node | OOD | WAT-SC | nan | 0.42 |
| Ovary | OOD | WAT-SC | nan | 0.36 |
| Pancreas | OOD | KIDNEY | nan | 0.55 |
| Placenta | OOD | LUNG | nan | 0.54 |
| Prostate | OOD | WAT-SC | nan | 0.42 |
| Rectum | OOD | WAT-SC | nan | 0.53 |
| Salivary gland | OOD | KIDNEY | nan | 0.53 |
| Small intestine | OOD | KIDNEY | nan | 0.28 |
| Smooth muscle | OOD | WAT-SC | nan | 0.86 |
| Spleen | OOD | LUNG | nan | 0.62 |
| Stomach | OOD | WAT-SC | nan | 0.53 |
| Testis | OOD | SKM-GN | nan | 0.61 |
| Thyroid | OOD | SKM-GN | nan | 0.49 |
| Tonsil | OOD | LUNG | nan | 0.77 |
| Urinary bladder | OOD | LUNG | nan | 0.73 |

## Conformal sets with MoTrPAC calibration (α = 0.1) — `conformal_transfer.csv`, `recalibration.csv`, `ood_sets.csv`

| model | conformal | n_mapped | coverage_mapped | frac_empty_mapped | avg_set_size_mapped | ood_frac_empty | ood_avg_set_size |
|---|---|---|---|---|---|---|---|
| k20 | marginal | 6 | 0.333 | 0.667 | 0.333 | 1.000 | 0.000 |
| k20 | mondrian | 6 | 0.333 | 0.667 | 0.333 | 1.000 | 0.000 |
| k20 | floored | 6 | 0.333 | 0.667 | 0.333 | 1.000 | 0.000 |
| k50 | marginal | 6 | 0.333 | 0.667 | 0.333 | 1.000 | 0.000 |
| k50 | mondrian | 6 | 0.333 | 0.667 | 0.333 | 1.000 | 0.000 |
| k50 | floored | 6 | 0.333 | 0.667 | 0.333 | 1.000 | 0.000 |
| full | marginal | 6 | 0.167 | 0.833 | 0.167 | 1.000 | 0.000 |
| full | mondrian | 6 | 0.000 | 1.000 | 0.000 | 1.000 | 0.000 |
| full | floored | 6 | 0.167 | 0.833 | 0.167 | 1.000 | 0.000 |

| model | n_recal | n_test_individuals | draws | coverage_recalibrated | coverage_source_cal_same_test | frac_empty_recalibrated | frac_empty_source_cal | set_size_recalibrated |
|---|---|---|---|---|---|---|---|---|
| k20 | 3 | 26 | 9 | 1.000 | 0.344 | 0.000 | 0.656 | 7.000 |
| k20 | 5 | 24 | 16 | 1.000 | 0.350 | 0.000 | 0.650 | 7.000 |
| k50 | 3 | 26 | 9 | 1.000 | 0.317 | 0.000 | 0.683 | 7.000 |
| k50 | 5 | 24 | 16 | 1.000 | 0.355 | 0.000 | 0.645 | 7.000 |
| full | 3 | 26 | 11 | 1.000 | 0.127 | 0.000 | 0.873 | 7.000 |
| full | 5 | 24 | 11 | 1.000 | 0.155 | 0.000 | 0.845 | 7.000 |

## What this does not show

- One sample per tissue: no within-tissue spread, no donor-level statement; recalibration 'individuals' are single tissue samples.

_Run time 0.1 min._