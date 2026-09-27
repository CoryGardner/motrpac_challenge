# Phase 4 — metabolite transfer, deep_mw

Built by `scripts/multiomic/04_metab_transfer.py` on 2026-09-27 07:35 UTC; numbers from the CSVs here.

- Source `deep`: 451 vials, 54 animals, 9 tissues, 340 RefMet-named metabolites present in every tissue (`../source_deep_tissues.csv`).
- Target `mw`: 840 samples, 70 individuals, 12 tissues, 190 RefMet keys; matched 44 (`feature_overlap.csv`).
- Tissue map: Plasma → PLASMA; Brain → CORTEX/HIPPOC/HYPOTH; Heart → HEART; Kidney → KIDNEY; Liver → LIVER; Lung → LUNG; Muscle (Quad) → SKM-GN/SKM-VL; Spleen → SPLEEN; other target tissues are OOD.
- Values: MoTrPAC sample-centred log2 (R package); target log2(peak area + 1) centred on the sample median; z-scored per metabolite within each dataset by the shared code path.

Matched RefMet names: **44** (`matched_features.csv`); source classes 9; mapped target samples 490 from 70 individuals.

## Accuracy — `accuracy_overall.csv`, `accuracy_by_tissue.csv`, `accuracy_by_stage.csv`

| model | accuracy_sample_weighted | accuracy_macro_over_target_tissues | n_samples_mapped | n_individuals_mapped | n_target_tissues_mapped | n_source_classes | chance | n_matched_features |
|---|---|---|---|---|---|---|---|---|
| full | 0.753 | 0.753 | 490 | 70 | 7 | 9 | 0.111 | 44 |
| k10 | 0.496 | 0.496 | 490 | 70 | 7 | 9 | 0.111 | 44 |
| k20 | 0.637 | 0.637 | 490 | 70 | 7 | 9 | 0.111 | 44 |

| target_tissue | rat_classes | n | n_individuals | full | k10 | k20 |
|---|---|---|---|---|---|---|
| Brain | HIPPOC | 70 | 70 | 1.00 | 1.00 | 1.00 |
| Heart | HEART | 70 | 70 | 0.80 | 0.36 | 0.96 |
| Kidney | KIDNEY | 70 | 70 | 0.96 | 0.66 | 0.06 |
| Liver | LIVER | 70 | 70 | 1.00 | 0.99 | 1.00 |
| Lung | LUNG | 70 | 70 | 0.80 | 0.13 | 0.51 |
| Muscle (Quad) | SKM-GN | 70 | 70 | 0.33 | 0.34 | 0.27 |
| Plasma | PLASMA | 70 | 70 | 0.39 | 0.00 | 0.66 |

Top prediction and main wrong call at k20:

| target_tissue | rat_classes | n | top_prediction | top_prediction_frac | main_wrong_call | main_wrong_call_frac |
|---|---|---|---|---|---|---|
| Bladder | OOD | 70 | KIDNEY | 0.63 | KIDNEY | 0.63 |
| Brain | HIPPOC | 70 | HIPPOC | 1.00 |  | 0.00 |
| Heart | HEART | 70 | HEART | 0.96 | SKM-GN | 0.04 |
| Kidney | KIDNEY | 70 | LIVER | 0.90 | LIVER | 0.90 |
| Liver | LIVER | 70 | LIVER | 1.00 |  | 0.00 |
| Lung | LUNG | 70 | LUNG | 0.51 | KIDNEY | 0.19 |
| Muscle (Quad) | SKM-GN | 70 | BAT | 0.40 | BAT | 0.40 |
| Pancreas | OOD | 70 | LUNG | 0.73 | LUNG | 0.73 |
| Plasma | PLASMA | 70 | PLASMA | 0.66 | SKM-GN | 0.17 |
| Spleen | OOD | 70 | BAT | 0.41 | BAT | 0.41 |
| Thymus | OOD | 70 | BAT | 0.96 | BAT | 0.96 |
| Tongue | OOD | 70 | KIDNEY | 0.34 | KIDNEY | 0.34 |

By stage (mapped samples):

| stage | full | k10 | k20 |
|---|---|---|---|
| 1 | 0.71 | 0.55 | 0.59 |
| 18 | 0.72 | 0.51 | 0.63 |
| 21 | 0.78 | 0.45 | 0.64 |
| 24 | 0.74 | 0.44 | 0.66 |
| 3 | 0.81 | 0.53 | 0.65 |

## Conformal transfer (α = 0.1, LAC, calibration on 16 held-out MoTrPAC animals) — `conformal_transfer.csv`, `coverage_by_tissue.csv`, `ood_sets.csv`, `recalibration.csv`

| stage | model | conformal | n_mapped | coverage_mapped | frac_empty_mapped | avg_set_size_mapped | ood_frac_empty | ood_avg_set_size | n_cal_animals |
|---|---|---|---|---|---|---|---|---|---|
| 1 | k10 | marginal | 98 | 0.143 | 0.786 | 0.214 | 1.000 | 0.000 | 16 |
| 1 | k10 | mondrian | 98 | 0.071 | 0.806 | 0.194 | 1.000 | 0.000 | 16 |
| 1 | k10 | floored | 98 | 0.163 | 0.714 | 0.286 | 1.000 | 0.000 | 16 |
| 18 | k10 | marginal | 98 | 0.163 | 0.816 | 0.184 | 1.000 | 0.000 | 16 |
| 18 | k10 | mondrian | 98 | 0.143 | 0.765 | 0.235 | 1.000 | 0.000 | 16 |
| 18 | k10 | floored | 98 | 0.163 | 0.745 | 0.255 | 1.000 | 0.000 | 16 |
| 21 | k10 | marginal | 98 | 0.143 | 0.816 | 0.184 | 1.000 | 0.000 | 16 |
| 21 | k10 | mondrian | 98 | 0.112 | 0.755 | 0.245 | 1.000 | 0.000 | 16 |
| 21 | k10 | floored | 98 | 0.153 | 0.714 | 0.286 | 1.000 | 0.000 | 16 |
| 24 | k10 | marginal | 98 | 0.173 | 0.796 | 0.204 | 1.000 | 0.000 | 16 |
| 24 | k10 | mondrian | 98 | 0.194 | 0.684 | 0.316 | 1.000 | 0.000 | 16 |
| 24 | k10 | floored | 98 | 0.204 | 0.673 | 0.327 | 1.000 | 0.000 | 16 |
| 3 | k10 | marginal | 98 | 0.153 | 0.827 | 0.173 | 0.986 | 0.014 | 16 |
| 3 | k10 | mondrian | 98 | 0.092 | 0.776 | 0.224 | 1.000 | 0.000 | 16 |
| 3 | k10 | floored | 98 | 0.153 | 0.714 | 0.286 | 0.986 | 0.014 | 16 |
| 1 | k20 | marginal | 98 | 0.122 | 0.878 | 0.122 | 1.000 | 0.000 | 16 |
| 1 | k20 | mondrian | 98 | 0.133 | 0.867 | 0.133 | 1.000 | 0.000 | 16 |
| 1 | k20 | floored | 98 | 0.143 | 0.857 | 0.143 | 1.000 | 0.000 | 16 |
| 18 | k20 | marginal | 98 | 0.214 | 0.786 | 0.214 | 1.000 | 0.000 | 16 |
| 18 | k20 | mondrian | 98 | 0.153 | 0.847 | 0.153 | 1.000 | 0.000 | 16 |
| 18 | k20 | floored | 98 | 0.214 | 0.786 | 0.214 | 1.000 | 0.000 | 16 |
| 21 | k20 | marginal | 98 | 0.214 | 0.786 | 0.214 | 1.000 | 0.000 | 16 |
| 21 | k20 | mondrian | 98 | 0.143 | 0.857 | 0.143 | 1.000 | 0.000 | 16 |
| 21 | k20 | floored | 98 | 0.235 | 0.765 | 0.235 | 1.000 | 0.000 | 16 |
| 24 | k20 | marginal | 98 | 0.173 | 0.827 | 0.173 | 1.000 | 0.000 | 16 |
| 24 | k20 | mondrian | 98 | 0.143 | 0.857 | 0.143 | 1.000 | 0.000 | 16 |
| 24 | k20 | floored | 98 | 0.194 | 0.806 | 0.194 | 1.000 | 0.000 | 16 |
| 3 | k20 | marginal | 98 | 0.153 | 0.847 | 0.153 | 1.000 | 0.000 | 16 |
| 3 | k20 | mondrian | 98 | 0.143 | 0.857 | 0.143 | 1.000 | 0.000 | 16 |
| 3 | k20 | floored | 98 | 0.163 | 0.837 | 0.163 | 1.000 | 0.000 | 16 |
| 1 | full | marginal | 98 | 0.000 | 1.000 | 0.000 | 1.000 | 0.000 | 16 |
| 1 | full | mondrian | 98 | 0.000 | 1.000 | 0.000 | 1.000 | 0.000 | 16 |
| 1 | full | floored | 98 | 0.000 | 1.000 | 0.000 | 1.000 | 0.000 | 16 |
| 18 | full | marginal | 98 | 0.000 | 1.000 | 0.000 | 1.000 | 0.000 | 16 |
| 18 | full | mondrian | 98 | 0.000 | 1.000 | 0.000 | 1.000 | 0.000 | 16 |
| 18 | full | floored | 98 | 0.000 | 1.000 | 0.000 | 1.000 | 0.000 | 16 |
| 21 | full | marginal | 98 | 0.000 | 1.000 | 0.000 | 1.000 | 0.000 | 16 |
| 21 | full | mondrian | 98 | 0.000 | 1.000 | 0.000 | 1.000 | 0.000 | 16 |
| 21 | full | floored | 98 | 0.000 | 1.000 | 0.000 | 1.000 | 0.000 | 16 |
| 24 | full | marginal | 98 | 0.000 | 1.000 | 0.000 | 1.000 | 0.000 | 16 |
| 24 | full | mondrian | 98 | 0.000 | 1.000 | 0.000 | 1.000 | 0.000 | 16 |
| 24 | full | floored | 98 | 0.000 | 1.000 | 0.000 | 1.000 | 0.000 | 16 |
| 3 | full | marginal | 98 | 0.000 | 1.000 | 0.000 | 1.000 | 0.000 | 16 |
| 3 | full | mondrian | 98 | 0.000 | 1.000 | 0.000 | 1.000 | 0.000 | 16 |
| 3 | full | floored | 98 | 0.000 | 1.000 | 0.000 | 1.000 | 0.000 | 16 |

| model | n_recal | n_test_individuals | draws | coverage_recalibrated | coverage_source_cal_same_test | frac_empty_recalibrated | frac_empty_source_cal | set_size_recalibrated |
|---|---|---|---|---|---|---|---|---|
| k10 | 3 | 11 | 20 | 0.931 | 0.155 | 0.000 | 0.823 | 4.701 |
| k10 | 5 | 9 | 20 | 0.921 | 0.155 | 0.000 | 0.824 | 4.629 |
| k20 | 3 | 11 | 20 | 0.924 | 0.154 | 0.001 | 0.846 | 2.344 |
| k20 | 5 | 9 | 20 | 0.920 | 0.146 | 0.000 | 0.854 | 2.233 |
| full | 3 | 11 | 20 | 0.912 | 0.000 | 0.011 | 1.000 | 1.549 |
| full | 5 | 9 | 20 | 0.928 | 0.000 | 0.006 | 1.000 | 1.688 |

Per tissue, k20 (marginal / Mondrian / floored, MoTrPAC calibration):

| target_tissue | floored | marginal | mondrian |
|---|---|---|---|
| Bladder |  |  |  |
| Brain | 0.34 | 0.34 | 0.01 |
| Heart | 0.00 | 0.00 | 0.00 |
| Kidney | 0.00 | 0.00 | 0.00 |
| Liver | 0.99 | 0.89 | 0.99 |
| Lung | 0.00 | 0.00 | 0.00 |
| Muscle (Quad) | 0.00 | 0.00 | 0.00 |
| Pancreas |  |  |  |
| Plasma | 0.00 | 0.00 | 0.00 |
| Spleen |  |  |  |
| Thymus |  |  |  |
| Tongue |  |  |  |

OOD tissues at k20: `ood_sets.csv`.

Panel at k20: `panel_k20.csv`.
