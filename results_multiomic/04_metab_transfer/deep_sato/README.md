# Phase 4 — metabolite transfer, deep_sato

Built by `scripts/multiomic/04_metab_transfer.py` on 2026-09-27 07:41 UTC; numbers from the CSVs here.

- Source `deep`: 451 vials, 54 animals, 9 tissues, 340 RefMet-named metabolites present in every tissue (`../source_deep_tissues.csv`).
- Target `sato`: 191 samples, 24 individuals, 8 tissues, 237 RefMet keys; matched 83 (`feature_overlap.csv`).
- Tissue map: BAT → BAT; eWAT → WAT-SC; iWAT → WAT-SC; HEART → HEART; HYPOTHALAMUS → HYPOTH; LIVER → LIVER; MUSCLE → SKM-GN/SKM-VL; SERUM → PLASMA; other target tissues are OOD.
- Values: MoTrPAC sample-centred log2 (R package); target log2(peak area + 1) centred on the sample median; z-scored per metabolite within each dataset by the shared code path.

## Invariance test (native Sato panel, sedentary → exercised) — `invariance_native_sedentary_to_exercised.csv`

| model | n_fit_animals | n_cal_animals | n_test_exercised_samples | n_test_animals | accuracy_exercised | coverage_exercised | avg_set_size | frac_empty | rna_reference_accuracy | rna_reference_coverage |
|---|---|---|---|---|---|---|---|---|---|---|
| k10 | 8 | 4 | 95 | 12 | 0.979 | 0.926 | 0.968 | 0.053 | 0.961 | 0.903 |
| k20 | 8 | 4 | 95 | 12 | 0.947 | 0.905 | 0.916 | 0.084 | 0.961 | 0.903 |
| k50 | 8 | 4 | 95 | 12 | 0.968 | 0.968 | 1.011 | 0.000 | 0.961 | 0.903 |
| full | 8 | 4 | 95 | 12 | 1.000 | 0.926 | 0.926 | 0.074 | 0.961 | 0.903 |

Matched RefMet names: **83** (`matched_features.csv`); source classes 9; mapped target samples 167 from 24 individuals.

## Accuracy — `accuracy_overall.csv`, `accuracy_by_tissue.csv`, `accuracy_by_stage.csv`

| model | accuracy_sample_weighted | accuracy_macro_over_target_tissues | n_samples_mapped | n_individuals_mapped | n_target_tissues_mapped | n_source_classes | chance | n_matched_features | acc_ci95_low_animal_boot | acc_ci95_high_animal_boot | n_boot |
|---|---|---|---|---|---|---|---|---|---|---|---|
| full | 0.527 | 0.524 | 167 | 24 | 7 | 9 | 0.111 | 83 | 0.500 | 0.555 | 1000 |
| k10 | 0.192 | 0.190 | 167 | 24 | 7 | 9 | 0.111 | 83 | 0.162 | 0.220 | 1000 |
| k20 | 0.329 | 0.327 | 167 | 24 | 7 | 9 | 0.111 | 83 | 0.269 | 0.388 | 1000 |
| k50 | 0.533 | 0.530 | 167 | 24 | 7 | 9 | 0.111 | 83 | 0.494 | 0.566 | 1000 |

| target_tissue | rat_classes | n | n_individuals | full | k10 | k20 | k50 |
|---|---|---|---|---|---|---|---|
| BAT | BAT | 24 | 24 | 1.00 | 0.12 | 0.67 | 0.92 |
| HEART | HEART | 24 | 24 | 1.00 | 0.42 | 0.75 | 1.00 |
| LIVER | LIVER | 24 | 24 | 0.67 | 0.79 | 0.38 | 0.96 |
| MUSCLE | SKM-GN | 24 | 24 | 0.00 | 0.00 | 0.00 | 0.00 |
| SERUM | PLASMA | 24 | 24 | 1.00 | 0.00 | 0.50 | 0.83 |
| eWAT | WAT-SC | 24 | 24 | 0.00 | 0.00 | 0.00 | 0.00 |
| iWAT | WAT-SC | 23 | 23 | 0.00 | 0.00 | 0.00 | 0.00 |

Top prediction and main wrong call at k20:

| target_tissue | rat_classes | n | top_prediction | top_prediction_frac | main_wrong_call | main_wrong_call_frac |
|---|---|---|---|---|---|---|
| BAT | BAT | 24 | BAT | 0.67 | HEART | 0.29 |
| HEART | HEART | 24 | HEART | 0.75 | SKM-GN | 0.25 |
| HYPOTHALAMUS | OOD | 24 | HIPPOC | 1.00 | HIPPOC | 1.00 |
| LIVER | LIVER | 24 | WAT-SC | 0.62 | WAT-SC | 0.62 |
| MUSCLE | SKM-GN | 24 | BAT | 0.71 | BAT | 0.71 |
| SERUM | PLASMA | 24 | SKM-GN | 0.50 | SKM-GN | 0.50 |
| eWAT | WAT-SC | 24 | BAT | 1.00 | BAT | 1.00 |
| iWAT | WAT-SC | 23 | BAT | 1.00 | BAT | 1.00 |

By stage (mapped samples):

| stage | full | k10 | k20 | k50 |
|---|---|---|---|---|
| Exercise | 0.53 | 0.17 | 0.43 | 0.57 |
| Sedentary | 0.52 | 0.21 | 0.23 | 0.50 |

## Conformal transfer (α = 0.1, LAC, calibration on 16 held-out MoTrPAC animals) — `conformal_transfer.csv`, `coverage_by_tissue.csv`, `ood_sets.csv`, `recalibration.csv`

| stage | model | conformal | n_mapped | coverage_mapped | frac_empty_mapped | avg_set_size_mapped | ood_frac_empty | ood_avg_set_size | n_cal_animals |
|---|---|---|---|---|---|---|---|---|---|
| Exercise | k10 | marginal | 83 | 0.000 | 1.000 | 0.000 | 1.000 | 0.000 | 16 |
| Exercise | k10 | mondrian | 83 | 0.096 | 0.687 | 0.313 | 0.500 | 0.500 | 16 |
| Exercise | k10 | floored | 83 | 0.096 | 0.687 | 0.313 | 0.500 | 0.500 | 16 |
| Sedentary | k10 | marginal | 84 | 0.000 | 1.000 | 0.000 | 1.000 | 0.000 | 16 |
| Sedentary | k10 | mondrian | 84 | 0.107 | 0.690 | 0.310 | 0.167 | 0.833 | 16 |
| Sedentary | k10 | floored | 84 | 0.107 | 0.690 | 0.310 | 0.167 | 0.833 | 16 |
| Exercise | k20 | marginal | 83 | 0.000 | 1.000 | 0.000 | 1.000 | 0.000 | 16 |
| Exercise | k20 | mondrian | 83 | 0.000 | 1.000 | 0.000 | 1.000 | 0.000 | 16 |
| Exercise | k20 | floored | 83 | 0.000 | 1.000 | 0.000 | 1.000 | 0.000 | 16 |
| Sedentary | k20 | marginal | 84 | 0.000 | 1.000 | 0.000 | 1.000 | 0.000 | 16 |
| Sedentary | k20 | mondrian | 84 | 0.012 | 0.988 | 0.012 | 1.000 | 0.000 | 16 |
| Sedentary | k20 | floored | 84 | 0.012 | 0.988 | 0.012 | 1.000 | 0.000 | 16 |
| Exercise | k50 | marginal | 83 | 0.000 | 1.000 | 0.000 | 1.000 | 0.000 | 16 |
| Exercise | k50 | mondrian | 83 | 0.000 | 1.000 | 0.000 | 1.000 | 0.000 | 16 |
| Exercise | k50 | floored | 83 | 0.000 | 1.000 | 0.000 | 1.000 | 0.000 | 16 |
| Sedentary | k50 | marginal | 84 | 0.000 | 1.000 | 0.000 | 1.000 | 0.000 | 16 |
| Sedentary | k50 | mondrian | 84 | 0.000 | 1.000 | 0.000 | 1.000 | 0.000 | 16 |
| Sedentary | k50 | floored | 84 | 0.000 | 1.000 | 0.000 | 1.000 | 0.000 | 16 |
| Exercise | full | marginal | 83 | 0.000 | 1.000 | 0.000 | 1.000 | 0.000 | 16 |
| Exercise | full | mondrian | 83 | 0.000 | 1.000 | 0.000 | 1.000 | 0.000 | 16 |
| Exercise | full | floored | 83 | 0.000 | 1.000 | 0.000 | 1.000 | 0.000 | 16 |
| Sedentary | full | marginal | 84 | 0.000 | 1.000 | 0.000 | 1.000 | 0.000 | 16 |
| Sedentary | full | mondrian | 84 | 0.000 | 1.000 | 0.000 | 1.000 | 0.000 | 16 |
| Sedentary | full | floored | 84 | 0.000 | 1.000 | 0.000 | 1.000 | 0.000 | 16 |

| model | n_recal | n_test_individuals | draws | coverage_recalibrated | coverage_source_cal_same_test | frac_empty_recalibrated | frac_empty_source_cal | set_size_recalibrated |
|---|---|---|---|---|---|---|---|---|
| k10 | 3 | 9 | 20 | 0.905 | 0.000 | 0.000 | 1.000 | 3.617 |
| k10 | 5 | 7 | 20 | 0.918 | 0.000 | 0.000 | 1.000 | 3.703 |
| k20 | 3 | 9 | 20 | 0.910 | 0.000 | 0.000 | 1.000 | 5.352 |
| k20 | 5 | 7 | 20 | 0.910 | 0.000 | 0.000 | 1.000 | 5.277 |
| k50 | 3 | 9 | 20 | 0.913 | 0.000 | 0.000 | 1.000 | 5.071 |
| k50 | 5 | 7 | 20 | 0.926 | 0.000 | 0.000 | 1.000 | 5.149 |
| full | 3 | 9 | 20 | 0.916 | 0.000 | 0.000 | 1.000 | 5.788 |
| full | 5 | 7 | 20 | 0.909 | 0.000 | 0.000 | 1.000 | 5.744 |

Per tissue, k20 (marginal / Mondrian / floored, MoTrPAC calibration):

| target_tissue | floored | marginal | mondrian |
|---|---|---|---|
| BAT | 0.00 | 0.00 | 0.00 |
| HEART | 0.04 | 0.00 | 0.04 |
| HYPOTHALAMUS |  |  |  |
| LIVER | 0.00 | 0.00 | 0.00 |
| MUSCLE | 0.00 | 0.00 | 0.00 |
| SERUM | 0.00 | 0.00 | 0.00 |
| eWAT | 0.00 | 0.00 | 0.00 |
| iWAT | 0.00 | 0.00 | 0.00 |

OOD tissues at k20: `ood_sets.csv`.

Panel at k20: `panel_k20.csv`.
