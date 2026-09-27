# Phase 4 — metabolite transfer, hilic_sato

Built by `scripts/multiomic/04_metab_transfer.py` on 2026-09-27 07:41 UTC; numbers from the CSVs here.

- Source `hilic`: 876 vials, 54 animals, 19 tissues, 129 RefMet-named metabolites present in every tissue (`../source_hilic_tissues.csv`).
- Target `sato`: 191 samples, 24 individuals, 8 tissues, 237 RefMet keys; matched 58 (`feature_overlap.csv`).
- Tissue map: BAT → BAT; eWAT → WAT-SC; iWAT → WAT-SC; HEART → HEART; HYPOTHALAMUS → HYPOTH; LIVER → LIVER; MUSCLE → SKM-GN/SKM-VL; SERUM → PLASMA; other target tissues are OOD.
- Values: MoTrPAC sample-centred log2 (R package); target log2(peak area + 1) centred on the sample median; z-scored per metabolite within each dataset by the shared code path.

## Invariance test (native Sato panel, sedentary → exercised) — `invariance_native_sedentary_to_exercised.csv`

| model | n_fit_animals | n_cal_animals | n_test_exercised_samples | n_test_animals | accuracy_exercised | coverage_exercised | avg_set_size | frac_empty | rna_reference_accuracy | rna_reference_coverage |
|---|---|---|---|---|---|---|---|---|---|---|
| k10 | 8 | 4 | 95 | 12 | 0.979 | 0.926 | 0.968 | 0.053 | 0.961 | 0.903 |
| k20 | 8 | 4 | 95 | 12 | 0.947 | 0.905 | 0.916 | 0.084 | 0.961 | 0.903 |
| k50 | 8 | 4 | 95 | 12 | 0.968 | 0.968 | 1.011 | 0.000 | 0.961 | 0.903 |
| full | 8 | 4 | 95 | 12 | 1.000 | 0.926 | 0.926 | 0.074 | 0.961 | 0.903 |

Matched RefMet names: **58** (`matched_features.csv`); source classes 19; mapped target samples 191 from 24 individuals.

## Accuracy — `accuracy_overall.csv`, `accuracy_by_tissue.csv`, `accuracy_by_stage.csv`

| model | accuracy_sample_weighted | accuracy_macro_over_target_tissues | n_samples_mapped | n_individuals_mapped | n_target_tissues_mapped | n_source_classes | chance | n_matched_features | acc_ci95_low_animal_boot | acc_ci95_high_animal_boot | n_boot |
|---|---|---|---|---|---|---|---|---|---|---|---|
| full | 0.445 | 0.443 | 191 | 24 | 8 | 19 | 0.053 | 58 | 0.401 | 0.484 | 1000 |
| k10 | 0.215 | 0.214 | 191 | 24 | 8 | 19 | 0.053 | 58 | 0.188 | 0.241 | 1000 |
| k20 | 0.277 | 0.276 | 191 | 24 | 8 | 19 | 0.053 | 58 | 0.257 | 0.300 | 1000 |
| k50 | 0.445 | 0.443 | 191 | 24 | 8 | 19 | 0.053 | 58 | 0.419 | 0.471 | 1000 |

| target_tissue | rat_classes | n | n_individuals | full | k10 | k20 | k50 |
|---|---|---|---|---|---|---|---|
| BAT | BAT | 24 | 24 | 0.67 | 0.00 | 0.00 | 0.46 |
| HEART | HEART | 24 | 24 | 0.08 | 0.00 | 0.00 | 0.04 |
| HYPOTHALAMUS | HYPOTH | 24 | 24 | 1.00 | 0.96 | 1.00 | 1.00 |
| LIVER | LIVER | 24 | 24 | 1.00 | 0.00 | 1.00 | 1.00 |
| MUSCLE | SKM-GN/SKM-VL | 24 | 24 | 0.50 | 0.00 | 0.21 | 1.00 |
| SERUM | PLASMA | 24 | 24 | 0.29 | 0.71 | 0.00 | 0.04 |
| eWAT | WAT-SC | 24 | 24 | 0.00 | 0.00 | 0.00 | 0.00 |
| iWAT | WAT-SC | 23 | 23 | 0.00 | 0.04 | 0.00 | 0.00 |

Top prediction and main wrong call at k20:

| target_tissue | rat_classes | n | top_prediction | top_prediction_frac | main_wrong_call | main_wrong_call_frac |
|---|---|---|---|---|---|---|
| BAT | BAT | 24 | KIDNEY | 1.00 | KIDNEY | 1.00 |
| HEART | HEART | 24 | VENACV | 0.46 | VENACV | 0.46 |
| HYPOTHALAMUS | HYPOTH | 24 | HYPOTH | 1.00 |  | 0.00 |
| LIVER | LIVER | 24 | LIVER | 1.00 |  | 0.00 |
| MUSCLE | SKM-GN/SKM-VL | 24 | OVARY | 0.38 | OVARY | 0.38 |
| SERUM | PLASMA | 24 | WAT-SC | 0.67 | WAT-SC | 0.67 |
| eWAT | WAT-SC | 24 | BAT | 0.67 | BAT | 0.67 |
| iWAT | WAT-SC | 23 | BAT | 0.74 | BAT | 0.74 |

By stage (mapped samples):

| stage | full | k10 | k20 | k50 |
|---|---|---|---|---|
| Exercise | 0.42 | 0.23 | 0.27 | 0.48 |
| Sedentary | 0.47 | 0.20 | 0.28 | 0.41 |

## Conformal transfer (α = 0.1, LAC, calibration on 16 held-out MoTrPAC animals) — `conformal_transfer.csv`, `coverage_by_tissue.csv`, `ood_sets.csv`, `recalibration.csv`

| stage | model | conformal | n_mapped | coverage_mapped | frac_empty_mapped | avg_set_size_mapped | ood_frac_empty | ood_avg_set_size | n_cal_animals |
|---|---|---|---|---|---|---|---|---|---|
| Exercise | k10 | marginal | 95 | 0.084 | 0.800 | 0.200 |  |  | 16 |
| Exercise | k10 | mondrian | 95 | 0.000 | 0.000 | 1.579 |  |  | 16 |
| Exercise | k10 | floored | 95 | 0.084 | 0.000 | 1.716 |  |  | 16 |
| Sedentary | k10 | marginal | 96 | 0.052 | 0.740 | 0.260 |  |  | 16 |
| Sedentary | k10 | mondrian | 96 | 0.000 | 0.000 | 1.562 |  |  | 16 |
| Sedentary | k10 | floored | 96 | 0.052 | 0.000 | 1.688 |  |  | 16 |
| Exercise | k20 | marginal | 95 | 0.053 | 0.947 | 0.053 |  |  | 16 |
| Exercise | k20 | mondrian | 95 | 0.126 | 0.000 | 1.516 |  |  | 16 |
| Exercise | k20 | floored | 95 | 0.126 | 0.000 | 1.516 |  |  | 16 |
| Sedentary | k20 | marginal | 96 | 0.052 | 0.948 | 0.052 |  |  | 16 |
| Sedentary | k20 | mondrian | 96 | 0.052 | 0.000 | 1.417 |  |  | 16 |
| Sedentary | k20 | floored | 96 | 0.052 | 0.000 | 1.417 |  |  | 16 |
| Exercise | k50 | marginal | 95 | 0.000 | 1.000 | 0.000 |  |  | 16 |
| Exercise | k50 | mondrian | 95 | 0.095 | 0.000 | 1.453 |  |  | 16 |
| Exercise | k50 | floored | 95 | 0.095 | 0.000 | 1.453 |  |  | 16 |
| Sedentary | k50 | marginal | 96 | 0.000 | 1.000 | 0.000 |  |  | 16 |
| Sedentary | k50 | mondrian | 96 | 0.125 | 0.000 | 1.542 |  |  | 16 |
| Sedentary | k50 | floored | 96 | 0.125 | 0.000 | 1.542 |  |  | 16 |
| Exercise | full | marginal | 95 | 0.000 | 1.000 | 0.000 |  |  | 16 |
| Exercise | full | mondrian | 95 | 0.126 | 0.000 | 1.505 |  |  | 16 |
| Exercise | full | floored | 95 | 0.126 | 0.000 | 1.505 |  |  | 16 |
| Sedentary | full | marginal | 96 | 0.000 | 1.000 | 0.000 |  |  | 16 |
| Sedentary | full | mondrian | 96 | 0.125 | 0.000 | 1.562 |  |  | 16 |
| Sedentary | full | floored | 96 | 0.125 | 0.000 | 1.562 |  |  | 16 |

| model | n_recal | n_test_individuals | draws | coverage_recalibrated | coverage_source_cal_same_test | frac_empty_recalibrated | frac_empty_source_cal | set_size_recalibrated |
|---|---|---|---|---|---|---|---|---|
| k10 | 3 | 9 | 20 | 0.946 | 0.054 | 0.000 | 0.738 | 14.797 |
| k10 | 5 | 7 | 20 | 0.920 | 0.050 | 0.000 | 0.747 | 13.671 |
| k20 | 3 | 9 | 20 | 0.925 | 0.049 | 0.000 | 0.951 | 12.718 |
| k20 | 5 | 7 | 20 | 0.899 | 0.054 | 0.000 | 0.946 | 11.865 |
| k50 | 3 | 9 | 20 | 0.933 | 0.000 | 0.000 | 1.000 | 4.601 |
| k50 | 5 | 7 | 20 | 0.913 | 0.000 | 0.000 | 1.000 | 4.256 |
| full | 3 | 9 | 20 | 0.935 | 0.000 | 0.000 | 1.000 | 4.631 |
| full | 5 | 7 | 20 | 0.905 | 0.000 | 0.000 | 1.000 | 4.132 |

Per tissue, k20 (marginal / Mondrian / floored, MoTrPAC calibration):

| target_tissue | floored | marginal | mondrian |
|---|---|---|---|
| BAT | 0.00 | 0.00 | 0.00 |
| HEART | 0.29 | 0.00 | 0.29 |
| HYPOTHALAMUS | 0.00 | 0.00 | 0.00 |
| LIVER | 0.42 | 0.42 | 0.42 |
| MUSCLE | 0.00 | 0.00 | 0.00 |
| SERUM | 0.00 | 0.00 | 0.00 |
| eWAT | 0.00 | 0.00 | 0.00 |
| iWAT | 0.00 | 0.00 | 0.00 |

OOD tissues at k20: `ood_sets.csv`.

Panel at k20: `panel_k20.csv`.
