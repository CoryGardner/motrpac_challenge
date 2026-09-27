# Run 2 A — verification

Built by `scripts/multiomic/08_verification.py` on 2026-09-27 09:00 UTC; every number from a CSV in this directory or the one it names.

## Per-tissue accuracy, k20, with n samples and n donors — `per_tissue_k20.csv`

| run | layer | jiang_tissue | rat_class | n_samples | n_donors | accuracy | top_prediction |
|---|---|---|---|---|---|---|---|
| jiang_relative_k20 | protein | Brain - Cortex | CORTEX | 2 | 2 | 1.000 | CORTEX |
| jiang_relative_k20 | protein | Heart - Atrial Appendage | HEART | 11 | 11 | 0.000 | KIDNEY |
| jiang_relative_k20 | protein | Heart - Left Ventricle | HEART | 7 | 7 | 0.000 | KIDNEY |
| jiang_relative_k20 | protein | Liver | LIVER | 5 | 5 | 0.800 | LIVER |
| jiang_relative_k20 | protein | Lung | LUNG | 8 | 8 | 0.875 | LUNG |
| jiang_relative_k20 | protein | Muscle - Skeletal | SKM-GN | 11 | 11 | 0.636 | SKM-GN |
| jiang_rawppm_k20 | protein | Brain - Cortex | CORTEX | 4 | 2 | 1.000 | CORTEX |
| jiang_rawppm_k20 | protein | Heart - Atrial Appendage | HEART | 23 | 11 | 0.348 | KIDNEY |
| jiang_rawppm_k20 | protein | Heart - Left Ventricle | HEART | 15 | 7 | 0.733 | HEART |
| jiang_rawppm_k20 | protein | Liver | LIVER | 11 | 5 | 0.818 | LIVER |
| jiang_rawppm_k20 | protein | Lung | LUNG | 17 | 8 | 0.882 | LUNG |
| jiang_rawppm_k20 | protein | Muscle - Skeletal | SKM-GN | 24 | 11 | 0.917 | SKM-GN |
| same42_RNA_k20 | RNA | Brain - Cortex | CORTEX | 2 | 2 | 1.000 | CORTEX |
| same42_RNA_k20 | RNA | Heart - Atrial Appendage;Heart - Left Ventricle | HEART | 17 | 11 | 1.000 | HEART |
| same42_RNA_k20 | RNA | Liver | LIVER | 5 | 5 | 1.000 | LIVER |
| same42_RNA_k20 | RNA | Lung | LUNG | 8 | 8 | 1.000 | LUNG |
| same42_RNA_k20 | RNA | Muscle - Skeletal | SKM-GN | 10 | 10 | 1.000 | SKM-GN |
| same42_protein_k20 | protein | Brain - Cortex | CORTEX | 2 | 2 | 1.000 | CORTEX |
| same42_protein_k20 | protein | Heart - Atrial Appendage;Heart - Left Ventricle | HEART | 17 | 11 | 0.000 | KIDNEY |
| same42_protein_k20 | protein | Liver | LIVER | 5 | 5 | 0.800 | LIVER |
| same42_protein_k20 | protein | Lung | LUNG | 8 | 8 | 0.875 | LUNG |
| same42_protein_k20 | protein | Muscle - Skeletal | SKM-GN | 10 | 10 | 0.700 | SKM-GN |
| same42_late_mean_k20 | late_mean | Brain - Cortex | CORTEX | 2 | 2 | 1.000 | CORTEX |
| same42_late_mean_k20 | late_mean | Heart - Atrial Appendage;Heart - Left Ventricle | HEART | 17 | 11 | 0.353 | KIDNEY |
| same42_late_mean_k20 | late_mean | Liver | LIVER | 5 | 5 | 1.000 | LIVER |
| same42_late_mean_k20 | late_mean | Lung | LUNG | 8 | 8 | 1.000 | LUNG |
| same42_late_mean_k20 | late_mean | Muscle - Skeletal | SKM-GN | 10 | 10 | 1.000 | SKM-GN |
| same42_stacked_LR_k20 | stacked_LR | Brain - Cortex | CORTEX | 2 | 2 | 1.000 | CORTEX |
| same42_stacked_LR_k20 | stacked_LR | Heart - Atrial Appendage;Heart - Left Ventricle | HEART | 17 | 11 | 0.353 | KIDNEY |
| same42_stacked_LR_k20 | stacked_LR | Liver | LIVER | 5 | 5 | 1.000 | LIVER |
| same42_stacked_LR_k20 | stacked_LR | Lung | LUNG | 8 | 8 | 1.000 | LUNG |
| same42_stacked_LR_k20 | stacked_LR | Muscle - Skeletal | SKM-GN | 10 | 10 | 1.000 | SKM-GN |

## k20 accuracies with donor-bootstrap intervals — `accuracy_ci.csv`

| run | layer | model | accuracy | ci95_low | ci95_high | n_samples | n_donors | bootstrap_unit | n_boot | source |
|---|---|---|---|---|---|---|---|---|---|---|
| jiang_relative | protein | k20 | 0.455 | 0.357 | 0.553 | 44 | 13 | donor | 1000 | results_multiomic/03_prot_transfer/accuracy_overall.csv |
| jiang_rawppm | protein | k20 | 0.734 | 0.640 | 0.827 | 94 | 13 | donor | 1000 | results_multiomic/03_prot_transfer/rawppm/accuracy_overall.csv |
| same42 | RNA | k20 | 1.000 | 1.000 | 1.000 | 42 | 12 | donor | 1000 | results_multiomic/05_fusion_transfer/fusion_transfer_summary.csv |
| same42 | protein | k20 | 0.476 | 0.386 | 0.571 | 42 | 12 | donor | 1000 | results_multiomic/05_fusion_transfer/fusion_transfer_summary.csv |
| same42 | late_mean | k20 | 0.738 | 0.611 | 0.857 | 42 | 12 | donor | 1000 | results_multiomic/05_fusion_transfer/fusion_transfer_summary.csv |
| same42 | stacked_LR | k20 | 0.738 | 0.611 | 0.857 | 42 | 12 | donor | 1000 | results_multiomic/05_fusion_transfer/fusion_transfer_summary.csv |

## Every recalibrated coverage with its set size and label space — `recalibration_set_sizes.csv`

| source | model | n_recal | n_test_individuals | draws | coverage_recalibrated | coverage_source_cal_same_test | set_size_recalibrated | n_classes_label_space | set_size_frac_of_classes | frac_empty_recalibrated | note |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Phase 3 protein → Jiang 2020 (cleaned relative) | k20 | 3 | 11 | 20 | 0.900 | 0.070 | 3.497 | 7 | 0.500 | 0.000 |  |
| Phase 3 protein → Jiang 2020 (cleaned relative) | k20 | 5 | 9 | 20 | 0.917 | 0.077 | 2.954 | 7 | 0.422 | 0.000 |  |
| Phase 3 protein → Jiang 2020 (cleaned relative) | k50 | 3 | 11 | 20 | 0.926 | 0.046 | 4.027 | 7 | 0.575 | 0.000 |  |
| Phase 3 protein → Jiang 2020 (cleaned relative) | k50 | 5 | 9 | 20 | 0.968 | 0.053 | 3.110 | 7 | 0.444 | 0.000 |  |
| Phase 3 protein → Jiang 2020 (cleaned relative) | full | 3 | 11 | 20 | 0.960 | 0.050 | 5.359 | 7 | 0.766 | 0.000 |  |
| Phase 3 protein → Jiang 2020 (cleaned relative) | full | 5 | 9 | 20 | 0.950 | 0.045 | 4.348 | 7 | 0.621 | 0.000 |  |
| Phase 3 protein → Jiang 2020 (raw ppm) | k20 | 3 | 11 | 20 | 0.921 | 0.081 | 1.634 | 7 | 0.233 | 0.000 |  |
| Phase 3 protein → Jiang 2020 (raw ppm) | k20 | 5 | 9 | 20 | 0.923 | 0.085 | 1.538 | 7 | 0.220 | 0.000 |  |
| Phase 3 protein → Jiang 2020 (raw ppm) | k50 | 3 | 11 | 20 | 0.916 | 0.063 | 1.483 | 7 | 0.212 | 0.001 |  |
| Phase 3 protein → Jiang 2020 (raw ppm) | k50 | 5 | 9 | 20 | 0.917 | 0.071 | 1.445 | 7 | 0.206 | 0.000 |  |
| Phase 3 protein → Jiang 2020 (raw ppm) | full | 3 | 11 | 20 | 0.904 | 0.000 | 1.419 | 7 | 0.203 | 0.011 |  |
| Phase 3 protein → Jiang 2020 (raw ppm) | full | 5 | 9 | 20 | 0.911 | 0.000 | 1.223 | 7 | 0.175 | 0.002 |  |
| Phase 3b protein → Wang 2019 | k20 | 3 | 26 | 9 | 1.000 | 0.344 | 7.000 | 7 | 1.000 | 0.000 | recalibration individuals are single tissue samples |
| Phase 3b protein → Wang 2019 | k20 | 5 | 24 | 16 | 1.000 | 0.350 | 7.000 | 7 | 1.000 | 0.000 | recalibration individuals are single tissue samples |
| Phase 3b protein → Wang 2019 | k50 | 3 | 26 | 9 | 1.000 | 0.317 | 7.000 | 7 | 1.000 | 0.000 | recalibration individuals are single tissue samples |
| Phase 3b protein → Wang 2019 | k50 | 5 | 24 | 16 | 1.000 | 0.355 | 7.000 | 7 | 1.000 | 0.000 | recalibration individuals are single tissue samples |
| Phase 3b protein → Wang 2019 | full | 3 | 26 | 11 | 1.000 | 0.127 | 7.000 | 7 | 1.000 | 0.000 | recalibration individuals are single tissue samples |
| Phase 3b protein → Wang 2019 | full | 5 | 24 | 11 | 1.000 | 0.155 | 7.000 | 7 | 1.000 | 0.000 | recalibration individuals are single tissue samples |
| Phase 3b protein → Geiger 2013 | k20 | 3 | 26 | 14 | 1.000 | 0.119 | 7.000 | 7 | 1.000 | 0.000 | recalibration individuals are single tissue samples |
| Phase 3b protein → Geiger 2013 | k20 | 5 | 24 | 19 | 1.000 | 0.129 | 7.000 | 7 | 1.000 | 0.000 | recalibration individuals are single tissue samples |
| Phase 3b protein → Geiger 2013 | k50 | 3 | 26 | 12 | 1.000 | 0.112 | 7.000 | 7 | 1.000 | 0.000 | recalibration individuals are single tissue samples |
| Phase 3b protein → Geiger 2013 | k50 | 5 | 24 | 16 | 1.000 | 0.128 | 7.000 | 7 | 1.000 | 0.000 | recalibration individuals are single tissue samples |
| Phase 3b protein → Geiger 2013 | full | 3 | 26 | 16 | 1.000 | 0.000 | 7.000 | 7 | 1.000 | 0.000 | recalibration individuals are single tissue samples |
| Phase 3b protein → Geiger 2013 | full | 5 | 24 | 19 | 1.000 | 0.000 | 7.000 | 7 | 1.000 | 0.000 | recalibration individuals are single tissue samples |
| Phase 4 hilic_sato | k10 | 3 | 9 | 20 | 0.946 | 0.054 | 14.797 | 19 | 0.779 | 0.000 |  |
| Phase 4 hilic_sato | k10 | 5 | 7 | 20 | 0.920 | 0.050 | 13.671 | 19 | 0.720 | 0.000 |  |
| Phase 4 hilic_sato | k20 | 3 | 9 | 20 | 0.925 | 0.049 | 12.718 | 19 | 0.669 | 0.000 |  |
| Phase 4 hilic_sato | k20 | 5 | 7 | 20 | 0.899 | 0.054 | 11.865 | 19 | 0.624 | 0.000 |  |
| Phase 4 hilic_sato | k50 | 3 | 9 | 20 | 0.933 | 0.000 | 4.601 | 19 | 0.242 | 0.000 |  |
| Phase 4 hilic_sato | k50 | 5 | 7 | 20 | 0.913 | 0.000 | 4.256 | 19 | 0.224 | 0.000 |  |
| Phase 4 hilic_sato | full | 3 | 9 | 20 | 0.935 | 0.000 | 4.631 | 19 | 0.244 | 0.000 |  |
| Phase 4 hilic_sato | full | 5 | 7 | 20 | 0.905 | 0.000 | 4.132 | 19 | 0.217 | 0.000 |  |
| Phase 4 deep_sato | k10 | 3 | 9 | 20 | 0.905 | 0.000 | 3.617 | 9 | 0.402 | 0.000 |  |
| Phase 4 deep_sato | k10 | 5 | 7 | 20 | 0.918 | 0.000 | 3.703 | 9 | 0.411 | 0.000 |  |
| Phase 4 deep_sato | k20 | 3 | 9 | 20 | 0.910 | 0.000 | 5.352 | 9 | 0.595 | 0.000 |  |
| Phase 4 deep_sato | k20 | 5 | 7 | 20 | 0.910 | 0.000 | 5.277 | 9 | 0.586 | 0.000 |  |
| Phase 4 deep_sato | k50 | 3 | 9 | 20 | 0.913 | 0.000 | 5.071 | 9 | 0.563 | 0.000 |  |
| Phase 4 deep_sato | k50 | 5 | 7 | 20 | 0.926 | 0.000 | 5.149 | 9 | 0.572 | 0.000 |  |
| Phase 4 deep_sato | full | 3 | 9 | 20 | 0.916 | 0.000 | 5.788 | 9 | 0.643 | 0.000 |  |
| Phase 4 deep_sato | full | 5 | 7 | 20 | 0.909 | 0.000 | 5.744 | 9 | 0.638 | 0.000 |  |
| Phase 4 deep_mw | k10 | 3 | 11 | 20 | 0.931 | 0.155 | 4.701 | 9 | 0.522 | 0.000 |  |
| Phase 4 deep_mw | k10 | 5 | 9 | 20 | 0.921 | 0.155 | 4.629 | 9 | 0.514 | 0.000 |  |
| Phase 4 deep_mw | k20 | 3 | 11 | 20 | 0.924 | 0.154 | 2.344 | 9 | 0.260 | 0.001 |  |
| Phase 4 deep_mw | k20 | 5 | 9 | 20 | 0.920 | 0.146 | 2.233 | 9 | 0.248 | 0.000 |  |
| Phase 4 deep_mw | full | 3 | 11 | 20 | 0.912 | 0.000 | 1.549 | 9 | 0.172 | 0.011 |  |
| Phase 4 deep_mw | full | 5 | 9 | 20 | 0.928 | 0.000 | 1.688 | 9 | 0.188 | 0.006 |  |
| frozen phase 13 RNA → GTEx | k20 | 3 | 859 | 20 | 0.954 | 0.364 | 11.697 | 19 | 0.616 | 0.009 | the RNA ladder rows quoted for comparison |
| frozen phase 13 RNA → GTEx | k20 | 5 | 857 | 20 | 0.933 | 0.364 | 6.282 | 19 | 0.331 | 0.014 | the RNA ladder rows quoted for comparison |
| frozen phase 13 RNA → GTEx | k50 | 3 | 859 | 20 | 0.959 | 0.363 | 11.624 | 19 | 0.612 | 0.008 | the RNA ladder rows quoted for comparison |
| frozen phase 13 RNA → GTEx | k50 | 5 | 857 | 20 | 0.936 | 0.363 | 4.036 | 19 | 0.212 | 0.000 | the RNA ladder rows quoted for comparison |
| frozen phase 13 RNA → GTEx | full | 3 | 859 | 20 | 0.951 | 0.062 | 11.178 | 19 | 0.588 | 0.019 | the RNA ladder rows quoted for comparison |
| frozen phase 13 RNA → GTEx | full | 5 | 857 | 20 | 0.866 | 0.062 | 2.045 | 19 | 0.108 | 0.071 | the RNA ladder rows quoted for comparison |
| frozen phase 12 RNA → BodyMap | k20 | 3 | 5 | 20 | 0.943 | 0.619 | 0.998 | 19 | 0.053 | 0.057 | the RNA ladder rows quoted for comparison |
| frozen phase 12 RNA → BodyMap | k20 | 5 | 3 | 20 | 0.914 | 0.620 | 0.920 | 19 | 0.048 | 0.086 | the RNA ladder rows quoted for comparison |
| frozen phase 12 RNA → BodyMap | k50 | 3 | 5 | 20 | 0.937 | 0.662 | 0.953 | 19 | 0.050 | 0.047 | the RNA ladder rows quoted for comparison |
| frozen phase 12 RNA → BodyMap | k50 | 5 | 3 | 20 | 0.920 | 0.658 | 0.934 | 19 | 0.049 | 0.066 | the RNA ladder rows quoted for comparison |
| frozen phase 12 RNA → BodyMap | full | 3 | 5 | 20 | 0.946 | 0.711 | 0.946 | 19 | 0.050 | 0.054 | the RNA ladder rows quoted for comparison |
| frozen phase 12 RNA → BodyMap | full | 5 | 3 | 20 | 0.919 | 0.710 | 0.919 | 19 | 0.048 | 0.081 | the RNA ladder rows quoted for comparison |

## Sato 2022 design behind the invariance test — `sato_design.csv`, `sato_invariance_design.csv`

| tissue | n_samples | n_mice | n_sedentary_mice | n_exercised_mice | n_sedentary_samples | n_exercised_samples | time_after_exercise_h | zeitgeber_time_points |
|---|---|---|---|---|---|---|---|---|
| BAT | 24 | 24 | 12 | 12 | 12 | 12 | 0 | 4;16 |
| eWAT | 24 | 24 | 12 | 12 | 12 | 12 | 0 | 4;16 |
| HEART | 24 | 24 | 12 | 12 | 12 | 12 | 0 | 4;16 |
| HYPOTHALAMUS | 24 | 24 | 12 | 12 | 12 | 12 | 0 | 4;16 |
| iWAT | 23 | 23 | 12 | 11 | 12 | 11 | 0 | 4;16 |
| LIVER | 24 | 24 | 12 | 12 | 12 | 12 | 0 | 4;16 |
| MUSCLE | 24 | 24 | 12 | 12 | 12 | 12 |  | 4;16 |
| SERUM | 24 | 24 | 12 | 12 | 12 | 12 | 0 | 4;16 |

| quantity | value |
|---|---|
| n_sedentary_mice_union | 12 |
| n_exercised_mice_union | 12 |
| n_sedentary_samples_total | 96 |
| n_exercised_samples_total | 95 |
| animal_ids_shared_across_tissues | True |
| exercise_design | acute: a single treadmill bout; the Metabolon sample metadata records TIME AFTER EXERCISE = 0 h and two Zeitgeber time points (ZT4 light, ZT16 dark) with 6 mice per treatment × time point per tissue; no training programme (Sato et al. 2022, Cell Metab 34:329, 'time-dependent signatures of metabolic homeostasis' after acute exercise) |
| invariance_fit_animals | 8 |
| invariance_cal_animals | 4 |
| invariance_test_exercised_samples | 95 |
| invariance_test_exercised_animals | 12 |
| invariance_k20_accuracy | 0.9473684210526316 |
| invariance_k20_coverage | 0.9052631578947368 |
| invariance_k20_set_size | 0.9157894736842104 |
| invariance_full_accuracy | 1.0 |
| invariance_full_coverage | 0.9263157894736842 |
| invariance_full_set_size | 0.9263157894736842 |
| invariance_n_classes | 8 |

## Flags — `flags.csv`

| item | flag | value | wording |
|---|---|---|---|
| Phase 3 protein → Jiang 2020 (cleaned relative) k50 n_recal=3 | recalibrated set holds ≥ half the label space | set size 4.03 of 7 classes | the number is restored, the information is not (GTEx pattern); never quote this coverage without the set size |
| Phase 3 protein → Jiang 2020 (cleaned relative) full n_recal=3 | recalibrated set holds ≥ half the label space | set size 5.36 of 7 classes | the number is restored, the information is not (GTEx pattern); never quote this coverage without the set size |
| Phase 3 protein → Jiang 2020 (cleaned relative) full n_recal=5 | recalibrated set holds ≥ half the label space | set size 4.35 of 7 classes | the number is restored, the information is not (GTEx pattern); never quote this coverage without the set size |
| Phase 3b protein → Wang 2019 k20 n_recal=3 | recalibrated set holds ≥ half the label space | set size 7.00 of 7 classes | the number is restored, the information is not (GTEx pattern); never quote this coverage without the set size |
| Phase 3b protein → Wang 2019 k20 n_recal=5 | recalibrated set holds ≥ half the label space | set size 7.00 of 7 classes | the number is restored, the information is not (GTEx pattern); never quote this coverage without the set size |
| Phase 3b protein → Wang 2019 k50 n_recal=3 | recalibrated set holds ≥ half the label space | set size 7.00 of 7 classes | the number is restored, the information is not (GTEx pattern); never quote this coverage without the set size |
| Phase 3b protein → Wang 2019 k50 n_recal=5 | recalibrated set holds ≥ half the label space | set size 7.00 of 7 classes | the number is restored, the information is not (GTEx pattern); never quote this coverage without the set size |
| Phase 3b protein → Wang 2019 full n_recal=3 | recalibrated set holds ≥ half the label space | set size 7.00 of 7 classes | the number is restored, the information is not (GTEx pattern); never quote this coverage without the set size |
| Phase 3b protein → Wang 2019 full n_recal=5 | recalibrated set holds ≥ half the label space | set size 7.00 of 7 classes | the number is restored, the information is not (GTEx pattern); never quote this coverage without the set size |
| Phase 3b protein → Geiger 2013 k20 n_recal=3 | recalibrated set holds ≥ half the label space | set size 7.00 of 7 classes | the number is restored, the information is not (GTEx pattern); never quote this coverage without the set size |
| Phase 3b protein → Geiger 2013 k20 n_recal=5 | recalibrated set holds ≥ half the label space | set size 7.00 of 7 classes | the number is restored, the information is not (GTEx pattern); never quote this coverage without the set size |
| Phase 3b protein → Geiger 2013 k50 n_recal=3 | recalibrated set holds ≥ half the label space | set size 7.00 of 7 classes | the number is restored, the information is not (GTEx pattern); never quote this coverage without the set size |
| Phase 3b protein → Geiger 2013 k50 n_recal=5 | recalibrated set holds ≥ half the label space | set size 7.00 of 7 classes | the number is restored, the information is not (GTEx pattern); never quote this coverage without the set size |
| Phase 3b protein → Geiger 2013 full n_recal=3 | recalibrated set holds ≥ half the label space | set size 7.00 of 7 classes | the number is restored, the information is not (GTEx pattern); never quote this coverage without the set size |
| Phase 3b protein → Geiger 2013 full n_recal=5 | recalibrated set holds ≥ half the label space | set size 7.00 of 7 classes | the number is restored, the information is not (GTEx pattern); never quote this coverage without the set size |
| Phase 4 hilic_sato k10 n_recal=3 | recalibrated set holds ≥ half the label space | set size 14.80 of 19 classes | the number is restored, the information is not (GTEx pattern); never quote this coverage without the set size |
| Phase 4 hilic_sato k10 n_recal=5 | recalibrated set holds ≥ half the label space | set size 13.67 of 19 classes | the number is restored, the information is not (GTEx pattern); never quote this coverage without the set size |
| Phase 4 hilic_sato k20 n_recal=3 | recalibrated set holds ≥ half the label space | set size 12.72 of 19 classes | the number is restored, the information is not (GTEx pattern); never quote this coverage without the set size |
| Phase 4 hilic_sato k20 n_recal=5 | recalibrated set holds ≥ half the label space | set size 11.87 of 19 classes | the number is restored, the information is not (GTEx pattern); never quote this coverage without the set size |
| Phase 4 deep_sato k20 n_recal=3 | recalibrated set holds ≥ half the label space | set size 5.35 of 9 classes | the number is restored, the information is not (GTEx pattern); never quote this coverage without the set size |
| Phase 4 deep_sato k20 n_recal=5 | recalibrated set holds ≥ half the label space | set size 5.28 of 9 classes | the number is restored, the information is not (GTEx pattern); never quote this coverage without the set size |
| Phase 4 deep_sato k50 n_recal=3 | recalibrated set holds ≥ half the label space | set size 5.07 of 9 classes | the number is restored, the information is not (GTEx pattern); never quote this coverage without the set size |
| Phase 4 deep_sato k50 n_recal=5 | recalibrated set holds ≥ half the label space | set size 5.15 of 9 classes | the number is restored, the information is not (GTEx pattern); never quote this coverage without the set size |
| Phase 4 deep_sato full n_recal=3 | recalibrated set holds ≥ half the label space | set size 5.79 of 9 classes | the number is restored, the information is not (GTEx pattern); never quote this coverage without the set size |
| Phase 4 deep_sato full n_recal=5 | recalibrated set holds ≥ half the label space | set size 5.74 of 9 classes | the number is restored, the information is not (GTEx pattern); never quote this coverage without the set size |
| Phase 4 deep_mw k10 n_recal=3 | recalibrated set holds ≥ half the label space | set size 4.70 of 9 classes | the number is restored, the information is not (GTEx pattern); never quote this coverage without the set size |
| Phase 4 deep_mw k10 n_recal=5 | recalibrated set holds ≥ half the label space | set size 4.63 of 9 classes | the number is restored, the information is not (GTEx pattern); never quote this coverage without the set size |
| frozen phase 13 RNA → GTEx k20 n_recal=3 | recalibrated set holds ≥ half the label space | set size 11.70 of 19 classes | the number is restored, the information is not (GTEx pattern); never quote this coverage without the set size |
| frozen phase 13 RNA → GTEx k50 n_recal=3 | recalibrated set holds ≥ half the label space | set size 11.62 of 19 classes | the number is restored, the information is not (GTEx pattern); never quote this coverage without the set size |
| frozen phase 13 RNA → GTEx full n_recal=3 | recalibrated set holds ≥ half the label space | set size 11.18 of 19 classes | the number is restored, the information is not (GTEx pattern); never quote this coverage without the set size |
| Jiang k20 n_recal=3 | pattern classification | coverage 0.900 at set size 3.50 of 7 | GTEx pattern (number restored, information partly: sets hold 50% of the label space), not the BodyMap pattern (one tissue per set) |
| Jiang k20 n_recal=5 | pattern classification | coverage 0.917 at set size 2.95 of 7 | GTEx pattern (number restored, information partly: sets hold 42% of the label space), not the BodyMap pattern (one tissue per set) |
| Sato invariance test | wording | acute single bout, TIME AFTER EXERCISE = 0 h | 'exercised', never 'trained', for Sato mice; the RNA analogue (controls → trained) is an 8-week training programme, so the two tests are not the same shift |
| RNA k20 on the same 42 Jiang samples | n per class | CORTEX 2 samples / 2 donors; HEART 17 samples / 11 donors; LIVER 5 samples / 5 donors; LUNG 8 samples / 8 donors; SKM-GN 10 samples / 10 donors | two classes have ≤ 5 samples; the 1.000 is over 42 samples from 12 donors |
