# Run 2 B — two extensions

Built by `scripts/multiomic/09_extensions.py` on 2026-09-27 09:00 UTC; pre-registered in `docs/PREREGISTRATION_MULTIOMIC.md` (run 2).

## B1 — RNA vs protein panels — `panel_overlap.csv`, `overlap_summary.csv`

| quantity | value |
|---|---|
| n_rna_candidates | 51 |
| n_rna_k20 | 20 |
| n_protein_candidates | 35 |
| n_protein_k20 | 20 |
| n_protein_stable_core_ge_0.8 | 16 |
| k20_intersection | 0 |
| rna_k20_x_protein_core | 0 |
| genes_in_both_candidate_lists | 0 |
| genes_in_both_with_rna_marker_in_prot7 | 0 |
| n_markers_agree | 0 |
| frac_markers_agree | nan |
| prereg_B1_threshold | 0.7 |
| prereg_B1_pass | False |
| genes_in_both |  |

| ensembl_gene | gene_symbol | in_rna_k20 | rna_stability_frequency | rna_marker_tissue | in_protein_k20 | protein_stability_frequency | protein_in_stable_core | protein_marker_tissue | in_both_layers | rna_marker_in_prot7 | markers_agree |
|---|---|---|---|---|---|---|---|---|---|---|---|
| ENSRNOG00000001091 | Hip1r | False |  | nan | True | 0.68 | False | CORTEX | False | None | None |
| ENSRNOG00000002343 | Uchl1 | False |  | nan | True | 1.00 | True | CORTEX | False | None | None |
| ENSRNOG00000005096 | Bzw2 | False |  | nan | True | 0.90 | True | CORTEX | False | None | None |
| ENSRNOG00000005796 | Ctnna1 | False |  | nan | True | 1.00 | True | LUNG | False | None | None |
| ENSRNOG00000007862 | Acat1 | False |  | nan | True | 1.00 | True | KIDNEY | False | None | None |
| ENSRNOG00000011483 | S100a9 | False |  | nan | True | 1.00 | True | LUNG | False | None | None |
| ENSRNOG00000016610 | Arhgap1 | False |  | nan | True | 0.96 | True | LUNG | False | None | None |
| ENSRNOG00000016690 | Idi1 | False |  | nan | True | 0.72 | False | CORTEX | False | None | None |
| ENSRNOG00000017209 | Tubb3 | False |  | nan | True | 1.00 | True | CORTEX | False | None | None |
| ENSRNOG00000017328 | Pter | False |  | nan | True | 0.80 | True | KIDNEY | False | None | None |
| ENSRNOG00000019688 | Diaph1 | False |  | nan | True | 1.00 | True | LUNG | False | None | None |
| ENSRNOG00000019689 | Vwf | False |  | nan | True | 0.84 | True | LUNG | False | None | None |
| ENSRNOG00000020235 | Hnrnpl | False |  | nan | True | 0.70 | False | LIVER | False | None | None |
| ENSRNOG00000020425 | Stim1 | False |  | nan | True | 0.70 | False | SKM-GN | False | None | None |
| ENSRNOG00000021735 | Akr1c15 | False |  | nan | True | 0.98 | True | LUNG | False | None | None |
| ENSRNOG00000022593 | Pdpr | False |  | nan | True | 1.00 | True | CORTEX | False | None | None |
| ENSRNOG00000026415 | Col14a1 | False |  | nan | True | 1.00 | True | WAT-SC | False | None | None |
| ENSRNOG00000028319 | Cryz | False |  | nan | True | 1.00 | True | KIDNEY | False | None | None |
| ENSRNOG00000033835 | Dnm1 | False |  | nan | True | 0.82 | True | CORTEX | False | None | None |
| ENSRNOG00000043094 | Oxct1 | False |  | nan | True | 1.00 | True | KIDNEY | False | None | None |

## B2 — the protein ladder inside MoTrPAC (RII) beside the RNA ladder — `protein_ladder.csv`, `ladder_side_by_side.csv`, `protein_ladder_recalibration.csv`

| split | arm | layer | accuracy | accuracy_sd | coverage | coverage_sd | set_size | frac_empty | n_test_samples | n_test_animals | n_classes |
|---|---|---|---|---|---|---|---|---|---|---|---|
| in_distribution | full | protein (RII) | 1.000 | 0.000 | 0.919 | 0.033 | 0.919 | 0.081 | 420 | 60 | 7 |
| in_distribution | k20 | protein (RII) | 1.000 | 0.000 | 0.914 | 0.026 | 0.914 | 0.086 | 420 | 60 | 7 |
| train_male_test_female | full | protein (RII) | 1.000 |  | 0.814 |  | 0.814 | 0.186 | 210 | 30 | 7 |
| train_male_test_female | k20 | protein (RII) | 1.000 |  | 1.000 |  | 1.000 | 0.000 | 210 | 30 | 7 |
| train_female_test_male | full | protein (RII) | 1.000 |  | 0.876 |  | 0.876 | 0.124 | 210 | 30 | 7 |
| train_female_test_male | k20 | protein (RII) | 1.000 |  | 0.876 |  | 0.876 | 0.124 | 210 | 30 | 7 |
| train_control_test_trained | full | protein (RII) | 1.000 |  | 0.970 |  | 0.970 | 0.030 | 336 | 48 | 7 |
| train_control_test_trained | k20 | protein (RII) | 1.000 |  | 0.958 |  | 0.958 | 0.042 | 336 | 48 | 7 |
| train_male_test_female | full | RNA (counts) | 0.913 |  | 0.764 |  | 0.780 | 0.220 | 449 | 25 | 19 |
| train_male_test_female | k20 | RNA (counts) | 0.893 |  | 0.788 |  | 0.817 | 0.183 | 449 | 25 | 19 |
| train_female_test_male | full | RNA (counts) | 0.942 |  | 0.824 |  | 0.824 | 0.176 | 450 | 25 | 19 |
| train_female_test_male | k20 | RNA (counts) | 0.900 |  | 0.833 |  | 0.833 | 0.167 | 450 | 25 | 19 |
| train_control_test_trained | full | RNA (counts) | 0.979 |  | 0.880 |  | 0.890 | 0.110 | 719 | 40 | 19 |
| train_control_test_trained | k20 | RNA (counts) | 0.961 |  | 0.903 |  | 0.919 | 0.081 | 719 | 40 | 19 |
| in_distribution | k20 | RNA (counts) | 0.976 |  | 0.900 |  | 0.909 | 0.091 | 899 | 50 | 19 |
| in_distribution | full | RNA (counts) |  |  | 0.908 |  | 0.910 | 0.090 | 899 | 50 | 19 |

| split | arm | n_recal_animals | repeats | coverage_target_recalibrated | coverage_target_recal_sd | coverage_target_source_cal_same_test | set_size_recalibrated | set_size_source_cal | n_recal_dropped_unseen_mean | n_classes |
|---|---|---|---|---|---|---|---|---|---|---|
| train_male_test_female | full | 3 | 10 | 0.923 | 0.027 | 0.813 | 0.923 | 0.813 | 0.000 | 7 |
| train_male_test_female | full | 5 | 10 | 0.937 | 0.039 | 0.811 | 0.937 | 0.811 | 0.000 | 7 |
| train_male_test_female | panel_k20 | 3 | 10 | 0.912 | 0.050 | 1.000 | 0.912 | 1.000 | 0.000 | 7 |
| train_male_test_female | panel_k20 | 5 | 10 | 0.934 | 0.041 | 1.000 | 0.934 | 1.000 | 0.000 | 7 |
| train_female_test_male | full | 3 | 10 | 0.915 | 0.037 | 0.876 | 0.915 | 0.876 | 0.000 | 7 |
| train_female_test_male | full | 5 | 10 | 0.927 | 0.019 | 0.877 | 0.927 | 0.877 | 0.000 | 7 |
| train_female_test_male | panel_k20 | 3 | 10 | 0.931 | 0.037 | 0.877 | 0.931 | 0.877 | 0.000 | 7 |
| train_female_test_male | panel_k20 | 5 | 10 | 0.930 | 0.036 | 0.875 | 0.930 | 0.875 | 0.000 | 7 |
| train_control_test_trained | full | 3 | 10 | 0.943 | 0.034 | 0.970 | 0.943 | 0.970 | 0.000 | 7 |
| train_control_test_trained | full | 5 | 10 | 0.911 | 0.022 | 0.969 | 0.911 | 0.969 | 0.000 | 7 |
| train_control_test_trained | panel_k20 | 3 | 10 | 0.931 | 0.041 | 0.958 | 0.931 | 0.958 | 0.000 | 7 |
| train_control_test_trained | panel_k20 | 5 | 10 | 0.915 | 0.033 | 0.958 | 0.915 | 0.958 | 0.000 | 7 |

Pre-registered rules: B2_i_held_out_animals PASS, B2_ii_controls_to_trained PASS, B2_iii_held_out_sex_either PASS
