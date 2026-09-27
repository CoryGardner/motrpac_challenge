# Numbers reconciliation — what the site shows, where it comes from, and what moved with the conformal-quantile fix

Generated 2026-09-27 04:45 UTC by `scripts/30_export_site_data.py --reconciliation` from `site/data/provenance.json`.
`results/` (post-fix, 2026-09-25) is the truth; the pre-fix values come from `docs/reconciliation/pre_quantile_fix_values.csv`, generated once from the pre-fix results of 2026-09-25 (not in the repository). Values are shown to 4 decimals; the JSON holds them unrounded.

## 1. Headline numbers (the home-page tiles and the transfer ladder)

| id | post-fix value | pre-fix value | changed | source | selector | column | agg |
|---|---|---|---|---|---|---|---|
| `tile_acc_k20` | 0.9758 |  |  | `results/05_panels/TRNSCRPT/panel_curve.csv` | `{"k": 20}` | `balanced_accuracy` | mean |
| `tile_acc_k20_sd` | 0.0080 |  |  | `results/05_panels/TRNSCRPT/panel_curve.csv` | `{"k": 20}` | `balanced_accuracy` | std |
| `acc_k50` | 0.9926 |  |  | `results/05_panels/TRNSCRPT/panel_curve.csv` | `{"k": 50}` | `balanced_accuracy` | mean |
| `acc_full` | 0.9947 |  |  | `results/04_baselines/TRNSCRPT/summary.csv` | `{"model": "logreg_l2"}` | `balanced_accuracy_mean` | value |
| `acc_fclassif_k20` | 0.3989 |  |  | `results/05_panels/TRNSCRPT/panel_curve_fclassif.csv` | `{"k": 20}` | `balanced_accuracy` | mean |
| `tile_bodymap_k20` | 1.0000 | 1.0000 | no | `results/12_bodymap/age_shift_accuracy.csv` | `{"stage_weeks": 21}` | `k20` | value |
| `tile_bodymap_cov_k20` | 0.6176 | 0.6176 | no | `results/12_bodymap/conformal_transfer.csv` | `{"stage_weeks": 21, "model": "k20", "conformal": "marginal"}` | `coverage_mapped` | value |
| `tile_bodymap_empty_k20` | 0.3824 | 0.3824 | no | `results/12_bodymap/conformal_transfer.csv` | `{"stage_weeks": 21, "model": "k20", "conformal": "marginal"}` | `frac_empty_mapped` | value |
| `bodymap_floored_k20` | 0.6912 | 0.6912 | no | `results/12_bodymap/conformal_transfer.csv` | `{"stage_weeks": 21, "model": "k20", "conformal": "floored"}` | `coverage_mapped` | value |
| `recal3_bodymap_k20` | 0.9426 | 0.9696 | yes | `results/12_bodymap/recalibration.csv` | `{"model": "k20", "n_recal": 3}` | `coverage_recalibrated` | value |
| `recal3size_bodymap_k20` | 0.9979 | 1.1058 | yes | `results/12_bodymap/recalibration.csv` | `{"model": "k20", "n_recal": 3}` | `set_size_recalibrated` | value |
| `tile_estimable` | 1 |  |  | `results/16_identifiability/estimable_pairs.csv` | `{"assay": "TRNSCRPT"}` | `n_pairs_estimable` | value |
| `tile_estimable_total` | 171 |  |  | `results/16_identifiability/estimable_pairs.csv` | `{"assay": "TRNSCRPT"}` | `n_pairs_total` | value |
| `cov_id_full_marginal_pooled` | 0.9077 | 0.9110 | yes | `results/06_conformal/TRNSCRPT/coverage.csv` | `{"calibration": "pooled", "conformal": "marginal", "method": "lac", "alpha": 0.1}` | `coverage` | mean |
| `cov_id_full_marginal_one_per_animal` | 0.9165 | 0.9622 | yes | `results/06_conformal/TRNSCRPT/coverage.csv` | `{"calibration": "one_per_animal", "conformal": "marginal", "method": "lac", "alpha": 0.1}` | `coverage` | mean |
| `cov_id_full_mondrian_pooled` | 0.9188 | 0.9600 | yes | `results/06_conformal/TRNSCRPT/coverage.csv` | `{"calibration": "pooled", "conformal": "mondrian", "method": "lac", "alpha": 0.1}` | `coverage` | mean |
| `cov_id_full_floored_pooled` | 0.9700 | 0.9855 | yes | `results/06_conformal/TRNSCRPT/coverage.csv` | `{"calibration": "pooled", "conformal": "floored", "method": "lac", "alpha": 0.1}` | `coverage` | mean |
| `cov_train_male_test_female_k20_marginal` | 0.8329 | 0.8376 | yes | `results/08_shift/TRNSCRPT/shift_table.csv` | `{"split": "train_male_test_female", "arm": "panel_k20"}` | `coverage_target_seen` | value |
| `cov_train_male_test_female_full_marginal` | 0.8071 | 0.8188 | yes | `results/08_shift/TRNSCRPT/shift_table.csv` | `{"split": "train_male_test_female", "arm": "full"}` | `coverage_target_seen` | value |
| `cov_train_female_test_male_k20_marginal` | 0.8824 | 0.8824 | no | `results/08_shift/TRNSCRPT/shift_table.csv` | `{"split": "train_female_test_male", "arm": "panel_k20"}` | `coverage_target_seen` | value |
| `acc_gtex_k20` | 0.6539 | 0.6539 | no | `results/13_gtex/accuracy_overall.csv` | `{"model": "k20"}` | `accuracy_sample_weighted` | value |
| `acc_gtex_k50` | 0.7815 | 0.7815 | no | `results/13_gtex/accuracy_overall.csv` | `{"model": "k50"}` | `accuracy_sample_weighted` | value |
| `acc_gtex_full` | 0.8547 | 0.8547 | no | `results/13_gtex/accuracy_overall.csv` | `{"model": "full"}` | `accuracy_sample_weighted` | value |
| `cov_gtex_k20_marginal` | 0.3638 | 0.3678 | yes | `results/13_gtex/conformal_transfer.csv` | `{"stage": "adult", "model": "k20", "conformal": "marginal"}` | `coverage_mapped` | value |
| `empty_gtex_k20_marginal` | 0.6165 | 0.6113 | yes | `results/13_gtex/conformal_transfer.csv` | `{"stage": "adult", "model": "k20", "conformal": "marginal"}` | `frac_empty_mapped` | value |
| `cov_gtex_full_marginal` | 0.0620 | 0.0636 | yes | `results/13_gtex/conformal_transfer.csv` | `{"stage": "adult", "model": "full", "conformal": "marginal"}` | `coverage_mapped` | value |
| `recal3_gtex_k20` | 0.9544 | 0.8795 | yes | `results/13_gtex/recalibration.csv` | `{"model": "k20", "n_recal": 3}` | `coverage_recalibrated` | value |
| `recal3size_gtex_k20` | 11.6966 | 4.1494 | yes | `results/13_gtex/recalibration.csv` | `{"model": "k20", "n_recal": 3}` | `set_size_recalibrated` | value |
| `gtex_recal_k20_n3_frac_inf` | 0.4500 |  |  | `results/31_site_regen/13_gtex/recal_thresholds.csv` | `{}` | `` | recomputed |
| `qc_technical` | 0.8742 |  |  | `results/16_identifiability/qc_only_summary.csv` | `{"features": "technical"}` | `bal_acc_mean` | value |
| `qc_composition` | 0.9516 |  |  | `results/16_identifiability/qc_only_summary.csv` | `{"features": "composition"}` | `bal_acc_mean` | value |
| `qc_all` | 0.9758 |  |  | `results/16_identifiability/qc_only_summary.csv` | `{"features": "all"}` | `bal_acc_mean` | value |

## 2. Every exported number that moved with the fix (113 of 203 comparable entries)

| id | post-fix | pre-fix | source | selector | column |
|---|---|---|---|---|---|
| `cov_id_full_marginal_pooled` | 0.9077 | 0.9110 | `results/06_conformal/TRNSCRPT/coverage.csv` | `{"calibration": "pooled", "conformal": "marginal", "method": "lac", "alpha": 0.1}` | `coverage` |
| `cov_id_full_marginal_pooled_sd` | 0.0334 | 0.0402 | `results/06_conformal/TRNSCRPT/coverage.csv` | `{"calibration": "pooled", "conformal": "marginal", "method": "lac", "alpha": 0.1}` | `coverage` |
| `empty_id_full_marginal_pooled` | 0.0901 | 0.0867 | `results/06_conformal/TRNSCRPT/coverage.csv` | `{"calibration": "pooled", "conformal": "marginal", "method": "lac", "alpha": 0.1}` | `frac_empty` |
| `size_id_full_marginal_pooled` | 0.9099 | 0.9133 | `results/06_conformal/TRNSCRPT/coverage.csv` | `{"calibration": "pooled", "conformal": "marginal", "method": "lac", "alpha": 0.1}` | `avg_set_size` |
| `cov_id_full_marginal_one_per_animal` | 0.9165 | 0.9622 | `results/06_conformal/TRNSCRPT/coverage.csv` | `{"calibration": "one_per_animal", "conformal": "marginal", "method": "lac", "alpha": 0.1}` | `coverage` |
| `cov_id_full_marginal_one_per_animal_sd` | 0.0635 | 0.0502 | `results/06_conformal/TRNSCRPT/coverage.csv` | `{"calibration": "one_per_animal", "conformal": "marginal", "method": "lac", "alpha": 0.1}` | `coverage` |
| `empty_id_full_marginal_one_per_animal` | 0.0802 | 0.0367 | `results/06_conformal/TRNSCRPT/coverage.csv` | `{"calibration": "one_per_animal", "conformal": "marginal", "method": "lac", "alpha": 0.1}` | `frac_empty` |
| `size_id_full_marginal_one_per_animal` | 0.9198 | 1.1805 | `results/06_conformal/TRNSCRPT/coverage.csv` | `{"calibration": "one_per_animal", "conformal": "marginal", "method": "lac", "alpha": 0.1}` | `avg_set_size` |
| `cov_id_full_mondrian_pooled` | 0.9188 | 0.9600 | `results/06_conformal/TRNSCRPT/coverage.csv` | `{"calibration": "pooled", "conformal": "mondrian", "method": "lac", "alpha": 0.1}` | `coverage` |
| `cov_id_full_mondrian_pooled_sd` | 0.0246 | 0.0082 | `results/06_conformal/TRNSCRPT/coverage.csv` | `{"calibration": "pooled", "conformal": "mondrian", "method": "lac", "alpha": 0.1}` | `coverage` |
| `empty_id_full_mondrian_pooled` | 0.0789 | 0.0334 | `results/06_conformal/TRNSCRPT/coverage.csv` | `{"calibration": "pooled", "conformal": "mondrian", "method": "lac", "alpha": 0.1}` | `frac_empty` |
| `size_id_full_mondrian_pooled` | 0.9211 | 0.9922 | `results/06_conformal/TRNSCRPT/coverage.csv` | `{"calibration": "pooled", "conformal": "mondrian", "method": "lac", "alpha": 0.1}` | `avg_set_size` |
| `cov_id_full_floored_pooled` | 0.9700 | 0.9855 | `results/06_conformal/TRNSCRPT/coverage.csv` | `{"calibration": "pooled", "conformal": "floored", "method": "lac", "alpha": 0.1}` | `coverage` |
| `cov_id_full_floored_pooled_sd` | 0.0155 | 0.0093 | `results/06_conformal/TRNSCRPT/coverage.csv` | `{"calibration": "pooled", "conformal": "floored", "method": "lac", "alpha": 0.1}` | `coverage` |
| `empty_id_full_floored_pooled` | 0.0278 | 0.0089 | `results/06_conformal/TRNSCRPT/coverage.csv` | `{"calibration": "pooled", "conformal": "floored", "method": "lac", "alpha": 0.1}` | `frac_empty` |
| `size_id_full_floored_pooled` | 0.9722 | 1.0178 | `results/06_conformal/TRNSCRPT/coverage.csv` | `{"calibration": "pooled", "conformal": "floored", "method": "lac", "alpha": 0.1}` | `avg_set_size` |
| `cov_train_control_test_trained_k20_marginal` | 0.9026 | 0.9040 | `results/08_shift/TRNSCRPT/shift_table.csv` | `{"split": "train_control_test_trained", "arm": "panel_k20"}` | `coverage_target_seen` |
| `size_train_control_test_trained_k20_marginal` | 0.9193 | 0.9221 | `results/08_shift/TRNSCRPT/shift_table.csv` | `{"split": "train_control_test_trained", "arm": "panel_k20"}` | `avg_set_size_target` |
| `recal3_train_control_test_trained_k20` | 0.9243 | 0.9322 | `results/08_shift/TRNSCRPT/shift_table.csv` | `{"split": "train_control_test_trained", "arm": "panel_k20"}` | `cov_target_recal_N3` |
| `cov_train_control_test_trained_k20_mondrian` | 1.0000 | 0.7107 | `results/08_shift/TRNSCRPT/shift_table.csv` | `{"split": "train_control_test_trained", "arm": "panel_k20"}` | `coverage_target_seen_mondrian` |
| `size_train_control_test_trained_k20_mondrian` | 19.0000 | 0.7177 | `results/08_shift/TRNSCRPT/shift_table.csv` | `{"split": "train_control_test_trained", "arm": "panel_k20"}` | `avg_set_size_target_mondrian` |
| `cov_train_control_test_trained_k20_floored` | 1.0000 | 0.9263 | `results/08_shift/TRNSCRPT/shift_table.csv` | `{"split": "train_control_test_trained", "arm": "panel_k20"}` | `coverage_target_seen_floored` |
| `size_train_control_test_trained_k20_floored` | 19.0000 | 0.9458 | `results/08_shift/TRNSCRPT/shift_table.csv` | `{"split": "train_control_test_trained", "arm": "panel_k20"}` | `avg_set_size_target_floored` |
| `cov_train_control_test_trained_full_marginal` | 0.8804 | 0.8846 | `results/08_shift/TRNSCRPT/shift_table.csv` | `{"split": "train_control_test_trained", "arm": "full"}` | `coverage_target_seen` |
| `size_train_control_test_trained_full_marginal` | 0.8901 | 0.8943 | `results/08_shift/TRNSCRPT/shift_table.csv` | `{"split": "train_control_test_trained", "arm": "full"}` | `avg_set_size_target` |
| `recal3_train_control_test_trained_full` | 0.9184 | 0.9346 | `results/08_shift/TRNSCRPT/shift_table.csv` | `{"split": "train_control_test_trained", "arm": "full"}` | `cov_target_recal_N3` |
| `cov_train_control_test_trained_full_mondrian` | 1.0000 | 0.7191 | `results/08_shift/TRNSCRPT/shift_table.csv` | `{"split": "train_control_test_trained", "arm": "full"}` | `coverage_target_seen_mondrian` |
| `size_train_control_test_trained_full_mondrian` | 19.0000 | 0.7218 | `results/08_shift/TRNSCRPT/shift_table.csv` | `{"split": "train_control_test_trained", "arm": "full"}` | `avg_set_size_target_mondrian` |
| `cov_train_control_test_trained_full_floored` | 1.0000 | 0.9110 | `results/08_shift/TRNSCRPT/shift_table.csv` | `{"split": "train_control_test_trained", "arm": "full"}` | `coverage_target_seen_floored` |
| `size_train_control_test_trained_full_floored` | 19.0000 | 0.9221 | `results/08_shift/TRNSCRPT/shift_table.csv` | `{"split": "train_control_test_trained", "arm": "full"}` | `avg_set_size_target_floored` |
| `cov_train_male_test_female_k20_marginal` | 0.8329 | 0.8376 | `results/08_shift/TRNSCRPT/shift_table.csv` | `{"split": "train_male_test_female", "arm": "panel_k20"}` | `coverage_target_seen` |
| `size_train_male_test_female_k20_marginal` | 0.8174 | 0.8218 | `results/08_shift/TRNSCRPT/shift_table.csv` | `{"split": "train_male_test_female", "arm": "panel_k20"}` | `avg_set_size_target` |
| `covsrc_train_male_test_female_k20` | 0.9306 | 0.9444 | `results/08_shift/TRNSCRPT/shift_table.csv` | `{"split": "train_male_test_female", "arm": "panel_k20"}` | `coverage_source_id` |
| `recal3_train_male_test_female_k20` | 0.9144 | 0.9388 | `results/08_shift/TRNSCRPT/shift_table.csv` | `{"split": "train_male_test_female", "arm": "panel_k20"}` | `cov_target_recal_N3` |
| `cov_train_male_test_female_k20_mondrian` | 1.0000 | 0.7035 | `results/08_shift/TRNSCRPT/shift_table.csv` | `{"split": "train_male_test_female", "arm": "panel_k20"}` | `coverage_target_seen_mondrian` |
| `size_train_male_test_female_k20_mondrian` | 18.0000 | 0.7016 | `results/08_shift/TRNSCRPT/shift_table.csv` | `{"split": "train_male_test_female", "arm": "panel_k20"}` | `avg_set_size_target_mondrian` |
| `cov_train_male_test_female_k20_floored` | 1.0000 | 0.9247 | `results/08_shift/TRNSCRPT/shift_table.csv` | `{"split": "train_male_test_female", "arm": "panel_k20"}` | `coverage_target_seen_floored` |
| `size_train_male_test_female_k20_floored` | 18.0000 | 0.9109 | `results/08_shift/TRNSCRPT/shift_table.csv` | `{"split": "train_male_test_female", "arm": "panel_k20"}` | `avg_set_size_target_floored` |
| `cov_train_male_test_female_full_marginal` | 0.8071 | 0.8188 | `results/08_shift/TRNSCRPT/shift_table.csv` | `{"split": "train_male_test_female", "arm": "full"}` | `coverage_target_seen` |
| `size_train_male_test_female_full_marginal` | 0.7795 | 0.7929 | `results/08_shift/TRNSCRPT/shift_table.csv` | `{"split": "train_male_test_female", "arm": "full"}` | `avg_set_size_target` |
| `recal3_train_male_test_female_full` | 0.9035 | 0.9307 | `results/08_shift/TRNSCRPT/shift_table.csv` | `{"split": "train_male_test_female", "arm": "full"}` | `cov_target_recal_N3` |
| `cov_train_male_test_female_full_mondrian` | 1.0000 | 0.6988 | `results/08_shift/TRNSCRPT/shift_table.csv` | `{"split": "train_male_test_female", "arm": "full"}` | `coverage_target_seen_mondrian` |
| `size_train_male_test_female_full_mondrian` | 18.0000 | 0.6927 | `results/08_shift/TRNSCRPT/shift_table.csv` | `{"split": "train_male_test_female", "arm": "full"}` | `avg_set_size_target_mondrian` |
| `cov_train_male_test_female_full_floored` | 1.0000 | 0.8941 | `results/08_shift/TRNSCRPT/shift_table.csv` | `{"split": "train_male_test_female", "arm": "full"}` | `coverage_target_seen_floored` |
| `size_train_male_test_female_full_floored` | 18.0000 | 0.8775 | `results/08_shift/TRNSCRPT/shift_table.csv` | `{"split": "train_male_test_female", "arm": "full"}` | `avg_set_size_target_floored` |
| `covsrc_train_female_test_male_k20` | 0.9306 | 0.9583 | `results/08_shift/TRNSCRPT/shift_table.csv` | `{"split": "train_female_test_male", "arm": "panel_k20"}` | `coverage_source_id` |
| `recal3_train_female_test_male_k20` | 0.9214 | 0.9278 | `results/08_shift/TRNSCRPT/shift_table.csv` | `{"split": "train_female_test_male", "arm": "panel_k20"}` | `cov_target_recal_N3` |
| `cov_train_female_test_male_k20_mondrian` | 1.0000 | 0.7976 | `results/08_shift/TRNSCRPT/shift_table.csv` | `{"split": "train_female_test_male", "arm": "panel_k20"}` | `coverage_target_seen_mondrian` |
| `size_train_female_test_male_k20_mondrian` | 18.0000 | 0.8800 | `results/08_shift/TRNSCRPT/shift_table.csv` | `{"split": "train_female_test_male", "arm": "panel_k20"}` | `avg_set_size_target_mondrian` |
| `cov_train_female_test_male_k20_floored` | 1.0000 | 0.9153 | `results/08_shift/TRNSCRPT/shift_table.csv` | `{"split": "train_female_test_male", "arm": "panel_k20"}` | `coverage_target_seen_floored` |
| `size_train_female_test_male_k20_floored` | 18.0000 | 0.9911 | `results/08_shift/TRNSCRPT/shift_table.csv` | `{"split": "train_female_test_male", "arm": "panel_k20"}` | `avg_set_size_target_floored` |
| `cov_train_female_test_male_full_marginal` | 0.8729 | 0.8965 | `results/08_shift/TRNSCRPT/shift_table.csv` | `{"split": "train_female_test_male", "arm": "full"}` | `coverage_target_seen` |
| `size_train_female_test_male_full_marginal` | 0.8244 | 0.8467 | `results/08_shift/TRNSCRPT/shift_table.csv` | `{"split": "train_female_test_male", "arm": "full"}` | `avg_set_size_target` |
| `recal3_train_female_test_male_full` | 0.9139 | 0.9404 | `results/08_shift/TRNSCRPT/shift_table.csv` | `{"split": "train_female_test_male", "arm": "full"}` | `cov_target_recal_N3` |
| `cov_train_female_test_male_full_mondrian` | 1.0000 | 0.7976 | `results/08_shift/TRNSCRPT/shift_table.csv` | `{"split": "train_female_test_male", "arm": "full"}` | `coverage_target_seen_mondrian` |
| `size_train_female_test_male_full_mondrian` | 18.0000 | 0.8089 | `results/08_shift/TRNSCRPT/shift_table.csv` | `{"split": "train_female_test_male", "arm": "full"}` | `avg_set_size_target_mondrian` |
| `cov_train_female_test_male_full_floored` | 1.0000 | 0.9365 | `results/08_shift/TRNSCRPT/shift_table.csv` | `{"split": "train_female_test_male", "arm": "full"}` | `coverage_target_seen_floored` |
| `size_train_female_test_male_full_floored` | 18.0000 | 0.9400 | `results/08_shift/TRNSCRPT/shift_table.csv` | `{"split": "train_female_test_male", "arm": "full"}` | `avg_set_size_target_floored` |
| `recal3_bodymap_k20` | 0.9426 | 0.9696 | `results/12_bodymap/recalibration.csv` | `{"model": "k20", "n_recal": 3}` | `coverage_recalibrated` |
| `recal3size_bodymap_k20` | 0.9979 | 1.1058 | `results/12_bodymap/recalibration.csv` | `{"model": "k20", "n_recal": 3}` | `set_size_recalibrated` |
| `cov_bodymap_k20_mondrian` | 0.3382 | 0.2794 | `results/12_bodymap/conformal_transfer.csv` | `{"stage_weeks": 21, "model": "k20", "conformal": "mondrian"}` | `coverage_mapped` |
| `empty_bodymap_k20_mondrian` | 0.0000 | 0.6324 | `results/12_bodymap/conformal_transfer.csv` | `{"stage_weeks": 21, "model": "k20", "conformal": "mondrian"}` | `frac_empty_mapped` |
| `size_bodymap_k20_mondrian` | 2.3676 | 0.3676 | `results/12_bodymap/conformal_transfer.csv` | `{"stage_weeks": 21, "model": "k20", "conformal": "mondrian"}` | `avg_set_size_mapped` |
| `empty_bodymap_k20_floored` | 0.0000 | 0.2353 | `results/12_bodymap/conformal_transfer.csv` | `{"stage_weeks": 21, "model": "k20", "conformal": "floored"}` | `frac_empty_mapped` |
| `size_bodymap_k20_floored` | 2.7206 | 0.7794 | `results/12_bodymap/conformal_transfer.csv` | `{"stage_weeks": 21, "model": "k20", "conformal": "floored"}` | `avg_set_size_mapped` |
| `recal3_bodymap_k50` | 0.9365 | 0.9424 | `results/12_bodymap/recalibration.csv` | `{"model": "k50", "n_recal": 3}` | `coverage_recalibrated` |
| `recal3size_bodymap_k50` | 0.9530 | 1.1077 | `results/12_bodymap/recalibration.csv` | `{"model": "k50", "n_recal": 3}` | `set_size_recalibrated` |
| `cov_bodymap_k50_mondrian` | 0.2353 | 0.1765 | `results/12_bodymap/conformal_transfer.csv` | `{"stage_weeks": 21, "model": "k50", "conformal": "mondrian"}` | `coverage_mapped` |
| `empty_bodymap_k50_mondrian` | 0.0000 | 0.7353 | `results/12_bodymap/conformal_transfer.csv` | `{"stage_weeks": 21, "model": "k50", "conformal": "mondrian"}` | `frac_empty_mapped` |
| `size_bodymap_k50_mondrian` | 2.2647 | 0.2647 | `results/12_bodymap/conformal_transfer.csv` | `{"stage_weeks": 21, "model": "k50", "conformal": "mondrian"}` | `avg_set_size_mapped` |
| `empty_bodymap_k50_floored` | 0.0000 | 0.2059 | `results/12_bodymap/conformal_transfer.csv` | `{"stage_weeks": 21, "model": "k50", "conformal": "floored"}` | `frac_empty_mapped` |
| `size_bodymap_k50_floored` | 2.7353 | 0.7941 | `results/12_bodymap/conformal_transfer.csv` | `{"stage_weeks": 21, "model": "k50", "conformal": "floored"}` | `avg_set_size_mapped` |
| `recal3_bodymap_full` | 0.9461 | 0.9707 | `results/12_bodymap/recalibration.csv` | `{"model": "full", "n_recal": 3}` | `coverage_recalibrated` |
| `recal3size_bodymap_full` | 0.9461 | 0.9707 | `results/12_bodymap/recalibration.csv` | `{"model": "full", "n_recal": 3}` | `set_size_recalibrated` |
| `empty_bodymap_full_mondrian` | 0.0000 | 0.7206 | `results/12_bodymap/conformal_transfer.csv` | `{"stage_weeks": 21, "model": "full", "conformal": "mondrian"}` | `frac_empty_mapped` |
| `size_bodymap_full_mondrian` | 2.2206 | 0.2794 | `results/12_bodymap/conformal_transfer.csv` | `{"stage_weeks": 21, "model": "full", "conformal": "mondrian"}` | `avg_set_size_mapped` |
| `empty_bodymap_full_floored` | 0.0000 | 0.2500 | `results/12_bodymap/conformal_transfer.csv` | `{"stage_weeks": 21, "model": "full", "conformal": "floored"}` | `frac_empty_mapped` |
| `size_bodymap_full_floored` | 2.6912 | 0.7500 | `results/12_bodymap/conformal_transfer.csv` | `{"stage_weeks": 21, "model": "full", "conformal": "floored"}` | `avg_set_size_mapped` |
| `cov_gtex_k20_marginal` | 0.3638 | 0.3678 | `results/13_gtex/conformal_transfer.csv` | `{"stage": "adult", "model": "k20", "conformal": "marginal"}` | `coverage_mapped` |
| `empty_gtex_k20_marginal` | 0.6165 | 0.6113 | `results/13_gtex/conformal_transfer.csv` | `{"stage": "adult", "model": "k20", "conformal": "marginal"}` | `frac_empty_mapped` |
| `size_gtex_k20_marginal` | 0.3835 | 0.3887 | `results/13_gtex/conformal_transfer.csv` | `{"stage": "adult", "model": "k20", "conformal": "marginal"}` | `avg_set_size_mapped` |
| `recal3_gtex_k20` | 0.9544 | 0.8795 | `results/13_gtex/recalibration.csv` | `{"model": "k20", "n_recal": 3}` | `coverage_recalibrated` |
| `recal3size_gtex_k20` | 11.6966 | 4.1494 | `results/13_gtex/recalibration.csv` | `{"model": "k20", "n_recal": 3}` | `set_size_recalibrated` |
| `cov_gtex_k20_mondrian` | 0.3292 | 0.2085 | `results/13_gtex/conformal_transfer.csv` | `{"stage": "adult", "model": "k20", "conformal": "mondrian"}` | `coverage_mapped` |
| `empty_gtex_k20_mondrian` | 0.0000 | 0.6503 | `results/13_gtex/conformal_transfer.csv` | `{"stage": "adult", "model": "k20", "conformal": "mondrian"}` | `frac_empty_mapped` |
| `size_gtex_k20_mondrian` | 2.3529 | 0.3529 | `results/13_gtex/conformal_transfer.csv` | `{"stage": "adult", "model": "k20", "conformal": "mondrian"}` | `avg_set_size_mapped` |
| `cov_gtex_k20_floored` | 0.5215 | 0.4044 | `results/13_gtex/conformal_transfer.csv` | `{"stage": "adult", "model": "k20", "conformal": "floored"}` | `coverage_mapped` |
| `empty_gtex_k20_floored` | 0.0000 | 0.4523 | `results/13_gtex/conformal_transfer.csv` | `{"stage": "adult", "model": "k20", "conformal": "floored"}` | `frac_empty_mapped` |
| `size_gtex_k20_floored` | 2.5481 | 0.5525 | `results/13_gtex/conformal_transfer.csv` | `{"stage": "adult", "model": "k20", "conformal": "floored"}` | `avg_set_size_mapped` |
| `cov_gtex_k50_marginal` | 0.3630 | 0.3678 | `results/13_gtex/conformal_transfer.csv` | `{"stage": "adult", "model": "k50", "conformal": "marginal"}` | `coverage_mapped` |
| `empty_gtex_k50_marginal` | 0.6314 | 0.6266 | `results/13_gtex/conformal_transfer.csv` | `{"stage": "adult", "model": "k50", "conformal": "marginal"}` | `frac_empty_mapped` |
| `size_gtex_k50_marginal` | 0.3686 | 0.3734 | `results/13_gtex/conformal_transfer.csv` | `{"stage": "adult", "model": "k50", "conformal": "marginal"}` | `avg_set_size_mapped` |
| `recal3_gtex_k50` | 0.9592 | 0.9015 | `results/13_gtex/recalibration.csv` | `{"model": "k50", "n_recal": 3}` | `coverage_recalibrated` |
| `recal3size_gtex_k50` | 11.6244 | 3.4557 | `results/13_gtex/recalibration.csv` | `{"model": "k50", "n_recal": 3}` | `set_size_recalibrated` |
| `cov_gtex_k50_mondrian` | 0.3018 | 0.1811 | `results/13_gtex/conformal_transfer.csv` | `{"stage": "adult", "model": "k50", "conformal": "mondrian"}` | `coverage_mapped` |
| `empty_gtex_k50_mondrian` | 0.0000 | 0.7835 | `results/13_gtex/conformal_transfer.csv` | `{"stage": "adult", "model": "k50", "conformal": "mondrian"}` | `frac_empty_mapped` |
| `size_gtex_k50_mondrian` | 2.2165 | 0.2165 | `results/13_gtex/conformal_transfer.csv` | `{"stage": "adult", "model": "k50", "conformal": "mondrian"}` | `avg_set_size_mapped` |
| `cov_gtex_k50_floored` | 0.5103 | 0.3932 | `results/13_gtex/conformal_transfer.csv` | `{"stage": "adult", "model": "k50", "conformal": "floored"}` | `coverage_mapped` |
| `empty_gtex_k50_floored` | 0.0000 | 0.5682 | `results/13_gtex/conformal_transfer.csv` | `{"stage": "adult", "model": "k50", "conformal": "floored"}` | `frac_empty_mapped` |
| `size_gtex_k50_floored` | 2.4282 | 0.4318 | `results/13_gtex/conformal_transfer.csv` | `{"stage": "adult", "model": "k50", "conformal": "floored"}` | `avg_set_size_mapped` |
| `cov_gtex_full_marginal` | 0.0620 | 0.0636 | `results/13_gtex/conformal_transfer.csv` | `{"stage": "adult", "model": "full", "conformal": "marginal"}` | `coverage_mapped` |
| `empty_gtex_full_marginal` | 0.9380 | 0.9364 | `results/13_gtex/conformal_transfer.csv` | `{"stage": "adult", "model": "full", "conformal": "marginal"}` | `frac_empty_mapped` |
| `size_gtex_full_marginal` | 0.0620 | 0.0636 | `results/13_gtex/conformal_transfer.csv` | `{"stage": "adult", "model": "full", "conformal": "marginal"}` | `avg_set_size_mapped` |
| `recal3_gtex_full` | 0.9507 | 0.8801 | `results/13_gtex/recalibration.csv` | `{"model": "full", "n_recal": 3}` | `coverage_recalibrated` |
| `recal3size_gtex_full` | 11.1775 | 1.5144 | `results/13_gtex/recalibration.csv` | `{"model": "full", "n_recal": 3}` | `set_size_recalibrated` |
| `cov_gtex_full_mondrian` | 0.1984 | 0.0777 | `results/13_gtex/conformal_transfer.csv` | `{"stage": "adult", "model": "full", "conformal": "mondrian"}` | `coverage_mapped` |
| `empty_gtex_full_mondrian` | 0.0000 | 0.9203 | `results/13_gtex/conformal_transfer.csv` | `{"stage": "adult", "model": "full", "conformal": "mondrian"}` | `frac_empty_mapped` |
| `size_gtex_full_mondrian` | 2.0797 | 0.0797 | `results/13_gtex/conformal_transfer.csv` | `{"stage": "adult", "model": "full", "conformal": "mondrian"}` | `avg_set_size_mapped` |
| `cov_gtex_full_floored` | 0.2519 | 0.1328 | `results/13_gtex/conformal_transfer.csv` | `{"stage": "adult", "model": "full", "conformal": "floored"}` | `coverage_mapped` |
| `empty_gtex_full_floored` | 0.0000 | 0.8652 | `results/13_gtex/conformal_transfer.csv` | `{"stage": "adult", "model": "full", "conformal": "floored"}` | `frac_empty_mapped` |
| `size_gtex_full_floored` | 2.1332 | 0.1348 | `results/13_gtex/conformal_transfer.csv` | `{"stage": "adult", "model": "full", "conformal": "floored"}` | `avg_set_size_mapped` |
| `skmgn_cov_marginal` | 0.4200 | 0.4400 | `results/06_conformal/TRNSCRPT/per_tissue_marginal_vs_mondrian_alpha0.1.csv` | `{"y_true": "SKM-GN"}` | `coverage_marginal` |
| `skmgn_cov_mondrian` | 0.8800 | 0.9800 | `results/06_conformal/TRNSCRPT/per_tissue_marginal_vs_mondrian_alpha0.1.csv` | `{"y_true": "SKM-GN"}` | `coverage_mondrian` |

## 3. Recomputed and pending entries

| id | value | files / reason |
|---|---|---|
| `bodymap_acc_2wk_k20_ci` | [0.676056338028169, 0.8461538461538461] | `results/31_site_regen/12_bodymap/scores_target_probs.csv`: 95 % cluster bootstrap over the 8 BodyMap animals at 2 weeks (k20) |
| `bodymap_acc_2wk_k50_ci` | [0.7887323943661971, 0.8615384615384616] | `results/31_site_regen/12_bodymap/scores_target_probs.csv`: 95 % cluster bootstrap over the 8 BodyMap animals at 2 weeks (k50) |
| `bodymap_acc_2wk_full_ci` | [0.9014084507042254, 0.9846153846153847] | `results/31_site_regen/12_bodymap/scores_target_probs.csv`: 95 % cluster bootstrap over the 8 BodyMap animals at 2 weeks (full) |
| `bodymap_acc_6wk_k20_ci` | [0.9264705882352942, 1.0] | `results/31_site_regen/12_bodymap/scores_target_probs.csv`: 95 % cluster bootstrap over the 8 BodyMap animals at 6 weeks (k20) |
| `bodymap_acc_6wk_k50_ci` | [0.9264705882352942, 1.0] | `results/31_site_regen/12_bodymap/scores_target_probs.csv`: 95 % cluster bootstrap over the 8 BodyMap animals at 6 weeks (k50) |
| `bodymap_acc_6wk_full_ci` | [0.6305833524471807, 1.0] | `results/31_site_regen/12_bodymap/scores_target_probs.csv`: 95 % cluster bootstrap over the 8 BodyMap animals at 6 weeks (full) |
| `bodymap_acc_21wk_k20_ci` | [0.6305833524471807, 1.0] | `results/31_site_regen/12_bodymap/scores_target_probs.csv`: 95 % cluster bootstrap over the 8 BodyMap animals at 21 weeks (k20) |
| `bodymap_acc_21wk_k50_ci` | [0.9558823529411765, 1.0] | `results/31_site_regen/12_bodymap/scores_target_probs.csv`: 95 % cluster bootstrap over the 8 BodyMap animals at 21 weeks (k50) |
| `bodymap_acc_21wk_full_ci` | [0.6305833524471807, 1.0] | `results/31_site_regen/12_bodymap/scores_target_probs.csv`: 95 % cluster bootstrap over the 8 BodyMap animals at 21 weeks (full) |
| `bodymap_acc_104wk_k20_ci` | [0.8461538461538461, 0.967741935483871] | `results/31_site_regen/12_bodymap/scores_target_probs.csv`: 95 % cluster bootstrap over the 8 BodyMap animals at 104 weeks (k20) |
| `bodymap_acc_104wk_k50_ci` | [0.8888888888888888, 0.96875] | `results/31_site_regen/12_bodymap/scores_target_probs.csv`: 95 % cluster bootstrap over the 8 BodyMap animals at 104 weeks (k50) |
| `bodymap_acc_104wk_full_ci` | [0.8970588235294118, 0.9841269841269841] | `results/31_site_regen/12_bodymap/scores_target_probs.csv`: 95 % cluster bootstrap over the 8 BodyMap animals at 104 weeks (full) |
| `gtex_recal_infinite_draws` | table | `results/31_site_regen/13_gtex/recal_thresholds.csv`: fraction of recalibration draws whose threshold is +∞ (too few mapped samples for a finite rank), from the per-draw thresholds |
| `n_trnscrpt_animals` | 50 | `results/04_baselines/TRNSCRPT/per_fold.csv`: train + test animals of one fold |
| `n_bodymap_adult_animals` | 8 | `results/12_bodymap/recalibration.csv`: recalibration + test individuals |
| `acc06_k20` | 0.9779 | `results/31_site_regen/06_conformal/TRNSCRPT/scores_test_probs.csv`: balanced accuracy of the phase-06 k20 models (18 fit animals per fold) on their 10 test animals, mean over 5 folds |
| `acc06_k50` | 0.9821 | `results/31_site_regen/06_conformal/TRNSCRPT/scores_test_probs.csv`: balanced accuracy of the phase-06 k50 models (18 fit animals per fold) on their 10 test animals, mean over 5 folds |
| `acc06_full` | 0.9916 | `results/31_site_regen/06_conformal/TRNSCRPT/scores_test_probs.csv`: balanced accuracy of the phase-06 full models (18 fit animals per fold) on their 10 test animals, mean over 5 folds |
| `cov_id_k20_marginal_pooled` | 0.8999 | `results/31_site_regen/06_conformal/TRNSCRPT/scores_calibration.csv`, `results/31_site_regen/06_conformal/TRNSCRPT/scores_test_probs.csv`: mean over 5 folds of the coverage recomputed from the phase-06 design scores (k20, marginal, pooled) |
| `cov_id_k20_marginal_one_per_animal` | 0.9242 | `results/31_site_regen/06_conformal/TRNSCRPT/scores_calibration.csv`, `results/31_site_regen/06_conformal/TRNSCRPT/scores_test_probs.csv`: mean over 5 folds of the coverage recomputed from the phase-06 design scores (k20, marginal, one_per_animal) |
| `cov_id_k20_mondrian_pooled` | 0.9166 | `results/31_site_regen/06_conformal/TRNSCRPT/scores_calibration.csv`, `results/31_site_regen/06_conformal/TRNSCRPT/scores_test_probs.csv`: mean over 5 folds of the coverage recomputed from the phase-06 design scores (k20, mondrian, pooled) |
| `cov_id_k20_mondrian_one_per_animal` | 0.9621 | `results/31_site_regen/06_conformal/TRNSCRPT/scores_calibration.csv`, `results/31_site_regen/06_conformal/TRNSCRPT/scores_test_probs.csv`: mean over 5 folds of the coverage recomputed from the phase-06 design scores (k20, mondrian, one_per_animal) |
| `cov_id_k20_floored_pooled` | 0.9744 | `results/31_site_regen/06_conformal/TRNSCRPT/scores_calibration.csv`, `results/31_site_regen/06_conformal/TRNSCRPT/scores_test_probs.csv`: mean over 5 folds of the coverage recomputed from the phase-06 design scores (k20, floored, pooled) |
| `cov_id_k20_floored_one_per_animal` | 0.9621 | `results/31_site_regen/06_conformal/TRNSCRPT/scores_calibration.csv`, `results/31_site_regen/06_conformal/TRNSCRPT/scores_test_probs.csv`: mean over 5 folds of the coverage recomputed from the phase-06 design scores (k20, floored, one_per_animal) |
| `cov_id_k50_marginal_pooled` | 0.9010 | `results/31_site_regen/06_conformal/TRNSCRPT/scores_calibration.csv`, `results/31_site_regen/06_conformal/TRNSCRPT/scores_test_probs.csv`: mean over 5 folds of the coverage recomputed from the phase-06 design scores (k50, marginal, pooled) |
| `cov_id_k50_marginal_one_per_animal` | 0.9065 | `results/31_site_regen/06_conformal/TRNSCRPT/scores_calibration.csv`, `results/31_site_regen/06_conformal/TRNSCRPT/scores_test_probs.csv`: mean over 5 folds of the coverage recomputed from the phase-06 design scores (k50, marginal, one_per_animal) |
| `cov_id_k50_mondrian_pooled` | 0.9166 | `results/31_site_regen/06_conformal/TRNSCRPT/scores_calibration.csv`, `results/31_site_regen/06_conformal/TRNSCRPT/scores_test_probs.csv`: mean over 5 folds of the coverage recomputed from the phase-06 design scores (k50, mondrian, pooled) |
| `cov_id_k50_mondrian_one_per_animal` | 0.9566 | `results/31_site_regen/06_conformal/TRNSCRPT/scores_calibration.csv`, `results/31_site_regen/06_conformal/TRNSCRPT/scores_test_probs.csv`: mean over 5 folds of the coverage recomputed from the phase-06 design scores (k50, mondrian, one_per_animal) |
| `cov_id_k50_floored_pooled` | 0.9822 | `results/31_site_regen/06_conformal/TRNSCRPT/scores_calibration.csv`, `results/31_site_regen/06_conformal/TRNSCRPT/scores_test_probs.csv`: mean over 5 folds of the coverage recomputed from the phase-06 design scores (k50, floored, pooled) |
| `cov_id_k50_floored_one_per_animal` | 0.9566 | `results/31_site_regen/06_conformal/TRNSCRPT/scores_calibration.csv`, `results/31_site_regen/06_conformal/TRNSCRPT/scores_test_probs.csv`: mean over 5 folds of the coverage recomputed from the phase-06 design scores (k50, floored, one_per_animal) |
| `cov_id_full_mondrian_one_per_animal` | 0.9643 | `results/31_site_regen/06_conformal/TRNSCRPT/scores_calibration.csv`, `results/31_site_regen/06_conformal/TRNSCRPT/scores_test_probs.csv`: mean over 5 folds of the coverage recomputed from the phase-06 design scores (full, mondrian, one_per_animal) |
| `cov_id_full_floored_one_per_animal` | 0.9643 | `results/31_site_regen/06_conformal/TRNSCRPT/scores_calibration.csv`, `results/31_site_regen/06_conformal/TRNSCRPT/scores_test_probs.csv`: mean over 5 folds of the coverage recomputed from the phase-06 design scores (full, floored, one_per_animal) |
| `empty_seen_train_control_test_trained_k20` | 0.0807 | `results/31_site_regen/08_shift_k20/scores_target_vials.csv`: empty-set rate over the target vials of seen classes (the unseen sex-specific tissue's vials excluded) |
| `empty_seen_train_control_test_trained_k50` | 0.1043 | `results/31_site_regen/08_shift_k50/scores_target_vials.csv`: empty-set rate over the target vials of seen classes (the unseen sex-specific tissue's vials excluded) |
| `empty_seen_train_control_test_trained_full` | 0.1099 | `results/31_site_regen/08_shift_k20/scores_target_vials.csv`: empty-set rate over the target vials of seen classes (the unseen sex-specific tissue's vials excluded) |
| `empty_seen_train_male_test_female_k20` | 0.1365 | `results/31_site_regen/08_shift_k20/scores_target_vials.csv`: empty-set rate over the target vials of seen classes (the unseen sex-specific tissue's vials excluded) |
| `empty_seen_train_male_test_female_k50` | 0.1624 | `results/31_site_regen/08_shift_k50/scores_target_vials.csv`: empty-set rate over the target vials of seen classes (the unseen sex-specific tissue's vials excluded) |
| `empty_seen_train_male_test_female_full` | 0.1765 | `results/31_site_regen/08_shift_k20/scores_target_vials.csv`: empty-set rate over the target vials of seen classes (the unseen sex-specific tissue's vials excluded) |
| `empty_seen_train_female_test_male_k20` | 0.1176 | `results/31_site_regen/08_shift_k20/scores_target_vials.csv`: empty-set rate over the target vials of seen classes (the unseen sex-specific tissue's vials excluded) |
| `empty_seen_train_female_test_male_k50` | 0.0776 | `results/31_site_regen/08_shift_k50/scores_target_vials.csv`: empty-set rate over the target vials of seen classes (the unseen sex-specific tissue's vials excluded) |
| `empty_seen_train_female_test_male_full` | 0.1271 | `results/31_site_regen/08_shift_k20/scores_target_vials.csv`: empty-set rate over the target vials of seen classes (the unseen sex-specific tissue's vials excluded) |
| `acc_bodymap_k20_ci` | [0.6305833524471807, 1.0] | `results/31_site_regen/12_bodymap/scores_target_probs.csv`: 95 % cluster bootstrap over the 8 adult animals |
| `acc_bodymap_k50_ci` | [0.9558823529411765, 1.0] | `results/31_site_regen/12_bodymap/scores_target_probs.csv`: 95 % cluster bootstrap over the 8 adult animals |
| `acc_bodymap_full_ci` | [0.6305833524471807, 1.0] | `results/31_site_regen/12_bodymap/scores_target_probs.csv`: 95 % cluster bootstrap over the 8 adult animals |
| `acc_gtex_k20_ci` | [0.6364354467303392, 0.670906708249313] | `results/31_site_regen/13_gtex/scores_target_probs.csv`: 95 % cluster bootstrap over the 862 donors |
| `acc_gtex_k50_ci` | [0.7667314563741295, 0.795568797011502] | `results/31_site_regen/13_gtex/scores_target_probs.csv`: 95 % cluster bootstrap over the 862 donors |
| `acc_gtex_full_ci` | [0.8418673798425484, 0.8681280359103697] | `results/31_site_regen/13_gtex/scores_target_probs.csv`: 95 % cluster bootstrap over the 862 donors |
| `gtex_heart_k20_to_skm_frac` | 0.9067 | `results/13_gtex/confusion_k20.csv`: fraction of GTEx heart samples called either skeletal muscle class (SKM-GN + SKM-VL) by the k20 panel |
| `gtex_recal_k20_n3_frac_inf` | 0.4500 | `results/31_site_regen/13_gtex/recal_thresholds.csv`: fraction of the 20 three-donor recalibration draws (k20) whose threshold is +∞ |
| `venacv_bat_calls` | 7 | `results/31_site_regen/06_conformal/TRNSCRPT/scores_test_probs.csv`: held-out vena cava vials the 20-gene model calls brown fat (phase-06 design, all folds) |
| `venacv_bat_calls_flagged` | 7 | `results/31_site_regen/06_conformal/TRNSCRPT/scores_test_probs.csv`, `results/15_time_course/design/flagged_vials.csv`: of those, vials the consortium flagged as brown-fat contaminated |
| `venacv_vials` | 50 | `results/31_site_regen/06_conformal/TRNSCRPT/scores_test_probs.csv`: vena cava vials scored in the phase-06 design |
| `fbd_acc_k20_1w` | 0.9444 | `results/31_site_regen/08_shift_k20/scores_target_vials.csv`: accuracy of the control-fit k20 model on the 1w animals' vials (train_control_test_trained split) |
| `fbd_cov_k20_1w` | 0.8833 | `results/31_site_regen/08_shift_k20/scores_target_vials.csv`: coverage of the α = 0.10 marginal sets on the 1w animals' vials |
| `fbd_acc_k20_2w` | 0.9611 | `results/31_site_regen/08_shift_k20/scores_target_vials.csv`: accuracy of the control-fit k20 model on the 2w animals' vials (train_control_test_trained split) |
| `fbd_cov_k20_2w` | 0.8944 | `results/31_site_regen/08_shift_k20/scores_target_vials.csv`: coverage of the α = 0.10 marginal sets on the 2w animals' vials |
| `fbd_acc_k20_4w` | 0.9665 | `results/31_site_regen/08_shift_k20/scores_target_vials.csv`: accuracy of the control-fit k20 model on the 4w animals' vials (train_control_test_trained split) |
| `fbd_cov_k20_4w` | 0.9162 | `results/31_site_regen/08_shift_k20/scores_target_vials.csv`: coverage of the α = 0.10 marginal sets on the 4w animals' vials |
| `fbd_acc_k20_8w` | 0.9722 | `results/31_site_regen/08_shift_k20/scores_target_vials.csv`: accuracy of the control-fit k20 model on the 8w animals' vials (train_control_test_trained split) |
| `fbd_cov_k20_8w` | 0.9167 | `results/31_site_regen/08_shift_k20/scores_target_vials.csv`: coverage of the α = 0.10 marginal sets on the 8w animals' vials |
| `fbd_acc_k50_1w` | 0.9611 | `results/31_site_regen/08_shift_k50/scores_target_vials.csv`: accuracy of the control-fit k50 model on the 1w animals' vials (train_control_test_trained split) |
| `fbd_cov_k50_1w` | 0.8611 | `results/31_site_regen/08_shift_k50/scores_target_vials.csv`: coverage of the α = 0.10 marginal sets on the 1w animals' vials |
| `fbd_acc_k50_2w` | 0.9667 | `results/31_site_regen/08_shift_k50/scores_target_vials.csv`: accuracy of the control-fit k50 model on the 2w animals' vials (train_control_test_trained split) |
| `fbd_cov_k50_2w` | 0.8778 | `results/31_site_regen/08_shift_k50/scores_target_vials.csv`: coverage of the α = 0.10 marginal sets on the 2w animals' vials |
| `fbd_acc_k50_4w` | 0.9665 | `results/31_site_regen/08_shift_k50/scores_target_vials.csv`: accuracy of the control-fit k50 model on the 4w animals' vials (train_control_test_trained split) |
| `fbd_cov_k50_4w` | 0.8883 | `results/31_site_regen/08_shift_k50/scores_target_vials.csv`: coverage of the α = 0.10 marginal sets on the 4w animals' vials |
| `fbd_acc_k50_8w` | 0.9778 | `results/31_site_regen/08_shift_k50/scores_target_vials.csv`: accuracy of the control-fit k50 model on the 8w animals' vials (train_control_test_trained split) |
| `fbd_cov_k50_8w` | 0.8889 | `results/31_site_regen/08_shift_k50/scores_target_vials.csv`: coverage of the α = 0.10 marginal sets on the 8w animals' vials |
| `fbd_acc_full_1w` | 0.9611 | `results/31_site_regen/08_shift_k20/scores_target_vials.csv`: accuracy of the control-fit full model on the 1w animals' vials (train_control_test_trained split) |
| `fbd_cov_full_1w` | 0.8667 | `results/31_site_regen/08_shift_k20/scores_target_vials.csv`: coverage of the α = 0.10 marginal sets on the 1w animals' vials |
| `fbd_acc_full_2w` | 0.9667 | `results/31_site_regen/08_shift_k20/scores_target_vials.csv`: accuracy of the control-fit full model on the 2w animals' vials (train_control_test_trained split) |
| `fbd_cov_full_2w` | 0.8778 | `results/31_site_regen/08_shift_k20/scores_target_vials.csv`: coverage of the α = 0.10 marginal sets on the 2w animals' vials |
| `fbd_acc_full_4w` | 0.9888 | `results/31_site_regen/08_shift_k20/scores_target_vials.csv`: accuracy of the control-fit full model on the 4w animals' vials (train_control_test_trained split) |
| `fbd_cov_full_4w` | 0.8883 | `results/31_site_regen/08_shift_k20/scores_target_vials.csv`: coverage of the α = 0.10 marginal sets on the 4w animals' vials |
| `fbd_acc_full_8w` | 1.0000 | `results/31_site_regen/08_shift_k20/scores_target_vials.csv`: accuracy of the control-fit full model on the 8w animals' vials (train_control_test_trained split) |
| `fbd_cov_full_8w` | 0.8889 | `results/31_site_regen/08_shift_k20/scores_target_vials.csv`: coverage of the α = 0.10 marginal sets on the 8w animals' vials |
