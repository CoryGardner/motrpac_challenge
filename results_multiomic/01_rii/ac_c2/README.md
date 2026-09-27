# Phase 1 — the RII rescue (prot-ac, release c2.0)

Built by `scripts/multiomic/01_rii_rescue.py` on 2026-09-27 07:44 UTC; every number below is read from a CSV in this directory.

## Data
- Peptide-level reporter-ion intensities from the portal quant-id `prot-ac` folders, 7 tissues, 420 vials, 60 animals (`join_summary.csv`, `tissue_summary.csv`, `plex_summary.csv`).
- Per plex: reference channel dropped, contaminants dropped, peptides summed per protein, proteins with < 2 quantified peptides dropped, each channel normalised to its total (log2 ppm).
- Union 2538 proteins; **114 quantified in every tissue** (the main matrix; NaN fraction 0.097, 35 with no missing value at all); outer NaN fraction 0.768.
- **Plex is nested in tissue here too** (one plex = 10 samples of one tissue + that tissue's reference pool; `plex_id` has one tissue per level), so everything below is within-study evidence.

## 1. Variance partition (PCA on the stacked matrix, phase-03 code) — `variance_partition.csv`, `variance_partition_ratio.csv`

| PC | explained_RII | R2_tissue_RII | R2_tissue_null95_RII | R2_plex_id_RII | R2_sex_RII | R2_pid_RII | R2_tissue_ratio_same_code | R2_tissue_ratio_frozen_phase03 |
|---|---|---|---|---|---|---|---|---|
| PC1 | 0.2915 | 0.7917 | 0.0332 | 0.9499 | 0.0075 | 0.0376 | 0.0009 | 0.0008 |
| PC2 | 0.1717 | 0.7480 | 0.0267 | 0.9637 | 0.0022 | 0.0388 | 0.0039 | 0.0039 |
| PC3 | 0.1266 | 0.8568 | 0.0299 | 0.9715 | 0.0039 | 0.0234 | 0.2184 | 0.2078 |

Pre-registration (a) asks for tissue R² of PC1 > 0.5: observed **0.792** → PASS. Median-normalised sensitivity: PC1 R² 0.979 (`variance_partition_median_norm.csv`).
R2_plex_id equals R2_tissue up to the within-tissue plex split because plex is nested in tissue; it cannot be separated here.

## 2. Phase-04 diagnostic on every animal-grouped fold — `diagnostic_accuracy.csv`, `diagnostic_accuracy_summary.csv`

| quantity | mean | sd | n_folds | n_test_animals_mean | n_missingness_features_outer | n_missingness_features_inner | frozen_ratio_fold0 |
|---|---|---|---|---|---|---|---|
| model_as_fitted_acc | 1.000 | 0.000 | 5 | 12.000 | 2503 | 79 | 0.976 |
| model_as_fitted_bal_acc | 1.000 | 0.000 | 5 | 12.000 | 2503 | 79 |  |
| missingness_outer_acc | 1.000 | 0.000 | 5 | 12.000 | 2503 | 79 | 1.000 |
| missingness_outer_bal_acc | 1.000 | 0.000 | 5 | 12.000 | 2503 | 79 |  |
| missingness_inner_acc | 0.952 | 0.000 | 5 | 12.000 | 2503 | 79 |  |
| missingness_inner_bal_acc | 0.952 | 0.000 | 5 | 12.000 | 2503 | 79 |  |
| per_tissue_means_removed_acc | 0.143 | 0.000 | 5 | 12.000 | 2503 | 79 | 0.193 |
| per_tissue_means_removed_bal_acc | 0.143 | 0.000 | 5 | 12.000 | 2503 | 79 |  |
| chance | 0.143 | 0.000 | 5 | 12.000 | 2503 | 79 | 0.143 |

Read: `model_as_fitted` is the pipeline's logreg_l2 on the RII matrix (5 folds, 12.0 test animals per fold). `missingness_outer` classifies tissue from the NaN pattern of the union matrix alone (which proteins were quantified in which tissue — an artefact of per-tissue searches, still present on the RII scale); `missingness_inner` does the same on the every-tissue matrix (within-tissue plex gaps only). `per_tissue_means_removed` erases every protein's tissue mean using the labels: on ratios it collapsed to chance because the classifier lived on normalisation offsets; on RII the tissue means ARE the fingerprint, so the collapse is expected — the difference between the two matrices is what the means are (section 1 and 4), not whether the classifier uses them.

## 3. Protein panel curve (RoundRobinSelector → logreg_l2, animal-grouped folds, pipeline unchanged) — `panel_curve_summary.csv`, `panel_null_k20.csv`

| k | bal_acc_mean | bal_acc_sd | acc_mean | macro_f1_mean | n_folds | n_test_animals_mean | null_q95_k20 |
|---|---|---|---|---|---|---|---|
| 1 | 0.362 | 0.022 | 0.362 | 0.334 | 5 | 12.000 |  |
| 2 | 0.750 | 0.039 | 0.750 | 0.731 | 5 | 12.000 |  |
| 3 | 0.845 | 0.024 | 0.845 | 0.842 | 5 | 12.000 |  |
| 5 | 0.914 | 0.028 | 0.914 | 0.913 | 5 | 12.000 |  |
| 8 | 0.998 | 0.005 | 0.998 | 0.998 | 5 | 12.000 |  |
| 10 | 1.000 | 0.000 | 1.000 | 1.000 | 5 | 12.000 |  |
| 15 | 0.998 | 0.005 | 0.998 | 0.998 | 5 | 12.000 |  |
| 20 | 0.998 | 0.005 | 0.998 | 0.998 | 5 | 12.000 | 0.178 |
| 30 | 1.000 | 0.000 | 1.000 | 1.000 | 5 | 12.000 |  |
| 50 | 1.000 | 0.000 | 1.000 | 1.000 | 5 | 12.000 |  |
| 100 | 1.000 | 0.000 | 1.000 | 1.000 | 5 | 12.000 |  |

Label-permutation null at k = 20 (10 permutations, quick fits): mean 0.150, 95th percentile 0.178, chance 0.143; observed 0.998.
Pre-registration (b) asks for mean balanced accuracy ≥ 0.95 at some k ≤ 100: best 1.000 at k = 10 → PASS (within-study; context, not a finding).

The k = 20 panel selected on all animals (`panel_k20_all_animals.csv`):

| protein_id | gene_symbol | ensembl_gene | marker_tissue | mean_z_in_marker | next_highest_tissue | effect_z | frac_nan |
|---|---|---|---|---|---|---|---|
| NP_001012183.1 | Cryz | ENSRNOG00000028319 | KIDNEY | 1.68 | LIVER | 1.22 | 0.29 |
| NP_001014157.1 | Wdr1 | ENSRNOG00000028498 | CORTEX | 1.17 | LUNG | 0.61 | 0.05 |
| NP_001014183.1 | Idh2 | ENSRNOG00000013949 | HEART | 1.76 | SKM-GN | 1.25 | 0.00 |
| NP_001028853.1 | Eif5a | ENSRNOG00000016478 | SKM-GN | 1.52 | LUNG | 1.14 | 0.00 |
| NP_001129343.1 | Ppa2 | ENSRNOG00000012091 | KIDNEY | 1.47 | HEART | 1.18 | 0.10 |
| NP_001170776.1 | Aldoa | ENSRNOG00000052802 | SKM-GN | 1.54 | CORTEX | 0.80 | 0.00 |
| NP_036702.1 | Glud1 | ENSRNOG00000057367 | LIVER | 1.74 | CORTEX | 0.89 | 0.00 |
| NP_058797.1 | Ppia | ENSRNOG00000027864 | LUNG | 1.26 | CORTEX | 0.59 | 0.07 |
| NP_071609.2 | Aldh9a1 | ENSRNOG00000004027 | LIVER | 1.81 | KIDNEY | 1.50 | 0.13 |
| NP_074042.1 | Aldh5a1 | ENSRNOG00000023538 | CORTEX | 2.19 | SKM-GN | 2.19 | 0.24 |
| NP_075211.2 | Tpi1 | ENSRNOG00000015290 | SKM-GN | 1.77 | CORTEX | 1.20 | 0.00 |
| NP_112319.2 | Aldh6a1 | ENSRNOG00000011419 | KIDNEY | 1.29 | LIVER | 0.28 | 0.05 |
| NP_112644.1 | Vdac2 | ENSRNOG00000013505 | CORTEX | 1.74 | LUNG | 1.51 | 0.21 |
| NP_113791.2 | Ywhae | nan | CORTEX | 1.23 | LUNG | 0.43 | 0.10 |
| NP_445743.2 | Pgk1 | ENSRNOG00000058249 | SKM-GN | 1.90 | CORTEX | 1.35 | 0.00 |
| NP_446028.1 | Prdx6 | ENSRNOG00000002896 | LUNG | 1.62 | HEART | 1.47 | 0.03 |
| NP_599191.1 | Atp5f1b | ENSRNOG00000002840 | CORTEX | 0.90 | LIVER | 0.04 | 0.00 |
| XP_006253736.1 | Auh | ENSRNOG00000011684 | CORTEX | 1.73 | HEART | 1.14 | 0.14 |
| XP_006254040.1 | LOC102549061 | nan | LUNG | 1.70 | CORTEX | 1.54 | 0.07 |
| NP_001013128.1 | Tf | ENSRNOG00000030625 | LUNG | 1.31 | WAT-SC | 0.29 | 0.10 |

**Missingness-free check** (`variance_partition_complete.csv`, `diagnostic_accuracy_complete.csv`, `panel_curve_complete_summary.csv`): on the 35 proteins quantified in every vial (no NaN, no imputation), tissue R² of PC1 = 0.795, PC2 = 0.922; logreg_l2 accuracy 1.000 ± 0.000, means-removed 0.143; panel k = 20 balanced accuracy 1.000 ± 0.000.

Confused pairs at k = 20 (pooled over folds; 1 errors):

| k | true | predicted | count | frac_of_true |
|---|---|---|---|---|
| 20 | LUNG | WAT-SC | 1 | 0.02 |

Confused pairs at k = 30 (pooled over folds; 0 errors):

_none_

## What this phase does not show

- Nothing about transfer: one plex holds one tissue, so tissue and plex are confounded exactly as in the audit; the tissue axis is real on this scale but its size cannot be separated from a per-plex processing offset without an external dataset (Phases 3 and 6).
- Within-study accuracy is a ceiling task and is reported as context.
- The normalisation (channel total over the kept proteins) makes every value relative to that channel's quantified proteome; a tissue that quantifies fewer proteins has larger shares. The median-normalised sensitivity (`variance_partition_median_norm.csv`) checks that the PC1 result does not depend on this choice.

_Run time 0.7 min._