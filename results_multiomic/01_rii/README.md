# Phase 1 — the RII rescue (prot-pr, release c1.0)

Built by `scripts/multiomic/01_rii_rescue.py` on 2026-09-27 07:05 UTC; every number below is read from a CSV in this directory.

## Data
- Peptide-level reporter-ion intensities from the portal quant-id `prot-pr` folders, 7 tissues, 420 vials, 60 animals (`join_summary.csv`, `tissue_summary.csv`, `plex_summary.csv`).
- Per plex: reference channel dropped, contaminants dropped, peptides summed per protein, proteins with < 2 quantified peptides dropped, each channel normalised to its total (log2 ppm).
- Union 17396 proteins; **3637 quantified in every tissue** (the main matrix; NaN fraction 0.035, 2393 with no missing value at all); outer NaN fraction 0.548.
- **Plex is nested in tissue here too** (one plex = 10 samples of one tissue + that tissue's reference pool; `plex_id` has one tissue per level), so everything below is within-study evidence.

## 1. Variance partition (PCA on the stacked matrix, phase-03 code) — `variance_partition.csv`, `variance_partition_ratio.csv`

| PC | explained_RII | R2_tissue_RII | R2_tissue_null95_RII | R2_plex_id_RII | R2_sex_RII | R2_pid_RII | R2_tissue_ratio_same_code | R2_tissue_ratio_frozen_phase03 |
|---|---|---|---|---|---|---|---|---|
| PC1 | 0.4501 | 0.9910 | 0.0294 | 0.9949 | 0.0004 | 0.0014 | 0.0009 | 0.0008 |
| PC2 | 0.1159 | 0.9983 | 0.0276 | 0.9993 | 0.0000 | 0.0003 | 0.0039 | 0.0039 |
| PC3 | 0.1125 | 0.9931 | 0.0276 | 0.9961 | 0.0002 | 0.0011 | 0.2184 | 0.2078 |

Pre-registration (a) asks for tissue R² of PC1 > 0.5: observed **0.991** → PASS. Median-normalised sensitivity: PC1 R² 0.997 (`variance_partition_median_norm.csv`).
R2_plex_id equals R2_tissue up to the within-tissue plex split because plex is nested in tissue; it cannot be separated here.

## 2. Phase-04 diagnostic on every animal-grouped fold — `diagnostic_accuracy.csv`, `diagnostic_accuracy_summary.csv`

| quantity | mean | sd | n_folds | n_test_animals_mean | n_missingness_features_outer | n_missingness_features_inner | frozen_ratio_fold0 |
|---|---|---|---|---|---|---|---|
| model_as_fitted_acc | 1.000 | 0.000 | 5 | 12.000 | 15003 | 1244 | 0.976 |
| model_as_fitted_bal_acc | 1.000 | 0.000 | 5 | 12.000 | 15003 | 1244 |  |
| missingness_outer_acc | 1.000 | 0.000 | 5 | 12.000 | 15003 | 1244 | 1.000 |
| missingness_outer_bal_acc | 1.000 | 0.000 | 5 | 12.000 | 15003 | 1244 |  |
| missingness_inner_acc | 1.000 | 0.000 | 5 | 12.000 | 15003 | 1244 |  |
| missingness_inner_bal_acc | 1.000 | 0.000 | 5 | 12.000 | 15003 | 1244 |  |
| per_tissue_means_removed_acc | 0.143 | 0.000 | 5 | 12.000 | 15003 | 1244 | 0.193 |
| per_tissue_means_removed_bal_acc | 0.143 | 0.000 | 5 | 12.000 | 15003 | 1244 |  |
| chance | 0.143 | 0.000 | 5 | 12.000 | 15003 | 1244 | 0.143 |

Read: `model_as_fitted` is the pipeline's logreg_l2 on the RII matrix (5 folds, 12.0 test animals per fold). `missingness_outer` classifies tissue from the NaN pattern of the union matrix alone (which proteins were quantified in which tissue — an artefact of per-tissue searches, still present on the RII scale); `missingness_inner` does the same on the every-tissue matrix (within-tissue plex gaps only). `per_tissue_means_removed` erases every protein's tissue mean using the labels: on ratios it collapsed to chance because the classifier lived on normalisation offsets; on RII the tissue means ARE the fingerprint, so the collapse is expected — the difference between the two matrices is what the means are (section 1 and 4), not whether the classifier uses them.

## 3. Protein panel curve (RoundRobinSelector → logreg_l2, animal-grouped folds, pipeline unchanged) — `panel_curve_summary.csv`, `panel_null_k20.csv`

| k | bal_acc_mean | bal_acc_sd | acc_mean | macro_f1_mean | n_folds | n_test_animals_mean | null_q95_k20 |
|---|---|---|---|---|---|---|---|
| 1 | 0.390 | 0.042 | 0.390 | 0.326 | 5 | 12.000 |  |
| 2 | 0.855 | 0.016 | 0.855 | 0.852 | 5 | 12.000 |  |
| 3 | 0.871 | 0.051 | 0.871 | 0.871 | 5 | 12.000 |  |
| 5 | 0.993 | 0.011 | 0.993 | 0.993 | 5 | 12.000 |  |
| 8 | 1.000 | 0.000 | 1.000 | 1.000 | 5 | 12.000 |  |
| 10 | 1.000 | 0.000 | 1.000 | 1.000 | 5 | 12.000 |  |
| 15 | 1.000 | 0.000 | 1.000 | 1.000 | 5 | 12.000 |  |
| 20 | 1.000 | 0.000 | 1.000 | 1.000 | 5 | 12.000 | 0.156 |
| 30 | 1.000 | 0.000 | 1.000 | 1.000 | 5 | 12.000 |  |
| 50 | 1.000 | 0.000 | 1.000 | 1.000 | 5 | 12.000 |  |
| 100 | 1.000 | 0.000 | 1.000 | 1.000 | 5 | 12.000 |  |

Label-permutation null at k = 20 (20 permutations, quick fits): mean 0.135, 95th percentile 0.156, chance 0.143; observed 1.000.
Pre-registration (b) asks for mean balanced accuracy ≥ 0.95 at some k ≤ 100: best 1.000 at k = 8 → PASS (within-study; context, not a finding).

The k = 20 panel selected on all animals (`panel_k20_all_animals.csv`):

| protein_id | gene_symbol | ensembl_gene | marker_tissue | mean_z_in_marker | next_highest_tissue | effect_z | frac_nan |
|---|---|---|---|---|---|---|---|
| NP_001007146.1 | Ctnna1 | ENSRNOG00000005796 | LUNG | 0.82 | LIVER | 0.11 | 0.00 |
| NP_001012183.1 | Cryz | ENSRNOG00000028319 | KIDNEY | 2.14 | LIVER | 1.88 | 0.00 |
| NP_001100863.2 | Diaph1 | ENSRNOG00000019688 | LUNG | 1.11 | KIDNEY | 0.47 | 0.00 |
| NP_001100900.1 | Pdpr | ENSRNOG00000022593 | CORTEX | 0.71 | LUNG | 0.16 | 0.00 |
| NP_001101217.2 | Arhgap1 | ENSRNOG00000016610 | LUNG | 1.02 | LIVER | 0.58 | 0.00 |
| NP_001101966.2 | Stim1 | ENSRNOG00000020425 | SKM-GN | 0.87 | CORTEX | 0.21 | 0.00 |
| NP_001103370.1 | Akr1c15 | ENSRNOG00000021735 | LUNG | 2.07 | WAT-SC | 1.43 | 0.00 |
| NP_001121052.1 | Oxct1 | ENSRNOG00000043094 | KIDNEY | 1.17 | HEART | 0.69 | 0.00 |
| NP_001124020.1 | Col14a1 | ENSRNOG00000026415 | WAT-SC | 2.21 | KIDNEY | 2.06 | 0.00 |
| NP_001128232.1 | Hnrnpl | ENSRNOG00000020235 | LIVER | 0.88 | LUNG | 0.08 | 0.00 |
| NP_001128235.1 | Hip1r | ENSRNOG00000001091 | CORTEX | 1.00 | LIVER | 0.55 | 0.00 |
| NP_058771.2 | Acat1 | ENSRNOG00000007862 | KIDNEY | 2.26 | CORTEX | 2.44 | 0.00 |
| NP_058933.2 | Uchl1 | ENSRNOG00000002343 | CORTEX | 2.23 | LIVER | 2.12 | 0.02 |
| NP_071560.2 | Pter | ENSRNOG00000017328 | KIDNEY | 2.05 | WAT-SC | 1.81 | 0.10 |
| NP_445991.2 | Idi1 | ENSRNOG00000016690 | CORTEX | 0.54 | LIVER | 0.01 | 0.12 |
| NP_446039.1 | S100a9 | ENSRNOG00000011483 | LUNG | 2.09 | LIVER | 1.61 | 0.05 |
| NP_446341.1 | Vwf | ENSRNOG00000019689 | LUNG | 2.05 | WAT-SC | 1.46 | 0.00 |
| NP_542420.1 | Dnm1 | ENSRNOG00000033835 | CORTEX | 2.20 | LUNG | 1.91 | 0.00 |
| NP_599229.1 | Bzw2 | ENSRNOG00000005096 | CORTEX | 0.92 | SKM-GN | 0.50 | 0.00 |
| NP_640347.2 | Tubb3 | ENSRNOG00000017209 | CORTEX | 2.37 | WAT-SC | 2.55 | 0.21 |

**Missingness-free check** (`variance_partition_complete.csv`, `diagnostic_accuracy_complete.csv`, `panel_curve_complete_summary.csv`): on the 2393 proteins quantified in every vial (no NaN, no imputation), tissue R² of PC1 = 0.990, PC2 = 0.996; logreg_l2 accuracy 1.000 ± 0.000, means-removed 0.143; panel k = 20 balanced accuracy 1.000 ± 0.000.

Confused pairs at k = 20 (pooled over folds; 0 errors):

_none_

Confused pairs at k = 30 (pooled over folds; 0 errors):

_none_

## 4. Cross-tissue RNA–protein correlation on the same animals — `rna_protein_correlation.csv`, `rna_protein_correlation_summary.csv`

Tissue means over the animals shared by RNA-seq and proteomics in each tissue (CORTEX:50;SKM-GN:50;HEART:50;KIDNEY:50;LUNG:50;LIVER:50;WAT-SC:50); Spearman across the 7 tissues per gene (protein means averaged over the proteins mapping to the gene; a tissue mean needs ≥ 3 quantified animals). n = 3701 genes.
- Spearman median **0.679** (IQR 0.429–0.821); 68.4 % of genes > 0.5; 7.0 % < 0.
- Mismatched-pair null (RNA of one gene vs the protein of another, 19994 pairs): median 0.286, 95th percentile 0.857; 16.6 % of genes exceed the null 95th percentile.
- Same marker tissue (highest tissue mean) at RNA and protein: 44.9 % of genes (chance 1/7 = 14.3 %).
- By protein cross-tissue range: median Spearman 0.429 (Q1) → 0.821 (Q4) (`rna_protein_correlation_by_protein_range.csv`).

### RNA panel genes at the protein level — `rna_protein_panel_genes.csv`, `rna_protein_panel_summary.csv`

51 RNA panel genes (stable core, stability-51, k=20 all-animal panel); 33 have a quantified protein, 2 in all 7 tissues; 23 have their RNA marker tissue among the 7 proteomics tissues.
- Pre-registration (c): of the **20** testable genes (marker tissue among the 7, protein quantified), the protein's highest tissue mean is the RNA marker tissue in **17** (85 %) → PASS (threshold 70 %).
- Median Spearman over panel genes with all-7 means 0.500 vs 0.679 over all genes; 0 panel genes above the null 95th percentile.

| gene_symbol | lists | marker_tissue_19 | protein_quantified_any_tissue | prot_marker_tissue | rna_marker_tissue | spearman | testable_c | protein_marker_equals_rna_marker19 |
|---|---|---|---|---|---|---|---|---|
| Cfhr1 | stability_k20 | LIVER | True | KIDNEY | LIVER | 0.57 | True | False |
| F9 | stability_k20 | LIVER | True | KIDNEY | LIVER | 0.43 | True | False |
| Gabrd | panel_k20_all_animals;stability_k20;stable_core_T5 | CORTEX | True | CORTEX | CORTEX |  | True | True |
| Lrrc10 | panel_k20_all_animals;stability_k20;stable_core_T5 | HEART | True | HEART | HEART |  | True | True |
| Umod | panel_k20_all_animals;stability_k20;stable_core_T5 | KIDNEY | True | KIDNEY | KIDNEY |  | True | True |
| Cldn18 | panel_k20_all_animals;stability_k20;stable_core_T5 | LUNG | True | LUNG | LUNG |  | True | True |
| Ephx4 | stability_k20 | CORTEX | True | CORTEX | CORTEX |  | True | True |
| Tbr1 | stability_k20 | CORTEX | True | CORTEX | CORTEX |  | True | True |
| Gnb3 | stability_k20 | HEART | True | HEART | HEART |  | True | True |
| Cfhr2 | panel_k20_all_animals;stability_k20 | LIVER | True | LIVER | LIVER |  | True | True |
| Cpn2 | stability_k20 | LIVER | True | LUNG | LIVER |  | True | False |
| Mbl2 | stability_k20 | LIVER | True | LIVER | LIVER |  | True | True |
| Slco1b2 | stability_k20 | LIVER | True | LIVER | LIVER |  | True | True |
| Sftpb | stability_k20 | LUNG | True | LUNG | LUNG |  | True | True |
| Ankrd2 | stability_k20 | SKM-GN | True | SKM-GN | SKM-GN |  | True | True |
| Mss51 | panel_k20_all_animals;stability_k20 | SKM-GN | True | SKM-GN | SKM-GN |  | True | True |
| RGD1565323 | stability_k20 | SKM-GN | True | SKM-GN | SKM-GN |  | True | True |
| Tnni1 | stability_k20 | SKM-GN | True | SKM-GN | SKM-GN |  | True | True |
| Trim29 | panel_k20_all_animals;stability_k20 | WAT-SC | True | WAT-SC | WAT-SC |  | True | True |
| Csn3 | stability_k20 | WAT-SC | True | WAT-SC | WAT-SC |  | True | True |
| Cyp21a1 | panel_k20_all_animals;stability_k20;stable_core_T5 | ADRNL | False | nan | nan |  | False | False |
| Fibcd1l1 | panel_k20_all_animals;stability_k20;stable_core_T5 | HIPPOC | True | CORTEX | CORTEX |  | False | False |
| Pmch | panel_k20_all_animals;stability_k20;stable_core_T5 | HYPOTH | True | CORTEX | CORTEX |  | False | False |
| Akr1c3 | panel_k20_all_animals;stability_k20;stable_core_T5 | OVARY | True | LIVER | WAT-SC |  | False | False |
| Ccl25 | panel_k20_all_animals;stability_k20;stable_core_T5 | SMLINT | False | nan | nan |  | False | False |
| Fcrl5 | panel_k20_all_animals;stability_k20;stable_core_T5 | SPLEEN | False | nan | nan |  | False | False |
| Qrfpr | panel_k20_all_animals;stability_k20 | ADRNL | False | nan | nan |  | False | False |
| Dbh | stability_k20 | ADRNL | True | CORTEX | HEART |  | False | False |
| Mc2r | stability_k20 | ADRNL | False | nan | nan |  | False | False |
| Otop1 | panel_k20_all_animals;stability_k20 | BAT | False | nan | nan |  | False | False |
| Ucp1 | stability_k20 | BAT | True | LUNG | KIDNEY |  | False | False |
| Grifin | stability_k20 | BAT | True | WAT-SC | WAT-SC |  | False | False |
| Ebf2 | stability_k20 | BAT | True | WAT-SC | WAT-SC |  | False | False |
| Hbq1b | panel_k20_all_animals;stability_k20 | BLOOD | True | LUNG | nan |  | False | False |
| Treml1 | stability_k20 | BLOOD | True | LUNG | LUNG |  | False | False |
| Hbb-b1 | stability_k20 | BLOOD | True | WAT-SC | LUNG |  | False | False |
| Abo3 | panel_k20_all_animals;stability_k20 | COLON | False | nan | nan |  | False | False |
| Hoxb13 | stability_k20 | COLON | False | nan | nan |  | False | False |
| St6galnac1 | stability_k20 | COLON | False | nan | nan |  | False | False |
| Slc17a6 | stability_k20 | HYPOTH | True | CORTEX | CORTEX |  | False | False |
| AABR07069371.1 | stability_k20 | LUNG | False | nan | nan |  | False | False |
| Lbx1 | stability_k20 | SKM-GN | False | nan | nan |  | False | False |
| Mybph | panel_k20_all_animals;stability_k20 | SKM-VL | True | SKM-GN | SKM-GN |  | False | False |
| Ccr3 | stability_k20 | SPLEEN | False | nan | nan |  | False | False |
| Pgk2 | panel_k20_all_animals;stability_k20 | TESTES | True | SKM-GN | nan |  | False | False |
| LOC100271845 | stability_k20 | TESTES | False | nan | nan |  | False | False |
| Ubqln3 | stability_k20 | TESTES | False | nan | nan |  | False | False |
| Ubqlnl | stability_k20 | TESTES | False | nan | nan |  | False | False |
| Smcp | stability_k20 | TESTES | False | nan | nan |  | False | False |
| Gdf10 | panel_k20_all_animals;stability_k20 | VENACV | False | nan | nan |  | False | False |
| Snorc | stability_k20 | WAT-SC | False | nan | nan |  | False | False |

## What this phase does not show

- Nothing about transfer: one plex holds one tissue, so tissue and plex are confounded exactly as in the audit; the tissue axis is real on this scale but its size cannot be separated from a per-plex processing offset without an external dataset (Phases 3 and 6).
- Within-study accuracy is a ceiling task and is reported as context.
- The normalisation (channel total over the kept proteins) makes every value relative to that channel's quantified proteome; a tissue that quantifies fewer proteins has larger shares. The median-normalised sensitivity (`variance_partition_median_norm.csv`) checks that the PC1 result does not depend on this choice.

_Run time 1.9 min._