# Phase 1 — the RII rescue (prot-ph, release c1.0)

Built by `scripts/multiomic/01_rii_rescue.py` on 2026-09-27 07:06 UTC; every number below is read from a CSV in this directory.

## Data
- Peptide-level reporter-ion intensities from the portal quant-id `prot-ph` folders, 7 tissues, 420 vials, 60 animals (`join_summary.csv`, `tissue_summary.csv`, `plex_summary.csv`).
- Per plex: reference channel dropped, contaminants dropped, peptides summed per protein, proteins with < 2 quantified peptides dropped, each channel normalised to its total (log2 ppm).
- Union 11762 proteins; **1236 quantified in every tissue** (the main matrix; NaN fraction 0.099, 417 with no missing value at all); outer NaN fraction 0.674.
- **Plex is nested in tissue here too** (one plex = 10 samples of one tissue + that tissue's reference pool; `plex_id` has one tissue per level), so everything below is within-study evidence.

## 1. Variance partition (PCA on the stacked matrix, phase-03 code) — `variance_partition.csv`, `variance_partition_ratio.csv`

| PC | explained_RII | R2_tissue_RII | R2_tissue_null95_RII | R2_plex_id_RII | R2_sex_RII | R2_pid_RII | R2_tissue_ratio_same_code | R2_tissue_ratio_frozen_phase03 |
|---|---|---|---|---|---|---|---|---|
| PC1 | 0.3412 | 0.9864 | 0.0257 | 0.9934 | 0.0006 | 0.0021 | 0.0009 | 0.0008 |
| PC2 | 0.1079 | 0.9947 | 0.0264 | 0.9983 | 0.0001 | 0.0005 | 0.0039 | 0.0039 |
| PC3 | 0.0887 | 0.9950 | 0.0277 | 0.9993 | 0.0000 | 0.0003 | 0.2184 | 0.2078 |

Pre-registration (a) asks for tissue R² of PC1 > 0.5: observed **0.986** → PASS. Median-normalised sensitivity: PC1 R² 0.989 (`variance_partition_median_norm.csv`).
R2_plex_id equals R2_tissue up to the within-tissue plex split because plex is nested in tissue; it cannot be separated here.

## 2. Phase-04 diagnostic on every animal-grouped fold — `diagnostic_accuracy.csv`, `diagnostic_accuracy_summary.csv`

| quantity | mean | sd | n_folds | n_test_animals_mean | n_missingness_features_outer | n_missingness_features_inner | frozen_ratio_fold0 |
|---|---|---|---|---|---|---|---|
| model_as_fitted_acc | 1.000 | 0.000 | 5 | 12.000 | 11345 | 819 | 0.976 |
| model_as_fitted_bal_acc | 1.000 | 0.000 | 5 | 12.000 | 11345 | 819 |  |
| missingness_outer_acc | 1.000 | 0.000 | 5 | 12.000 | 11345 | 819 | 1.000 |
| missingness_outer_bal_acc | 1.000 | 0.000 | 5 | 12.000 | 11345 | 819 |  |
| missingness_inner_acc | 1.000 | 0.000 | 5 | 12.000 | 11345 | 819 |  |
| missingness_inner_bal_acc | 1.000 | 0.000 | 5 | 12.000 | 11345 | 819 |  |
| per_tissue_means_removed_acc | 0.143 | 0.000 | 5 | 12.000 | 11345 | 819 | 0.193 |
| per_tissue_means_removed_bal_acc | 0.143 | 0.000 | 5 | 12.000 | 11345 | 819 |  |
| chance | 0.143 | 0.000 | 5 | 12.000 | 11345 | 819 | 0.143 |

Read: `model_as_fitted` is the pipeline's logreg_l2 on the RII matrix (5 folds, 12.0 test animals per fold). `missingness_outer` classifies tissue from the NaN pattern of the union matrix alone (which proteins were quantified in which tissue — an artefact of per-tissue searches, still present on the RII scale); `missingness_inner` does the same on the every-tissue matrix (within-tissue plex gaps only). `per_tissue_means_removed` erases every protein's tissue mean using the labels: on ratios it collapsed to chance because the classifier lived on normalisation offsets; on RII the tissue means ARE the fingerprint, so the collapse is expected — the difference between the two matrices is what the means are (section 1 and 4), not whether the classifier uses them.

## 3. Protein panel curve (RoundRobinSelector → logreg_l2, animal-grouped folds, pipeline unchanged) — `panel_curve_summary.csv`, `panel_null_k20.csv`

| k | bal_acc_mean | bal_acc_sd | acc_mean | macro_f1_mean | n_folds | n_test_animals_mean | null_q95_k20 |
|---|---|---|---|---|---|---|---|
| 1 | 0.436 | 0.022 | 0.436 | 0.405 | 5 | 12.000 |  |
| 2 | 0.752 | 0.032 | 0.752 | 0.733 | 5 | 12.000 |  |
| 3 | 0.840 | 0.034 | 0.840 | 0.839 | 5 | 12.000 |  |
| 5 | 0.995 | 0.007 | 0.995 | 0.995 | 5 | 12.000 |  |
| 8 | 1.000 | 0.000 | 1.000 | 1.000 | 5 | 12.000 |  |
| 10 | 1.000 | 0.000 | 1.000 | 1.000 | 5 | 12.000 |  |
| 15 | 1.000 | 0.000 | 1.000 | 1.000 | 5 | 12.000 |  |
| 20 | 1.000 | 0.000 | 1.000 | 1.000 | 5 | 12.000 | 0.169 |
| 30 | 1.000 | 0.000 | 1.000 | 1.000 | 5 | 12.000 |  |
| 50 | 1.000 | 0.000 | 1.000 | 1.000 | 5 | 12.000 |  |
| 100 | 1.000 | 0.000 | 1.000 | 1.000 | 5 | 12.000 |  |

Label-permutation null at k = 20 (10 permutations, quick fits): mean 0.132, 95th percentile 0.169, chance 0.143; observed 1.000.
Pre-registration (b) asks for mean balanced accuracy ≥ 0.95 at some k ≤ 100: best 1.000 at k = 8 → PASS (within-study; context, not a finding).

The k = 20 panel selected on all animals (`panel_k20_all_animals.csv`):

| protein_id | gene_symbol | ensembl_gene | marker_tissue | mean_z_in_marker | next_highest_tissue | effect_z | frac_nan |
|---|---|---|---|---|---|---|---|
| NP_001004245.1 | Esam | ENSRNOG00000033217 | LUNG | 1.92 | WAT-SC | 1.48 | 0.10 |
| NP_001011991.1 | Ndrg1 | ENSRNOG00000007393 | KIDNEY | 1.99 | WAT-SC | 1.63 | 0.00 |
| NP_001012025.1 | Ubxn4 | ENSRNOG00000003625 | LIVER | 2.19 | KIDNEY | 2.18 | 0.50 |
| NP_001015022.1 | Tbc1d10a | ENSRNOG00000006394 | KIDNEY | 1.99 | LUNG | 1.32 | 0.10 |
| NP_001035266.2 | Mon2 | ENSRNOG00000004185 | LIVER | 2.10 | KIDNEY | 2.17 | 0.12 |
| NP_001094371.1 | Snta1 | ENSRNOG00000016062 | HEART | 1.98 | SKM-GN | 1.28 | 0.00 |
| NP_001099282.1 | Slc43a2 | ENSRNOG00000003835 | KIDNEY | 2.05 | WAT-SC | 1.92 | 0.05 |
| NP_001128428.1 | Ahdc1 | ENSRNOG00000042855 | LUNG | 0.72 | WAT-SC | 0.14 | 0.00 |
| NP_001164010.1 | Stt3b | ENSRNOG00000011922 | LIVER | 2.04 | KIDNEY | 1.38 | 0.00 |
| NP_001178622.1 | Arhgap21 | ENSRNOG00000008659 | HEART | 0.71 | LUNG | 0.05 | 0.00 |
| NP_036658.2 | Chgb | ENSRNOG00000021269 | CORTEX | 2.16 | HEART | 1.87 | 0.29 |
| NP_037270.1 | Itpr3 | ENSRNOG00000052795 | LUNG | 1.99 | KIDNEY | 1.71 | 0.21 |
| NP_037349.1 | Afdn | ENSRNOG00000023753 | LIVER | 0.84 | LUNG | 0.09 | 0.00 |
| NP_058891.1 | Gap43 | ENSRNOG00000001528 | CORTEX | 2.23 | WAT-SC | 2.21 | 0.10 |
| NP_073160.1 | Scg2 | ENSRNOG00000015055 | CORTEX | 2.17 | HEART | 1.88 | 0.29 |
| NP_112313.1 | Mif | ENSRNOG00000006589 | WAT-SC | 1.94 | HEART | 1.75 | 0.14 |
| NP_112413.2 | Mdh2 | ENSRNOG00000001440 | HEART | 2.05 | LIVER | 2.02 | 0.29 |
| NP_001020878.1 | Tfeb | ENSRNOG00000014666 | WAT-SC | 2.00 | LUNG | 1.52 | 0.10 |
| NP_001101857.2 | Sucla2 | ENSRNOG00000017481 | HEART | 2.08 | KIDNEY | 1.40 | 0.21 |
| NP_001094160.1 | Myh14 | ENSRNOG00000020014 | LUNG | 1.86 | HEART | 1.37 | 0.43 |

**Missingness-free check** (`variance_partition_complete.csv`, `diagnostic_accuracy_complete.csv`, `panel_curve_complete_summary.csv`): on the 417 proteins quantified in every vial (no NaN, no imputation), tissue R² of PC1 = 0.990, PC2 = 0.996; logreg_l2 accuracy 1.000 ± 0.000, means-removed 0.143; panel k = 20 balanced accuracy 1.000 ± 0.000.

Confused pairs at k = 20 (pooled over folds; 0 errors):

_none_

Confused pairs at k = 30 (pooled over folds; 0 errors):

_none_

## 4. Cross-tissue RNA–protein correlation on the same animals — `rna_protein_correlation.csv`, `rna_protein_correlation_summary.csv`

Tissue means over the animals shared by RNA-seq and proteomics in each tissue (CORTEX:50;SKM-GN:50;HEART:50;KIDNEY:50;LUNG:50;LIVER:50;WAT-SC:50); Spearman across the 7 tissues per gene (protein means averaged over the proteins mapping to the gene; a tissue mean needs ≥ 3 quantified animals). n = 1234 genes.
- Spearman median **0.429** (IQR 0.143–0.714); 43.8 % of genes > 0.5; 16.9 % < 0.
- Mismatched-pair null (RNA of one gene vs the protein of another, 19985 pairs): median 0.143, 95th percentile 0.786; 15.9 % of genes exceed the null 95th percentile.
- Same marker tissue (highest tissue mean) at RNA and protein: 32.8 % of genes (chance 1/7 = 14.3 %).
- By protein cross-tissue range: median Spearman 0.286 (Q1) → 0.643 (Q4) (`rna_protein_correlation_by_protein_range.csv`).

### RNA panel genes at the protein level — `rna_protein_panel_genes.csv`, `rna_protein_panel_summary.csv`

51 RNA panel genes (stable core, stability-51, k=20 all-animal panel); 17 have a quantified protein, 1 in all 7 tissues; 23 have their RNA marker tissue among the 7 proteomics tissues.
- Pre-registration (c): of the **10** testable genes (marker tissue among the 7, protein quantified), the protein's highest tissue mean is the RNA marker tissue in **8** (80 %) → PASS (threshold 70 %).
- Median Spearman over panel genes with all-7 means 0.107 vs 0.429 over all genes; 0 panel genes above the null 95th percentile.

| gene_symbol | lists | marker_tissue_19 | protein_quantified_any_tissue | prot_marker_tissue | rna_marker_tissue | spearman | testable_c | protein_marker_equals_rna_marker19 |
|---|---|---|---|---|---|---|---|---|
| Cldn18 | panel_k20_all_animals;stability_k20;stable_core_T5 | LUNG | True | LUNG | LUNG |  | True | True |
| Tbr1 | stability_k20 | CORTEX | True | CORTEX | CORTEX |  | True | True |
| Gnb3 | stability_k20 | HEART | True | CORTEX | HEART |  | True | False |
| F9 | stability_k20 | LIVER | True | KIDNEY | LIVER |  | True | False |
| Slco1b2 | stability_k20 | LIVER | True | LIVER | LIVER |  | True | True |
| Ankrd2 | stability_k20 | SKM-GN | True | SKM-GN | SKM-GN |  | True | True |
| RGD1565323 | stability_k20 | SKM-GN | True | SKM-GN | SKM-GN |  | True | True |
| Tnni1 | stability_k20 | SKM-GN | True | SKM-GN | SKM-GN |  | True | True |
| Trim29 | panel_k20_all_animals;stability_k20 | WAT-SC | True | WAT-SC | WAT-SC |  | True | True |
| Csn3 | stability_k20 | WAT-SC | True | WAT-SC | WAT-SC |  | True | True |
| Hbb-b1 | stability_k20 | BLOOD | True | KIDNEY | LUNG | 0.11 | False | False |
| Cyp21a1 | panel_k20_all_animals;stability_k20;stable_core_T5 | ADRNL | False | nan | nan |  | False | False |
| Gabrd | panel_k20_all_animals;stability_k20;stable_core_T5 | CORTEX | False | nan | nan |  | False | False |
| Lrrc10 | panel_k20_all_animals;stability_k20;stable_core_T5 | HEART | False | nan | nan |  | False | False |
| Fibcd1l1 | panel_k20_all_animals;stability_k20;stable_core_T5 | HIPPOC | False | nan | nan |  | False | False |
| Pmch | panel_k20_all_animals;stability_k20;stable_core_T5 | HYPOTH | False | nan | nan |  | False | False |
| Umod | panel_k20_all_animals;stability_k20;stable_core_T5 | KIDNEY | False | nan | nan |  | False | False |
| Akr1c3 | panel_k20_all_animals;stability_k20;stable_core_T5 | OVARY | False | nan | nan |  | False | False |
| Ccl25 | panel_k20_all_animals;stability_k20;stable_core_T5 | SMLINT | False | nan | nan |  | False | False |
| Fcrl5 | panel_k20_all_animals;stability_k20;stable_core_T5 | SPLEEN | False | nan | nan |  | False | False |
| Qrfpr | panel_k20_all_animals;stability_k20 | ADRNL | False | nan | nan |  | False | False |
| Dbh | stability_k20 | ADRNL | False | nan | nan |  | False | False |
| Mc2r | stability_k20 | ADRNL | False | nan | nan |  | False | False |
| Otop1 | panel_k20_all_animals;stability_k20 | BAT | False | nan | nan |  | False | False |
| Ucp1 | stability_k20 | BAT | False | nan | nan |  | False | False |
| Grifin | stability_k20 | BAT | True | WAT-SC | WAT-SC |  | False | False |
| Ebf2 | stability_k20 | BAT | True | HEART | WAT-SC |  | False | False |
| Hbq1b | panel_k20_all_animals;stability_k20 | BLOOD | False | nan | nan |  | False | False |
| Treml1 | stability_k20 | BLOOD | True | LUNG | LUNG |  | False | False |
| Abo3 | panel_k20_all_animals;stability_k20 | COLON | False | nan | nan |  | False | False |
| Hoxb13 | stability_k20 | COLON | False | nan | nan |  | False | False |
| St6galnac1 | stability_k20 | COLON | False | nan | nan |  | False | False |
| Ephx4 | stability_k20 | CORTEX | False | nan | nan |  | False | False |
| Slc17a6 | stability_k20 | HYPOTH | True | CORTEX | CORTEX |  | False | False |
| Cfhr2 | panel_k20_all_animals;stability_k20 | LIVER | False | nan | nan |  | False | False |
| Cfhr1 | stability_k20 | LIVER | False | nan | nan |  | False | False |
| Cpn2 | stability_k20 | LIVER | False | nan | nan |  | False | False |
| Mbl2 | stability_k20 | LIVER | False | nan | nan |  | False | False |
| Sftpb | stability_k20 | LUNG | False | nan | nan |  | False | False |
| AABR07069371.1 | stability_k20 | LUNG | False | nan | nan |  | False | False |
| Mss51 | panel_k20_all_animals;stability_k20 | SKM-GN | False | nan | nan |  | False | False |
| Lbx1 | stability_k20 | SKM-GN | False | nan | nan |  | False | False |
| Mybph | panel_k20_all_animals;stability_k20 | SKM-VL | True | SKM-GN | SKM-GN |  | False | False |
| Ccr3 | stability_k20 | SPLEEN | False | nan | nan |  | False | False |
| Pgk2 | panel_k20_all_animals;stability_k20 | TESTES | True | HEART | nan |  | False | False |
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

_Run time 0.8 min._