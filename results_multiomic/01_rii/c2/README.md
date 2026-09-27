# Phase 1 — the RII rescue (prot-pr, release c2.0)

Built by `scripts/multiomic/01_rii_rescue.py` on 2026-09-27 07:46 UTC; every number below is read from a CSV in this directory.

## Data
- Peptide-level reporter-ion intensities from the portal quant-id `prot-pr` folders, 7 tissues, 420 vials, 60 animals (`join_summary.csv`, `tissue_summary.csv`, `plex_summary.csv`).
- Per plex: reference channel dropped, contaminants dropped, peptides summed per protein, proteins with < 2 quantified peptides dropped, each channel normalised to its total (log2 ppm).
- Union 15802 proteins; **3827 quantified in every tissue** (the main matrix; NaN fraction 0.038, 2513 with no missing value at all); outer NaN fraction 0.521.
- **Plex is nested in tissue here too** (one plex = 10 samples of one tissue + that tissue's reference pool; `plex_id` has one tissue per level), so everything below is within-study evidence.

## 1. Variance partition (PCA on the stacked matrix, phase-03 code) — `variance_partition.csv`, `variance_partition_ratio.csv`

| PC | explained_RII | R2_tissue_RII | R2_tissue_null95_RII | R2_plex_id_RII | R2_sex_RII | R2_pid_RII | R2_tissue_ratio_same_code | R2_tissue_ratio_frozen_phase03 |
|---|---|---|---|---|---|---|---|---|
| PC1 | 0.3692 | 0.9888 | 0.0285 | 0.9919 | 0.0008 | 0.0024 | 0.0009 | 0.0008 |
| PC2 | 0.1315 | 0.9943 | 0.0281 | 0.9965 | 0.0000 | 0.0009 | 0.0039 | 0.0039 |
| PC3 | 0.1280 | 0.9967 | 0.0242 | 0.9984 | 0.0003 | 0.0007 | 0.2184 | 0.2078 |

Pre-registration (a) asks for tissue R² of PC1 > 0.5: observed **0.989** → PASS. Median-normalised sensitivity: PC1 R² 0.997 (`variance_partition_median_norm.csv`).
R2_plex_id equals R2_tissue up to the within-tissue plex split because plex is nested in tissue; it cannot be separated here.

## 2. Phase-04 diagnostic on every animal-grouped fold — `diagnostic_accuracy.csv`, `diagnostic_accuracy_summary.csv`

| quantity | mean | sd | n_folds | n_test_animals_mean | n_missingness_features_outer | n_missingness_features_inner | frozen_ratio_fold0 |
|---|---|---|---|---|---|---|---|
| model_as_fitted_acc | 1.000 | 0.000 | 5 | 12.000 | 13289 | 1314 | 0.976 |
| model_as_fitted_bal_acc | 1.000 | 0.000 | 5 | 12.000 | 13289 | 1314 |  |
| missingness_outer_acc | 1.000 | 0.000 | 5 | 12.000 | 13289 | 1314 | 1.000 |
| missingness_outer_bal_acc | 1.000 | 0.000 | 5 | 12.000 | 13289 | 1314 |  |
| missingness_inner_acc | 1.000 | 0.000 | 5 | 12.000 | 13289 | 1314 |  |
| missingness_inner_bal_acc | 1.000 | 0.000 | 5 | 12.000 | 13289 | 1314 |  |
| per_tissue_means_removed_acc | 0.143 | 0.000 | 5 | 12.000 | 13289 | 1314 | 0.193 |
| per_tissue_means_removed_bal_acc | 0.143 | 0.000 | 5 | 12.000 | 13289 | 1314 |  |
| chance | 0.143 | 0.000 | 5 | 12.000 | 13289 | 1314 | 0.143 |

Read: `model_as_fitted` is the pipeline's logreg_l2 on the RII matrix (5 folds, 12.0 test animals per fold). `missingness_outer` classifies tissue from the NaN pattern of the union matrix alone (which proteins were quantified in which tissue — an artefact of per-tissue searches, still present on the RII scale); `missingness_inner` does the same on the every-tissue matrix (within-tissue plex gaps only). `per_tissue_means_removed` erases every protein's tissue mean using the labels: on ratios it collapsed to chance because the classifier lived on normalisation offsets; on RII the tissue means ARE the fingerprint, so the collapse is expected — the difference between the two matrices is what the means are (section 1 and 4), not whether the classifier uses them.

## 3. Protein panel curve (RoundRobinSelector → logreg_l2, animal-grouped folds, pipeline unchanged) — `panel_curve_summary.csv`, `panel_null_k20.csv`

| k | bal_acc_mean | bal_acc_sd | acc_mean | macro_f1_mean | n_folds | n_test_animals_mean | null_q95_k20 |
|---|---|---|---|---|---|---|---|
| 1 | 0.350 | 0.014 | 0.350 | 0.303 | 5 | 12.000 |  |
| 2 | 0.462 | 0.031 | 0.462 | 0.434 | 5 | 12.000 |  |
| 3 | 0.800 | 0.028 | 0.800 | 0.792 | 5 | 12.000 |  |
| 5 | 0.993 | 0.007 | 0.993 | 0.993 | 5 | 12.000 |  |
| 8 | 1.000 | 0.000 | 1.000 | 1.000 | 5 | 12.000 |  |
| 10 | 1.000 | 0.000 | 1.000 | 1.000 | 5 | 12.000 |  |
| 15 | 1.000 | 0.000 | 1.000 | 1.000 | 5 | 12.000 |  |
| 20 | 1.000 | 0.000 | 1.000 | 1.000 | 5 | 12.000 | 0.165 |
| 30 | 1.000 | 0.000 | 1.000 | 1.000 | 5 | 12.000 |  |
| 50 | 1.000 | 0.000 | 1.000 | 1.000 | 5 | 12.000 |  |
| 100 | 1.000 | 0.000 | 1.000 | 1.000 | 5 | 12.000 |  |

Label-permutation null at k = 20 (10 permutations, quick fits): mean 0.141, 95th percentile 0.165, chance 0.143; observed 1.000.
Pre-registration (b) asks for mean balanced accuracy ≥ 0.95 at some k ≤ 100: best 1.000 at k = 8 → PASS (within-study; context, not a finding).

The k = 20 panel selected on all animals (`panel_k20_all_animals.csv`):

| protein_id | gene_symbol | ensembl_gene | marker_tissue | mean_z_in_marker | next_highest_tissue | effect_z | frac_nan |
|---|---|---|---|---|---|---|---|
| NP_001004225.1 | Gatd3a | ENSRNOG00000001211 | HEART | 2.02 | LIVER | 1.75 | 0.00 |
| NP_001007146.1 | Ctnna1 | ENSRNOG00000005796 | LUNG | 0.83 | KIDNEY | 0.18 | 0.00 |
| NP_001012183.1 | Cryz | ENSRNOG00000028319 | KIDNEY | 2.14 | LIVER | 1.50 | 0.00 |
| NP_001013148.1 | Sptbn1 | ENSRNOG00000005434 | CORTEX | 1.00 | LUNG | 0.13 | 0.00 |
| NP_001100242.2 | Ndufa7 | ENSRNOG00000006939 | HEART | 2.01 | KIDNEY | 1.58 | 0.00 |
| NP_001100900.1 | Pdpr | ENSRNOG00000022593 | HEART | 1.02 | CORTEX | 0.41 | 0.00 |
| NP_001128232.1 | Hnrnpl | ENSRNOG00000020235 | LUNG | 0.83 | LIVER | 0.16 | 0.00 |
| NP_036621.1 | Acaa1a | ENSRNOG00000032908 | LIVER | 2.10 | LUNG | 1.76 | 0.00 |
| NP_036640.1 | Atp2b2 | ENSRNOG00000030269 | CORTEX | 2.27 | LUNG | 2.26 | 0.17 |
| NP_036795.1 | Vamp2 | ENSRNOG00000006989 | CORTEX | 2.27 | LIVER | 2.10 | 0.02 |
| NP_037228.1 | Hba-a1 | ENSRNOG00000029886 | LUNG | 1.10 | WAT-SC | 0.11 | 0.00 |
| NP_077325.2 | Ak1 | nan | SKM-GN | 0.91 | HEART | 0.09 | 0.00 |
| NP_113876.2 | Sncg | nan | WAT-SC | 2.16 | KIDNEY | 1.89 | 0.02 |
| NP_446039.2 | S100a9 | nan | LUNG | 2.06 | LIVER | 1.57 | 0.00 |
| NP_446274.2 | S100a8 | ENSRNOG00000011557 | LUNG | 2.01 | LIVER | 1.50 | 0.00 |
| NP_446341.1 | Vwf | ENSRNOG00000019689 | LUNG | 2.06 | WAT-SC | 1.63 | 0.00 |
| NP_612537.1 | Ero1a | ENSRNOG00000006462 | LIVER | 0.83 | CORTEX | 0.01 | 0.00 |
| NP_640347.2 | Tubb3 | ENSRNOG00000017209 | CORTEX | 2.32 | HEART | 2.51 | 0.21 |
| XP_006243777.1 | Acy1 | ENSRNOG00000011189 | KIDNEY | 2.08 | LIVER | 1.52 | 0.00 |
| XP_006256952.1 | Txlng | ENSRNOG00000004971 | WAT-SC | 0.59 | LIVER | 0.00 | 0.12 |

**Missingness-free check** (`variance_partition_complete.csv`, `diagnostic_accuracy_complete.csv`, `panel_curve_complete_summary.csv`): on the 2513 proteins quantified in every vial (no NaN, no imputation), tissue R² of PC1 = 0.988, PC2 = 0.992; logreg_l2 accuracy 1.000 ± 0.000, means-removed 0.143; panel k = 20 balanced accuracy 1.000 ± 0.000.

Confused pairs at k = 20 (pooled over folds; 0 errors):

_none_

Confused pairs at k = 30 (pooled over folds; 0 errors):

_none_

## 4. Cross-tissue RNA–protein correlation on the same animals — `rna_protein_correlation.csv`, `rna_protein_correlation_summary.csv`

Tissue means over the animals shared by RNA-seq and proteomics in each tissue (CORTEX:50;SKM-GN:50;HEART:50;KIDNEY:50;LUNG:50;LIVER:50;WAT-SC:50); Spearman across the 7 tissues per gene (protein means averaged over the proteins mapping to the gene; a tissue mean needs ≥ 3 quantified animals). n = 3478 genes.
- Spearman median **0.714** (IQR 0.429–0.857); 70.5 % of genes > 0.5; 6.4 % < 0.
- Mismatched-pair null (RNA of one gene vs the protein of another, 19993 pairs): median 0.250, 95th percentile 0.857; 21.8 % of genes exceed the null 95th percentile.
- Same marker tissue (highest tissue mean) at RNA and protein: 47.8 % of genes (chance 1/7 = 14.3 %).
- By protein cross-tissue range: median Spearman 0.429 (Q1) → 0.857 (Q4) (`rna_protein_correlation_by_protein_range.csv`).

### RNA panel genes at the protein level — `rna_protein_panel_genes.csv`, `rna_protein_panel_summary.csv`

51 RNA panel genes (stable core, stability-51, k=20 all-animal panel); 23 have a quantified protein, 2 in all 7 tissues; 23 have their RNA marker tissue among the 7 proteomics tissues.
- Pre-registration (c): of the **15** testable genes (marker tissue among the 7, protein quantified), the protein's highest tissue mean is the RNA marker tissue in **11** (73 %) → PASS (threshold 70 %).
- Median Spearman over panel genes with all-7 means 0.536 vs 0.714 over all genes; 0 panel genes above the null 95th percentile.

| gene_symbol | lists | marker_tissue_19 | protein_quantified_any_tissue | prot_marker_tissue | rna_marker_tissue | spearman | testable_c | protein_marker_equals_rna_marker19 |
|---|---|---|---|---|---|---|---|---|
| Cfhr1 | stability_k20 | LIVER | True | KIDNEY | LIVER | 0.54 | True | False |
| F9 | stability_k20 | LIVER | True | LUNG | LIVER | 0.54 | True | False |
| Gabrd | panel_k20_all_animals;stability_k20;stable_core_T5 | CORTEX | True | CORTEX | CORTEX |  | True | True |
| Lrrc10 | panel_k20_all_animals;stability_k20;stable_core_T5 | HEART | True | HEART | HEART |  | True | True |
| Cldn18 | panel_k20_all_animals;stability_k20;stable_core_T5 | LUNG | True | LUNG | LUNG |  | True | True |
| Tbr1 | stability_k20 | CORTEX | True | CORTEX | CORTEX |  | True | True |
| Gnb3 | stability_k20 | HEART | True | HEART | HEART |  | True | True |
| Cfhr2 | panel_k20_all_animals;stability_k20 | LIVER | True | WAT-SC | LIVER |  | True | False |
| Cpn2 | stability_k20 | LIVER | True | LUNG | LIVER |  | True | False |
| Mbl2 | stability_k20 | LIVER | True | LIVER | LIVER |  | True | True |
| Ankrd2 | stability_k20 | SKM-GN | True | SKM-GN | SKM-GN |  | True | True |
| Mss51 | panel_k20_all_animals;stability_k20 | SKM-GN | True | SKM-GN | SKM-GN |  | True | True |
| RGD1565323 | stability_k20 | SKM-GN | True | SKM-GN | SKM-GN |  | True | True |
| Trim29 | panel_k20_all_animals;stability_k20 | WAT-SC | True | WAT-SC | WAT-SC |  | True | True |
| Csn3 | stability_k20 | WAT-SC | True | WAT-SC | WAT-SC |  | True | True |
| Cyp21a1 | panel_k20_all_animals;stability_k20;stable_core_T5 | ADRNL | False | nan | nan |  | False | False |
| Fibcd1l1 | panel_k20_all_animals;stability_k20;stable_core_T5 | HIPPOC | True | CORTEX | CORTEX |  | False | False |
| Pmch | panel_k20_all_animals;stability_k20;stable_core_T5 | HYPOTH | False | nan | nan |  | False | False |
| Umod | panel_k20_all_animals;stability_k20;stable_core_T5 | KIDNEY | False | nan | nan |  | False | False |
| Akr1c3 | panel_k20_all_animals;stability_k20;stable_core_T5 | OVARY | True | WAT-SC | WAT-SC |  | False | False |
| Ccl25 | panel_k20_all_animals;stability_k20;stable_core_T5 | SMLINT | False | nan | nan |  | False | False |
| Fcrl5 | panel_k20_all_animals;stability_k20;stable_core_T5 | SPLEEN | False | nan | nan |  | False | False |
| Qrfpr | panel_k20_all_animals;stability_k20 | ADRNL | False | nan | nan |  | False | False |
| Dbh | stability_k20 | ADRNL | False | nan | nan |  | False | False |
| Mc2r | stability_k20 | ADRNL | False | nan | nan |  | False | False |
| Otop1 | panel_k20_all_animals;stability_k20 | BAT | False | nan | nan |  | False | False |
| Ucp1 | stability_k20 | BAT | True | LUNG | KIDNEY |  | False | False |
| Grifin | stability_k20 | BAT | True | WAT-SC | WAT-SC |  | False | False |
| Ebf2 | stability_k20 | BAT | True | WAT-SC | WAT-SC |  | False | False |
| Hbq1b | panel_k20_all_animals;stability_k20 | BLOOD | False | nan | nan |  | False | False |
| Treml1 | stability_k20 | BLOOD | True | LUNG | LUNG |  | False | False |
| Hbb-b1 | stability_k20 | BLOOD | False | nan | nan |  | False | False |
| Abo3 | panel_k20_all_animals;stability_k20 | COLON | False | nan | nan |  | False | False |
| Hoxb13 | stability_k20 | COLON | False | nan | nan |  | False | False |
| St6galnac1 | stability_k20 | COLON | False | nan | nan |  | False | False |
| Ephx4 | stability_k20 | CORTEX | False | nan | nan |  | False | False |
| Slc17a6 | stability_k20 | HYPOTH | True | CORTEX | CORTEX |  | False | False |
| Slco1b2 | stability_k20 | LIVER | False | nan | nan |  | False | False |
| Sftpb | stability_k20 | LUNG | False | nan | nan |  | False | False |
| AABR07069371.1 | stability_k20 | LUNG | False | nan | nan |  | False | False |
| Tnni1 | stability_k20 | SKM-GN | False | nan | nan |  | False | False |
| Lbx1 | stability_k20 | SKM-GN | False | nan | nan |  | False | False |
| Mybph | panel_k20_all_animals;stability_k20 | SKM-VL | False | nan | nan |  | False | False |
| Ccr3 | stability_k20 | SPLEEN | False | nan | nan |  | False | False |
| Pgk2 | panel_k20_all_animals;stability_k20 | TESTES | True | KIDNEY | nan |  | False | False |
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

_Run time 2.3 min._