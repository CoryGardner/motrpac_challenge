# Phase 3 — protein fingerprint transfer to Jiang 2020 (authors’ cleaned relative abundance)

Built by `scripts/multiomic/03_prot_transfer.py` on 2026-09-27 07:17 UTC; every number is read from a CSV in this directory.

## Data and matching — `gene_overlap.csv`
- Source: MoTrPAC RII, 3637 proteins quantified in every tissue → 3570 genes; 3167 with a 1:1 human ortholog; **2731 matched** to the target after its completeness filter (≥ 80% non-missing; 6958 of 12627 target genes).
- Target: 201 samples, 14 donors, 32 tissues; **44 samples from 13 donors map to CORTEX;HEART;LIVER;LUNG;SKM-GN**; KIDNEY;WAT-SC have no human counterpart in the atlas; 26 human tissues are out-of-distribution for the 7-class protein fingerprint.
- Human atlas caveats: adult donors of both sexes, post-mortem GTEx tissue, one pooled all-tissue reference per TMT run, different search engine and FDR — species, age, death and processing are all mixed into the shift.

## Accuracy per target tissue (super-class scoring; OOD rows have no accuracy) — `accuracy_by_tissue.csv`, `accuracy_overall.csv`

| jiang_tissue | rat_class | n | n_donors | full | k20 | k50 |
|---|---|---|---|---|---|---|
| Brain - Cortex | CORTEX | 2 | 2 | 1.00 | 1.00 | 1.00 |
| Heart - Atrial Appendage | HEART | 11 | 11 | 0.00 | 0.00 | 0.00 |
| Heart - Left Ventricle | HEART | 7 | 7 | 0.00 | 0.00 | 0.00 |
| Liver | LIVER | 5 | 5 | 1.00 | 0.80 | 0.20 |
| Lung | LUNG | 8 | 8 | 1.00 | 0.88 | 1.00 |
| Muscle - Skeletal | SKM-GN | 11 | 11 | 0.73 | 0.64 | 0.91 |

Overall, sample-weighted over mapped samples: full 0.523, k20 0.455, k50 0.477 (chance 1/7 = 0.143).

Top prediction and main wrong call at k20:

| jiang_tissue | rat_class | n | top_prediction | top_prediction_frac | main_wrong_call | main_wrong_call_frac |
|---|---|---|---|---|---|---|
| Adrenal Gland | OOD | 7 | LIVER | 0.71 | LIVER | 0.71 |
| Artery - Aorta | OOD | 10 | WAT-SC | 1.00 | WAT-SC | 1.00 |
| Artery - Coronary | OOD | 2 | WAT-SC | 1.00 | WAT-SC | 1.00 |
| Artery - Tibial | OOD | 4 | WAT-SC | 1.00 | WAT-SC | 1.00 |
| Brain - Cerebellum | OOD | 3 | CORTEX | 1.00 | CORTEX | 1.00 |
| Brain - Cortex | CORTEX | 2 | CORTEX | 1.00 |  | 0.00 |
| Breast - Mammary Tissue | OOD | 7 | WAT-SC | 0.86 | WAT-SC | 0.86 |
| Colon - Sigmoid | OOD | 11 | HEART | 0.45 | HEART | 0.45 |
| Colon - Transverse | OOD | 4 | WAT-SC | 1.00 | WAT-SC | 1.00 |
| Esophagus - Gastroesophageal Junction | OOD | 10 | CORTEX | 0.40 | CORTEX | 0.40 |
| Esophagus - Mucosa | OOD | 4 | LUNG | 1.00 | LUNG | 1.00 |
| Esophagus - Muscularis | OOD | 8 | CORTEX | 0.88 | CORTEX | 0.88 |
| Heart - Atrial Appendage | HEART | 11 | KIDNEY | 1.00 | KIDNEY | 1.00 |
| Heart - Left Ventricle | HEART | 7 | KIDNEY | 1.00 | KIDNEY | 1.00 |
| Liver | LIVER | 5 | LIVER | 0.80 | KIDNEY | 0.20 |
| Lung | LUNG | 8 | LUNG | 0.88 | WAT-SC | 0.12 |
| Minor Salivary Gland | OOD | 1 | LUNG | 1.00 | LUNG | 1.00 |
| Muscle - Skeletal | SKM-GN | 11 | SKM-GN | 0.64 | CORTEX | 0.18 |
| Nerve - Tibial | OOD | 6 | WAT-SC | 0.50 | WAT-SC | 0.50 |
| Ovary | OOD | 4 | WAT-SC | 0.50 | WAT-SC | 0.50 |
| Pancreas | OOD | 8 | LIVER | 0.50 | LIVER | 0.50 |
| Pituitary | OOD | 2 | CORTEX | 1.00 | CORTEX | 1.00 |
| Prostate | OOD | 2 | LUNG | 0.50 | LUNG | 0.50 |
| Skin - Not Sun Exposed (Suprapubic) | OOD | 8 | WAT-SC | 0.62 | WAT-SC | 0.62 |
| Skin - Sun Exposed (Lower leg) | OOD | 7 | WAT-SC | 0.57 | WAT-SC | 0.57 |
| Small Intestine - Terminal Ileum | OOD | 3 | LIVER | 0.67 | LIVER | 0.67 |
| Spleen | OOD | 8 | LUNG | 1.00 | LUNG | 1.00 |
| Stomach | OOD | 11 | KIDNEY | 0.91 | KIDNEY | 0.91 |
| Testis | OOD | 3 | SKM-GN | 1.00 | SKM-GN | 1.00 |
| Thyroid | OOD | 11 | WAT-SC | 1.00 | WAT-SC | 1.00 |
| Uterus | OOD | 6 | WAT-SC | 1.00 | WAT-SC | 1.00 |
| Vagina | OOD | 7 | WAT-SC | 0.71 | WAT-SC | 0.71 |

## Conformal transfer (α = 0.1, LAC, calibration on 18 held-out MoTrPAC animals) — `conformal_transfer.csv`, `coverage_by_tissue.csv`, `ood_sets.csv`, `recalibration.csv`

| model | conformal | n_mapped | coverage_mapped | frac_empty_mapped | avg_set_size_mapped | ood_frac_empty | ood_avg_set_size | n_cal_animals |
|---|---|---|---|---|---|---|---|---|
| k20 | marginal | 44 | 0.068 | 0.932 | 0.068 | 0.975 | 0.025 | 18 |
| k20 | mondrian | 44 | 0.068 | 0.932 | 0.068 | 0.936 | 0.064 | 18 |
| k20 | floored | 44 | 0.068 | 0.932 | 0.068 | 0.930 | 0.070 | 18 |
| k50 | marginal | 44 | 0.045 | 0.955 | 0.045 | 0.981 | 0.019 | 18 |
| k50 | mondrian | 44 | 0.045 | 0.955 | 0.045 | 0.981 | 0.019 | 18 |
| k50 | floored | 44 | 0.045 | 0.955 | 0.045 | 0.981 | 0.019 | 18 |
| full | marginal | 44 | 0.045 | 0.955 | 0.045 | 1.000 | 0.000 | 18 |
| full | mondrian | 44 | 0.045 | 0.955 | 0.045 | 1.000 | 0.000 | 18 |
| full | floored | 44 | 0.045 | 0.955 | 0.045 | 1.000 | 0.000 | 18 |

Per tissue, k20, marginal vs Mondrian vs floored (rat calibration):

| rat_class_or_ood_tissue | floored | marginal | mondrian |
|---|---|---|---|
| Adrenal Gland |  |  |  |
| Artery - Aorta |  |  |  |
| Artery - Coronary |  |  |  |
| Artery - Tibial |  |  |  |
| Brain - Cerebellum |  |  |  |
| Breast - Mammary Tissue |  |  |  |
| CORTEX | 1.00 | 1.00 | 1.00 |
| Colon - Sigmoid |  |  |  |
| Colon - Transverse |  |  |  |
| Esophagus - Gastroesophageal Junction |  |  |  |
| Esophagus - Mucosa |  |  |  |
| Esophagus - Muscularis |  |  |  |
| HEART | 0.00 | 0.00 | 0.00 |
| LIVER | 0.00 | 0.00 | 0.00 |
| LUNG | 0.12 | 0.12 | 0.12 |
| Minor Salivary Gland |  |  |  |
| Nerve - Tibial |  |  |  |
| Ovary |  |  |  |
| Pancreas |  |  |  |
| Pituitary |  |  |  |
| Prostate |  |  |  |
| SKM-GN | 0.00 | 0.00 | 0.00 |
| Skin - Not Sun Exposed (Suprapubic) |  |  |  |
| Skin - Sun Exposed (Lower leg) |  |  |  |
| Small Intestine - Terminal Ileum |  |  |  |
| Spleen |  |  |  |
| Stomach |  |  |  |
| Testis |  |  |  |
| Thyroid |  |  |  |
| Uterus |  |  |  |
| Vagina |  |  |  |

What the sets hold for the out-of-distribution tissues at k20 (`ood_sets.csv`):

| model | conformal | organ | n | frac_empty | avg_set_size | max_set_size | most_frequent_members |
|---|---|---|---|---|---|---|---|
| k20 | marginal | Adrenal Gland | 7 | 1.00 | 0.00 | 0 |  |
| k20 | marginal | Artery - Aorta | 10 | 1.00 | 0.00 | 0 |  |
| k20 | marginal | Artery - Coronary | 2 | 1.00 | 0.00 | 0 |  |
| k20 | marginal | Artery - Tibial | 4 | 0.75 | 0.25 | 1 | WAT-SC (1) |
| k20 | marginal | Brain - Cerebellum | 3 | 0.00 | 1.00 | 1 | CORTEX (3) |
| k20 | marginal | Breast - Mammary Tissue | 7 | 1.00 | 0.00 | 0 |  |
| k20 | marginal | Colon - Sigmoid | 11 | 1.00 | 0.00 | 0 |  |
| k20 | marginal | Colon - Transverse | 4 | 1.00 | 0.00 | 0 |  |
| k20 | marginal | Esophagus - Gastroesophageal Junction | 10 | 1.00 | 0.00 | 0 |  |
| k20 | marginal | Esophagus - Mucosa | 4 | 1.00 | 0.00 | 0 |  |
| k20 | marginal | Esophagus - Muscularis | 8 | 1.00 | 0.00 | 0 |  |
| k20 | marginal | Minor Salivary Gland | 1 | 1.00 | 0.00 | 0 |  |
| k20 | marginal | Nerve - Tibial | 6 | 1.00 | 0.00 | 0 |  |
| k20 | marginal | Ovary | 4 | 1.00 | 0.00 | 0 |  |
| k20 | marginal | Pancreas | 8 | 1.00 | 0.00 | 0 |  |
| k20 | marginal | Pituitary | 2 | 1.00 | 0.00 | 0 |  |
| k20 | marginal | Prostate | 2 | 1.00 | 0.00 | 0 |  |
| k20 | marginal | Skin - Not Sun Exposed (Suprapubic) | 8 | 1.00 | 0.00 | 0 |  |
| k20 | marginal | Skin - Sun Exposed (Lower leg) | 7 | 1.00 | 0.00 | 0 |  |
| k20 | marginal | Small Intestine - Terminal Ileum | 3 | 1.00 | 0.00 | 0 |  |
| k20 | marginal | Spleen | 8 | 1.00 | 0.00 | 0 |  |
| k20 | marginal | Stomach | 11 | 1.00 | 0.00 | 0 |  |
| k20 | marginal | Testis | 3 | 1.00 | 0.00 | 0 |  |
| k20 | marginal | Thyroid | 11 | 1.00 | 0.00 | 0 |  |
| k20 | marginal | Uterus | 6 | 1.00 | 0.00 | 0 |  |
| k20 | marginal | Vagina | 7 | 1.00 | 0.00 | 0 |  |
| k20 | mondrian | Adrenal Gland | 7 | 1.00 | 0.00 | 0 |  |
| k20 | mondrian | Artery - Aorta | 10 | 1.00 | 0.00 | 0 |  |
| k20 | mondrian | Artery - Coronary | 2 | 0.50 | 0.50 | 1 | WAT-SC (1) |
| k20 | mondrian | Artery - Tibial | 4 | 0.00 | 1.00 | 1 | WAT-SC (4) |
| k20 | mondrian | Brain - Cerebellum | 3 | 0.33 | 0.67 | 1 | CORTEX (2) |
| k20 | mondrian | Breast - Mammary Tissue | 7 | 0.57 | 0.43 | 1 | WAT-SC (3) |
| k20 | mondrian | Colon - Sigmoid | 11 | 1.00 | 0.00 | 0 |  |
| k20 | mondrian | Colon - Transverse | 4 | 1.00 | 0.00 | 0 |  |
| k20 | mondrian | Esophagus - Gastroesophageal Junction | 10 | 1.00 | 0.00 | 0 |  |
| k20 | mondrian | Esophagus - Mucosa | 4 | 1.00 | 0.00 | 0 |  |
| k20 | mondrian | Esophagus - Muscularis | 8 | 1.00 | 0.00 | 0 |  |
| k20 | mondrian | Minor Salivary Gland | 1 | 1.00 | 0.00 | 0 |  |
| k20 | mondrian | Nerve - Tibial | 6 | 1.00 | 0.00 | 0 |  |
| k20 | mondrian | Ovary | 4 | 1.00 | 0.00 | 0 |  |
| k20 | mondrian | Pancreas | 8 | 1.00 | 0.00 | 0 |  |
| k20 | mondrian | Pituitary | 2 | 1.00 | 0.00 | 0 |  |
| k20 | mondrian | Prostate | 2 | 1.00 | 0.00 | 0 |  |
| k20 | mondrian | Skin - Not Sun Exposed (Suprapubic) | 8 | 1.00 | 0.00 | 0 |  |
| k20 | mondrian | Skin - Sun Exposed (Lower leg) | 7 | 1.00 | 0.00 | 0 |  |
| k20 | mondrian | Small Intestine - Terminal Ileum | 3 | 1.00 | 0.00 | 0 |  |
| k20 | mondrian | Spleen | 8 | 1.00 | 0.00 | 0 |  |
| k20 | mondrian | Stomach | 11 | 1.00 | 0.00 | 0 |  |
| k20 | mondrian | Testis | 3 | 1.00 | 0.00 | 0 |  |
| k20 | mondrian | Thyroid | 11 | 1.00 | 0.00 | 0 |  |
| k20 | mondrian | Uterus | 6 | 1.00 | 0.00 | 0 |  |
| k20 | mondrian | Vagina | 7 | 1.00 | 0.00 | 0 |  |
| k20 | floored | Adrenal Gland | 7 | 1.00 | 0.00 | 0 |  |
| k20 | floored | Artery - Aorta | 10 | 1.00 | 0.00 | 0 |  |
| k20 | floored | Artery - Coronary | 2 | 0.50 | 0.50 | 1 | WAT-SC (1) |
| k20 | floored | Artery - Tibial | 4 | 0.00 | 1.00 | 1 | WAT-SC (4) |
| k20 | floored | Brain - Cerebellum | 3 | 0.00 | 1.00 | 1 | CORTEX (3) |
| k20 | floored | Breast - Mammary Tissue | 7 | 0.57 | 0.43 | 1 | WAT-SC (3) |
| k20 | floored | Colon - Sigmoid | 11 | 1.00 | 0.00 | 0 |  |
| k20 | floored | Colon - Transverse | 4 | 1.00 | 0.00 | 0 |  |

_(first 60 of 78 rows)_

Recalibration on Jiang donors (pooled tissues), tested on the remaining donors:

| model | n_recal | n_test_individuals | draws | coverage_recalibrated | coverage_source_cal_same_test | frac_empty_recalibrated | frac_empty_source_cal | set_size_recalibrated |
|---|---|---|---|---|---|---|---|---|
| k20 | 3 | 11 | 20 | 0.900 | 0.070 | 0.000 | 0.930 | 3.497 |
| k20 | 5 | 9 | 20 | 0.917 | 0.077 | 0.000 | 0.923 | 2.954 |
| k50 | 3 | 11 | 20 | 0.926 | 0.046 | 0.000 | 0.954 | 4.027 |
| k50 | 5 | 9 | 20 | 0.968 | 0.053 | 0.000 | 0.947 | 3.110 |
| full | 3 | 11 | 20 | 0.960 | 0.050 | 0.000 | 0.950 | 5.359 |
| full | 5 | 9 | 20 | 0.950 | 0.045 | 0.000 | 0.955 | 4.348 |

## The ladder: protein beside RNA — `ladder_protein_vs_rna.csv`

| layer | model | accuracy_sample_weighted | n_samples | n_individuals | n_target_classes | coverage_source_cal_marginal | frac_empty_source_cal | avg_set_size_source_cal | coverage_recal_3 | coverage_recal_5 | set_size_recal_5 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| protein (RII → Jiang 2020) | k20 | 0.455 | 44 | 13 | 5 | 0.068 | 0.932 | 0.068 | 0.900 | 0.917 | 2.954 |
| protein (RII → Jiang 2020) | k50 | 0.477 | 44 | 13 | 5 | 0.045 | 0.955 | 0.045 | 0.926 | 0.968 | 3.110 |
| protein (RII → Jiang 2020) | full | 0.523 | 44 | 13 | 5 | 0.045 | 0.955 | 0.045 | 0.960 | 0.950 | 4.348 |
| RNA (counts → GTEx v8, frozen phase 13) | full | 0.855 | 2485 | 862 | 17 | 0.062 | 0.938 | 0.062 | 0.951 | 0.866 | 2.045 |
| RNA (counts → GTEx v8, frozen phase 13) | k20 | 0.654 | 2485 | 862 | 17 | 0.364 | 0.616 | 0.384 | 0.954 | 0.933 | 6.282 |
| RNA (counts → GTEx v8, frozen phase 13) | k50 | 0.781 | 2485 | 862 | 17 | 0.363 | 0.631 | 0.369 | 0.959 | 0.936 | 4.036 |

The two rows are different targets (Jiang protein: 5 mapped tissues, few donors; GTEx RNA: 17 tissues, 862 donors); Phase 5 puts both layers on the same Jiang samples.

## Reverse direction: panels selected on Jiang, scored on MoTrPAC RII — `reverse_direction.csv`, `reverse_direction_overall.csv`

| model | motrpac_tissue | in_jiang_classes | n_vials | n_animals | accuracy | top_prediction | top_prediction_frac |
|---|---|---|---|---|---|---|---|
| k20 | CORTEX | True | 60 | 60 | 0.00 | HEART | 1.00 |
| k20 | HEART | True | 60 | 60 | 0.88 | HEART | 0.88 |
| k20 | KIDNEY | False | 60 | 60 |  | HEART | 1.00 |
| k20 | LIVER | True | 60 | 60 | 1.00 | LIVER | 1.00 |
| k20 | LUNG | True | 60 | 60 | 0.00 | HEART | 1.00 |
| k20 | SKM-GN | True | 60 | 60 | 1.00 | SKM-GN | 1.00 |
| k20 | WAT-SC | False | 60 | 60 |  | SKM-GN | 0.70 |
| k50 | CORTEX | True | 60 | 60 | 0.00 | HEART | 1.00 |
| k50 | HEART | True | 60 | 60 | 0.83 | HEART | 0.83 |
| k50 | KIDNEY | False | 60 | 60 |  | HEART | 1.00 |
| k50 | LIVER | True | 60 | 60 | 1.00 | LIVER | 1.00 |
| k50 | LUNG | True | 60 | 60 | 0.03 | HEART | 0.67 |
| k50 | SKM-GN | True | 60 | 60 | 1.00 | SKM-GN | 1.00 |
| k50 | WAT-SC | False | 60 | 60 |  | LUNG | 0.58 |
| full | CORTEX | True | 60 | 60 | 0.00 | HEART | 1.00 |
| full | HEART | True | 60 | 60 | 0.83 | HEART | 0.83 |
| full | KIDNEY | False | 60 | 60 |  | HEART | 1.00 |
| full | LIVER | True | 60 | 60 | 1.00 | LIVER | 1.00 |
| full | LUNG | True | 60 | 60 | 0.10 | HEART | 0.90 |
| full | SKM-GN | True | 60 | 60 | 1.00 | SKM-GN | 1.00 |
| full | WAT-SC | False | 60 | 60 |  | HEART | 0.92 |

## Per-gene check of the k20 panel — `panel_gene_check.csv`, `panel_k20.csv`

| feature_ID | gene_symbol | marker_tissue | target_organ | source_effect_z | target_effect_z | target_top_organ | fails_in_target | weakened | human_gene |
|---|---|---|---|---|---|---|---|---|---|
| ENSRNOG00000017209 | Tubb3 | CORTEX | CORTEX | 2.55 | 3.57 | CORTEX | False | False | ENSG00000258947 |
| ENSRNOG00000021881 | Metap2 | WAT-SC | WAT-SC | -2.85 |  | nan | nan | nan | ENSG00000111142 |
| ENSRNOG00000028319 | Cryz | KIDNEY | KIDNEY | 1.88 |  | nan | nan | nan | ENSG00000116791 |
| ENSRNOG00000005796 | Ctnna1 | SKM-GN | SKM-GN | -3.05 | -2.19 | HEART | True | False | ENSG00000044115 |
| ENSRNOG00000019688 | Diaph1 | HEART | HEART | -3.34 | -2.40 | LIVER | True | False | ENSG00000131504 |
| ENSRNOG00000033835 | Dnm1 | CORTEX | CORTEX | 1.91 | 4.77 | CORTEX | False | False | ENSG00000106976 |
| ENSRNOG00000005096 | Bzw2 | LIVER | LIVER | -3.19 | -2.90 | SKM-GN | True | False | ENSG00000136261 |
| ENSRNOG00000016690 | Idi1 | HEART | HEART | -2.78 | -3.01 | LIVER | True | False | ENSG00000067064 |
| ENSRNOG00000007862 | Acat1 | KIDNEY | KIDNEY | 2.44 |  | nan | nan | nan | ENSG00000075239 |
| ENSRNOG00000016610 | Arhgap1 | SKM-GN | SKM-GN | -3.26 | -2.38 | HEART | True | False | ENSG00000175220 |
| ENSRNOG00000043094 | Oxct1 | LIVER | LIVER | -3.42 | -3.81 | HEART | True | False | ENSG00000083720 |
| ENSRNOG00000022593 | Pdpr | LIVER | LIVER | -3.06 | -4.37 | HEART | True | False | ENSG00000090857 |
| ENSRNOG00000011557 | S100a8 | LUNG | LUNG | 1.78 | 0.76 | LUNG | False | True | ENSG00000143546 |
| ENSRNOG00000011483 | S100a9 | LUNG | LUNG | 1.61 | 0.69 | LUNG | False | True | ENSG00000163220 |
| ENSRNOG00000020425 | Stim1 | HEART | HEART | -3.13 | -3.29 | SKM-GN | True | False | ENSG00000167323 |
| ENSRNOG00000002343 | Uchl1 | CORTEX | CORTEX | 2.12 | 3.16 | CORTEX | False | False | ENSG00000154277 |
| ENSRNOG00000026415 | Col14a1 | WAT-SC | WAT-SC | 2.06 |  | nan | nan | nan | ENSG00000187955 |
| ENSRNOG00000019689 | Vwf | LUNG | LUNG | 1.46 | 0.67 | LUNG | False | True | ENSG00000110799 |
| ENSRNOG00000001091 | Hip1r | SKM-GN | SKM-GN | -3.24 | -1.97 | CORTEX | True | False | ENSG00000130787 |
| ENSRNOG00000017328 | Pter | KIDNEY | KIDNEY | 1.81 |  | nan | nan | nan | ENSG00000165983 |

15 of 20 panel genes have their tissue in the target; 9 fail (effect ≤ 0 in the target: Ctnna1, Diaph1, Bzw2, Idi1, Arhgap1, Oxct1, Pdpr, Stim1, Hip1r).

## What this phase does not show

- Only 5 of the 7 protein classes have a human counterpart here (no kidney, no adipose in Jiang 2020); the accuracy is over those five, with 2 cortex samples.
- Species, adult age, post-mortem interval, TMT reference design and search pipeline shift together; the design does not separate them.
- Coverage numbers are on ≤ 50 mapped samples from ≤ 14 donors; the recalibration draws share donors across draws.

_Run time 0.0 min._