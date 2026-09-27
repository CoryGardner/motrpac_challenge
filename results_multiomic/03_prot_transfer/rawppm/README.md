# Phase 3 — protein fingerprint transfer to Jiang 2020 (raw reporter intensities re-normalised to channel-total log2 ppm)

Built by `scripts/multiomic/03_prot_transfer.py` on 2026-09-27 07:24 UTC; every number is read from a CSV in this directory.

## Data and matching — `gene_overlap.csv`
- Source: MoTrPAC RII, 3637 proteins quantified in every tissue → 3570 genes; 3167 with a 1:1 human ortholog; **2414 matched** to the target after its completeness filter (≥ 80% non-missing; 4737 of 12627 target genes).
- Target: 420 samples, 14 donors, 32 tissues; **94 samples from 13 donors map to CORTEX;HEART;LIVER;LUNG;SKM-GN**; KIDNEY;WAT-SC have no human counterpart in the atlas; 26 human tissues are out-of-distribution for the 7-class protein fingerprint.
- Human atlas caveats: adult donors of both sexes, post-mortem GTEx tissue, one pooled all-tissue reference per TMT run, different search engine and FDR — species, age, death and processing are all mixed into the shift.

## Accuracy per target tissue (super-class scoring; OOD rows have no accuracy) — `accuracy_by_tissue.csv`, `accuracy_overall.csv`

| jiang_tissue | rat_class | n | n_donors | full | k20 | k50 |
|---|---|---|---|---|---|---|
| Brain - Cortex | CORTEX | 4 | 2 | 1.00 | 1.00 | 1.00 |
| Heart - Atrial Appendage | HEART | 23 | 11 | 0.39 | 0.35 | 0.39 |
| Heart - Left Ventricle | HEART | 15 | 7 | 0.87 | 0.73 | 0.67 |
| Liver | LIVER | 11 | 5 | 1.00 | 0.82 | 0.09 |
| Lung | LUNG | 17 | 8 | 0.94 | 0.88 | 0.94 |
| Muscle - Skeletal | SKM-GN | 24 | 11 | 0.92 | 0.92 | 0.88 |

Overall, sample-weighted over mapped samples: full 0.798, k20 0.734, k50 0.649 (chance 1/7 = 0.143).

Top prediction and main wrong call at k20:

| jiang_tissue | rat_class | n | top_prediction | top_prediction_frac | main_wrong_call | main_wrong_call_frac |
|---|---|---|---|---|---|---|
| Adrenal Gland | OOD | 14 | KIDNEY | 0.50 | KIDNEY | 0.50 |
| Artery - Aorta | OOD | 21 | WAT-SC | 0.86 | WAT-SC | 0.86 |
| Artery - Coronary | OOD | 4 | WAT-SC | 1.00 | WAT-SC | 1.00 |
| Artery - Tibial | OOD | 8 | WAT-SC | 1.00 | WAT-SC | 1.00 |
| Brain - Cerebellum | OOD | 6 | CORTEX | 1.00 | CORTEX | 1.00 |
| Brain - Cortex | CORTEX | 4 | CORTEX | 1.00 |  | 0.00 |
| Breast - Mammary Tissue | OOD | 14 | WAT-SC | 0.86 | WAT-SC | 0.86 |
| Colon - Sigmoid | OOD | 23 | HEART | 0.30 | HEART | 0.30 |
| Colon - Transverse | OOD | 9 | WAT-SC | 1.00 | WAT-SC | 1.00 |
| Esophagus - Gastroesophageal Junction | OOD | 21 | WAT-SC | 0.29 | WAT-SC | 0.29 |
| Esophagus - Mucosa | OOD | 8 | LUNG | 1.00 | LUNG | 1.00 |
| Esophagus - Muscularis | OOD | 18 | CORTEX | 0.33 | CORTEX | 0.33 |
| Heart - Atrial Appendage | HEART | 23 | KIDNEY | 0.65 | KIDNEY | 0.65 |
| Heart - Left Ventricle | HEART | 15 | HEART | 0.73 | KIDNEY | 0.27 |
| Liver | LIVER | 11 | LIVER | 0.82 | KIDNEY | 0.18 |
| Lung | LUNG | 17 | LUNG | 0.88 | WAT-SC | 0.06 |
| Minor Salivary Gland | OOD | 2 | LUNG | 1.00 | LUNG | 1.00 |
| Muscle - Skeletal | SKM-GN | 24 | SKM-GN | 0.92 | HEART | 0.08 |
| Nerve - Tibial | OOD | 12 | CORTEX | 0.58 | CORTEX | 0.58 |
| Ovary | OOD | 8 | WAT-SC | 0.62 | WAT-SC | 0.62 |
| Pancreas | OOD | 17 | KIDNEY | 0.94 | KIDNEY | 0.94 |
| Pituitary | OOD | 4 | CORTEX | 1.00 | CORTEX | 1.00 |
| Prostate | OOD | 4 | LIVER | 0.50 | LIVER | 0.50 |
| Skin - Not Sun Exposed (Suprapubic) | OOD | 17 | WAT-SC | 0.65 | WAT-SC | 0.65 |
| Skin - Sun Exposed (Lower leg) | OOD | 14 | WAT-SC | 0.29 | WAT-SC | 0.29 |
| Small Intestine - Terminal Ileum | OOD | 6 | LIVER | 0.50 | LIVER | 0.50 |
| Spleen | OOD | 17 | LUNG | 1.00 | LUNG | 1.00 |
| Stomach | OOD | 23 | KIDNEY | 0.91 | KIDNEY | 0.91 |
| Testis | OOD | 6 | SKM-GN | 0.83 | SKM-GN | 0.83 |
| Thyroid | OOD | 23 | WAT-SC | 0.48 | WAT-SC | 0.48 |
| Uterus | OOD | 13 | WAT-SC | 0.62 | WAT-SC | 0.62 |
| Vagina | OOD | 14 | LUNG | 0.64 | LUNG | 0.64 |

## Conformal transfer (α = 0.1, LAC, calibration on 18 held-out MoTrPAC animals) — `conformal_transfer.csv`, `coverage_by_tissue.csv`, `ood_sets.csv`, `recalibration.csv`

| model | conformal | n_mapped | coverage_mapped | frac_empty_mapped | avg_set_size_mapped | ood_frac_empty | ood_avg_set_size | n_cal_animals |
|---|---|---|---|---|---|---|---|---|
| k20 | marginal | 94 | 0.085 | 0.915 | 0.085 | 0.966 | 0.034 | 18 |
| k20 | mondrian | 94 | 0.106 | 0.894 | 0.106 | 0.951 | 0.049 | 18 |
| k20 | floored | 94 | 0.106 | 0.894 | 0.106 | 0.951 | 0.049 | 18 |
| k50 | marginal | 94 | 0.064 | 0.936 | 0.064 | 0.988 | 0.012 | 18 |
| k50 | mondrian | 94 | 0.053 | 0.947 | 0.053 | 0.997 | 0.003 | 18 |
| k50 | floored | 94 | 0.064 | 0.936 | 0.064 | 0.988 | 0.012 | 18 |
| full | marginal | 94 | 0.000 | 1.000 | 0.000 | 1.000 | 0.000 | 18 |
| full | mondrian | 94 | 0.000 | 1.000 | 0.000 | 1.000 | 0.000 | 18 |
| full | floored | 94 | 0.000 | 1.000 | 0.000 | 1.000 | 0.000 | 18 |

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
| LUNG | 0.00 | 0.00 | 0.00 |
| Minor Salivary Gland |  |  |  |
| Nerve - Tibial |  |  |  |
| Ovary |  |  |  |
| Pancreas |  |  |  |
| Pituitary |  |  |  |
| Prostate |  |  |  |
| SKM-GN | 0.25 | 0.17 | 0.25 |
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
| k20 | marginal | Adrenal Gland | 14 | 1.00 | 0.00 | 0 |  |
| k20 | marginal | Artery - Aorta | 21 | 1.00 | 0.00 | 0 |  |
| k20 | marginal | Artery - Coronary | 4 | 1.00 | 0.00 | 0 |  |
| k20 | marginal | Artery - Tibial | 8 | 0.88 | 0.12 | 1 | WAT-SC (1) |
| k20 | marginal | Brain - Cerebellum | 6 | 0.17 | 0.83 | 1 | CORTEX (5) |
| k20 | marginal | Breast - Mammary Tissue | 14 | 1.00 | 0.00 | 0 |  |
| k20 | marginal | Colon - Sigmoid | 23 | 1.00 | 0.00 | 0 |  |
| k20 | marginal | Colon - Transverse | 9 | 1.00 | 0.00 | 0 |  |
| k20 | marginal | Esophagus - Gastroesophageal Junction | 21 | 1.00 | 0.00 | 0 |  |
| k20 | marginal | Esophagus - Mucosa | 8 | 1.00 | 0.00 | 0 |  |
| k20 | marginal | Esophagus - Muscularis | 18 | 1.00 | 0.00 | 0 |  |
| k20 | marginal | Minor Salivary Gland | 2 | 0.50 | 0.50 | 1 | LUNG (1) |
| k20 | marginal | Nerve - Tibial | 12 | 1.00 | 0.00 | 0 |  |
| k20 | marginal | Ovary | 8 | 1.00 | 0.00 | 0 |  |
| k20 | marginal | Pancreas | 17 | 1.00 | 0.00 | 0 |  |
| k20 | marginal | Pituitary | 4 | 1.00 | 0.00 | 0 |  |
| k20 | marginal | Prostate | 4 | 1.00 | 0.00 | 0 |  |
| k20 | marginal | Skin - Not Sun Exposed (Suprapubic) | 17 | 1.00 | 0.00 | 0 |  |
| k20 | marginal | Skin - Sun Exposed (Lower leg) | 14 | 1.00 | 0.00 | 0 |  |
| k20 | marginal | Small Intestine - Terminal Ileum | 6 | 0.83 | 0.17 | 1 | LUNG (1) |
| k20 | marginal | Spleen | 17 | 0.82 | 0.18 | 1 | LUNG (3) |
| k20 | marginal | Stomach | 23 | 1.00 | 0.00 | 0 |  |
| k20 | marginal | Testis | 6 | 1.00 | 0.00 | 0 |  |
| k20 | marginal | Thyroid | 23 | 1.00 | 0.00 | 0 |  |
| k20 | marginal | Uterus | 13 | 1.00 | 0.00 | 0 |  |
| k20 | marginal | Vagina | 14 | 1.00 | 0.00 | 0 |  |
| k20 | mondrian | Adrenal Gland | 14 | 1.00 | 0.00 | 0 |  |
| k20 | mondrian | Artery - Aorta | 21 | 1.00 | 0.00 | 0 |  |
| k20 | mondrian | Artery - Coronary | 4 | 1.00 | 0.00 | 0 |  |
| k20 | mondrian | Artery - Tibial | 8 | 0.62 | 0.38 | 1 | WAT-SC (3) |
| k20 | mondrian | Brain - Cerebellum | 6 | 0.17 | 0.83 | 1 | CORTEX (5) |
| k20 | mondrian | Breast - Mammary Tissue | 14 | 0.79 | 0.21 | 1 | WAT-SC (3) |
| k20 | mondrian | Colon - Sigmoid | 23 | 1.00 | 0.00 | 0 |  |
| k20 | mondrian | Colon - Transverse | 9 | 1.00 | 0.00 | 0 |  |
| k20 | mondrian | Esophagus - Gastroesophageal Junction | 21 | 1.00 | 0.00 | 0 |  |
| k20 | mondrian | Esophagus - Mucosa | 8 | 1.00 | 0.00 | 0 |  |
| k20 | mondrian | Esophagus - Muscularis | 18 | 1.00 | 0.00 | 0 |  |
| k20 | mondrian | Minor Salivary Gland | 2 | 0.50 | 0.50 | 1 | LUNG (1) |
| k20 | mondrian | Nerve - Tibial | 12 | 1.00 | 0.00 | 0 |  |
| k20 | mondrian | Ovary | 8 | 1.00 | 0.00 | 0 |  |
| k20 | mondrian | Pancreas | 17 | 1.00 | 0.00 | 0 |  |
| k20 | mondrian | Pituitary | 4 | 1.00 | 0.00 | 0 |  |
| k20 | mondrian | Prostate | 4 | 1.00 | 0.00 | 0 |  |
| k20 | mondrian | Skin - Not Sun Exposed (Suprapubic) | 17 | 1.00 | 0.00 | 0 |  |
| k20 | mondrian | Skin - Sun Exposed (Lower leg) | 14 | 1.00 | 0.00 | 0 |  |
| k20 | mondrian | Small Intestine - Terminal Ileum | 6 | 0.83 | 0.17 | 1 | LUNG (1) |
| k20 | mondrian | Spleen | 17 | 0.82 | 0.18 | 1 | LUNG (3) |
| k20 | mondrian | Stomach | 23 | 1.00 | 0.00 | 0 |  |
| k20 | mondrian | Testis | 6 | 1.00 | 0.00 | 0 |  |
| k20 | mondrian | Thyroid | 23 | 1.00 | 0.00 | 0 |  |
| k20 | mondrian | Uterus | 13 | 1.00 | 0.00 | 0 |  |
| k20 | mondrian | Vagina | 14 | 1.00 | 0.00 | 0 |  |
| k20 | floored | Adrenal Gland | 14 | 1.00 | 0.00 | 0 |  |
| k20 | floored | Artery - Aorta | 21 | 1.00 | 0.00 | 0 |  |
| k20 | floored | Artery - Coronary | 4 | 1.00 | 0.00 | 0 |  |
| k20 | floored | Artery - Tibial | 8 | 0.62 | 0.38 | 1 | WAT-SC (3) |
| k20 | floored | Brain - Cerebellum | 6 | 0.17 | 0.83 | 1 | CORTEX (5) |
| k20 | floored | Breast - Mammary Tissue | 14 | 0.79 | 0.21 | 1 | WAT-SC (3) |
| k20 | floored | Colon - Sigmoid | 23 | 1.00 | 0.00 | 0 |  |
| k20 | floored | Colon - Transverse | 9 | 1.00 | 0.00 | 0 |  |

_(first 60 of 78 rows)_

Recalibration on Jiang donors (pooled tissues), tested on the remaining donors:

| model | n_recal | n_test_individuals | draws | coverage_recalibrated | coverage_source_cal_same_test | frac_empty_recalibrated | frac_empty_source_cal | set_size_recalibrated |
|---|---|---|---|---|---|---|---|---|
| k20 | 3 | 11 | 20 | 0.921 | 0.081 | 0.000 | 0.919 | 1.634 |
| k20 | 5 | 9 | 20 | 0.923 | 0.085 | 0.000 | 0.915 | 1.538 |
| k50 | 3 | 11 | 20 | 0.916 | 0.063 | 0.001 | 0.937 | 1.483 |
| k50 | 5 | 9 | 20 | 0.917 | 0.071 | 0.000 | 0.929 | 1.445 |
| full | 3 | 11 | 20 | 0.904 | 0.000 | 0.011 | 1.000 | 1.419 |
| full | 5 | 9 | 20 | 0.911 | 0.000 | 0.002 | 1.000 | 1.223 |

## The ladder: protein beside RNA — `ladder_protein_vs_rna.csv`

| layer | model | accuracy_sample_weighted | n_samples | n_individuals | n_target_classes | coverage_source_cal_marginal | frac_empty_source_cal | avg_set_size_source_cal | coverage_recal_3 | coverage_recal_5 | set_size_recal_5 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| protein (RII → Jiang 2020) | k20 | 0.734 | 94 | 13 | 5 | 0.085 | 0.915 | 0.085 | 0.921 | 0.923 | 1.538 |
| protein (RII → Jiang 2020) | k50 | 0.649 | 94 | 13 | 5 | 0.064 | 0.936 | 0.064 | 0.916 | 0.917 | 1.445 |
| protein (RII → Jiang 2020) | full | 0.798 | 94 | 13 | 5 | 0.000 | 1.000 | 0.000 | 0.904 | 0.911 | 1.223 |
| RNA (counts → GTEx v8, frozen phase 13) | full | 0.855 | 2485 | 862 | 17 | 0.062 | 0.938 | 0.062 | 0.951 | 0.866 | 2.045 |
| RNA (counts → GTEx v8, frozen phase 13) | k20 | 0.654 | 2485 | 862 | 17 | 0.364 | 0.616 | 0.384 | 0.954 | 0.933 | 6.282 |
| RNA (counts → GTEx v8, frozen phase 13) | k50 | 0.781 | 2485 | 862 | 17 | 0.363 | 0.631 | 0.369 | 0.959 | 0.936 | 4.036 |

The two rows are different targets (Jiang protein: 5 mapped tissues, few donors; GTEx RNA: 17 tissues, 862 donors); Phase 5 puts both layers on the same Jiang samples.

## Reverse direction: panels selected on Jiang, scored on MoTrPAC RII — `reverse_direction.csv`, `reverse_direction_overall.csv`

| model | motrpac_tissue | in_jiang_classes | n_vials | n_animals | accuracy | top_prediction | top_prediction_frac |
|---|---|---|---|---|---|---|---|
| k20 | CORTEX | True | 60 | 60 | 0.00 | HEART | 1.00 |
| k20 | HEART | True | 60 | 60 | 0.82 | HEART | 0.82 |
| k20 | KIDNEY | False | 60 | 60 |  | HEART | 1.00 |
| k20 | LIVER | True | 60 | 60 | 1.00 | LIVER | 1.00 |
| k20 | LUNG | True | 60 | 60 | 1.00 | LUNG | 1.00 |
| k20 | SKM-GN | True | 60 | 60 | 1.00 | SKM-GN | 1.00 |
| k20 | WAT-SC | False | 60 | 60 |  | LUNG | 0.58 |
| k50 | CORTEX | True | 60 | 60 | 0.00 | HEART | 1.00 |
| k50 | HEART | True | 60 | 60 | 0.38 | SKM-GN | 0.62 |
| k50 | KIDNEY | False | 60 | 60 |  | HEART | 1.00 |
| k50 | LIVER | True | 60 | 60 | 1.00 | LIVER | 1.00 |
| k50 | LUNG | True | 60 | 60 | 1.00 | LUNG | 1.00 |
| k50 | SKM-GN | True | 60 | 60 | 1.00 | SKM-GN | 1.00 |
| k50 | WAT-SC | False | 60 | 60 |  | LUNG | 0.82 |
| full | CORTEX | True | 60 | 60 | 0.00 | HEART | 1.00 |
| full | HEART | True | 60 | 60 | 0.00 | SKM-GN | 1.00 |
| full | KIDNEY | False | 60 | 60 |  | HEART | 1.00 |
| full | LIVER | True | 60 | 60 | 1.00 | LIVER | 1.00 |
| full | LUNG | True | 60 | 60 | 1.00 | LUNG | 1.00 |
| full | SKM-GN | True | 60 | 60 | 1.00 | SKM-GN | 1.00 |
| full | WAT-SC | False | 60 | 60 |  | LUNG | 1.00 |

## Per-gene check of the k20 panel — `panel_gene_check.csv`, `panel_k20.csv`

| feature_ID | gene_symbol | marker_tissue | target_organ | source_effect_z | target_effect_z | target_top_organ | fails_in_target | weakened | human_gene |
|---|---|---|---|---|---|---|---|---|---|
| ENSRNOG00000017209 | Tubb3 | CORTEX | CORTEX | 2.55 | 3.34 | CORTEX | False | False | ENSG00000258947 |
| ENSRNOG00000021881 | Metap2 | WAT-SC | WAT-SC | -2.85 |  | nan | nan | nan | ENSG00000111142 |
| ENSRNOG00000028319 | Cryz | KIDNEY | KIDNEY | 1.88 |  | nan | nan | nan | ENSG00000116791 |
| ENSRNOG00000005796 | Ctnna1 | SKM-GN | SKM-GN | -3.05 | -2.61 | LUNG | True | False | ENSG00000044115 |
| ENSRNOG00000019688 | Diaph1 | HEART | HEART | -3.34 | -2.16 | LIVER | True | False | ENSG00000131504 |
| ENSRNOG00000033835 | Dnm1 | CORTEX | CORTEX | 1.91 | 4.63 | CORTEX | False | False | ENSG00000106976 |
| ENSRNOG00000005096 | Bzw2 | LIVER | LIVER | -3.19 | -1.83 | SKM-GN | True | False | ENSG00000136261 |
| ENSRNOG00000016690 | Idi1 | HEART | HEART | -2.78 | -2.43 | LIVER | True | False | ENSG00000067064 |
| ENSRNOG00000007862 | Acat1 | KIDNEY | KIDNEY | 2.44 |  | nan | nan | nan | ENSG00000075239 |
| ENSRNOG00000016610 | Arhgap1 | SKM-GN | SKM-GN | -3.26 | -2.53 | HEART | True | False | ENSG00000175220 |
| ENSRNOG00000043094 | Oxct1 | LIVER | LIVER | -3.42 | -3.66 | HEART | True | False | ENSG00000083720 |
| ENSRNOG00000022593 | Pdpr | LIVER | LIVER | -3.06 | -2.79 | HEART | True | False | ENSG00000090857 |
| ENSRNOG00000011557 | S100a8 | LUNG | LUNG | 1.78 | 0.75 | LUNG | False | True | ENSG00000143546 |
| ENSRNOG00000011483 | S100a9 | LUNG | LUNG | 1.61 | 0.81 | LUNG | False | False | ENSG00000163220 |
| ENSRNOG00000020425 | Stim1 | HEART | HEART | -3.13 | -2.41 | LIVER | True | False | ENSG00000167323 |
| ENSRNOG00000002343 | Uchl1 | CORTEX | CORTEX | 2.12 | 3.35 | CORTEX | False | False | ENSG00000154277 |
| ENSRNOG00000026415 | Col14a1 | WAT-SC | WAT-SC | 2.06 |  | nan | nan | nan | ENSG00000187955 |
| ENSRNOG00000019689 | Vwf | LUNG | LUNG | 1.46 | 1.42 | LUNG | False | False | ENSG00000110799 |
| ENSRNOG00000001091 | Hip1r | SKM-GN | SKM-GN | -3.24 | -2.28 | CORTEX | True | False | ENSG00000130787 |
| ENSRNOG00000017328 | Pter | KIDNEY | KIDNEY | 1.81 |  | nan | nan | nan | ENSG00000165983 |

15 of 20 panel genes have their tissue in the target; 9 fail (effect ≤ 0 in the target: Ctnna1, Diaph1, Bzw2, Idi1, Arhgap1, Oxct1, Pdpr, Stim1, Hip1r).

## What this phase does not show

- Only 5 of the 7 protein classes have a human counterpart here (no kidney, no adipose in Jiang 2020); the accuracy is over those five, with 2 cortex samples.
- Species, adult age, post-mortem interval, TMT reference design and search pipeline shift together; the design does not separate them.
- Coverage numbers are on ≤ 50 mapped samples from ≤ 14 donors; the recalibration draws share donors across draws.

_Run time 0.0 min._