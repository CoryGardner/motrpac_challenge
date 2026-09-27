*Draft content for a new site page (Markdown only; `site/` unchanged). Each figure names the CSV it would be drawn from.*

### Title: The proteome and metabolome carry the fingerprint too — on the right scale

**Lead.** The submission left proteomics and metabolomics out of the cross-tissue fingerprint. On the portal's reporter-ion intensities tissue explains R² 0.991 of the first principal component (distributed ratios: 0.0009); a 20-protein panel then transfers to a human TMT atlas and a 20-metabolite panel to a mouse metabolome atlas, both with the coverage collapse and recalibration recovery already seen for RNA.

**Figure 1 — Same proteins, two scales.** Tissue R² of PC1–PC3 on the distributed ratio matrix vs the reporter-ion (log2 ppm) matrix, with the label-permutation null. Source: `results_multiomic/01_rii/variance_partition.csv`, `variance_partition_ratio.csv`, `variance_partition_complete.csv` (missingness-free repeat).

**Figure 2 — The protein ladder beside the RNA ladder.** Accuracy and 90 %-set coverage (MoTrPAC calibration) for k20 / k50 / full: protein → Jiang 2020 (44 samples, 13 donors; primary scale and raw-ppm scale) and RNA → GTEx (frozen). Source: `results_multiomic/03_prot_transfer/ladder_protein_vs_rna.csv`, `rawppm/accuracy_overall.csv`, `coverage_ci.csv` (donor-bootstrap intervals).

**Figure 3 — Recalibration on a few donors.** Coverage vs number of target individuals used for recalibration (0, 3, 5) with set size, for protein (Jiang) and metabolites (Sato, MW ST003188). Source: `results_multiomic/03_prot_transfer/recalibration.csv`, `results_multiomic/04_metab_transfer/<leg>/recalibration.csv`.

**Figure 4 — Where the protein panel is wrong.** Confusion of Jiang tissues (rows, mapped and out-of-distribution) against the 7 rat classes at k20; out-of-distribution tissues get empty sets in 0.975 of cases. Source: `results_multiomic/03_prot_transfer/confusion_k20.csv`, `ood_sets.csv`, `accuracy_by_tissue.csv`.

**Figure 5 — Metabolites: 70 mice, 12 organs, 5 ages.** Per-organ accuracy of the deep-platform panel on the mouse aging atlas (k20 0.637, full 0.753, chance 0.111), flat across ages; and the Sato sedentary → exercised invariance test (native panel: accuracy 1.000, coverage 0.926). Source: `results_multiomic/04_metab_transfer/deep_mw/accuracy_by_tissue.csv`, `accuracy_by_stage.csv`, `hilic_sato/invariance_native_sedentary_to_exercised.csv`.

**Figure 6 — Fusion under shift.** Accuracy, coverage and empty-set fraction of RNA, protein, late-mean fusion and stacked fusion on the same 42 Jiang samples (k20 and full), with donor-bootstrap intervals. Source: `results_multiomic/05_fusion_transfer/fusion_transfer_summary.csv`.

**Figure 7 — The design that is not nested.** Cramér's V (tissue × batch variable) and the fraction of estimable tissue pairs for every MoTrPAC layer and every external dataset; Jiang 2020 at V 0.324 and 424/528. Source: `results_multiomic/06_external_identifiability/design_comparison.csv`.

**Figure 8 — RNA markers at the protein level.** Distribution of the cross-tissue RNA–protein Spearman over 3701 genes with the mismatched-pair null, and the marker-tissue agreement of the RNA panel genes (17 of 20). Source: `results_multiomic/01_rii/rna_protein_correlation.csv`, `rna_protein_panel_genes.csv`.

**Caveats box.** Within-study accuracies are context: plex is nested in tissue on the reporter-ion scale too, and the NaN pattern alone identifies the tissue. Human atlases are adult, post-mortem and differently processed; kidney and adipose have no human protein target. RefMet name matches across platforms are name matches. Every number on this page is read from the CSV named beside it.
