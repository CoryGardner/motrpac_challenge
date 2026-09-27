# Multiomic overnight report — can the proteome and metabolome carry a tissue fingerprint that transfers?

Branch `multiomic-overnight`. Pre-registration: `docs/PREREGISTRATION_MULTIOMIC.md`. Log: `docs/MULTIOMIC_LOG.md`.
Every number in this file is read from a CSV under `results_multiomic/` by `scripts/multiomic/build_report.py`; the file is named beside each number. Within-study accuracy is context, never a finding (the audit: batch is nested in tissue).

## Significant findings so far (2026-09-27 07:43 UTC; phases with results: 0, 1, 2, 3, 4, 5, 6, 7)

1. **Protein fingerprint transfers to a human proteome atlas, with the RNA coverage pattern.** A 20-protein panel selected on MoTrPAC RII names the tissue of 0.455 of 44 Jiang 2020 samples (donor-bootstrap 95 % CI 0.36–0.55; 13 donors, 5 mapped tissues; k50 0.477, full 0.523; chance 0.143; on the raw reporter scale re-normalised like the RII, k20 rises to 0.734 on 94 samples); marginal conformal coverage with MoTrPAC calibration is 0.068 (CI 0.00–0.14; empty sets 0.932) against a nominal 0.90, and recalibration on 5 donors restores 0.917 at set size 2.954. Human adults, post-mortem, different TMT design — all shifts at once. — `results_multiomic/03_prot_transfer/accuracy_overall.csv`, `conformal_transfer.csv`, `recalibration.csv`
2. **An external TMT design with tissue crossed with plex exists.** In Jiang 2020 each TMT run holds up to 8 tissues from several donors (56 runs, 448 samples, Cramér's V run × tissue 0.324, 424 of 528 tissue pairs estimable within a run) versus MoTrPAC PROT (V = 1.000, 0 of 21). Nesting is a choice of design, not a property of the assay. — `results_multiomic/06_external_identifiability/design_comparison.csv`
3. **Metabolite fingerprint transfer (deep_mw).** A k20 RefMet-named metabolite panel selected on MoTrPAC (9 tissues, 44 matched names) names the tissue of 0.637 of 490 external samples (animal-bootstrap 95 % CI 0.61–0.66; 70 mice; full model 0.753; chance 0.111); coverage with MoTrPAC calibration 0.153, recalibrated on 5 mice 0.920. — `results_multiomic/04_metab_transfer/deep_mw/accuracy_overall.csv`, `conformal_transfer.csv`, `recalibration.csv`
4. **RII tissue-axis recovery.** On reporter-ion intensities normalised to the channel total, tissue explains R² = 0.991 of PC1 (PC2 0.998; label-permutation null 95th pct 0.029) versus 0.0009 on the distributed ratio matrix with the same code; n = 420 vials, 60 animals, 3637 proteins in every tissue; identical on the 2393 proteins with no missing value (R² 0.990). Within-study: plex is nested in tissue. — `results_multiomic/01_rii/variance_partition.csv`, `variance_partition_ratio.csv`, `variance_partition_complete.csv`
5. **RNA panel markers hold at the protein level.** Of the 20 RNA panel genes whose marker tissue is one of the 7 proteomics tissues and whose protein is quantified, 17 (85 %) have the same marker tissue at the protein level (chance 1/7); over all 3701 genes with both layers, the cross-tissue RNA–protein Spearman has median 0.679 (IQR 0.429–0.821; mismatched-pair null median 0.286, 95th pct 0.857); n = 50 shared animals per tissue. — `results_multiomic/01_rii/rna_protein_panel_summary.csv`, `rna_protein_correlation_summary.csv`
6. **Metabolite fingerprint transfer (hilic_sato).** A k20 RefMet-named metabolite panel selected on MoTrPAC (19 tissues, 58 matched names) names the tissue of 0.277 of 191 external samples (animal-bootstrap 95 % CI 0.26–0.30; 24 mice; full model 0.445; chance 0.053); coverage with MoTrPAC calibration 0.052, recalibrated on 5 mice 0.899; a native panel fit on sedentary mice keeps accuracy 1.000 / coverage 0.926 on exercised mice (RNA analogue 0.961 / 0.903). — `results_multiomic/04_metab_transfer/hilic_sato/accuracy_overall.csv`, `conformal_transfer.csv`, `recalibration.csv`
7. **Metabolite fingerprint transfer (deep_sato).** A k20 RefMet-named metabolite panel selected on MoTrPAC (9 tissues, 83 matched names) names the tissue of 0.329 of 167 external samples (animal-bootstrap 95 % CI 0.27–0.39; 24 mice; full model 0.527; chance 0.111); coverage with MoTrPAC calibration 0.000, recalibrated on 5 mice 0.910; a native panel fit on sedentary mice keeps accuracy 1.000 / coverage 0.926 on exercised mice (RNA analogue 0.961 / 0.903). — `results_multiomic/04_metab_transfer/deep_sato/accuracy_overall.csv`, `conformal_transfer.csv`, `recalibration.csv`
8. **Fusion under the species shift (same 42 Jiang samples, 12 donors, 5 tissues).** At k20 accuracy / MoTrPAC-calibrated coverage / empty sets: RNA 1.000 / 0.452 / 0.548, protein 0.476 / 0.071 / 0.929, late-mean fusion 0.738 / 0.048 / 0.952, stacked 0.738 / 0.048 / 0.952; pre-registered robustness rule: late mean FAIL, stacked FAIL. The two layers' confusion structures on Jiang correlate at 0.259 (soft off-diagonal, k20). — `results_multiomic/05_fusion_transfer/fusion_transfer_summary.csv`, `confusion_structure.csv`

## Phase 0 — setup and pre-registration

- question · fix predictions (a)–(f) before any data are analysed.
- data · none.
- design · `docs/PREREGISTRATION_MULTIOMIC.md` (six predictions with pass rules, fixed analysis choices, drop policy).
- result · written and committed before Phase 1 started.
- what it does not show · nothing; it is the contract.

## Phase 1 — the RII rescue (local data only)

- question · does MoTrPAC proteomics carry a tissue axis on a scale where cross-tissue comparison is defined (pre-registration a–c)?
- data · portal quant-id `prot-pr` reporter-ion intensities, 7 tissues, 420 vials, 60 animals; peptides summed per protein per plex, < 2 peptides dropped, channel-total normalisation, log2 ppm; 17396 proteins in the union, 3637 in every tissue (NaN 3.5 %), 2393 with no missing value (`join_summary.csv`).
- design · phase-03 PCA/R² code, phase-04 diagnostic on all 5 animal-grouped folds, `tfp.models.panel_curve` unchanged (RoundRobinSelector → logreg_l2), label-permutation nulls; RNA–protein Spearman of tissue means over the same animals with a mismatched-pair null.
- result · **(a) PASS**: tissue R² of PC1 = 0.991 (null95 0.029; ratio matrix 0.0009; median-normalised 0.997; complete proteins 0.990). **(b) PASS** (context): balanced accuracy ≥ 0.95 from k = 5; k = 20 gives 1.000 ± 0.000 over 5 folds (12 test animals each; permutation null 95th pct 0.156); the same on the complete-protein matrix (1.000). Diagnostic: missingness alone still classifies tissue (1.000 outer, 1.000 inner), per-tissue-mean removal collapses to 0.143 (chance 0.143; ratios n/a). **(c) PASS**: 17 of 20 testable RNA panel genes (85 %) keep their marker tissue at the protein level. Cross-tissue RNA–protein Spearman over 3701 genes: median 0.679, 68.4 % above 0.5, 16.6 % above the mismatched-pair null 95th percentile (0.857); same marker tissue in 44.9 % of genes (chance 14.3 %); correlation rises with the protein's cross-tissue range (median 0.429 → 0.821 by quartile).
- what it does not show · transfer. One plex is one tissue plus that tissue's reference pool, so tissue and plex are confounded exactly as in the audit; the axis is real on this scale but a per-plex processing offset cannot be excluded without an external dataset. Accuracy 1.000 is a ceiling-task number and is reported as context only. The mismatched-pair null median of the RNA–protein correlation is far above zero, meaning a large part of any gene's cross-tissue agreement is a shared tissue structure (e.g. muscle/heart vs brain), not gene-specific.
- files · `results_multiomic/01_rii/README.md` (built from the CSVs), matrices `rii_inner_log2ppm.parquet`, `rii_outer_log2ppm.parquet`, `rii_meta.csv`.
- secondary, phospho (`prot-ph`, `results_multiomic/01_rii/ph/`): 1236 phosphoprotein groups in every tissue; tissue R² of PC1 = 0.986; k = 20 balanced accuracy 1.000.

## Phase 2 — data discovery

- question · which external proteome and metabolome atlases are reachable without a login, and can they be matched to MoTrPAC's features?
- data · 19 files, 237 MB downloaded (`results_multiomic/02_discovery/download_log.csv`: URL, bytes, sha256, time); 15 attempts logged (`attempts.csv`).
- design · prompt order: Jiang 2020 → Wang 2019 → mouse atlas → rat atlas search → Sato 2022 → Metabolomics Workbench → MetaboLights; 20 min / 20 GB per attempt; processed tables only.
- result · **Jiang 2020 human proteome map** (protein (TMT, human); 201 samples, 32 tissues; 3097 of 3570 RII genes by 1:1 ortholog); **Wang 2019 human tissue atlas** (protein (label-free iBAQ, human); 29 samples, 29 tissues; 3123 of 3570 RII genes by 1:1 ortholog); **Geiger 2013 mouse tissue proteome** (protein (SILAC H/L ratios to one SILAC-mouse standard; mouse); 28 samples, 28 tissues; 2800 RII genes by gene symbol (rat↔mouse, caveat: symbol match, not orthology table)); **Sato 2022 atlas of exercise metabolism** (metabolite (Metabolon HD4 untargeted; mouse); 191 samples, 8 tissues; RefMet-matched names: 100 of 1159 queried, 7 in MoTrPAC (any platform), 7 in HILIC+; crude name overlap 191); **Metabolomics Workbench ST003188** (metabolite (targeted RP-negative triple-quad; mouse); 840 samples, 12 tissues; 109 RefMet names in MoTrPAC (any platform), 46 in HILIC+). Negative: no rat multi-tissue proteome exists in PRIDE/ProteomeXchange; Sato 2022 is not on Metabolomics Workbench; Wang 2019's PRIDE bundle (28 GB) is over the box (`datasets.csv`, `attempts.csv`).
- what it does not show · nothing about transfer; the tissue maps (`tissue_maps.csv`) carry imperfect matches (atrial appendage → HEART, serum → PLASMA, epididymal fat → WAT-SC) that the later phases inherit.
- files · `results_multiomic/02_discovery/README.md`.

## Phase 3 — protein fingerprint transfer

- question · does a protein panel selected on MoTrPAC RII name the tissue of an independently processed human proteome, and what happens to conformal coverage (pre-registration d)?
- data · Jiang et al. 2020 (Cell 183:269) cleaned relative protein abundances: 201 TMT samples, 14 GTEx donors, 32 tissues; 2731 genes matched to the 3570 RII genes through 1:1 orthologs; 44 samples from 13 donors in the 5 mapped classes (CORTEX;HEART;LIVER;LUNG;SKM-GN); KIDNEY and WAT-SC have no target; 26 human tissues are OOD (`results_multiomic/03_prot_transfer/gene_overlap.csv`).
- design · `tfp.transfer` unchanged: panels on all MoTrPAC animals, z-scores within dataset, super-class scoring, conformal sets calibrated on 30 % held-out MoTrPAC animals (α = 0.10), recalibration on 3 and 5 donors (20 draws), reverse direction, per-gene check; sensitivity run on raw reporter intensities re-normalised like the RII (`rawppm/`).
- result · accuracy over the 44 mapped samples: full **0.523** (donor-bootstrap 95 % CI 0.39–0.63), k20 **0.455** (donor-bootstrap 95 % CI 0.36–0.55), k50 **0.477** (donor-bootstrap 95 % CI 0.40–0.56) (chance 0.143). Marginal coverage with MoTrPAC calibration at k20: **0.068** (empty sets 0.932, set size 0.068); floored Mondrian 0.068. Recalibration on 3 donors → 0.900 (set size 3.497), on 5 donors → 0.917 (set size 2.954). OOD tissues at k20 (marginal): empty-set fraction 0.975. Reverse direction (panel selected on Jiang, scored on MoTrPAC RII): full 0.587, k20 0.577, k50 0.573. Pre-registration (d) asks for accuracy ≥ 0.429 and coverage < 0.90 with MoTrPAC calibration: PASS at k20 (accuracy 0.455, coverage 0.068); RNA ladder for comparison in `ladder_protein_vs_rna.csv`.
- what it does not show · which of species, age, post-mortem state and TMT design drives the shift; kidney and adipose are untested; n is small (13 donors) and the 2 cortex samples cannot support a per-tissue statement. Panel genes failing in the target at k20: Ctnna1, Diaph1, Bzw2, Idi1, Arhgap1, Oxct1, Pdpr, Stim1, Hip1r (`panel_gene_check.csv`).
- files · `results_multiomic/03_prot_transfer/README.md`.
- sensitivity, same scale both sides (`results_multiomic/03_prot_transfer/rawppm/`): Jiang raw reporter intensities re-normalised like the RII (channel total → log2 ppm, technical replicates kept), 2414 genes, 94 mapped samples / 13 donors: accuracy k20 **0.734**, k50 0.649, full 0.798; coverage with MoTrPAC calibration 0.085 (empty 0.915), recalibrated on 5 donors 0.923 (set size 1.538).
- secondary atlas (Wang 2019, 29 human tissues, one donor each, label-free; `results_multiomic/03_prot_transfer/wang2019/`): 3037 genes; 6 mapped tissue samples, 23 OOD; accuracy k20 0.500, k50 0.667, full 0.667; coverage with MoTrPAC calibration at k20 0.333 (empty 0.667); OOD empty-set fraction 1.000.
- secondary atlas (Geiger 2013, 28 mouse tissues, pooled mice, SILAC ratios, gene-symbol match; `results_multiomic/03_prot_transfer/geiger2013/`): 2151 genes; 9 mapped tissue samples, 20 OOD; accuracy k20 0.667, k50 0.889, full 1.000; coverage with MoTrPAC calibration at k20 0.111 (empty 0.889); OOD empty-set fraction 1.000.

## Phase 4 — metabolite fingerprint transfer

- question · does a metabolite panel selected on MoTrPAC transfer to an external rodent tissue metabolome (pre-registration e), and does it survive exercise?
- data · `hilic_mw`: source 876 vials / 54 animals / 19 tissues / 129 metabolites → target 840 samples / 70 individuals / 12 tissues, **25 RefMet-matched** metabolites (STOPPED: < 30); `hilic_sato`: source 876 vials / 54 animals / 19 tissues / 129 metabolites → target 191 samples / 24 individuals / 8 tissues, **58 RefMet-matched** metabolites; `deep_sato`: source 451 vials / 54 animals / 9 tissues / 340 metabolites → target 191 samples / 24 individuals / 8 tissues, **83 RefMet-matched** metabolites; `deep_mw`: source 451 vials / 54 animals / 9 tissues / 340 metabolites → target 840 samples / 70 individuals / 12 tissues, **44 RefMet-matched** metabolites (`results_multiomic/04_metab_transfer/legs_summary.csv`, `<leg>/feature_overlap.csv`).
- design · `tfp.transfer` unchanged (panels on all MoTrPAC animals, z-scores within dataset, super-classes, α = 0.10 conformal sets calibrated on held-out MoTrPAC animals, recalibration on 3 and 5 target animals, every split grouped on the mouse); MW ST003188 scored per age (primary 3 months); Sato scored per treatment (primary Sedentary) plus the native sedentary → exercised invariance test.
- result · `hilic_sato`: accuracy over 191 mapped samples (24 individuals) k20 **0.277** (animal-bootstrap 95 % CI 0.26–0.30), full 0.445 (0.40–0.48) (chance 0.053); coverage with MoTrPAC calibration at k20 (Sedentary) 0.052 (empty 0.948); recalibrated on 3 / 5 animals 0.925 / 0.899 (set size 11.865); native sedentary→exercised: accuracy k20 0.947 coverage 0.905, full 1.000 / 0.926 (RNA controls→trained 0.961 / 0.903); `deep_sato`: accuracy over 167 mapped samples (24 individuals) k20 **0.329** (animal-bootstrap 95 % CI 0.27–0.39), full 0.527 (0.50–0.55) (chance 0.111); coverage with MoTrPAC calibration at k20 (Sedentary) 0.000 (empty 1.000); recalibrated on 3 / 5 animals 0.910 / 0.910 (set size 5.277); native sedentary→exercised: accuracy k20 0.947 coverage 0.905, full 1.000 / 0.926 (RNA controls→trained 0.961 / 0.903); `deep_mw`: accuracy over 490 mapped samples (70 individuals) k20 **0.637** (animal-bootstrap 95 % CI 0.61–0.66), full 0.753 (0.72–0.78) (chance 0.111); coverage with MoTrPAC calibration at k20 (3) 0.153 (empty 0.847); recalibrated on 3 / 5 animals 0.924 / 0.920 (set size 2.233). Pre-registration (e) (accuracy ≥ 2 × chance with ≥ 30 matched metabolites): hilic_sato PASS, deep_sato PASS, deep_mw PASS.
- what it does not show · platform and chemistry differ (MoTrPAC HILIC+ / six platforms vs RP-negative triple-quad vs Metabolon HD4), so a RefMet name match is a name match, not an identical analyte measurement; serum vs plasma, whole brain vs three rat brain regions, quadriceps vs gastrocnemius are imperfect maps; sample centring on different metabolite sets adds an offset the z-scoring only partly removes.
- files · `results_multiomic/04_metab_transfer/<leg>/README.md`.

## Phase 5 — fusion judged by transfer

- question · is a two-layer fingerprint more robust under the species shift than either layer alone — accuracy AND coverage AND empty-set fraction on the same human samples — and do the layers confuse the same pairs?
- data · Jiang 2020: 182 samples with matched RNA and protein, 42 from 12 donors in the 5 mapped classes; MoTrPAC 7-tissue source, RNA 13800 / protein 2731 matched genes (`results_multiomic/05_fusion_transfer/overlap.csv`).
- design · both layers fit on the 7 proteomics tissues (one label space); late fusion = mean probability; stacked LR fit on animal-grouped out-of-fold MoTrPAC probabilities; one calibration set (30 % of animals, held out of both layers) for every model; every Jiang statement on donors.
- result · at k20 on the same 42 samples (accuracy with donor-bootstrap 95 % CI): RNA accuracy 1.000 (1.00–1.00) / coverage 0.452 / empty 0.548; protein 0.476 (0.39–0.57) / 0.071 / 0.929; late mean 0.738 (0.61–0.86) / 0.048 / 0.952; stacked LR 0.738 (0.61–0.86) / 0.048 / 0.952. Full models: RNA 1.000 / 0.024, protein 0.524 / 0.048, late mean 1.000 / 0.048, stacked 1.000 / 0.048. Pre-registered robustness rule at k20 (≥ RNA on accuracy and coverage, ≤ on empty sets): late mean FAIL, stacked FAIL. Confusion structure on Jiang (k20): off-diagonal Pearson hard n/a, soft 0.259; per-tissue accuracy correlation across layers n/a; locally (MoTrPAC out-of-fold, k20) soft 0.618. The 19-class RNA fingerprint on the same Jiang RNA: k20 0.526, full 0.816 over 76 mapped samples in 13 tissues (`rna19_on_jiang.csv`).
- what it does not show · a precise fusion benefit: 42 samples, 12 donors, two classes with ≤ 5 samples; coverage under MoTrPAC calibration is near zero for every model, so 'coverage ≥' is a comparison of collapses. The protein arm uses the cleaned relative scale (conservative, see Phase 3).
- files · `results_multiomic/05_fusion_transfer/README.md`.

## Phase 6 — identifiability of the external designs

- question · does any external dataset have tissue crossed with its processing batch (pre-registration f), and what does that say about the audit's conclusion?
- data · sample metadata of Jiang 2020 (448 samples), Wang 2019, Sato 2022, MW ST003188; MoTrPAC layers from `results_frozen/16_identifiability/` (`results_multiomic/06_external_identifiability/design_comparison.csv`, `nesting_all.csv`, `estimable_pairs.csv`).
- design · `tfp.batch.nesting_table` and `estimable_pairs` unchanged: levels per tissue, tissues per level, Cramér's V with tissue, tissue pairs sharing a level of every batch variable.
- result · **(f) PASS**: MoTrPAC IMMUNO — plate_id: 20 levels, up to 4 tissues per level, Cramér's V 0.612, 16 of 136 pairs estimable; Jiang2020 — tmt_run: 56 levels, up to 8 tissues per level, Cramér's V 0.324, 424 of 528 pairs estimable; Sato2022 — round: 6 levels, up to 7 tissues per level, Cramér's V 0.447, 21 of 28 pairs estimable. Nested designs: MoTrPAC TRNSCRPT (RNA_extr_plate_ID: V 1.000, 1/171); MoTrPAC METHYL (DNA_extr_plate_ID: V 1.000, 0/28); MoTrPAC ATAC (Nuclei_extr_date: V 0.918, 0/28); MoTrPAC PROT (plex_id: V 1.000, 0/21); MoTrPAC PHOSPHO (plex_id: V 1.000, 0/21); MoTrPAC ACETYL (plex_id: V 1.000, 0/1); MoTrPAC UBIQ (plex_id: V 1.000, 0/1); Wang2019 (ms_experiment: V 1.000, 0/465); MW_ST003188 (batch: V 1.000, 0/66). MoTrPAC IMMUNO (Luminex plates holding several tissues) was already the audit's one crossed layer; Sato's ROUND is a per-tissue-table label and is listed, not claimed.
- what it does not show · that MoTrPAC's fingerprint is or is not batch; only that the confound is a property of MoTrPAC's sample-to-plex allocation, not of TMT proteomics. Sato's ROUND labels repeat across tissue tables like MoTrPAC's S1–S6 and are not evidence of shared runs; Geiger has no per-sample batch variable.
- files · `results_multiomic/06_external_identifiability/README.md`.

## Phase 7 — synthesis

_pending_

## What this changes in the submission, if anything

Quoted sentences are from the submission as found (`README.md`, `site/limitations.html`); the replacement is what the overnight results support. **Nothing in the submission was edited.** Numbers: the CSV named in each item.

1. `README.md` (Scope): *"proteomics and metabolomics evaluated for within-tissue use only."*
   → *"proteomics evaluated cross-tissue on the portal reporter-ion scale (tissue R² of PC1 0.991 vs 0.0009 on the distributed ratios; `results_multiomic/01_rii/variance_partition.csv`) and transferred to a human proteome atlas (k20 accuracy 0.455 on 44 samples / 13 donors, coverage 0.068 under MoTrPAC calibration; `results_multiomic/03_prot_transfer/`); metabolomics transferred to two external mouse metabolomes by RefMet name (deep-platform panel: k20 0.637 on 490 samples / 70 mice; `results_multiomic/04_metab_transfer/legs_summary.csv`). Within-study, tissue is still nested in plex for both."*

2. `site/limitations.html`: *"Proteomics is within-tissue only (ratios to per-tissue reference pools; missingness identifies the tissue), so no cross-tissue protein fingerprint is claimed."*
   → *"The distributed proteomics (ratios to per-tissue reference pools) carries no cross-tissue axis; the portal reporter-ion intensities do (tissue R² of PC1 0.991, identical on the 2393 proteins with no missing value), and a protein panel selected on them transfers to Jiang 2020 with the RNA pattern: accuracy above chance, coverage collapse (0.068) under MoTrPAC calibration, recovery to 0.917 after recalibration on 5 donors. Missingness still identifies the tissue on that scale (1.000 accuracy from the NaN pattern alone), and plex is nested in tissue."*

3. `README.md` (Roadmap): *"a cross-tissue proteomics fingerprint;"*
   → done on the RII scale (branch `multiomic-overnight`, `results_multiomic/01_rii/`, `03_prot_transfer/`); remaining: a rat multi-tissue proteome as the BodyMap-equivalent (none exists in PRIDE/ProteomeXchange, `results_multiomic/02_discovery/attempts.csv`).

4. `README.md` (numbers table): *"Sedentary vs 8-week-trained within tissue: mean best single-omic AUROC / fusion beats single / attributable to training | 0.994 / 0 of 7 / 4 of 7"*
   → add a row: *"Fusion judged by transfer (Jiang 2020, 42 samples with matched RNA + protein, 12 donors): k20 accuracy RNA 1.000 / protein 0.476 / late fusion 0.738; coverage 0.452 / 0.071 / 0.048 — fusion is not more robust than RNA alone under the species shift"* (`results_multiomic/05_fusion_transfer/fusion_transfer_summary.csv`).

5. `site/identifiability.html` / README ("Estimable tissue pairs within study (RNA-seq) | 1 of 171"): add the external counterexample —
   *"In Jiang 2020's TMT design each run holds up to 8 tissues (56 runs; Cramér's V run × tissue 0.324; 424 of 528 tissue pairs estimable within a run) versus MoTrPAC PROT (V 1.000, 0 of 21): nesting is a choice of design, not a property of the assay."* (`results_multiomic/06_external_identifiability/design_comparison.csv`).

Not changed: every RNA number on the site and in the README; the audit's conclusion for MoTrPAC itself (still nested); the human-transfer numbers (GTEx).

## For the talk

1. Put MoTrPAC's proteomics back on the reporter-ion scale and tissue becomes its dominant axis (R² 0.991 of PC1 against 0.0009 on the distributed ratios; permutation null 0.029), so the layer we had to leave out of the fingerprint was a normalisation choice, not a property of the proteome (`results_multiomic/01_rii/variance_partition.csv`).
2. A 20-protein panel selected on those intensities names the tissue of 0.45 of 44 human TMT samples from 13 GTEx donors (chance 0.14; 0.73 when both sides are processed the same way), and its 90 % prediction sets cover only 0.07 with MoTrPAC calibration but 0.92 after recalibrating on five donors — exactly the RNA pattern (`results_multiomic/03_prot_transfer/`).
3. The audit's confound is a design choice, not the assay: in Jiang 2020 each TMT run holds up to 8 tissues and 424 of 528 tissue pairs are estimable within a run, versus 0 of 21 in MoTrPAC's one-tissue-per-plex design (`results_multiomic/06_external_identifiability/design_comparison.csv`).

## Not done / not possible

- **Metabolite transfer HILIC+ → MW ST003188 (`hilic_mw`)**: stopped under the pre-registered rule — only 25 RefMet names shared between the 129 HILIC+ metabolites named in every MoTrPAC tissue and the atlas's 190 (`results_multiomic/04_metab_transfer/hilic_mw/feature_overlap.csv`). The deep-platform source covered the atlas instead.
- **Wang 2019 as a full target**: the PRIDE MaxQuant bundle (28.2 GB) exceeds the 20 GB box and was not downloaded; the Europe PMC EV tables (gene-level intensities, one donor per tissue) were used, so Wang contributes 6 mapped tissue samples and no donor-level statement (`results_multiomic/03_prot_transfer/wang2019/`).
- **Geiger 2013**: matched by gene symbol across species (no rat–mouse orthology table in the repo) and one pooled sample per tissue; treated as a secondary check only (`results_multiomic/03_prot_transfer/geiger2013/`).
- **A rat multi-tissue proteome (the BodyMap-equivalent)**: none with processed tables exists on PRIDE or ProteomeXchange (searches logged in `results_multiomic/02_discovery/attempts.csv`); the reverse-direction test therefore uses a human atlas.
- **MetaboLights** was reachable but not searched: two matchable rodent metabolomes were already in hand within the 2 h discovery box.
- **Sato 2022 batch variable**: ROUND labels repeat across the per-tissue Metabolon tables and could not be verified as shared runs; listed in Phase 6, not claimed.
- **PXD082651** (2026 mouse lifespan multi-tissue atlas of non-canonical peptides) not pursued (non-canonical peptide focus).
- **Portal release c2.0 RII, acetyl and ubiquityl RII**: not run in the main line (c1.0 `prot-pr` and `prot-ph` only) unless a later log entry says otherwise.
- **Figures**: described in the site draft with their source CSVs; no image files were rendered and `site/` was not touched.
- **Kidney and adipose** have no human protein target in Jiang 2020; they are covered only by the mouse atlas (Geiger, n = 1 per tissue).
- **Download log** (`results_multiomic/02_discovery/download_log.csv`): 19 files, 237 MB — geiger2013 1 files, 10 MB; jiang2020 8 files, 169 MB; mw_ST003188 3 files, 6 MB; sato2022 7 files, 53 MB. 15 attempts recorded (`attempts.csv`), including the failures (raw-only PRIDE project, bot-walled publisher page, over-box bundle, no rat atlas, Sato absent from Metabolomics Workbench).

## Draft for site/_drafts/multiomic.html (Markdown only; the site is untouched)

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

RUN COMPLETE 2026-09-27T07:43:49Z
