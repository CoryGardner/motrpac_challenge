# Multiomic overnight report — can the proteome and metabolome carry a tissue fingerprint that transfers?

Branch `multiomic-overnight`. Pre-registration: `docs/PREREGISTRATION_MULTIOMIC.md`. Log: `docs/MULTIOMIC_LOG.md`.
Every number in this file is read from a CSV under `results_multiomic/` by `scripts/multiomic/build_report.py`; the file is named beside each number. Within-study accuracy is context, never a finding (the audit: batch is nested in tissue).

## Significant findings so far (2026-09-27 07:26 UTC; phases with results: 0, 1, 2, 3)

1. **RII tissue-axis recovery.** On reporter-ion intensities normalised to the channel total, tissue explains R² = 0.991 of PC1 (PC2 0.998; label-permutation null 95th pct 0.029) versus 0.0009 on the distributed ratio matrix with the same code; n = 420 vials, 60 animals, 3637 proteins in every tissue; identical on the 2393 proteins with no missing value (R² 0.990). Within-study: plex is nested in tissue. — `results_multiomic/01_rii/variance_partition.csv`, `variance_partition_ratio.csv`, `variance_partition_complete.csv`
2. **RNA panel markers hold at the protein level.** Of the 20 RNA panel genes whose marker tissue is one of the 7 proteomics tissues and whose protein is quantified, 17 (85 %) have the same marker tissue at the protein level (chance 1/7); over all 3701 genes with both layers, the cross-tissue RNA–protein Spearman has median 0.679 (IQR 0.429–0.821; mismatched-pair null median 0.286, 95th pct 0.857); n = 50 shared animals per tissue. — `results_multiomic/01_rii/rna_protein_panel_summary.csv`, `rna_protein_correlation_summary.csv`
3. **Protein fingerprint transfers to a human proteome atlas, with the RNA coverage pattern.** A 20-protein panel selected on MoTrPAC RII names the tissue of 0.455 of 44 Jiang 2020 samples (13 donors, 5 mapped tissues; k50 0.477, full 0.523; chance 0.143); marginal conformal coverage with MoTrPAC calibration is 0.068 (empty sets 0.932) against a nominal 0.90, and recalibration on 5 donors restores 0.917 at set size 2.954. Human adults, post-mortem, different TMT design — all shifts at once. — `results_multiomic/03_prot_transfer/accuracy_overall.csv`, `conformal_transfer.csv`, `recalibration.csv`

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
- result · accuracy over the 44 mapped samples: full **0.523**, k20 **0.455**, k50 **0.477** (chance 0.143). Marginal coverage with MoTrPAC calibration at k20: **0.068** (empty sets 0.932, set size 0.068); floored Mondrian 0.068. Recalibration on 3 donors → 0.900 (set size 3.497), on 5 donors → 0.917 (set size 2.954). OOD tissues at k20 (marginal): empty-set fraction 0.975. Reverse direction (panel selected on Jiang, scored on MoTrPAC RII): full 0.587, k20 0.577, k50 0.573. Pre-registration (d) asks for accuracy ≥ 0.429 and coverage < 0.90 with MoTrPAC calibration: PASS at k20 (accuracy 0.455, coverage 0.068); RNA ladder for comparison in `ladder_protein_vs_rna.csv`.
- what it does not show · which of species, age, post-mortem state and TMT design drives the shift; kidney and adipose are untested; n is small (13 donors) and the 2 cortex samples cannot support a per-tissue statement. Panel genes failing in the target at k20: Ctnna1, Diaph1, Bzw2, Idi1, Arhgap1, Oxct1, Pdpr, Stim1, Hip1r (`panel_gene_check.csv`).
- files · `results_multiomic/03_prot_transfer/README.md`.

## Phase 4 — metabolite fingerprint transfer

_pending_

## Phase 5 — fusion judged by transfer

_pending_

## Phase 6 — identifiability of the external designs

_pending_

## Phase 7 — synthesis

_pending_

