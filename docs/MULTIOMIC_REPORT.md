# Multiomic overnight report — can the proteome and metabolome carry a tissue fingerprint that transfers?

Branch `multiomic-overnight`. Pre-registration: `docs/PREREGISTRATION_MULTIOMIC.md`. Log: `docs/MULTIOMIC_LOG.md`.
Every number in this file is read from a CSV under `results_multiomic/` by `scripts/multiomic/build_report.py`; the file is named beside each number. Within-study accuracy is context, never a finding (the audit: batch is nested in tissue).

## Significant findings so far (2026-09-27 07:09 UTC; phases with results: 0, 1)

1. **RII tissue-axis recovery.** On reporter-ion intensities normalised to the channel total, tissue explains R² = 0.991 of PC1 (PC2 0.998; label-permutation null 95th pct 0.029) versus 0.0009 on the distributed ratio matrix with the same code; n = 420 vials, 60 animals, 3637 proteins in every tissue; identical on the 2393 proteins with no missing value (R² 0.990). Within-study: plex is nested in tissue. — `results_multiomic/01_rii/variance_partition.csv`, `variance_partition_ratio.csv`, `variance_partition_complete.csv`
2. **RNA panel markers hold at the protein level.** Of the 20 RNA panel genes whose marker tissue is one of the 7 proteomics tissues and whose protein is quantified, 17 (85 %) have the same marker tissue at the protein level (chance 1/7); over all 3701 genes with both layers, the cross-tissue RNA–protein Spearman has median 0.679 (IQR 0.429–0.821; mismatched-pair null median 0.286, 95th pct 0.857); n = 50 shared animals per tissue. — `results_multiomic/01_rii/rna_protein_panel_summary.csv`, `rna_protein_correlation_summary.csv`

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

_pending_

## Phase 3 — protein fingerprint transfer

_pending_

## Phase 4 — metabolite fingerprint transfer

_pending_

## Phase 5 — fusion judged by transfer

_pending_

## Phase 6 — identifiability of the external designs

_pending_

## Phase 7 — synthesis

_pending_

