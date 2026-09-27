# Pre-registration — overnight multiomic run (branch `multiomic-overnight`)

Written 2026-09-27 (Phase 0), before any RII value, external table or metabolite matrix was analysed. The prompt is
`prompts/overnight.md`; the log is `docs/MULTIOMIC_LOG.md`; the report is `docs/MULTIOMIC_REPORT.md`. Everything below
is fixed now and is not changed after the data are seen. Where a prediction fails, the failure is reported as such.

## Question

Can the proteome and the metabolome carry a tissue fingerprint that transfers? The submission's RNA fingerprint transfers
to rat BodyMap and to GTEx with a coverage collapse under rat calibration; the proteome as distributed (TMT ratios to a
per-tissue reference) has no cross-tissue axis (tissue R² of PC1 = 0.0008, PC2 = 0.0039 on the ratio matrix,
`results_frozen/03_eda/variance_partition_PROT.csv`); the metabolome has no external check.

## Predictions (fixed before the data are touched)

| # | Prediction | Pass rule | Where it is tested |
|---|---|---|---|
| (a) | RII-normalised MoTrPAC proteomics recovers tissue as the dominant axis | tissue R² of PC1 on the stacked RII matrix > 0.5 (ratios: 0.0008 on PC1, 0.0039 on PC2) | Phase 1, `results_multiomic/01_rii/variance_partition.csv` |
| (b) | A protein panel selected on RII classifies the 7 proteomics tissues | mean balanced accuracy ≥ 0.95 under the pipeline's animal-grouped 5-fold CV at some k ≤ 100 of the standard grid; the missingness-only baseline and the per-tissue-mean-removed classifier are reported beside it | Phase 1, `results_multiomic/01_rii/panel_curve.csv`, `diagnostic_accuracy.csv` |
| (c) | Proteins of the RNA panel genes that are quantified show tissue specificity at the protein level with the same marker tissue | among RNA panel genes whose marker tissue is one of the 7 proteomics tissues and whose protein is quantified, the protein's highest tissue mean is the RNA marker tissue in ≥ 70 % of cases | Phase 1, `results_multiomic/01_rii/rna_protein_panel_genes.csv` |
| (d) | An RII-selected protein panel transfers to a human proteome atlas | accuracy well above chance (≥ 3 × 1/n_classes, sample-weighted over mapped tissues, every split grouped on donor) AND marginal conformal coverage at α = 0.10 with MoTrPAC calibration < 0.90 (the RNA pattern: coverage collapses under shift), with recalibration on 3 and 5 donors restoring coverage toward 0.90 | Phase 3, `results_multiomic/03_prot_transfer/` |
| (e) | A metabolite panel selected on MoTrPAC HILIC+ transfers to an external rodent tissue metabolome | accuracy above chance (≥ 2 × 1/n_classes, sample-weighted over mapped tissues, grouped on animal) with ≥ 30 RefMet-matched named metabolites; the leg stops if the overlap is < 30 | Phase 4, `results_multiomic/04_metab_transfer/` |
| (f) | At least one external atlas has tissue crossed with its batch variable | for at least one obtained dataset and at least one batch variable (TMT plex, MS run/batch, acquisition date, plate), Cramér's V with tissue < 1 and > 0 estimable tissue pairs | Phase 6, `results_multiomic/06_external_identifiability/` |

Secondary, exploratory (no pass rule, reported as observed): the cross-tissue RNA–protein Spearman correlation per gene on the
same animals (Phase 1); fusion under shift (Phase 5): late fusion is called "more robust" only if its accuracy AND coverage
are both at least as high as the better single layer on the same target samples and its empty-set fraction is no higher.

## Fixed analysis choices

- Unit of every split: the animal (`pid`) in MoTrPAC, the donor or animal in every external dataset. `tfp.splits.grouped_kfold`,
  5 folds, seed 20260925 (`tfp.config.SEED`), unchanged.
- RII processing (Phase 1): peptide-level reporter-ion intensities from the c1.0 quant-id `prot-pr` folders; contaminants
  dropped; per plex (tissue × S1–S6, from the vial-metadata files) the reference channel is dropped, peptides are summed
  to their `protein_id`, proteins with < 2 quantified peptides in the plex are dropped, each channel is divided by its
  column total (× 10⁶ → parts per million of the channel's protein signal) and log2-transformed; plexes are combined
  within a tissue and tissues stacked. Main matrix: proteins present in every tissue (inner join); the outer join is kept
  for the missingness-only baseline.
- Panel selection: `tfp.models.RoundRobinSelector` inside `tfp.models.panel_curve`, classifier `logreg_l2`, standard grid
  1–100, variance prefilter 5000, inner grouped grid search over C — the pipeline as it stands.
- Transfer (Phases 3–5): `tfp.transfer.PanelModels` and `conformal_transfer` unchanged; panels selected on all MoTrPAC
  animals; z-scoring within each dataset; super-classes where a target tissue maps to several of ours; α = 0.10; LAC;
  recalibration on 3 and 5 target individuals, 20 draws.
- Metabolite matching: RefMet name first, then HMDB/KEGG/InChIKey, from `data/raw/metab_feature_id_map.csv`.
- Nulls: label-permutation nulls (200 permutations unless stated) wherever a within-study number is compared with a
  baseline; the RNA–protein correlation gets a mismatched-pair null (each RNA gene against a random other protein).
- Within-study accuracy is context, never a finding (the audit: plex is nested in tissue in the RII data too).

## Drop policy

Phases 3, 4, 5 depend on Phase 2. A dataset that cannot be matched (identifiers, units, tissues) is recorded as a
negative result. Time boxes are those of the prompt. Phases 1, 6 and 7 must exist.
