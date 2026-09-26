# GTEX_TRANSFER — the human leg (teammate's side)

Purpose: test whether a compact rat tissue fingerprint survives the jump to human tissue.
Two human sources, in order of availability:

1. **GTEx** (public, bulk RNA-seq, post-mortem). Tissues that match the MoTrPAC rat set:
   Whole Blood, Muscle - Skeletal, Adipose - Subcutaneous, Adipose - Visceral (Omentum),
   Heart - Left Ventricle, Liver, Kidney - Cortex, Lung, Brain - Cortex, Brain - Hippocampus,
   Brain - Hypothalamus, Colon - Transverse, Small Intestine - Terminal Ileum, Spleen,
   Adrenal Gland, Ovary, Testis, Artery - Aorta (closest to vena cava, imperfect).
2. **MoTrPAC human** (blood, vastus lateralis, subcutaneous adipose; acute endurance/resistance
   bouts; preprints March 2026). Access terms not verified — check https://motrpac-data.org/
   and ask organizers. If available, it is the better test because it is the same consortium's
   pipelines and the exercise-state question becomes possible.

## GTEx download (done 2026-09-17: v8 gene TPM GCT, 1.63 GB, and the v8 sample attributes; 2026-09-18: v8 gene read-count GCT, 0.92 GB; all under `data/external/gtex/`; `make gtex` builds the subsets if missing and runs the transfer, the representation comparison and the units test)

- Portal: https://gtexportal.org/home/downloads/adult-gtex (registration-free for open-access
  files). Get the **gene TPM** matrix (GCT, gzipped, ~1.5 GB for v8/v10) and the
  **sample attributes** file (`*SampleAttributesDS.txt`) which maps SAMPID → SMTSD (tissue).
- `scripts/11_gtex_prepare.py` subsets to the tissues above (≤ 150 donors per tissue, seeded) and writes
  `data/external/gtex_tpm_subset.csv` (log2 TPM+1), `gtex_meta.csv`, `gtex_gene_symbols.csv` and, from the
  read-count GCT (`--reads`), `gtex_cpm_subset.csv`: log2 CPM+1 with the library size = the sample's read
  total over every gene in the GCT, i.e. the unit of the MoTrPAC and BodyMap count matrices. Existing
  outputs are not rebuilt without `--force`. File names change per GTEx release; pass them explicitly.
- Sample IDs `GTEX-XXXX-...` share a donor prefix `GTEX-XXXX`; group splits on the donor,
  same rule as `pid` in rat.

## Transfer protocol (implemented in `scripts/13_gtex_transfer.py` on the shared `motrpac.transfer` code path; results in `results/13_gtex/`; representations and the transferability-aware selector in `scripts/14_transfer_representations.py`, `results/14_transfer/`)

Decisions: Muscle - Skeletal is scored as the {SKM-GN, SKM-VL} super-class (GTEx muscle is gastrocnemius, but the rat
muscles are the confusable pair); Artery - Aorta → VENACV is kept but read as a caveat; orthologs are 1:1 pairs from
`rat_to_human_gene.csv` matched to GTEx by human Ensembl id with symbol fallback (`motrpac.transfer.match_gtex_orthologs`,
shared by scripts 13 and 14 so that both quote the same gene set); panels are re-selected in the ortholog space on all
50 rat animals; every GTEx split is grouped on the donor.

**Accuracy convention.** Every GTEx (and BodyMap) accuracy in `results/SUMMARY.md`, `../../docs/findings/QUESTIONS_FOR_ORGANIZERS.md`
and `results/ABSTRACT.md` is *sample-weighted*: the mean over all mapped samples, each sample counting once. With
60–150 samples per tissue it is within 0.01 of the macro mean over tissues, which is kept in
`results/13_gtex/accuracy_overall.csv` and in the `accuracy_macro_over_tissues` column of
`results/14_transfer*/target_summary.csv` for anyone who prefers it.

**Units test (script 14, `results/14_transfer_cpm/`).** The within-sample rank and top-scoring-pair representations
transfer within species (BodyMap) but not to GTEx. Because the portal matrix is TPM (length-normalized) while MoTrPAC
and BodyMap are CPM, the ordering of genes within a sample differs by gene length between the two units. The test
re-scores GTEx on log2 CPM from the read counts with the same panels; if ranks and pairs recover to the z-score
representation the failure was units, if they stay where they were it is biology (the gene order differs between
species). The answer is in the "Units test" paragraph of the representations block of `results/SUMMARY.md`.

**Conformal variants on the target.** `conformal_transfer` reports marginal, Mondrian and marginal-floor Mondrian sets
(each class quantile floored at the marginal one, so the sets contain both the marginal and the Mondrian sets);
`results/13_gtex/coverage_by_tissue.csv` has all three per GTEx tissue.

## Transfer protocol (as planned)

1. Fit the rat panel (Ensembl rat gene IDs) on all rat animals; map to human via
   `data/raw/rat_to_human_gene.csv`; keep 1:1 orthologs only, report losses.
2. Re-standardize within species (per-gene z-score using each species' own data) — this is
   the minimal, honest correction; anything cleverer must be justified.
3. Score GTEx samples with the rat-trained classifier restricted to the panel genes; report
   accuracy per human tissue and a confusion matrix; then **re-calibrate conformal sets on a
   GTEx calibration split** (donor-grouped) and report coverage on GTEx test donors, versus
   coverage when using the rat calibration (expected to break — that's the point).
4. Compare with a GTEx-only panel of the same size selected under donor-grouped CV: how much
   is lost by using the rat panel?

## Known confounds to state up front

- Post-mortem ischemic time in GTEx; RIN differences; platform (poly-A vs total RNA).
- Blood cell composition differs between rat whole blood and GTEx whole blood.
- "Muscle - Skeletal" in GTEx is gastrocnemius, so it maps to SKM-GN, not SKM-VL.
