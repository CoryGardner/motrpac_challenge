# Competition compliance

Stanford Bioinformatics Center / MoTrPAC Hackathon, 25–27 September 2026. Track: **Molecular Tissue
Fingerprints — Can molecular signatures identify a tissue reliably?**

| Rule | How this submission meets it | Where to check |
|---|---|---|
| At least one MoTrPAC component | The dataset: MoTrPAC 6-month rat endurance training, portal release c1.0 (= `MotrpacRatTraining6moData` 2.0.0): 899 transcriptome vials of 50 animals in 19 tissues for the fingerprint; proteomics, phosphoproteomics, acetyl- and ubiquitylome, metabolomics, immunoassay, ATAC and RRBS metadata for the identifiability audit. | `docs/DATA_GUIDE.md`, `results/02_inventory/`, `scripts/16_identifiability.py` |
| At least one CFDE component | (1) **GTEx** (a CFDE Data Coordinating Center): GTEx v8 open-access expression, 2,485 samples, 862 donors, 17 tissues, used for the cross-species transfer (`scripts/11_gtex_prepare.py`, `scripts/13_gtex_transfer.py`, `results/13_gtex/`). (2) **Metabolomics Workbench** (a CFDE DCC): holds the MoTrPAC rat metabolomics as project PR001020; the metabolomics layer is part of the fusion and identifiability results (`results/07_fusion/`, `results/16_identifiability/layers.json`). The multiomic follow-up adds two more CFDE-linked sources: (3) **Metabolomics Workbench** study ST003188, "A metabolic atlas of mouse aging" (doi:10.21228/M88J0W, CC BY 4.0), the external target of the metabolite transfer (`scripts/multiomic/04_metab_transfer.py`, `results_multiomic/04_metab_transfer/deep_mw/`); (4) the **GTEx proteome**, Jiang et al., *Cell* 183:269 (2020), doi:10.1016/j.cell.2020.08.036, PRIDE PXD016999: TMT proteomics of GTEx donors' tissues, the external target of the protein transfer (`scripts/multiomic/03_prot_transfer.py`, `results_multiomic/03_prot_transfer/`). | `site/about.html`, `site/transfer.html`, `site/multiomic.html` |
| Open-source repository | This repository, MIT (`LICENSE`); the site is static HTML/JS in `site/` and deploys unchanged to GitHub Pages (`.github/workflows/pages.yml`). The full raw datasets are not redistributed (the small subsets that are committed are listed with their licences in `NOTICE.md`); the export and download instructions are (`R/`, `docs/GTEX_TRANSFER.md`, `docs/DATA_GUIDE.md` §9, `docs/REPRODUCIBILITY.md`). | `README.md`, `docs/REPRODUCIBILITY.md` |
| Data used as licensed | MoTrPAC: the `MotrpacRatTraining6moData` package is MIT-licensed and its data carry the consortium's data-use terms for public release c1.0. Rat BodyMap (GEO GSE53960): the Bioconductor package `bodymapRat` 1.28.0, **CC BY 4.0**. GTEx v8 open-access: the GTEx data-use policy (https://gtexportal.org/home/license). The repository contains small subsets and transformed values of these data (a four-sample BodyMap count fixture; log2 CPM of at most 100 genes and per-sample predictions for MoTrPAC, BodyMap and GTEx) and derived tables of the multiomic sources (Jiang 2020, Wang 2019, Geiger 2013, Sato 2022, Metabolomics Workbench ST003188, MoTrPAC reporter-ion files); `NOTICE.md` lists each item with its source, licence and what was changed. No licence-restricted data (DisGeNET, OMIM, PhosphoSitePlus, LINCS vectors) are committed. | `NOTICE.md`, `LICENSE`, `site/about.html` |
| Reproducibility and integrity | Pinned environments: `requirements.txt` (exact versions of the reference run), `environment.yml`, `environment-r.yml`; `R/install_deps.R` pins the data package to GitHub commit `f831a4f`. `results_frozen/` is a committed snapshot of the 179 derived result files the site reads, with a SHA-256 `MANIFEST.json` (`make verify-frozen`). CI (`.github/workflows/tests.yml`) runs the tests, checks the snapshot against its manifest, re-exports the site data from the snapshot and compares the export with the committed `site/data/`, then checks the site's links — all without the raw data. | `results_frozen/MANIFEST.json`, `.github/workflows/tests.yml`, `CHANGELOG.md` |
| Judging criteria | Scientific impact: the three-part answer (accuracy, guarantee, identifiability). Technical quality: animal-grouped splits, in-fold pipelines, tuned baselines, conformal sets with the textbook quantile, provenance for every number, tests. Presentation: the site (Check samples, the Reference atlas, The science), the five-minute tour, `figures/summary_figure.png`. | `docs/EVALUATION_RULES.md`, `tests/`, `site/` |

## What the site guarantees about its numbers

- Every number shown comes from `results/` (or, when a complete run is not present, from the committed snapshot
  `results_frozen/`) through `scripts/30_export_site_data.py`; `site/data/provenance.json` records the file, row
  selector, column and aggregation of every headline number and the source of every table.
- `tests/test_site_data.py` checks those entries against the results files; `tests/test_site_conformal.js` checks
  the browser's conformal port against Python-computed sets; `tests/test_frozen_results.py` checks the snapshot.
- A number that cannot be sourced from the results renders as "pending" with the reason rather than being typed in;
  after the final export none remains on the site.
- `docs/NUMBERS_RECONCILIATION.md` lists every headline number with its pre-fix value where the 2026-09-25
  conformal-quantile fix moved it.

## Team

The Rat PAC: Samuel Montalvo, Manasa Rapuru, Erol Evangelista and Cory Gardner (in no particular order); see
`site/about.html` and `CITATION.cff`.

AI tools assisted with code and writing; all results were verified by the team.

## Track outputs → where they are

| output the track names | what it is | where |
|---|---|---|
| A classifier | the 20-gene logistic regression with calibrated 90 % prediction sets, applied to your samples | [Check samples](https://corygardner.github.io/motrpac_challenge/) (calls, sets, claim checks, CSV and report); `site/data/panel_model.json` (the model the browser runs); [`src/tfp/models.py`](https://github.com/CoryGardner/motrpac_challenge/blob/main/src/tfp/models.py) |
| Minimal tissue-signature panel | the 20 genes and the 10-gene stable core | [Panel page](https://corygardner.github.io/motrpac_challenge/fingerprint.html); `site/data/panel_card.csv`, `site/data/panel_card.json` |
| Feature-selection workflow | the class-aware round-robin selector, fitted inside animal-grouped folds | [Methods: the selector](https://corygardner.github.io/motrpac_challenge/methods.html#selector); [Reference atlas: panel builder](https://corygardner.github.io/motrpac_challenge/explore.html#panel-builder); [`RoundRobinSelector`](https://github.com/CoryGardner/motrpac_challenge/blob/main/src/tfp/models.py#L58) in `src/tfp/models.py` |
| Interactive model-explanation tool | per sample: probabilities against the calibrated threshold, per-gene contributions ("why X, not Y"), gene values against the reference tissues, a reference map | [Check samples](https://corygardner.github.io/motrpac_challenge/) (the sample drawer); [Reference atlas](https://corygardner.github.io/motrpac_challenge/explore.html): [tissue card](https://corygardner.github.io/motrpac_challenge/explore.html#tissue-card), [gene explorer](https://corygardner.github.io/motrpac_challenge/explore.html#gene-explorer), [panel builder](https://corygardner.github.io/motrpac_challenge/explore.html#panel-builder) |

The same mapping is the "What this submission delivers" section of The science page (`site/science.html#deliverables`); the product itself is `site/index.html` (Check samples).
