# Competition compliance

Stanford Bioinformatics Center / MoTrPAC Hackathon, 25–27 September 2026. Track: **Molecular Tissue
Fingerprints — Can molecular signatures identify a tissue reliably?**

| Rule | How this submission meets it | Where to check |
|---|---|---|
| At least one MoTrPAC component | The dataset: MoTrPAC 6-month rat endurance training, portal release c1.0 (= `MotrpacRatTraining6moData` 2.0.0): 899 transcriptome vials of 50 animals in 19 tissues for the fingerprint; proteomics, phosphoproteomics, acetyl- and ubiquitylome, metabolomics, immunoassay, ATAC and RRBS metadata for the identifiability audit. | `docs/DATA_GUIDE.md`, `results/02_inventory/`, `scripts/16_identifiability.py` |
| At least one CFDE component | (1) **GTEx** (a CFDE Data Coordinating Center): GTEx v8 open-access expression, 2,485 samples, 862 donors, 17 tissues, used for the cross-species transfer (`scripts/11_gtex_prepare.py`, `scripts/13_gtex_transfer.py`, `results/13_gtex/`). (2) **Metabolomics Workbench** (a CFDE DCC): holds the MoTrPAC rat metabolomics as project PR001020; the metabolomics layer is part of the fusion and identifiability results (`results/07_fusion/`, `results/16_identifiability/layers.json`). | `site/about.html`, `site/transfer.html` |
| Open-source repository | This repository, MIT (`LICENSE`); the site is static HTML/JS in `site/` and deploys unchanged to GitHub Pages (`.github/workflows/pages.yml`). Raw data are not redistributed; the export and download instructions are (`R/`, `docs/GTEX_TRANSFER.md`, `docs/DATA_GUIDE.md` §9). | `README.md` |
| Data used as licensed | MoTrPAC: the `MotrpacRatTraining6moData` package is MIT-licensed and its data carry the consortium's data-use terms for public release c1.0. Rat BodyMap (GEO GSE53960): the Bioconductor package `bodymapRat` 1.28.0, **CC BY 4.0**. GTEx v8 open-access: the GTEx data-use policy (https://gtexportal.org/home/license). No restricted data are committed; `LICENSE` states each source's terms. | `LICENSE`, `site/about.html` |
| Reproducibility and integrity | Pinned environments: `requirements.txt` (exact versions of the reference run), `environment.yml`, `environment-r.yml`; `R/install_deps.R` pins the data package to GitHub commit `f831a4f`. `results_frozen/` is a committed snapshot of the 150 derived result files the site reads, with a SHA-256 `MANIFEST.json` (`make verify-frozen`). CI (`.github/workflows/tests.yml`) runs the tests, checks the snapshot against its manifest, re-exports the site data from the snapshot and compares the export with the committed `site/data/`, then checks the site's links — all without the raw data. | `results_frozen/MANIFEST.json`, `.github/workflows/tests.yml`, `CHANGELOG.md` |
| Judging criteria | Scientific impact: the three-part answer (accuracy, guarantee, identifiability). Technical quality: animal-grouped splits, in-fold pipelines, tuned baselines, conformal sets with the textbook quantile, provenance for every number, tests. Presentation: the site, the explorer, the five-minute tour, `figures/summary_figure.png`. | `docs/EVALUATION_RULES.md`, `tests/`, `site/` |

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
