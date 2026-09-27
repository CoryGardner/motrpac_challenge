<img src="site/assets/brand/badge-192.png" width="96" height="96" align="right" alt="The Rat PAC badge: Molecular Tissue Fingerprints, MoTrPAC Hackathon 2026">

# Molecular Tissue Fingerprints

**A 20-gene panel, its conformal guarantee, and what makes it credible as biology.**

A leakage-safe, calibrated evaluation of compact tissue fingerprints in the MoTrPAC rat multi-tissue transcriptomes,
replicated in an independent laboratory's rats (rat BodyMap) and across species (human GTEx), with a direct measurement
of processing effects on the consortium's bridging standards, an exercise-specific follow-up, and a static site on which
every number carries provenance and whose Explorer scores new samples with the panel in the browser.

[![tests](https://github.com/CoryGardner/motrpac_challenge/actions/workflows/tests.yml/badge.svg)](https://github.com/CoryGardner/motrpac_challenge/actions/workflows/tests.yml)
**Live site:** https://corygardner.github.io/motrpac_challenge/ · **Release:** tag `hackathon-submission-v7` (version 1.6.0) · **Licence:** MIT

![The home page: the question, the tiles, the one-picture diagram and the transfer ladder](figures/home.png)

| | |
|---|---|
| **Team** | The Rat PAC — Erol Evangelista, Cory Gardner, Samuel Montalvo, Manasa Rapuru |
| **Event** | Stanford Bioinformatics Center / MoTrPAC Hackathon, 25–27 September 2026, track *Molecular Tissue Fingerprints* |
| **Intended users** | judges and reviewers; MoTrPAC and CFDE analysts who need a tested harness for tissue classifiers, conformal prediction sets and shift tests; anyone reusing the panels or the evaluation rules |
| **Why it matters** | a signature that transfers to independently processed data, with a calibrated guarantee, is the check that a molecular signature reads biology rather than the processing design of the study it was learned in |

## The answer

**Twenty genes are enough.** Selected inside each animal-grouped fold by a class-aware round-robin rule, a 20-gene
panel identifies 19 rat tissues at 0.976 ± 0.008 balanced accuracy (50 genes 0.993, all genes 0.995). The panel fit on
all animals names every mapped adult organ correctly in another laboratory's rats (rat BodyMap: 9 of 11 organs have a
MoTrPAC counterpart; 68 samples from 8 animals).

**The guarantee travels honestly.** Its 90 % conformal guarantee holds within the study and on trained animals (fit on
the 10 sedentary controls alone, the panel names the tissue of all 40 trained animals at 0.961 with coverage 0.903).
Beyond the study it abstains rather than errs: calibrated on MoTrPAC it covers 0.618 of the BodyMap adults and 0.364 of
human GTEx samples, and the shortfall is empty sets, not confident error. Three animals from the new laboratory restore
it (0.943 coverage at one tissue per set); three human donors restore the number, not the information (0.954 at 11.7
tissues per set), which marks the species boundary honestly.

**It is biology, not processing.** Like every large multi-tissue design, this study processed each tissue as a unit, so
within-study accuracy alone cannot say how much of a fingerprint is biology. Two things can: the external replicate,
and MoTrPAC's bridging reference pools, on which batch measured directly is about 1.6 % of the variance that separates
tissues. Training itself barely moves the fingerprint: it is a within-tissue, minor-axis signal, smaller than the tissue
contrast on every panel gene (see the Exercise page).

## Key results

Every number below is read from `results/` by `scripts/30_export_site_data.py` and carries a provenance entry in
`site/data/provenance.json`; the table is the output of `python scripts/30_export_site_data.py --readme-table`.

| result | value | source |
|---|---|---|
| 20-gene panel, balanced accuracy (19 tissues, 5 animal-grouped folds) | 0.976 ± 0.008 | `results/05_panels/TRNSCRPT/panel_curve.csv` |
| 50-gene panel / all genes | 0.993 / 0.995 | `results/05_panels/TRNSCRPT/panel_curve.csv, results/04_baselines/TRNSCRPT/summary.csv` |
| F-test selector at k = 20 (why the selector matters) | 0.399 | `results/05_panels/TRNSCRPT/panel_curve_fclassif.csv` |
| Coverage of 90 % sets in-distribution, all genes (pooled / one vial per animal) | 0.908 / 0.916 | `results/06_conformal/TRNSCRPT/coverage.csv` |
| Coverage of 90 % sets in-distribution, 20-gene panel (pooled / one vial per animal) | 0.900 / 0.924 | `results/31_site_regen/06_conformal/TRNSCRPT/scores_*.csv (recomputed)` |
| Trained animals, panel fit on the sedentary controls only: accuracy k20 / coverage | 0.961 / 0.903 | `results/08_shift/TRNSCRPT/shift_table.csv` |
| BodyMap adults (another lab): accuracy k20 / coverage / empty sets | 1.000 / 0.618 / 0.382 | `results/12_bodymap/` |
| BodyMap recalibrated on 3 animals: coverage at set size | 0.943 at 1.00 | `results/12_bodymap/recalibration.csv` |
| GTEx (human): accuracy k20 / k50 / full | 0.654 / 0.781 / 0.855 | `results/13_gtex/accuracy_overall.csv` |
| GTEx coverage k20 / empty; recalibrated on 3 donors: coverage at set size | 0.364 / 0.616; 0.954 at 11.70 | `results/13_gtex/` |
| Estimable tissue pairs within study (RNA-seq) | 1 of 171 (ovary and testes) | `results/16_identifiability/estimable_pairs.csv` |
| Sedentary vs 8-week-trained within tissue: mean best single-omic AUROC / fusion beats single / attributable to training | 0.994 / 0 of 7 / 4 of 7 | `results/07_fusion/` |
| Trained animals, 8-week group only (cohort-matched with the controls): accuracy k20 / coverage | 0.972 / 0.917 | `results/31_site_regen/08_shift_k20/scores_target_vials.csv` |
| Batch measured on a bridging reference pool run on 6 plates (Σ V_batch / Σ V_tissue, all genes) | 0.016 | `results/16_identifiability/bridge_variance.csv` |
| QC covariates alone, balanced accuracy: technical / composition / all | 0.874 / 0.952 / 0.976 | `results/16_identifiability/qc_only_summary.csv` |

![Summary figure: the transfer ladder, the ten-gene stable core, and batch nested in tissue](figures/summary_figure.png)


## Research question

**Objective.** Can a compact molecular signature identify a tissue *reliably*, where reliably means: accurate on
held-out animals, calibrated (a 90 % prediction set really covers 90 %), and transferable to independently processed data?

**Scope.** MoTrPAC 6-month rat endurance-training study, portal release c1.0 (rn6), RNA-seq of 19 tissues from the same
50 animals (899 vials, 21,193 genes after the stacked filter; both sexes; sedentary controls and 1, 2, 4 and 8 weeks
of training); proteomics and metabolomics evaluated for within-tissue use only. External: rat BodyMap (316 samples,
11 organs, 4 ages) and GTEx v8 (2,485 samples, 862 donors, 17 tissues, 14,569 one-to-one orthologs present).

**Success criteria.** Balanced accuracy ≥ 0.95 at k ≤ 20 under animal-grouped cross-validation against tuned simple
baselines; in-distribution coverage within 0.02 of the nominal 0.90; every mapped organ of the external replicate
named correctly; every number reproducible from `results/` by tests. All four are met.

## Workflow

```mermaid
flowchart TD
  subgraph data["Data (not in the repository)"]
    M["MoTrPAC rat endurance training, release c1.0<br/>R package MotrpacRatTraining6moData 2.0.0<br/>R/export_motrpac.R → data/raw/ (1.9 GB)"]
    B["Rat BodyMap GSE53960<br/>bodymapRat 1.28.0 · R/export_bodymap.R"]
    G["GTEx v8 open access (3 files)<br/>scripts/11_gtex_prepare.py"]
  end
  M --> F["Animal-grouped folds — tfp.splits<br/>5 folds · 50 animals · every vial of an animal on one side"]
  F --> P["In-fold pipeline — scripts/04, 05<br/>median impute → variance prefilter (5,000) → z-score →<br/>round-robin selector (k genes) → logistic regression (L2)"]
  P --> C["Split-conformal calibration on held-out animals — scripts/06<br/>18 fit / 22 calibrate / 10 test animals per fold · α = 0.10<br/>certificate: smallest k with error ≤ α at confidence 1 − δ"]
  C --> L["Shift ladder<br/>training state and held-out sex — scripts/08<br/>another laboratory (BodyMap) — scripts/12<br/>another species (GTEx) — scripts/13, 14"]
  B --> L
  G --> L
  M --> I["Identifiability audit — scripts/16<br/>batch nesting · estimable tissue pairs · QC-only baseline<br/>batch measured on bridging reference pools"]
  M --> X["Exercise follow-up — scripts/07, 15, 05b<br/>within-tissue separability and its covariate check<br/>design dates, physiology, the fingerprint by training duration<br/>the panel genes' training response"]
  L --> E["Site export with provenance — scripts/30<br/>results/ or results_frozen/ → site/data/*.json + provenance.json"]
  I --> E
  X --> E
  E --> S["Static site — site/ · GitHub Pages"]
  E --> T["Tests and CI<br/>pytest · JS conformal port (2,700 cases) · anchors · snapshot checksums"]
```

Every data-dependent step (imputation, prefilter, scaling, selection, tuning) runs inside the training fold; the
calibration animals are disjoint from the fit and test animals; shift experiments and negative results are
first-class. The rules every phase follows are frozen in `docs/EVALUATION_RULES.md`.

## Quick start (no data needed)

```bash
git clone https://github.com/CoryGardner/motrpac_challenge && cd motrpac_challenge
conda env create -f environment.yml && conda activate tfp      # or: pip install -r requirements.txt
make test          # Python tests (the results-dependent ones read results_frozen/), the JS conformal test, snapshot checksums
make site          # export site/data from results_frozen/, check the sanity anchors, run the site tests
python -m http.server -d site 8000                              # open http://localhost:8000
make smoke         # synthetic data → phases 02, 04, 05, 06 in --quick mode → results_smoke/
```

Expected: pytest ends with every test passed (one comparison against a live `results/` run is skipped when none is present); the JS test prints
`ok: 2713 assertions, 2700 fixture cases (540 with infinite threshold)`; `make site` ends with `all anchors ok`;
`make smoke` ends with `wrote results_smoke/06_conformal/TRNSCRPT` and every banner says the data are synthetic.
A fresh export from the snapshot is byte-identical to the committed `site/data/` apart from its `_meta` block
(`python tools/compare_site_data.py site/data <fresh export dir>`); CI runs exactly this comparison on every push.

## Setup

| environment | file | contents |
|---|---|---|
| Python 3.12 | `environment.yml`, `requirements.txt` | numpy 2.5.3, pandas 3.0.6, scipy 1.18.1, scikit-learn 1.9.1, matplotlib 3.11.2, statsmodels 0.15.0, pyarrow 25.0.0, joblib, pytest 9.1.1 (the versions of the reference run; `pyproject.toml` keeps the lower bounds) |
| R 4.5.3 | `environment-r.yml`, `R/install_deps.R` | data.table 1.18.6.1, jsonlite 2.0.0, remotes, BiocManager (Bioconductor 3.22), SummarizedExperiment, ExperimentHub; `R/install_deps.R` installs `MotrpacRatTraining6moData` at commit `f831a4f` (DESCRIPTION version 2.0.0) and `bodymapRat` 1.28.0 |
| Node 20 | `tools/package.json` | playwright 1.47.2 (only for the render pass and `make figures`; needs a system Chrome) |
| browser | `site/` | Plotly.js 2.35.2 from jsDelivr with a vendored fallback (`site/vendor/`) |

Hardware: the full analysis run took about 45 minutes on a 20-core workstation (phase timings below); `make site`
from the snapshot takes well under a minute; the GTEx preparation reads two 1.6 GB files.

## Data: sources, access and provenance

Nothing under `data/` is in the repository. What the repository holds is derived: `results_frozen/` (the result tables
the site reads, with a sha256 manifest) and `site/data/` (aggregates, per-sample class probabilities and the log2 CPM
of about 100 genes).

| dataset | version | obtained | how | licence |
|---|---|---|---|---|
| MoTrPAC rat 6-month endurance training, portal release c1.0 (rn6) | R package `MotrpacRatTraining6moData` 2.0.0, GitHub commit `f831a4fe421ec11687640452484a8247137aa74a` (tag v2.1.0; the DESCRIPTION reads 2.0.0 at both v2.0.0 and v2.1.0) | 2026-09-17 | `make export` (`R/install_deps.R`, `R/export_motrpac.R` → `data/raw/`, 1.9 GB) | package MIT; consortium data-use terms |
| MoTrPAC portal subset for phase 16 (optional) | Data Hub release c1.0: the RNA-seq QC table and the 19 per-tissue RSEM count files (≈ 165 MB) | — | `MOTRPAC_PORTAL=<dir holding rat-training-06/> make identifiability`; `make portal-check` lists the files (`docs/DATA_GUIDE.md` §9) | consortium data-use terms |
| Rat BodyMap, GEO GSE53960 | Bioconductor `bodymapRat` 1.28.0 (ExperimentHub) | 2026-09-17 | `make bodymap` (`R/export_bodymap.R` → `data/external/bodymap_*.csv`) | CC BY 4.0 |
| GTEx v8 open access (release 2017-06-05, RNASeQCv1.1.9) | gene TPM, gene reads, sample attributes | 2026-09-17 (TPM, attributes), 2026-09-18 (reads) | download the three files below into `data/external/gtex/`, check the checksums, `make gtex` | GTEx data-use policy |
| Rat–human orthologs | `RAT_TO_HUMAN_GENE` as shipped in the MoTrPAC package (one-to-one pairs) | with the package | `data/raw/rat_to_human_gene.csv` | as the package |

GTEx files (sha256):

```
https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz     ec783825ebebb8ba8525b73cb7ff9577162bc0d2e8a7fdc408283b7639c306fb
https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz   148ab1c84609608a00772b0e1431a1f87d5bf40b63071997510181d0a684a110
https://storage.googleapis.com/adult-gtex/annotations/v8/metadata-files/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt  74f6ab4c34ed2648d708a0ae6e6dff324f6c86ea723ae7d1c37d76f5221148f0
```

Citations: MoTrPAC Study Group, *Nature* 629, 174–183 (2024); Yu et al., *Nat Commun* 5, 3230 (2014); GTEx
Consortium, *Science* 369, 1318–1330 (2020). The data acknowledgement sentence MoTrPAC asks for is in
`docs/DATA_GUIDE.md`.

## Inputs and outputs

| what | where | contract |
|---|---|---|
| phenotypes | `data/raw/pheno.csv` | one row per vial: `viallabel` (11 digits), `bid` = `viallabel[:5]` (one per animal), `pid` (the animal, the unit of every split), `sex`, `group` ∈ {control, 1w, 2w, 4w, 8w}, `tissue`; identifiers are strings |
| sample tables | `data/raw/{norm,counts}/<ASSAY>__<TISSUE>.csv` | features × samples: four leading columns `feature, feature_ID, tissue, assay`, then one column per `viallabel` |
| processing metadata | `data/raw/meta/<ASSAY>.csv` | one row per vial with plate, batch, flowcell, TMT plex and channel |
| external | `data/external/bodymap_{counts,meta}.csv`, `gtex_{tpm,cpm}_subset.csv`, `gtex_meta.csv` | BodyMap: `feature_ID` (ENSRNOG) × 316 samples, meta with organ, sex, stage_weeks, animal_id; GTEx: samples × genes (log2(x + 1)), meta with SAMPID, donor, SMTSD, rat_tissue |
| results | `results/<phase>/*.csv` (git-ignored) and `results_frozen/` (committed) | e.g. `05_panels/TRNSCRPT/panel_curve.csv`: k, fold, balanced_accuracy, n_train_animals, n_test_animals; `06_conformal/TRNSCRPT/coverage.csv`: calibration, conformal, alpha, coverage, frac_empty, avg_set_size; `08_shift/TRNSCRPT/shift_table.csv`: split, arm, accuracy_all, coverage_target_seen, …; `31_site_regen/*/scores_*.csv`: per-sample class probabilities and calibration scores |
| site data | `site/data/*.json` | `headline.json` (tiles, ladder rows, extras, design constants), `exercise.json`, `provenance.json` (one entry per number: file, row selector, column, aggregation, tolerance), `manifest.json` (the files read, with sha256; versions; data access dates), `samples_*.json`, `conformal_fixtures.json` |
| figures | `figures/summary_figure.png`, `figures/home.png` | `make figures` (the first needs no data; the second renders the site) |

The full layout and the identifier conventions are in `docs/DATA_GUIDE.md`.

## Methods and provenance

| step | where | what |
|---|---|---|
| splitting | `src/tfp/splits.py` (frozen) | `StratifiedGroupKFold` on `pid`; leave-one-sex-out; train on controls, test on trained; fit/calibration split of the training animals |
| models | `src/tfp/models.py` | `RoundRobinSelector` (the best remaining gene of each tissue in turn, by one-vs-rest effect size); L2 logistic regression, C tuned by an inner animal-grouped grid search; tuned baselines (nearest centroid, L1/L2 logistic regression, random forest) on the same folds |
| conformal sets | `src/tfp/conformal.py` | split conformal (LAC and APS): the ⌈(n+1)(1−α)⌉-th smallest calibration score, +∞ when that rank exceeds n; Mondrian and floored variants; Clopper–Pearson certificate with fixed-sequence testing; the same rule is ported to the browser in `site/assets/conformal.js` and tested against 2,700 Python-computed cases |
| transfer | `src/tfp/transfer.py` | z-scores within dataset; super-classes for muscle and brain; recalibration on 3 or 5 target individuals over 20 draws |
| identifiability | `src/tfp/batch.py`, `scripts/16_identifiability.py` | Cramér's V, estimable tissue pairs, a QC-only multinomial classifier, and the between-plate variance of the reference-standard RNA pools relative to the variance of the tissue means |
| exercise follow-up | `scripts/07_*`, `scripts/15_*`, `scripts/05_panel_training_response.py` | sedentary vs trained within tissue (13 arms, permutation null, covariate-only classifiers), the design dates and physiology, the fingerprint by training duration, the panel genes' training response |
| provenance | `scripts/30_export_site_data.py` | no number on the site is typed by hand; `provenance.json` records the file, row, column and aggregation of each; `tests/test_site_data.py` re-reads every entry from the results; sanity anchors stop the export when a headline number moves |

External code: scikit-learn, pandas, numpy, scipy, matplotlib (Python); `bodymapRat`, `MotrpacRatTraining6moData`
(R); Plotly.js 2.35.2 (MIT, `site/vendor/LICENSE.plotly.txt`). AI tools assisted with code and writing; all results
were verified by the team.

## Validation

| check | command | expected |
|---|---|---|
| unit, leakage and I/O tests; snapshot; site provenance | `make test` | every test passed; `verify-frozen` prints `ok` |
| JS conformal port against Python | `node tests/test_site_conformal.js` | `ok: 2713 assertions, 2700 fixture cases (540 with infinite threshold)` |
| the browser scoring tool against the pipeline's BodyMap scores, and its example recalibration | `node tests/test_score_tool.js` | `ok: 976 assertions over 316 BodyMap samples; max \|Δp\| = 2.59e-4` … `ok: 1075 assertions in all` |
| every site number against the result files | `PYTHONPATH=src pytest -q tests/test_site_data.py` | passed (runs from `results/` or `results_frozen/`) |
| regeneration reproduces the published tables cell for cell | `PYTHONPATH=src pytest -q tests/test_regen_scores.py` | passed |
| sanity anchors of the headline numbers | `make site-data` | `all anchors ok` |
| expected outputs | export from `results_frozen/` and `python tools/compare_site_data.py site/data <dir>` | `identical apart from _meta` |
| no-data smoke run | `make smoke` | four phases write `results_smoke/`; banners say the data are synthetic |
| render pass | `make screenshots` (Chrome + playwright) | `all page renders clean` in both themes at 1440, 1024 and 390 px |

Known failure modes and limits (details on the site's Limitations page):

- with fewer than 9 calibration scores in a class at α = 0.10, the Mondrian threshold is +∞ and the set holds every tissue;
- 45 % of the three-donor GTEx recalibration draws have no finite threshold, so the recalibrated GTEx coverage must be read with its set size;
- juvenile BodyMap testes and spleen are misnamed by single developmental markers (Pgk2, Hbq1b); the panel is an adult fingerprint;
- 7 of the 20 panel genes have no BodyMap organ and are untested within species; Artery–Aorta → vena cava is an imperfect human mapping;
- the transcript gene filter runs before the split (the one known deviation from the frozen rules);
- sex is confounded with arrival cohort, and only the 8-week animals share sacrifice dates with the controls;
- without the portal subset, phase 16 uses the package metadata for the QC baseline and reports the bridging measurement as pending;
- BodyMap animal identifiers are inferred from the replicate index (GEO carries none).

## Reuse

- **Licence:** MIT for code, site and derived tables (`LICENSE`); the data keep their own terms (table above).
- **Cite:** `CITATION.cff` (version 1.6.0, tag `hackathon-submission-v7`) and the three data papers.
- **Score your samples:** the Explorer's *Score your own samples* tool takes a CSV of log2 CPM for the 20 panel genes
  (a template is provided) and returns tissue calls and 90 % prediction sets in the browser, with optional recalibration
  on labelled samples; the panel card (`site/data/panel_card.csv`, `.json`) lists the genes with their mean expression per
  tissue, and `site/data/panel_model.json` holds the coefficients (`scripts/34_panel_model.py`, `make panel-model`).
- **Extend:** the library is `src/tfp/`; a new phase is a numbered script with a Makefile target; if the site shows its
  numbers, add provenance entries in `scripts/30_export_site_data.py` and re-freeze (`make freeze-results`).
- **Extras:** `extras/` holds analyses outside the submission path (discordance, the report builder, the GEO fallback,
  the time-course investigations, the frozen replication notebooks); not covered by the tests.
- **Roadmap:** rerun on portal release c2.0 (rn7); a cross-tissue proteomics fingerprint; MoTrPAC human tissues when
  released; exercise-response panels under the same rules; further independently processed rat cohorts.
- **Follow-up on the branch `multiomic-overnight` (not merged; a human decision):** the proteome and metabolome carry the
  fingerprint too, on the right scale. On the portal's reporter-ion intensities tissue explains R² 0.991 of PC1 against
  0.0009 on the distributed ratios (`results_multiomic/01_rii/variance_partition.csv`, `variance_partition_ratio.csv`); a
  20-protein panel selected there names the tissue of 0.455 of 44 human TMT samples from 13 GTEx donors (donor-bootstrap
  95 % CI 0.36–0.55; 0.734 of 94 samples when both sides are processed the same way; `results_multiomic/03_prot_transfer/
  accuracy_overall.csv`, `rawppm/accuracy_overall.csv`), with coverage 0.068 under MoTrPAC calibration and 0.917 after
  recalibrating on five donors at 2.95 of 7 classes per set (`recalibration.csv`); and in Jiang 2020's TMT design 424 of
  528 tissue pairs share a run (`results_multiomic/06_external_identifiability/design_comparison.csv`), the counterexample
  the audit lacked. Page: `site/multiomic.html`; report: `docs/MULTIOMIC_REPORT.md`; pre-registration:
  `docs/PREREGISTRATION_MULTIOMIC.md`. Within-study numbers there are context (plex is nested in tissue on that scale too).
- **History:** `CHANGELOG.md`.

## Example runs

No data: the quick start above, plus `make figures` (Chrome + playwright) and
`python scripts/30_export_site_data.py --readme-table` for the key-results table.

Full run:

```bash
conda env create -f environment-r.yml && conda activate tfp-r && Rscript R/install_deps.R
make export                                     # the MoTrPAC package → data/raw (1.9 GB, ~15 min)
Rscript R/export_bodymap.R data/external 100    # rat BodyMap counts
# download the three GTEx files into data/external/gtex/ and check the sha256 (Data section)
conda activate tfp
make check all external                         # phases 02–08, 11–14: about 45 min on 20 cores
MOTRPAC_PORTAL=/path/to/portal make identifiability   # phase 16 (bridge needs the portal subset)
make time-course panel-training                 # phase 15 and the panel genes' training response
make regen-scores                               # per-sample score reruns of phases 06, 08, 12, 13
make site freeze-results figures                # export, tests, refresh results_frozen/, figures
```

Phase timings of the reference run (2026-09-18, 20 cores): baselines 196 s, panels 389 s, conformal 525 s, fusion
1,559 s, shift 89 s, bodymap 155 s, gtex 541 s.

## Repository map

```
site/                 the static site: index, exercise, explore, transfer, fingerprint, identifiability, methods, limitations, about
  data/               JSON exported from the results, with provenance.json and manifest.json
  assets/             theme.css, charts.js, site.js, conformal.js (the port), ladder.js, overview.js, pages/, brand/
src/tfp/              the library: config, io, splits (frozen), models, conformal, transfer, batch, plots, report, cli
scripts/              numbered phases 02–16 (07, 15, 05b feed the Exercise page), 30 site export, 32 summary figure, 33 freeze
results_frozen/       the result files the site reads, with MANIFEST.json (sha256)
tests/                leakage, I/O, conformal, batch, config, regeneration, snapshot, site provenance; the JS conformal test
docs/                 EVALUATION_RULES (frozen), DATA_GUIDE, EXTERNAL_VALIDATION, GTEX_TRANSFER, ABSTRACT_SUBMISSION,
                      COMPETITION_COMPLIANCE, NUMBERS_RECONCILIATION + reconciliation/ (what moved when the conformal rule was corrected before release)
tools/                linkcheck.py, screenshot.js, compare_site_data.py, make_logo_assets.py
R/                    install_deps.R, export_motrpac.R, export_bodymap.R
figures/              summary_figure.png, home.png
branding/             the team logo source
extras/               analyses outside the submission path (see extras/README.md)
results/, data/       git-ignored: pipeline outputs, the data export and the external downloads
```
