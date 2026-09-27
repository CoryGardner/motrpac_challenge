# Changelog

All notable changes to this repository. The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/);
dates are UTC. Each section is a `hackathon-submission-v<n>` git tag; the version numbers are those of `CITATION.cff`
and `pyproject.toml`.

## [2.0.4] — 2026-09-27 (tag `hackathon-submission-v10.2.2`)

### Fixed
- Check samples: the explanation sentence names a gene only if its contribution is at least max(0.1, 0.1 × the largest),
  written as the chart writes it ("Umod (+2.9 sd)"), or says that no single gene dominates (the swapped kidney no longer
  cites Trim29 at z = +0.02).
- The example no longer recalibrates on its two deliberately swapped labels (at α = 0.05 that had made every set
  ambiguous and passed the swapped kidney as Consistent); the "Start here" note and "Calibrate to my lab" say so. Tests:
  both swaps stay Mismatch after recalibration at α = 0.05, 0.10 and 0.20.
- Set labels follow α (table header, summary tile, status key, drawer, charts, report); the CSV column is `set`.
- A forced scaling mode is described as the mode used, with the automatic rule's choice beside it; forcing within-set
  scaling on fewer than 8 samples or 3 tissues is a caveat.

### Changed
- Results table: sample IDs on one line, set and runner-up cells wrap, the "missing" column only when needed; the
  accepted-names table wraps its last column. The header is not sticky on phones (below 700 px).
- The science overview's recalibration tile reads "observed coverage in another lab after recalibrating on three of its
  animals".
- The vena cava case in both directions (`results_product/40_product/venacv_all.csv`, `venacv_summary.csv`, `pv_venacv_*`):
  of the held-out vena cava vials, the consortium flagged 10; the model calls 7 of them brown fat and none of the 40
  unflagged ones; one sentence on Check samples.
- README: the team repository beside the public mirror (top block, quick start), the site's rounding in the flag table;
  `CITATION.cff` repository-code is the team repository. `LICENSE` lists the multiomic sources and the BodyMap test
  fixture's CC BY 4.0 terms. `src/tfp/rii.py` defaults to `data/quant-id/rat-training-06` (the environment variable
  still overrides). The Pages workflow publishes only from the personal repository (Pages is not enabled on the team
  repository yet). `figures/home.png` re-rendered. The BodyMap fixture archive is written byte-reproducibly.

## [2.0.3] — 2026-09-27 (tag `hackathon-submission-v10.2.1`)

### Changed (audit fixes)
- Check samples: the results table puts the most actionable rows first (Mismatch, Can't confirm, not in reference,
  Unknown), a "Start here" note names the swapped pair of the example and, when most sets are empty, explains why and
  offers a one-click recalibration; a "missing" column marks samples with blank panel-gene cells.
- Input checks: a table with none of the panel genes is refused with a message; fewer than half the genes, and blank or
  non-numeric panel-gene cells (naming the samples), warn; "1 sample" grammar. A template of the 20 genes can be downloaded.
- The flag-quality table wraps its headers instead of clipping at 1440 px.
- "90 % set" is introduced as a conformal prediction set on Check samples, with a link to a new Glossary (The science
  page, also in the science sub-navigation).
- Naming: the Explorer is the Reference atlas everywhere; its leftover scoring stub is removed (the old
  `explore.html#score-tool` link forwards to Check samples); stale "Explorer" and "Home page" references fixed in the
  README, the compliance doc, Methods and the Makefile; the compliance doc's snapshot count is 179.
- Meaning: the browser model, the panel card and the per-fold panels are named as three fits of one method (14 of the
  browser model's 20 genes are on the panel card) on Check samples and the Panel page; the trust strip's "fails safe" is
  qualified (mixed uploads, not reference scaling); the science heading, the ladder subtitle, the Transfer target zone and
  the abstract no longer say the guarantee is restored by recalibration (the abstract now quotes the 20-gene panel's
  one-vial coverage and the five-donor GTEx result); the batch tile names its pool and gives the range and the panel-gene
  value; the Identifiability QC chart plots balanced accuracy, as its text and label say (it plotted accuracy).
- README: Validation rows for the Check samples tests, quick-start output, repository map, full-run recipe
  (`make panel-model product-validation`), the GTEx three-donor row label, and a note that `results/` paths are
  committed under `results_frozen/`.

## [2.0.2] — 2026-09-27 (tag `hackathon-submission-v10.2`)

### Fixed
- Check samples read an empty cell, a whitespace-only cell and every cell past the end of a short row as 0 ("not
  expressed") instead of missing (`Number("")` is 0), and did not list the gene as missing. One cell parser
  (`parseCell` in `site/assets/check-core.js`) now reads "", whitespace, NA, NaN, null, "-" and anything non-numeric as
  missing in both input formats; short rows are named in a warning; empty cells of a count matrix are left out of the
  library size and counted; a 20-gene table with values impossible for log2(CPM + 1) (< −0.5 or > 20) warns that it
  looks linear or TPM. Tests: blank and NA give identical probabilities, calls and sets in all 1,600 single-cell cases of
  the example, and the blanked gene is named each time; truncated-row, linear-scale and empty-count-cell cases.

### Added
- Composition sensitivity of within-set scaling: 200 random subsets (fixed seed; 8–80 samples, ≥ 3 organs) of the
  BodyMap 21-week adults, scored within-set and with reference scaling (`results_product/40_product/composition*.csv`);
  the 5th–95th percentile ranges on Check samples ("What to expect") and the Limitations page. Reference scaling is
  asserted per-sample (identical probabilities whatever else is uploaded). The automatic rule is unchanged.
- Per-draw recalibration coverage replayed from the frozen thresholds (`recal_draws.csv`, `recal_draws_summary.csv`;
  no rerun; the means reproduce `recalibration.csv`), exported as `pv_recal_*`.

### Changed
- "Restores the guarantee" after recalibration now reads "restores observed coverage", with the per-draw range, on The
  science, Transfer and Identifiability pages and in the README; for GTEx the text leads with five donors (every draw
  finite) and gives the three-donor split. Methods and Limitations explain why (exchangeability, the ≥ 9-animal
  arithmetic, inferred BodyMap animal IDs); "Calibrate to my lab" says the recalibrated coverage is observed, not
  guaranteed, and pooled within animals.
- LICENSE: the copyright holder is the Stanford Bioinformatics Center. The organizers' repository history is merged.

## [2.0.1] — 2026-09-27 (tag `hackathon-submission-v10.1`)

### Added
- Multiomic page, section 1: "The proteomics problem in one picture", two per-vial PCA panels (as distributed vs rebuilt
  from the reporter-ion intensities) with a colour-by tissue / sex / plex toggle and a shared legend; the R² bars move into
  "Show R² by component". `scripts/multiomic/01b_pca_scores.py` rebuilds the two matrices of
  `variance_partition.csv` and `variance_partition_ratio.csv` with the phase-1 code and seed and writes
  `results_multiomic/01_rii/pca_scores_{rii,ratio}.csv`, gated on the published tissue R² (4 decimals; it matched, as
  did the sex R² and explained variance). The export re-checks the gate; `mo_pca_*` provenance entries; a new test.
- Slide `results_multiomic/figures/fig0_pca_two_scales.png` (`figures_png.py --only fig0`).

## [2.0.0] — 2026-09-27 (tag `hackathon-submission-v10`)

### Added
- **Check samples** (`site/index.html`, `site/assets/pages/check.js`), the product: "Is this sample the tissue you think
  it is?" Input as a 20-gene log2 CPM table or a full raw-count matrix (either orientation, Ensembl IDs or symbols,
  parsed in a Web Worker with progress, converted with the pipeline's `io.log_cpm`), with an optional claimed-tissue
  column mapped through synonyms to the 19 classes. Per sample: call, calibrated 90 % set, status (Confident / Ambiguous
  / Unknown), claim status (Consistent / Mismatch / Can't confirm / Not in reference); an α slider; a sortable,
  keyboard-operable table with a flagged-only filter; a reference map (PCA of the 20-gene z-space fitted on the MoTrPAC
  vials); a sample drawer with probabilities against the threshold, "why X, not Y" per-gene contributions and gene values
  against the reference tissues' median and IQR, and a generated sentence; "Calibrate to my lab"; CSV (spreadsheet
  formula injection neutralised) and printable HTML report downloads; species and scaling guards with banners.
- `site/assets/check-core.js` (pure functions) and `site/assets/check-worker.js`; `tests/test_check_core.js` (status
  rules, count → log2 CPM parity with Python on real BodyMap counts in both orientations, projection parity, contributions
  + intercept = logit, synonyms, species and scaling guards, CSV safety, the example's swap detection).
- `scripts/40_product_validation.py` → `results_product/40_product/` (scaling on BodyMap in three modes, flag rates at
  α ∈ {0.05, 0.10, 0.20} in-study and on BodyMap with 1,000 simulated swaps, the vena cava → brown fat cases, the reference
  PCA, parity fixtures in `tests/fixtures/`); `scripts/41_export_product_data.py` → `site/data/product.json` with `pv_*`
  provenance entries; `tests/test_product_data.py`. `make product-validation`, `make product-data`.

### Changed
- Information architecture: `index.html` is the product; the former Home is `site/science.html` ("The science"). Primary
  nav: Check samples · Reference atlas · The science · About; the science pages carry a secondary row (Overview, Panel,
  Transfer, Identifiability, Exercise, Multiomic, Methods, Limitations). Every existing URL still works; the Explorer's
  scoring tool moved to Check samples (the `explore.html#score-tool` anchor links there).
- README opens with "Use it"; the track outputs map to the product.
- CI: the site-data comparison now re-runs the multiomic and product exports on the fresh copy (it had failed since
  v8 because `multiomic.json` was only in the committed data); the JS job runs `test_check_core.js`. `make test` no longer
  reports a failing JS test as "skipped". `make site-data` also re-runs the product export.

## [1.8.1] — 2026-09-27 (tag `hackathon-submission-v9.1`)

### Added
- Home: "What this submission delivers", four cards under the hero tiles mapping the track's four outputs (classifier,
  minimal panel, feature-selection workflow, interactive explanation tool) to the pages, data files and code.
- The same mapping as "Track outputs → where they are" in the README (after "The answer") and in
  `docs/COMPETITION_COMPLIANCE.md`.

### Changed
- Home word budget in `tools/screenshot.js` raised from 1200 to 1300 for the new section.

## [1.8.0] — 2026-09-27 (tag `hackathon-submission-v9`)

### Added
- Multiomic page: an "In plain words" box under the banner; a "How it was done" section before section 1 (why 7 of 19
  tissues, how the protein matrix is built, how the round-robin selector picks proteins high or low in a tissue, the
  20-protein transfer panel with its direction and rat and human effects); the per-tissue Jiang table taken out of its
  collapsed block and shown on both scales with n samples and n donors; in the metabolite section the within-MoTrPAC
  tissue R² of the first metabolite PC, the missing batch variable, a per-organ table for the mouse aging atlas at k20,
  and why metabolites are not the core fingerprint; a "Data sources" block with verified DOIs and accessions.
- `scripts/multiomic/export_site_data.py` exports the new values (`how`, `metab_within`, `deep_mw_by_tissue_k20`,
  `metabolites_stopped.source_metabolites`) with provenance entries `mo_how_*`, `mo_metab_*`, `mo_deep_mw_by_tissue_k20`,
  `mo_hilic_mw_src_metabolites`; two new checks in `tests/test_site_data.py`.
- Citations of Jiang 2020, Wang 2019, Geiger 2013, Sato 2022, Metabolomics Workbench ST003188 and the MoTrPAC
  reporter-ion files in the README data table (with download URLs and sha256), `site/about.html`,
  `docs/COMPETITION_COMPLIANCE.md` (CFDE: Metabolomics Workbench and the GTEx proteome) and `CITATION.cff`.
- Home: a seventh five-minute-tour item for the Multiomic page. Methods: a "Multiomic follow-up" subsection.

### Changed
- Multiomic source note: the work is merged into main; points to the report, pre-registration and log. The two
  hand-typed metabolite numbers on the page are now read from `multiomic.json`.
- `make site-data` re-runs the multiomic export after the main export, which rewrites `provenance.json`.
- Home tiles: the interval line may wrap, which removes the 15 px horizontal overflow at 768 px.
- Methods word budget in `tools/screenshot.js` raised from 1000 to 1200 for the new subsection.

## [1.7.0] — 2026-09-27 (tag `hackathon-submission-v8`)

### Added
- The multiomic follow-up (merged from branch `multiomic-overnight`; pre-registration `docs/PREREGISTRATION_MULTIOMIC.md`,
  report `docs/MULTIOMIC_REPORT.md`, results under `results_multiomic/`): MoTrPAC proteomics on the portal's reporter-ion
  scale, protein transfer to the Jiang 2020 GTEx proteome, metabolite transfer to two mouse atlases, fusion judged by
  transfer, and the identifiability of the external designs. New page `site/multiomic.html` with `site/data/multiomic.json`,
  provenance entries `mo_*` and checks in `tests/test_site_data.py`; slide figures in `results_multiomic/figures/`.

### Changed
- Navigation: "Multiomic" after Exercise; tighter nav spacing so the header stays on one row from 1024 to 1440 px.
- README scope, roadmap and follow-up paragraph, and the Limitations page: proteomics is within-tissue only as distributed
  (ratios), not on the reporter-ion scale.

## [1.6.0] — 2026-09-27 (tag `hackathon-submission-v7`)

### Changed
- Home tiles: 0.976 accuracy, 1.000 mapped organs, 0.943 coverage after recalibrating on three BodyMap animals
  (new `tile_bodymap_recal_k20`, from `results/12_bodymap/recalibration.csv`), 1.6 % batch on the bridging pools
  (relabelled, with the 1.6–5.3 % range). `headline.json` carries `home_tiles`; the 0.618 coverage tile stays in
  `tiles` for the Identifiability page.
- README key results: the in-distribution coverage row is labelled all genes and the 20-gene row (0.900 / 0.924) is
  added, from the generator.
- Ladder value labels sit above the whiskers; the abstain chart's labels on hatched segments have a background.
- Scoring tool counter reads "k labelled (minimum N)". `CITATION.cff` release description fixed.

## [1.5.0] — 2026-09-27 (tag `hackathon-submission-v6`)

### Changed
- Home §3 and panel c of `figures/summary_figure.png` show batch per bridging reference pool (six pools: gastrocnemius
  99 and 88 on six plates at both sites, liver and hippocampus 99 and 88 at one site) as Σ V_batch / Σ V_tissue over
  all genes, 1.6–5.3 %; the per-gene chart stays on the Identifiability page. `nesting.json` carries `bridge.pools`
  with a provenance entry per pool; the summary figure's panel titles are one row.
- The scoring tool's recalibration is demonstrable: "Load the example" fills the box with every 21-week rat BodyMap
  sample (80 rows) and a `true_tissue` column pre-filled for 12 of them; the recalibrate button names the number of
  labelled samples a finite threshold needs at the chosen α (⌈1/α⌉ − 1) and stays disabled until then, with a live
  counter; the banner explains abstentions with the MoTrPAC calibration (0.618 coverage on the BodyMap) and what
  "three animals" meant in the study (about 25 labelled samples). `panel_model.json` carries these numbers with
  provenance. `tests/test_score_tool.js` checks the example: 12 labels give a finite threshold, coverage 0.911 on the
  68 unlabelled mapped-organ samples with singleton sets wherever a single-tissue organ is covered, and thymus and
  uterus abstain on most samples.
- Cosmetic: the header badge is 72 px (56 px on phones); the theme toggle stays on the nav row at 1024 px; the "1 of 171" tile on the Identifiability page is a
  normal-width tile; the Home ladder title is one line with the rest in the subtitle; `make screenshots` renders
  1440, 1024 and 390 px.

## [1.4.0] — 2026-09-27 (tag `hackathon-submission-v5`)

### Added
- "Score your own samples" on the Explorer: paste or upload log2 CPM for the 20 panel genes (template provided) and
  read each sample's tissue call and 90 % conformal set, download the results, and recalibrate on labelled samples.
  Everything runs in the browser (`site/assets/score.js`); `scripts/34_panel_model.py` exports the transfer model
  (`site/data/panel_model.json`) and validates that the 20-gene arithmetic reproduces the pipeline; the panel card
  (`site/data/panel_card.csv`, `.json`) lists the genes with their mean expression per tissue.
  `tests/test_score_tool.js` checks that the 316 rat BodyMap samples come out with the pipeline's calls and sets.
- Home page rebuilt around the answer: the lede, four tiles (including the 1.6 % bridge measurement), the one-picture
  diagram, the ladder, three points, and three sections (the panel curve, the empty-set stack, the per-gene bridge
  chart); tour stop 6 for the scoring tool.
- `scripts/15_time_course_design.py` writes the consortium's flagged vials; the Panel page reports how many of the
  vena cava calls that read as brown adipose were flagged for brown-fat contamination.

### Changed
- Exercise page reordered: the training response first (VO2max, body fat), training visible in every omic layer, the
  fingerprint invariant to it, training a minor axis within a tissue, the design note, then the covariate check.
- Trained-animal results harmonised: the headline is the panel fit on sedentary controls scoring all 40 trained
  animals (0.961, coverage 0.903); the cohort-matched 8-week group (0.972 / 0.917) is labelled as such.
- Summary figure: one-line titles with subtitles; panel c is the per-gene bridge chart.
- Stat tiles show the interval on its own line; the overview boxes are sized to their content; the About versions
  and the identifiability nesting heatmap sit behind toggles; the Methods rules section and the compliance section
  are removed; the footer is one line.
- README: the answer in three paragraphs, the home figure, the 8-week-only row, the scoring tool.
- `Makefile`: `panel-model`; `make test` runs `tests/test_score_tool.js`; CI likewise.

## [1.3.0] — 2026-09-27 (tag `hackathon-submission-v4`)

### Added
- `results_frozen/`: a committed snapshot of the 171 derived result files the site reads, with a SHA-256
  `MANIFEST.json`; `scripts/33_freeze_results.py` writes it and `--verify` checks it (`make freeze-results`,
  `make verify-frozen`). `tfp.config.results_root()` resolves `results/` when a complete run is present, else the
  snapshot, for the export and the tests.
- Pinned environments: `requirements.txt` (the exact versions of the reference run), `environment.yml` (conda,
  Python 3.12.14, Node 20.20.2), `environment-r.yml` (R 4.5.3, Bioconductor 3.22). `R/install_deps.R` pins
  `MotrpacRatTraining6moData` to GitHub commit `f831a4f` (tag v2.1.0; the DESCRIPTION reads 2.0.0 at both v2.0.0
  and v2.1.0) with a tarball fallback, and installs `bodymapRat`; `R/export_bodymap.R` writes
  `data/external/bodymap_provenance.json`.
- The Exercise page (`site/exercise.html`): the training response of the panel genes in their marker tissues
  (`scripts/05_panel_training_response.py`, `make panel-training`).
- Home-page overview; the team logo and the derived badge, favicon and social-preview images (`make brand`).
- `extras/`: analyses outside the submission (discordance, the report notebooks, the GEO fallback for BodyMap, the
  quantile-fix diff, time-course investigations); not covered by the tests.
- `Makefile` targets: `help`, `smoke` (phases 02, 04, 05, 06 on synthetic data, no real data needed), `test`,
  `identifiability` with `--bridge` and `MOTRPAC_PORTAL`, `portal-check`, `regen-scores`, `freeze-results`,
  `verify-frozen`, `figures`, `brand`, `time-course`, `panel-training`.
- `scripts/16_identifiability.py --portal` / `--check-portal`: the portal inputs are read from `MOTRPAC_PORTAL`
  (`docs/DATA_GUIDE.md` §9 lists the 20 files).
- `tools/compare_site_data.py`: compares two `site/data` exports.
- `CHANGELOG.md`; `site/vendor/LICENSE.plotly.txt` (the MIT licence of the vendored Plotly.js 2.35.2).
- Tests: results-root resolver, snapshot integrity, manifest contents, integer counts, orphan files.

### Changed
- README rewritten to the submission scaffold; a tone and concision pass over the documentation, the site and the
  scripts' docstrings; `docs/GTEX_TRANSFER.md` carries the GTEx v8 download URLs, sizes and SHA-256 checksums.
- `docs/COMPETITION_COMPLIANCE.md` names the licences (rat BodyMap: CC BY 4.0), the pinned environments and the
  snapshot and CI as evidence of integrity, and carries the AI-use statement; `LICENSE` names the team and every
  data source's terms; `CITATION.cff` 1.3.0 with software references for the two data packages.
- The site export (`scripts/30_export_site_data.py`) takes `--results` and `--out`, writes canonical provenance
  paths and a manifest of the files it actually read, and builds the numbers reconciliation from the committed
  `docs/reconciliation/pre_quantile_fix_values.csv`.
- CI (`.github/workflows/tests.yml`) installs the pinned requirements and, after the tests, verifies the snapshot,
  re-exports the site data from it, compares the export with the committed `site/data/` and checks the site's links
  — all without the raw data.
- Team name spelled "The Rat PAC" everywhere.

### Removed
- From version control: notebook outputs and runs (52 MB), duplicate expected tables, the build's planning
  documents and log, orphan site data.

## [1.2.0] — 2026-09-27 02:11 UTC (tag `hackathon-submission-v3`, commit `6adc938`)

### Added
- Transfer ladder: a training-state rung (fit on controls, tested on trained animals) with animal-bootstrap
  intervals and a data-driven title.
- BodyMap age curve: n and 95 % intervals per point.
- Three new tests (53 pass).

### Changed
- About page and `CITATION.cff`: the team block; version 1.2.0.
- QC-only baseline reported as balanced accuracy, like the gene model; anchors, README and abstract regenerated.
- Dark mode: the sequential colour scale is monotonic from the surface; heatmap label contrast is automatic.
- Explorer: α and variant controls disabled in recalibrated mode; the stable-core table fits; annotation moved;
  limitations wording.
- The one estimable RNA-seq tissue pair (ovary and testes) is named wherever "1 of 171" appears.
- Site footer and README point to github.com/CoryGardner/motrpac_challenge and the GitHub Pages URL; the deployment
  is recorded in the docs.

## [1.1.0] — 2026-09-27 01:14 UTC (tag `hackathon-submission-v2`, commit `d960d96`)

`CITATION.cff` still read 1.0.0 at this tag; 1.1.0 is assigned here retrospectively.

### Added
- `scripts/08_shift_tests.py --save-scores` (per-vial sets, probabilities, calibration scores); reruns at k = 20 and
  k = 50 under `results/31_site_regen/08_shift_k{20,50}` reproduce `results/08` cell for cell
  (`tests/test_regen_scores.py`).
- Cluster-bootstrap (or exact animal-level) intervals on the BodyMap and GTEx rungs and tiles; k = 50 held-out-sex
  rows.

### Changed
- Exporter: the in-distribution rung comes from the phase-06 models' own accuracy; held-out-sex coverage, empty and
  wrong rates are over seen-class vials.
- Wording audit: accuracy "degrades gracefully"; "every mapped adult organ (9 of 11, super-class scoring)"; each
  tissue sits inside one plate, library batch and flowcell; the panel is re-selected inside each fold; abstract
  ≤ 200 words.

### Removed
- `site/beyond.html`, the training-vs-batch verdict section, unused JSON exports and references to parallel work.

## [1.0.0] — 2026-09-27 00:43 UTC (tag `hackathon-submission-v1`, commit `33b3733`)

The submission build on top of the pipeline as received (below), 2026-09-26 23:15 to 2026-09-27 00:43 UTC.

### Added
- Phase 16 (`scripts/16_identifiability.py`): the identifiability recompute — batch nesting, estimable tissue
  pairs, the QC-only baseline — and batch measured directly on the bridging reference pools (`--bridge`).
- `--save-scores` in phases 06, 12 and 13: per-sample probabilities, calibration scores and per-draw thresholds.
- The site data export with provenance (`scripts/30_export_site_data.py`: 235 provenance entries, fixtures with
  infinite thresholds, the numbers-reconciliation document) and the browser's conformal port with its tests
  (`site/assets/conformal.js`, `tests/test_site_conformal.js`).
- The site: design system, chart template and shared components; home (question, answer, four tiles, the transfer
  ladder, the five-minute tour); explorer (tissue card with client-side conformal sets, gene explorer, panel
  builder, animals-needed calculator); transfer (ladder, empty-set stack, recalibration cost, age and developmental
  markers, GTEx per tissue, representations, panel survival); fingerprint, identifiability, methods, limitations
  and about pages; the full-width pipeline diagram.
- Submission documents: README, abstract, competition compliance, `CITATION.cff`, the licence note, the GitHub
  workflows (tests without raw data, Pages deployment), `figures/summary_figure.png`.

### Changed
- The library renamed `src/motrpac` → `src/tfp` (the old name collided with the MoTrPAC R packages).
- Review pass: the sourced F-test number on Home, both model calls on the Explorer card, the ladder caption,
  hand-typed numbers wired to JSON, theme-toggle state, held-out-sex denominators, hide-the-answer, the
  regeneration test pins every table.

## Before the tags

Commit `5236735` (2026-09-26 23:14 UTC) is the pipeline as it stood at the start of the submission build: the
library (then `src/motrpac`), the numbered phase scripts 00–14, tests, notebooks and docs. Its inputs and
milestones, as recorded in the documents: the export of `MotrpacRatTraining6moData` 2.0.0 and of the rat BodyMap
through `bodymapRat` (2026-09-17); the GTEx v8 downloads (TPM and sample attributes 2026-09-17, read counts
2026-09-18); the full run of 2026-09-18; the conformal-quantile fix of 2026-09-25 (the textbook rule, the
⌈(n + 1)(1 − α)⌉-th smallest score, and a full set when that rank exceeds n), after which phases 06, 08 and 12–14
were rerun — its effect on every headline number is in `docs/NUMBERS_RECONCILIATION.md`. `results/` and `data/`
are not versioned; the results are reproducible with `make all`, `make external` and `make identifiability`.

[2.0.4]: https://github.com/CoryGardner/motrpac_challenge/compare/hackathon-submission-v10.2.1...hackathon-submission-v10.2.2
[2.0.3]: https://github.com/CoryGardner/motrpac_challenge/compare/hackathon-submission-v10.2...hackathon-submission-v10.2.1
[2.0.2]: https://github.com/CoryGardner/motrpac_challenge/compare/hackathon-submission-v10.1...hackathon-submission-v10.2
[2.0.1]: https://github.com/CoryGardner/motrpac_challenge/compare/hackathon-submission-v10...hackathon-submission-v10.1
[2.0.0]: https://github.com/CoryGardner/motrpac_challenge/compare/hackathon-submission-v9.1...hackathon-submission-v10
[1.8.1]: https://github.com/CoryGardner/motrpac_challenge/compare/hackathon-submission-v9...hackathon-submission-v9.1
[1.8.0]: https://github.com/CoryGardner/motrpac_challenge/compare/hackathon-submission-v8...hackathon-submission-v9
[1.7.0]: https://github.com/CoryGardner/motrpac_challenge/compare/hackathon-submission-v7...hackathon-submission-v8
[1.6.0]: https://github.com/CoryGardner/motrpac_challenge/compare/hackathon-submission-v6...hackathon-submission-v7
[1.5.0]: https://github.com/CoryGardner/motrpac_challenge/compare/hackathon-submission-v5...hackathon-submission-v6
[1.4.0]: https://github.com/CoryGardner/motrpac_challenge/compare/hackathon-submission-v4...hackathon-submission-v5
[1.3.0]: https://github.com/CoryGardner/motrpac_challenge/compare/hackathon-submission-v3...hackathon-submission-v4
[1.2.0]: https://github.com/CoryGardner/motrpac_challenge/compare/hackathon-submission-v2...hackathon-submission-v3
[1.1.0]: https://github.com/CoryGardner/motrpac_challenge/compare/hackathon-submission-v1...hackathon-submission-v2
[1.0.0]: https://github.com/CoryGardner/motrpac_challenge/releases/tag/hackathon-submission-v1
