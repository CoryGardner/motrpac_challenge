# Changelog

All notable changes to this repository. The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/);
dates are UTC. Each section is a `hackathon-submission-v<n>` git tag; the version numbers are those of `CITATION.cff`
and `pyproject.toml`.

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

[1.3.0]: https://github.com/CoryGardner/motrpac_challenge/compare/hackathon-submission-v3...hackathon-submission-v4
[1.2.0]: https://github.com/CoryGardner/motrpac_challenge/compare/hackathon-submission-v2...hackathon-submission-v3
[1.1.0]: https://github.com/CoryGardner/motrpac_challenge/compare/hackathon-submission-v1...hackathon-submission-v2
[1.0.0]: https://github.com/CoryGardner/motrpac_challenge/releases/tag/hackathon-submission-v1
