# Submission site, explorer and repo — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan
> task-by-task (the user supplied the spec and said there are no check-ins; execution is native, in
> this session, with subagents only for independent verification). Steps use checkbox (`- [ ]`) syntax.

**Goal:** Turn the verified post-fix `results/` of this pipeline into a static, judge-facing site
(`site/`), an interactive explorer with client-side conformal sets, and a deploy-ready open-source repo.

**Architecture:** One export script (`scripts/30_export_site_data.py`) is the only bridge between
`results/` (plus two small recomputations: phase 16 identifiability and the `--save-scores`
regeneration of phases 06/12/13) and `site/data/*.json`; every JSON value carries a provenance entry
naming the results file, row and column it came from. The site is plain HTML + CSS + ES-module JS
with Plotly (CDN + vendored fallback); `site/assets/conformal.js` is a faithful port of
`conformal.py` and is tested against Python-generated fixtures under Node.

**Tech Stack:** Python 3.12 (`motrpac-py` conda env: pandas 3.0, numpy 2.5, sklearn 1.9), Plotly.js
2.35.2 (cartesian bundle vendored), Node 20 + Playwright (system Chrome, `channel: "chrome"`) for
screenshots and the JS test, pytest, GitHub Actions.

**Spec:** `docs/superpowers/specs/2026-09-26-submission-site-spec.md` (the user's brief, verbatim).

## Global Constraints

- Repo root is `code/pipeline/`; all paths below are relative to it. It becomes a git repo in Task 0.
- `results/` is truth and is never modified. Regenerated outputs go to NEW directories:
  `results/16_identifiability/` (new analysis) and `results/31_site_regen/<phase>/` (reruns with
  `--save-scores`), produced with `MOTRPAC_NO_REPORT=1` so `results/REPORT.md` is not appended.
- Stale documents (`docs/findings/FINDINGS_REPORT.md` at workspace root, `results/SUMMARY.md`,
  `results/ABSTRACT.md`) get a banner only. Exception granted by the spec ("add a banner at the top
  of each") — this is the one write into `results/`.
- `src/motrpac/splits.py` is FROZEN: add functions only. `docs/EVALUATION_RULES.md` is FROZEN.
- No number typed by hand into HTML. Numbers reach pages only through `site/data/*.json`. A number
  without a source renders as "pending" with a reason string from the JSON.
- Every accuracy shown carries per-fold spread (sd) and n (animals, samples). Animal = `pid`.
- Every JSON ≤ 3 MB; `site/` ≤ 15 MB; floats quantised (probabilities 6 dp, expression 2 dp,
  aggregates as stored). No gene × sample matrix beyond the ≤ 100 panel genes.
- Design tokens exactly as in the spec §2; palette encodes ORGAN SYSTEM only; tissue identity by
  label/tooltip. One y-axis per chart. Title = claim; caption has "Source:" and "What it does not
  show:". Data-table toggle on every chart.
- Site works from `python -m http.server -d site 8000` and unchanged on GitHub Pages: relative
  paths, `.nojekyll`, no build step, no server routing.
- Pinned outside resources: Plotly.js 2.35.2; GTEx v8 (2017-06-05 release); BodyMap GSE53960; MoTrPAC
  c1.0 = `MotrpacRatTraining6moData` 2.0.0; human package not used here.
- Commit in logical chunks; tag `hackathon-submission-v1` at the end. Never push.
- Log to `docs/BUILD_LOG.md` roughly every 30 minutes of work.

## Decisions fixed during orientation (do not re-derive)

- Phases present: 02, 03(+inner, metab_core9), 04, 05, 06 (two runs), 07, 08, 09, 12, 13, 14,
  14_cpm, 15_time_course. Absent: 16–26 except 15. No pre-registration file exists anywhere in the
  workspace (grep "prereg|pre-registration" finds only unrelated starter READMEs).
- Post-fix in-distribution conformal coverage (TRNSCRPT full model, LAC, α = 0.10, mean of 5 folds):
  pooled vials 0.9077 (set size 0.910); one vial per animal (22 cal animals) 0.9165 (set size 0.920);
  Mondrian pooled 0.9188; floored pooled 0.9700. The 0.962 in old docs is pre-fix one-vial-per-animal.
  Source: `results/06_conformal/TRNSCRPT/coverage_marginal_vs_mondrian.csv` and `coverage.csv`.
- `notebooks/expected_values.csv`: `published_value` is PRE-fix; `reference_value` is the post-fix
  value to reproduce. Verification of regenerated aggregates uses `reference_value` and `tolerance`.
- Estimable pair definition (RNA-seq, recomputed from `data/raw/meta/TRNSCRPT.csv`, 899 study vials):
  a pair of tissues is estimable when the two tissues share a level of EVERY processing variable
  (RNA extraction plate, library batch, flowcell) so that a within-batch contrast exists. Plate:
  2 pairs share (Adrenal–Lung, Ovary–Testes); library batch: 2 (Hippocampus–Lung, Ovary–Testes);
  flowcell: 42; all three: 1 = Ovary–Testes, which is female-only vs male-only (the sex contrast).
  17 plates, 17 library batches, 4 flowcells, 2 sites.
- The QC-only baseline of notebook section 14 reads the consortium QC table at
  `../../data/quant-id/rat-training-06/c1.0/transcriptomics/qa-qc/motrpac_pass1b-06_transcript-rna-seq_qa-qc-metrics.csv`
  (935 rows: 899 study + 36 reference-standard vials, `Sample_category` = study/ref). Phase 16
  reproduces it on the phase-04 folds; expected: technical acc 0.8732, composition 0.9488, all 0.9755.
- Reference-standard ("ref") vials exist only in that QC table and the portal count files, not in
  `data/raw/`. The bridge-sample variance number (~1.7 %) belongs to the parallel audit and is NOT
  in this copy → render "pending" unless Task 6b's optional recomputation is done and labelled.
- Explorer data for MoTrPAC vials come from the phase-06 design (fit 18 / calibrate 22 / test 10
  animals per fold): test-vial probabilities and calibration scores from the same model, for the
  full model and the panel family at k = 20 and k = 50 (extended in `--save-scores`). These are
  out-of-fold for every vial. Phase-05 (40-animal fits) probabilities are optional extras.
- Recalibration thresholds: `transfer.conformal_transfer` gains an optional `collect` list that
  records `(model, n_recal, draw, chosen_ids, q_t, n_cal_scores)`; the export takes draw 0 per
  (target, model, n). If that changes any aggregate, stop and report.
- Site page list and order: index, explore, transfer, fingerprint, identifiability, beyond, methods,
  limitations, about. Explorer α slider exactness comes from the JS port; Python fixtures are made
  from the SAME quantised JSON the site ships, so equality is exact.
- Library rename: `src/motrpac/` → `src/tfp/` ("tissue fingerprint pipeline"). Imports, tests,
  scripts, investigations, Makefile, notebook builder and the two `.ipynb` (string patch of
  `src/motrpac` → `src/tfp` in cell sources only; outputs untouched) are updated; notebooks are NOT
  re-executed. No compatibility shim.

## Review Focus

1. α values on the slider that make `(n+1)(1−α)` land within 1e-9 of an integer (e.g. n = 22,
   α = 0.10 → 20.7; n = 19, α = 0.05 → 19.0): the JS port must apply the same `−1e-9` tolerance and
   return the 19th score, not +∞. Pinned by Task 7's fixture at α ∈ {0.05, 0.10, 0.20} on n = 22 and on
   the one-vial-per-animal set (n = 22) and by a unit case with n = 19, α = 0.05.
2. A class with zero calibration scores (OVARY/TESTES absent in a one-per-animal calibration of the
   opposite sex; BodyMap calibration missing a class): Mondrian must fall back to the marginal
   threshold, floored must be max(fallback, marginal). Pinned by Task 7's unit case.
3. A sample whose organ has no mapped class (BodyMap thymus/uterus, GTEx none): the card must show
   the set, "no mapped tissue" instead of correct/incorrect, and never crash. Pinned by Task 9b's
   explorer test list (manual render + a JS unit check in Task 7 on an `expected: null` fixture row).
4. Absent optional data (no phase 17/22–26, no bridge measurement, no k50 held-out-sex row): every
   page must render a labelled "pending" block rather than an empty chart or a JS error. Pinned by
   the provenance test (Task 8) asserting every `pending` entry has a `reason`, and by the render
   pass console check (Task 13).
5. Offline demo: CDN Plotly unreachable → the vendored bundle loads and every chart still renders.
   Pinned by Task 13 (render pass with `--block-cdn` flag in the screenshot script).

---

### Task 0: Repository baseline, build log, plan commit

**Files:**
- Create: `.git/` (git init), `docs/BUILD_LOG.md`
- Modify: `.gitignore` (add `site/_screenshots/`, `tools/node_modules/`, `results/31_site_regen/` is
  under `results/` already ignored)

- [ ] Step 1: `git init`, `git add -A`, commit "Pipeline as received before the submission build (2026-09-26)".
- [ ] Step 2: Write `docs/BUILD_LOG.md` with a header and the first entry (orientation done, decisions
  above). Commit with the spec and plan: "docs: spec, plan and build log for the submission build".

### Task 1: Phase-16 identifiability recompute (`scripts/16_identifiability.py`)

**Files:**
- Create: `scripts/16_identifiability.py`, `tests/test_identifiability.py`
- Output dir: `results/16_identifiability/`

**Interfaces:**
- Produces `results/16_identifiability/nesting_TRNSCRPT.csv`, `nesting_PROT.csv`, `nesting_PHOSPHO.csv`,
  `nesting_ACETYL.csv`, `nesting_UBIQ.csv`, `nesting_METAB.csv`, `nesting_IMMUNO.csv`, `nesting_ATAC.csv`,
  `nesting_METHYL.csv` (one row per batch variable: `assay, variable, n_levels, n_tissues,
  median_levels_per_tissue, max_levels_per_tissue, max_tissues_per_level, n_levels_shared,
  tissues_in_one_level, cramers_v, n_pairs_sharing_level, n_pairs_total`);
  `estimable_pairs.csv` (`assay, n_pairs_total, n_pairs_estimable, estimable_pairs, variables_used`);
  `qc_only_per_fold.csv` and `qc_only_summary.csv` (`features, n_features, acc_mean, acc_sd,
  bal_acc_mean, bal_acc_sd, n_folds, n_test_animals_mean`); `batch_counts.csv` (`assay, variable,
  n_levels`); `bridge_variance.csv` ONLY if Task 6b is done.
- Functions: `nesting_table(meta: pd.DataFrame, tissue_col: str, variables: list[str]) -> pd.DataFrame`,
  `estimable_pairs(meta, tissue_col, variables) -> tuple[list[tuple[str,str]], int]`,
  `cramers_v(a: pd.Series, b: pd.Series) -> float`.
- Metadata sources per assay: TRNSCRPT/METHYL/ATAC = `data/raw/meta/<assay>.csv` (columns
  `RNA_extr_plate_ID|DNA_extr_plate_ID`, `Lib_batch_ID`, `Seq_flowcell_ID`, `Seq_flowcell_lane`,
  `GET_site`, `*_extr_date`, `Lib_prep_date`, `Seq_date`); PROT/PHOSPHO/ACETYL/UBIQ = `tmt_plex`,
  `tmt11_channel`; METAB/IMMUNO = `dataset` (platform) from `*_VIALS.csv`, plate_id for IMMUNO.
  Study vials only (`viallabel` starts with 9). Tissue column: TRNSCRPT `Tissue` mapped to codes
  through `data/raw/codes/tissue_abbrev.csv` where possible, else the pipeline meta join via pheno.

- [ ] Step 1: Write `tests/test_identifiability.py` with a synthetic meta (4 tissues × 2 batches nested;
  one shared batch between tissues C and D) asserting `estimable_pairs` returns `[("C","D")]` of 6 and
  `cramers_v` = 1.0 for a perfectly nested variable and ≈ 0 for an independent one.
- [ ] Step 2: Run `PYTHONPATH=src python -m pytest tests/test_identifiability.py -q` → fails (import).
- [ ] Step 3: Implement the script (functions importable; `main()` writes the CSVs; QC-only baseline
  copied from `notebooks/_build/sections/s14_qc_baseline.py` with the same feature lists, model,
  folds `grouped_kfold(om.meta, "tissue", 5, C.SEED)` on the phase-04 stacked counts matrix).
- [ ] Step 4: Tests pass; run the script (`MOTRPAC_NO_REPORT=1`); check qc_only_summary against the
  expected values (acc technical 0.8732 ± 0.005, composition 0.9488, all 0.9755; 17 plates, 4 flowcells;
  Ovary–Testes the only estimable RNA pair). Any disagreement is reported, not tuned.
- [ ] Step 5: Commit "feat(16): identifiability recompute — batch nesting, estimable pairs, QC-only baseline".

### Task 2: `--save-scores` in phase 06 and regeneration

**Files:**
- Modify: `scripts/06_conformal_certify.py` (add flag; in `panel_family` add `proba(k, X)`; in the
  fold loop when `rep == 0` and flag set, write per-fold test probabilities and calibration scores)
- Output: `results/31_site_regen/06_conformal/TRNSCRPT/` with the standard files plus
  `scores_test_probs.csv` (`fold, model ∈ {full,k20,k50}, viallabel, pid, tissue, sex, group, p_<class>×19`),
  `scores_calibration.csv` (`fold, model, calibration ∈ {pooled,one_per_animal}, viallabel, pid, tissue,
  score_lac`), `classes.json`.

- [ ] Step 1: Add the flag and writers; `python scripts/06_conformal_certify.py --help` shows it.
- [ ] Step 2: Run the exact Makefile command with `--out results/31_site_regen/06_conformal/TRNSCRPT
  --save-scores` and `MOTRPAC_NO_REPORT=1` (≈ 9 min, background).
- [ ] Step 3: Verify: `coverage.csv` of the regen equals `results/06_conformal/TRNSCRPT/coverage.csv`
  cell-for-cell within 1e-6 (report any diff); expected keys s06.* within tolerance; recompute the
  coverage from `scores_*.csv` with `conformal.py` and check it equals the regen `coverage.csv`.
- [ ] Step 4: Commit "feat(06): --save-scores exports test probabilities and calibration scores".

### Task 3: `--save-scores` in phases 12 and 13, per-draw recalibration thresholds

**Files:**
- Modify: `src/motrpac/transfer.py` (`conformal_transfer(..., collect: list | None = None)` appends
  `{"model", "n_recal", "draw", "chosen", "q_t", "n_cal_scores"}`; `calibrate_models` also stores
  `"scores"` and `"y_idx"` in each entry), `scripts/12_bodymap_validate.py`, `scripts/13_gtex_transfer.py`
- Output: `results/31_site_regen/12_bodymap/` and `13_gtex/` with the standard files plus
  `scores_target_probs.csv` (`sample, organ, stage|tissue, group_id, model, p_<class>×19`),
  `scores_calibration.csv` (`model, viallabel, pid, tissue, score_lac`), `recal_thresholds.csv`
  (`model, n_recal, draw, chosen, q_t, n_cal_scores`), `classes.json`.

- [ ] Step 1: Implement; run 12 (`--out results/31_site_regen/12_bodymap --save-scores`) and 13
  (same pattern), `MOTRPAC_NO_REPORT=1`, background.
- [ ] Step 2: Verify regen `conformal_transfer.csv`, `recalibration.csv`, `accuracy_*` equal the
  `results/` versions (1e-6) and the s09/s10 expected reference values within tolerance; recompute
  `coverage_mapped` for adults from the exported scores and compare.
- [ ] Step 3: Commit "feat(12,13): --save-scores exports target probabilities, calibration scores, per-draw thresholds".

### Task 4: Export script `scripts/30_export_site_data.py` — manifest and aggregates

**Files:**
- Create: `scripts/30_export_site_data.py`, `site/data/` outputs, `tests/test_site_data.py`

**Interfaces (JSON contracts, all files carry `"_meta": {"generated", "git_hash", "sources": [...]}`):**
- `manifest.json`: `{generated, git_hash, phases: {name: {files: [{path, bytes, sha256}]}}, absent_phases:
  [...], site_data_files: [...]}`.
- `provenance.json`: list of `{id, value, file, row, column, note}` for every headline number and
  every aggregate table cell used on the home page tiles; `pending` entries have `{id, value: null,
  reason}`.
- `headline.json`: the stat tiles and the transfer ladder:
  `tiles: [{id, value, label, sub, source}]`, `ladder: [{rung, label, model, variant, accuracy,
  accuracy_sd, coverage, coverage_sd, empty, n_samples, n_individuals, source, pending?}]`.
- `panel_curve.json`, `stable_core.json`, `confusion_motrpac.json` (k20 and full/logreg_l2),
  `bodymap.json` (age accuracy, conformal, recalibration, coverage by organ, ood sets, panel survival,
  gene check), `gtex.json` (accuracy by tissue, overall, conformal, recalibration, confusion k20/k50/full,
  gene check, panel survival, native), `representations.json` (14 + 14_cpm), `certificate.json`
  (distribution, alpha-delta grid, sizing, validity), `qc_baseline.json`, `nesting.json`,
  `batch_verdict.json` (07), `fusion.json`, `discordance.json`, `shift.json` (08 table + recalibration),
  `beyond.json` (status of phases 17, 22–26).
- Per-sample: `samples_motrpac.json` `{classes, samples: [{id, pid, tissue, sex, group, fold,
  p: {full: [...19], k20: [...], k50: [...]}}], calibration: {fold: {pooled: {full: [[class_idx,
  score],...], k20, k50}, one_per_animal: {...}}}}`; `samples_bodymap.json` `{classes, organ_map,
  samples: [{id, organ, age_weeks, animal, p: {...}}], calibration: {full: [[class_idx, score]], k20,
  k50}, recal_thresholds: {model: {n: {draw0: q or "inf", chosen: [...]}}}}`; `samples_gtex.json`
  (same shape, `tissue`, `donor`).
- `genes.json`: `{genes: [{id, symbol, marker_tissue, freq, effect_size, r_pct_mrna, regulated, qc_flag,
  fails_bodymap, fails_gtex, in_k20, in_k50, in_core, human_gene}], sets: {k20: [...], k50: [...],
  core: [...], developmental: [...]}}`; `expr_motrpac.json` `{genes: [ids], samples: [{id, tissue, sex,
  group}], values: [[...per gene]] }` (2 dp), `expr_bodymap.json`, `expr_gtex.json`.
- `conformal_fixtures.json`: Python-computed sets for 40 fixed samples × α ∈ {0.05, 0.10, 0.20} × 3
  variants × 3 models, from the quantised JSON values.

- [ ] Step 1: Write `tests/test_site_data.py`: (a) every `provenance.json` entry with a file resolves
  and the value at (row, column) matches within tolerance 1e-6 (or the entry is `pending` with a
  reason); (b) the four home tiles have provenance ids; (c) every JSON < 3 MB; (d) sample counts
  899 / 316 / 2485 and 19 classes; (e) `ladder` rungs have `n_samples` and `n_individuals`.
- [ ] Step 2: Implement the exporter in sections (`export_manifest`, `export_aggregates`,
  `export_samples`, `export_genes`, `export_fixtures`), each idempotent, with a `--skip-expr` flag
  for fast reruns. Quantise with `round(x, 6)`; `float("inf")` → the string `"inf"`.
- [ ] Step 3: Run; check every sanity anchor of the spec against the JSON (a `--check-anchors` mode
  prints a table: anchor, exported, ok?). Any mismatch is investigated, not papered over.
- [ ] Step 4: Tests pass; commit "feat(30): site data export with provenance manifest".

### Task 5: Numbers reconciliation and banners

**Files:**
- Create: `docs/NUMBERS_RECONCILIATION.md`; Modify: banners at top of `results/SUMMARY.md`,
  `results/ABSTRACT.md`, `../../docs/findings/FINDINGS_REPORT.md` (one line each, dated).

- [ ] Step 1: Generate the table from `provenance.json` plus the pre-fix values read from
  `../../backup/pipeline_history/results_pre_quantile_fix_2026-09-25/` (same file/row/column) —
  written by `scripts/30_export_site_data.py --reconciliation`.
- [ ] Step 2: Add the three banners. Commit "docs: numbers reconciliation; banners on stale documents".

### Task 6: Design system and shared components

**Files:**
- Create: `site/assets/theme.css`, `site/assets/charts.js`, `site/assets/site.js` (nav, footer,
  theme toggle, figure block, table toggle, stat tile, callout helpers, data loader),
  `site/vendor/plotly-cartesian-2.35.2.min.js`, `site/.nojekyll`, `site/assets/logo.svg` (simple mark).

**Interfaces:**
- `charts.js` exports `template()`, `palette()`, `organSystem(tissue) -> {index, name, color}`,
  `figure(el, traces, layout, opts)` → renders with the template, adds table toggle, source caption;
  `bars`, `pairedBars`, `heatmap`, `strip`, `lineBand` helpers; theme change re-renders all figures.
- `site.js` exports `loadJSON(path)`, `mountChrome(page)` (nav + footer from `manifest.json` git hash),
  `statTile`, `callout(kind, html)`, `pendingBlock(reason)`, `fmt(value, digits)`.
- Dark mode: `prefers-color-scheme` + toggle stored in `localStorage("tfp-theme")` (try/catch), and a
  `?theme=dark|light` query override used by the screenshot script.

- [ ] Step 1: Load the `dataviz` skill; write `theme.css` with every token from the spec.
- [ ] Step 2: Write `charts.js` and `site.js`; a scratch page renders one bar chart and one heatmap in
  both themes; check in the browser via Playwright screenshot.
- [ ] Step 3: Download and vendor the Plotly cartesian bundle; CDN first, `onerror` fallback.
- [ ] Step 4: Commit "site: design system, chart template and shared components".

### Task 7: Client-side conformal (`site/assets/conformal.js`) with Node test

**Files:**
- Create: `site/assets/conformal.js`, `tests/test_site_conformal.js`

**Interfaces:**
- `conformalQuantile(scores: number[], alpha) -> number (may be Infinity)`;
  `perClassQuantiles(scores, yIdx, alpha, nClasses, fallback, floor) -> number[]`;
  `predictSet(p: number[], q: number) -> boolean[]`; `predictSetConditional(p, qs) -> boolean[]`;
  `setsFor(calibration, p, variant, alpha) -> {set: boolean[], q, qs}`.
- Rule: `k = ceil((n+1)(1-alpha) - 1e-9)`; `n == 0 || k > n → Infinity`; else sorted ascending
  `[max(k,1)-1]`. LAC set: `p[j] >= 1 - q`.

- [ ] Step 1: Write `tests/test_site_conformal.js` reading `site/data/conformal_fixtures.json` and
  asserting exact equality of every set, plus unit cases: n = 19, α = 0.05 (finite, rank 19); n = 22,
  α = 0.05 (Infinity); empty class → fallback; floor rule; `p >= 1 - Infinity` is true for all.
- [ ] Step 2: `node tests/test_site_conformal.js` fails (module missing). Implement. Passes.
- [ ] Step 3: Commit "site: client-side conformal port with fixture test".

### Task 8: Home page (`site/index.html`, `site/assets/pages/home.js`)

- [ ] H1, subtitle, four tiles from `headline.json`, ladder figure with model/variant controls
  (rungs with `pending` render as hatched placeholders with the reason in the tooltip and caption),
  three sections (Accuracy: panel curve mini; Guarantee: empty-set stack mini; Identifiability:
  nesting mini), five-minute tour list (also exported to `docs/` for the final report).
- [ ] Screenshot both themes; commit "site: home page".

### Task 9: Explorer (`site/explore.html`, `site/assets/pages/explore.js`)

- [ ] (a) Tissue card: picker (source → sub-picker: MoTrPAC tissue → vial; BodyMap age → organ →
  sample incl. Thymus/Uterus; GTEx tissue → sample), random button, hide-answer toggle, α slider
  (0.05–0.30 step 0.01), model, variant, calibration (source | recalibrated draw 0 when present;
  else summary-only note). Output chips, badge Confident/Ambiguous/Abstains + one sentence, top-5
  bars, panel-gene z strip vs tissue reference profiles (from `expr_*.json`: z within dataset; the
  reference profile = mean z per tissue).
- [ ] (b) Gene explorer: picker, three panels (strip + median; BodyMap organ × age; GTEx tissue),
  annotation card; preload Pgk2 and Gnb3 (two featured buttons).
- [ ] (c) Panel builder: k slider over the grid; genes at k (from phase-12/13 all-animal panel order
  = `panel_survival` k50 list order? NO — use `selected_by_fold.csv` fold-0 selections per k, labelled
  "fold 0 of 5"), accuracy ± sd from `panel_curve.json`, 19-tissue coverage strip (marker tissue of
  each selected gene from `genes.json`).
- [ ] (d) Animals-needed calculator: α, δ sliders → ceil(ln δ / ln(1−α)); four reference points from
  `certificate.json.sizing`.
- [ ] Screenshot; commit "site: explorer".

### Task 10: Transfer page (`site/transfer.html`)

- [ ] Ladder (full controls), empty-set stack, recalibration-cost scatter, BodyMap age curve + small
  multiples of developmental markers (from `expr_bodymap.json`; panel markers vs confirmatory),
  GTEx per-tissue accuracy + row-normalised confusion with k toggle, representation result (z vs
  rank vs pairs, CPM test, ovary exception), panel-survival table, the "three donors" callout with
  the 45 % infinite-draw fact (from expected_values `s10.frac_recal_draws_infinite_k20_n3` = 0.45,
  sourced to the notebook, or recomputed from `recal_thresholds.csv` draws).
- [ ] Screenshot; commit.

### Task 11: Fingerprint, Identifiability, Beyond pages

- [ ] `fingerprint.html`: panel curve both selectors (log x, sd bands), stable core table, per-tissue
  accuracy and confusion (k20 and full), hard-tissues section from `stable_core`/`genes.json` effect
  sizes (BAT Otop1 2.09, SKM-GN Ankrd2 2.44, SKM-VL Mybph 1.26, VENACV Gdf10 1.99 log2 CPM above next
  tissue — all read from `stability_k20_annotated.csv` via `genes.json`), LIVER redundancy (six
  candidates with ovr_score 54–64).
- [ ] `identifiability.html`: framing, nesting heatmap per layer (from `nesting.json`), estimable
  pairs per layer, QC-only (technical vs composition separate) with caveat, bridge measurement
  (pending unless Task 6b), two failure modes, batch verdict 3-of-7 table with the pending-revision
  note, ATAC note, resolution = transfer.
- [ ] `beyond.html`: `beyond.json` drives "in progress" sections (no pre-registration file exists →
  state that plainly); the page structure lists what will appear when phases 17 / 22–26 land.
- [ ] Screenshots; commit.

### Task 12: Methods, Limitations, About

- [ ] `methods.html` with an inline SVG pipeline diagram (theme-aware via currentColor and tokens),
  splits, selector, quantile rule written out, certificate, transfer protocol, evaluation rules link,
  notebook self-check scheme, software versions (from `manifest.json.versions`).
- [ ] `limitations.html`, `about.html` (team placeholder clearly marked; sources + licences; CFDE
  components: GTEx and Metabolomics Workbench PR001020; reproduce; cite; licence).
- [ ] Screenshots; commit.

### Task 13: Verification pass

**Files:** `tools/screenshot.js`, `tools/package.json`, `tools/linkcheck.py`, `site/_screenshots/` (ignored)

- [ ] Playwright script: every page × {1440×900, 390×844} × {light, dark}, console errors captured
  (fail on any), `--block-cdn` mode to prove the vendored fallback works. Look at every screenshot;
  fix collisions/overflow/empty charts.
- [ ] Link check (internal resolve, external well-formed), sizes (`du`, largest JSON), load time
  (Playwright `performance.timing` under 2 s on localhost), keyboard focus (all controls are native
  elements or have tabindex), alt text on all `<img>`.
- [ ] Commit "site: verification pass fixes and tooling".

### Task 14: Repo polish

- [ ] Rename `src/motrpac` → `src/tfp` (git mv), update imports/strings (`from motrpac` → `from tfp`,
  `import motrpac` → `import tfp`, `src/motrpac` → `src/tfp`, `motrpac.` namespace references in
  scripts/tests/investigations/notebook builder/notebooks), delete `__pycache__`, run pytest, run
  `notebooks/_build/build.py --only s01_inventory --out <scratch> --execute` to prove the drift check
  passes after the rename. Commit.
- [ ] `README.md` (submission README with screenshot `figures/home.png`, key-results table generated
  by `scripts/30_export_site_data.py --readme-table`), `LICENSE` (keep MIT; add data-licence note),
  `CITATION.cff`, `docs/COMPETITION_COMPLIANCE.md`, `.github/workflows/pages.yml`, `tests.yml`,
  Makefile `site`, `site-test`, `screenshots` targets, `docs/ABSTRACT_SUBMISSION.md` (numbers via the
  export's `--abstract`), `figures/summary_figure.png` (matplotlib, 300 dpi, three panels, tokens).
- [ ] Update `CLAUDE.md` for the rename and new targets. Commit chunks; tag `hackathon-submission-v1`.

### Task 15: Final report (in chat) and BUILD_LOG close-out

- [ ] Sections 1–7 of the spec §7, with the numbers taken from `provenance.json`.

### Optional Task 6b (only if time remains after Task 13): bridge-sample variance

- [ ] From the QC table's 36 reference-standard vials and the portal per-tissue RNA count files,
  compute log2 CPM of ref vials, the between-plate variance of the ref profile over the panel genes,
  and express it as a fraction of the between-tissue variance of study vials over the same genes.
  Label it "recomputed here, definition: …" — never as the audit's 1.7 %.
