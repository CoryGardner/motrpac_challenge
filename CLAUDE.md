# CLAUDE.md — pipeline rules (code/pipeline, formerly "the starter kit")

A tested, leakage-safe evaluation pipeline for the MoTrPAC 6-month rat endurance-training data
(R package `MotrpacRatTraining6moData` 2.0.0 = portal release c1.0, exported to `data/raw/`).
It produced the pre-hackathon findings in `../../docs/findings/` (tissue fingerprints with conformal
guarantees, shift tests, BodyMap/GTEx transfer, fusion, discordance); current numbers: `results/SUMMARY.md`.

Read `docs/EVALUATION_RULES.md` and `docs/DATA_GUIDE.md` before writing analysis code here.
Run with the `motrpac-py` conda env: `PYTHONPATH=src python -m pytest -q`, `make <phase> PY=python`.

## Non-negotiable evaluation rules (short form)

1. **Split by animal (`pid`), never by sample.** Every sample from one animal is on one side of
   every split. Use `motrpac.splits`. Tests in `tests/test_splits.py` must pass.
2. **Everything data-dependent happens inside the training fold**: scaling, imputation,
   feature selection, hyperparameter tuning. Use sklearn `Pipeline`s; never pre-select
   features on the full dataset.
3. **Every model is compared against a tuned simple baseline** (nearest centroid, L1/L2 logistic
   regression, random forest) under the *same* splits. Report the baseline first.
4. **Report uncertainty**: per-fold results, mean ± sd, and the number of animals in each fold.
   Small n is the norm here (~50–60 animals per tissue).
5. **Conformal calibration sets come from held-out animals**, disjoint from fitting and test animals.
   Report coverage *and* average set size; report per-tissue (conditional) coverage too.
6. **Shift experiments are first-class**: train on one sex / time point / training state and
   test on the other. Report coverage drop, not just accuracy drop.
7. **Say what didn't work.** Negative results go in `results/REPORT.md` with the same care.

## Data conventions

- Exported CSVs live under `data/raw/` (see `docs/DATA_GUIDE.md` for the exact layout).
  Sample-level files are features × samples with 4 leading columns
  (`feature`, `feature_ID`, `tissue`, `assay`) then one column per `viallabel`.
- `viallabel` (11-digit string) identifies a vial; `bid = viallabel[:5]` is the collection-event ID,
  which in this study is one per animal (147 bids ↔ 147 pids); `pid` (8-digit) identifies the animal.
  Join assays by `bid`/`pid`, never by `viallabel`. Reference-standard vials start with 8, study vials with 9.
- Keep IDs as strings everywhere (`dtype=str`); leading zeros and int overflow have bitten people.
- Tissue codes use hyphens in data (`SKM-GN`, `WAT-SC`) but not in R object/file names (`SKMGN`, `WATSC`).
  `motrpac.config.tissue_token()` / `tissue_from_token()` convert.
- Group labels: `control`, `1w`, `2w`, `4w`, `8w`. Sex: `male`, `female`.
- Use `motrpac.io.load_norm()` / `stack_tissues()` / `align_by_animal()`; don't hand-roll readers.

## Coding conventions

- Python 3.10+, numpy/pandas/scikit-learn/scipy/matplotlib only in the core path;
  anything else is optional and must be import-guarded.
- Scripts are numbered, idempotent, and accept `--quick` for fast runs; they write to
  `results/<phase>/` and append a short section to `results/REPORT.md` via `motrpac.report`.
- Seed everything with `motrpac.config.SEED`. Set `n_jobs=-1` where sklearn allows it.
- No notebooks in the core path; results must be reproducible from a script before they go in the report.
- Figures: matplotlib, PNG at 150 dpi, one idea per figure, axis labels with units/scales.
- Run `pytest -q` after touching `src/motrpac/`.

## Things not to do

- Don't fit anything on the full dataset "just to look" and then reuse those features/models.
- Don't evaluate tissue classification with plain `KFold` or `train_test_split`.
- Don't compare proteomics values across tissues without checking `docs/DATA_GUIDE.md` on
  normalization (per-tissue reference channels may make cross-tissue proteomics non-comparable).
- Don't silently drop samples or features; log counts before/after every filter.
- Don't add dependencies without a fallback path; the hackathon laptop may be offline for stretches.
- Don't invent results. If a phase can't run (missing data, missing R), say so in the report.

## Status

The pre-hackathon run is complete (`make all` + `bodymap gtex transfer`; clean run 2026-09-18, ~1 h).
**Conformal quantile fixed 2026-09-25.**
- `conformal.conformal_quantile` is now the textbook rule, the ceil((n+1)(1−α))-th smallest score, and +inf, a full set, when that rank exceeds n.
- Tests for the three α = 0.10 regimes are in `tests/test_conformal.py`.
- Phases 06, 08 and 12–14 were rerun and `SUMMARY.md` was rebuilt.
- What moved is in `results/QUANTILE_FIX_CHANGES.{md,csv}`; `scripts/quantile_fix_diff.py` makes the diff.
- The pre-fix results are in `../../backup/pipeline_history/results_pre_quantile_fix_2026-09-25/`.

**Known deviation, documented but not fixed** (`../../docs/findings/FINDINGS_REPORT.md` §14, 13(b)): the transcript gene filter runs before the split.

**Replication notebooks.**
- `notebooks/01_replication.ipynb` and `02_transfer.ipynb` are self-contained: the library is pasted in verbatim.
- `notebooks/_build/build.py` builds them from `notebooks/_build/sections/`. Rebuild after changing `src/motrpac/`.
Organizer questions raised by the data: `../../docs/findings/QUESTIONS_FOR_ORGANIZERS.md`.
