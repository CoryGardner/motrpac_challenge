# Conventions for replication-notebook section files

Repo: `/home/cory/projects/MoTrPAC`. Pipeline: `code/pipeline/` (called PIPE below). Python env:
`~/miniconda3/envs/motrpac-py/bin/python` (never base python). Kernel: `motrpac-py`.

## What we are building
Two self-contained Jupyter notebooks that replicate the pipeline's pre-hackathon findings as a readable
narrative: `01_replication.ipynb` (sections 1–8, 11–14) and `02_transfer.ipynb` (sections 9–10).
`_build/build.py` assembles them from **section files** in `_build/sections/sNN_name.py`, after a
header (`s00_header.py`) and the pipeline library pasted in verbatim. Read both files first.

## Hard rules (from the user — binding)
1. **Self-contained:** never `import motrpac` and never `sys.path` into `src/`. The library is already in
   the namespace (see below). Script logic you need (from `PIPE/scripts/*.py`) is **copied into the
   section**, verbatim where possible, trimmed only of argparse/CLI, `report.add_section`, and branches
   the notebook does not use. Every copied block starts with a provenance comment:
   `# from scripts/08_shift_tests.py::eval_split (verbatim)` or `(trimmed: dropped --quick branch)`.
2. **No numbers in markdown.** Every number the reader sees is printed by a code cell — computed in this
   run, or read from a `results/` CSV in this run. Markdown says "the table above", "far above chance",
   etc. Parameters of the design (k = 20, α = 0.10, 5 folds, 100 splits) may be stated in markdown;
   result values may not.
3. **Never write to `PIPE/results/`, `PIPE/data/`, or the repo's `data/`.** All outputs go under `OUT`
   (`notebooks/_outputs/`), e.g. `OUT / "s05"`. Note: the library's `C.RESULTS_DIR` points at
   `PIPE/results`; library plot functions take an explicit path — pass one under OUT.
4. **Do not modify** anything outside your own section files and your own `_build/expected/sNN.csv`.
   Not `src/`, not `scripts/`, not `docs/EVALUATION_RULES.md`, not `build.py`, not `s00_header.py`.
   If you think those need a change, say so in your final report.
5. **No cell may run over 5 minutes** in the default mode (RECOMPUTE=False). Do not execute any
   RECOMPUTE=True branch expected to take > 3 min; write it carefully, check it parses, and say in your
   report that it was not executed. A background pipeline rerun may be using CPU; keep `n_jobs` as the
   scripts use it.
6. **Do NOT touch the pre-split gene filter** (the transcript gene filter runs before the split — a
   known deviation). Reproduce the pipeline as it is.
7. Never invent numbers. If something does not reproduce, print both values and report it; do not tune
   it until it matches.

## Section shape (every section)
```
# %% [markdown]
# ## N. Title
# **Question.** … **Why it matters.** … **What to look for.** …
# %%
section("N title")          # first line of the first code cell: starts the section clock
... commented code ...
# %% [markdown]
# **What this shows.** … **What it does not show.** …
```
Several (b)/(c) pairs inside one section are fine (sub-questions). Keep cells short enough to read
(roughly ≤ 60 lines; library-style helper cells may be longer). Figures: show inline (`plt.show()`),
or write a PNG under OUT with the library's `plots.*` functions and `show_png(path)`.

## RECOMPUTE
`RECOMPUTE` (bool) is defined in the header. Expensive legs — the 100-split certificate, stability
selection (50 bootstraps), fusion tasks, batch-covariate permutation nulls, BodyMap / GTEx /
representations — follow this pattern:
```python
if RECOMPUTE:
    ... full computation (copied script logic) ... ; tbl.to_csv(OUT / "s07" / "certificate_validity.csv")
else:
    tbl = pd.read_csv(res("06_conformal", "TRNSCRPT", "certificate_validity.csv"))
```
Both branches must produce the same variable with the same columns. Everything cheap is recomputed live
in both modes (and where the pipeline wrote the same number to results/, print both side by side:
"recomputed vs published").

## Names available in the namespace (from the header + library cells)
- Header: `np, pd, plt, matplotlib, json, os, re, time, warnings, itertools, Path, SimpleNamespace,
  dataclass, field, Iterable, Iterator, Image, Markdown, display`; `PIPE, REPO, RES, NB, OUT,
  RECOMPUTE`; `section(name)`, `record(key, value, section_id, note="")`, `res(*parts) -> Path`
  (a results file, raises if missing), `show_png(path)`, `cached(key, fn)` (load a big matrix once per
  run — use shared keys: `"pheno"` → `io.load_pheno()`; `"TRNSCRPT_counts"` → the exact
  `io.stack_tissues(...)` call that scripts/04_fingerprint_baselines.py makes for TRNSCRPT).
- Library namespaces (verbatim copies of `src/motrpac/*.py`): `C` (config), `io`, `splits`, `models`,
  `cp` (conformal), `discordance`, `transfer`, `plots`. Every top-level name of those modules is also
  defined bare (e.g. `OmicsMatrix`, `grouped_kfold`). `cli.py` and `report.py` are **not** included:
  copy `cli.resolve_source` / `cli.parse_tissues` inline if you need them.
- Do not rebind any of the names above. Prefix section-level variables that later sections might
  collide with (e.g. `s05_curve`), or keep them obviously local.

## Self-check values
For each headline number, call `record("sNN.short_key", value, "NN")` in the cell that prints it, and
add a row to `_build/expected/sNN.csv` with columns
`key,section,published_value,reference_value,tolerance,source,note`:
- `published_value`: the value in `PIPE/results/SUMMARY.md` / the published results CSV **before the
  2026-09-25 conformal-quantile fix**. For conformal / shift / BodyMap / GTEx / transfer numbers, the
  pre-fix files are in `REPO/backup/pipeline_history/results_pre_quantile_fix_2026-09-25/`.
- `reference_value`: the value the notebook should reproduce now (post-fix results files where the
  quantile fix moved it; equal to published otherwise; for the NEW section 14, this run's value).
- `tolerance`: 0 for counts and values loaded from CSVs (use 1e-9); 0.005 for accuracies/coverages
  recomputed with the same seed; wider only where the pipeline itself is stochastic — say why in `note`.
  Never widen a tolerance to make a mismatch pass; report the mismatch instead.
- `source`: the results file (relative to PIPE) or SUMMARY line the reference comes from.

## Testing
```
cd PIPE/notebooks
~/miniconda3/envs/motrpac-py/bin/python _build/build.py --only s04_baselines,s05_panels \
    --out <your scratch dir>/test.ipynb --execute
```
(`--only` prepends the header and library.) Read the executed notebook's outputs (nbformat) to check
every cell ran and printed what you expected. Time each section (the `_TIMES` dict / `timing_table()`).

## Your final report (keep it under ~400 words)
Files written; cells per section; measured runtime per section (RECOMPUTE=False); which headline
numbers reproduced (recomputed vs published) and which did not, with both values and your best
explanation; code the notebook needed that does not exist in the pipeline; anything you think needs
the lead's decision.
