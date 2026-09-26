# EVALUATION_RULES — how every number in this project gets produced

These rules are the project's contribution. If a result is not produced this way it does not
go in the report.

## 1. Units of splitting

- The unit is the **animal (`pid`)**. `motrpac.splits.grouped_kfold()` wraps
  `StratifiedGroupKFold` (stratified on the label, grouped on `pid`).
- Default: 5 outer folds. With ~50–60 animals per tissue that is ~10–12 test animals per fold —
  report the count.
- Shift splits are deterministic: `leave_one_sex_out`, `leave_one_group_out` (time point),
  `train_controls_test_trained`.
- `assert_no_group_leak()` is called inside every evaluation loop and in tests.

## 2. What happens inside a fold

Order, all fit on the training fold only:
1. Impute (median) → 2. scale (z-score) → 3. feature selection → 4. classifier.
Implemented as one sklearn `Pipeline` so it is impossible to leak by accident.
Hyperparameters are tuned with an **inner** grouped CV (`GridSearchCV(cv=GroupKFold)`),
also on the training fold only.

## 3. Baselines (always reported first)

| Name | Model | Why |
|---|---|---|
| `centroid` | nearest (shrunken) centroid on top-k F-score features | the "one marker gene per class" floor |
| `logreg_l1` | L1 logistic regression (saga), C tuned | sparse, interpretable, hard to beat on this data |
| `logreg_l2` | L2 logistic regression, C tuned | dense linear reference |
| `rf` | random forest, 500 trees, `max_features` tuned | the nonlinear reference; what fusion must beat |

A "fancier" model earns a line in the report only if it beats the best tuned baseline on the
same outer folds by more than the fold-to-fold sd.

## 4. Metrics

- Multiclass: balanced accuracy, macro-F1, log-loss; confusion matrix (which pairs get confused).
- Binary: AUROC, AUPRC, balanced accuracy at 0.5.
- Always per fold; summarize as mean ± sd; keep the per-fold table in `results/`.

## 5. Compact panels

- Panel sizes on a grid: 1, 2, 3, 5, 8, 10, 15, 20, 30, 50, 100 features.
- Selection happens inside the fold (`SelectKBest(f_classif)` or L1 coefficients).
- **Stability**: run selection on bootstrap resamples of *animals*; report each feature's
  selection frequency. A panel is only "the panel" if its members are stable.
- Report the full accuracy-vs-k curve, not just the best k.

## 6. Conformal prediction sets

- Split conformal. Per outer fold, the training animals are split again into **fit** and
  **calibration** animals (default 70/30 by `pid`). Test animals are never used for either.
- Nonconformity scores: LAC (`1 − p̂_y`) by default; APS available.
- Report at α ∈ {0.05, 0.10, 0.20}: empirical coverage, average set size, fraction of singleton
  sets, and **per-tissue coverage** (conditional coverage is where things go wrong).
- Sanity check: coverage on in-distribution test animals should be ≈ 1 − α. If it is far above,
  the sets are trivially large; if far below, something leaked or shifted.

## 7. Certified panel size (LTT-style)

The question "what is the smallest panel whose error is ≤ α with confidence 1 − δ" is a
multiple-testing problem over panel sizes. `motrpac.conformal.certify_panel_size()`:
1. For each k on the grid, fit the panel on fit animals, measure errors on calibration animals.
2. Compute a Clopper–Pearson upper confidence bound on the error at level δ.
3. **Fixed-sequence testing from large k to small k**: walk down the grid and stop at the first
   k whose bound exceeds α; the previous k is certified. (Fixed-sequence testing controls the
   family-wise error without a Bonferroni penalty, as long as the order is fixed in advance.)
4. Report the certified k, its bound, and the calibration n. With n_cal ≈ 15 animals × ~18
   tissues the bounds are loose; say so.

## 8. Shift experiments (the interesting part)

For each of: held-out sex, held-out time point (8w; also 1w), controls → trained:
- Train + calibrate on the source; test on the target.
- Report: accuracy drop, **coverage drop at fixed α**, and set-size change.
- A coverage drop is the headline; the accuracy drop is context.

## 9. Fusion

- Early fusion: z-scored blocks concatenated (block-wise scaling so one omic can't dominate by
  feature count). Late fusion: average of per-omic calibrated probabilities.
- Same outer folds as the single-omic baselines; same animals; report the single-omic winner
  alongside each fusion row.
- Restrict to animals with all blocks present (report how many were dropped by that).

## 10. Reporting

- `motrpac.report.add_section()` appends to `results/REPORT.md` with a timestamp, the git
  hash (if any), the script name, and the exact parameters.
- Negative results are written up the same way as positive ones.
- Every table in the report has n (animals) next to the metric.
