# Part 4: physiology anchor

Script: `investigations/time_course/tc4_physiology.py`. Run it from `code/pipeline` with `PYTHONPATH=src`. It takes about 50 s and is fully seeded.
Phenotype source: `data/raw/pheno.csv` (c1.0 PHENO, dotted names). The c4.0 phenotype package was not used.

## 1. Which measurement is pre and which is post [measured]

- **Pre / post order.** `_1` is the pre-training test and `_2` the post-training test. The script asserts `days_vo2_1 < training.day1_days < days_vo2_2` for every trained animal.
  - control and 8w: VO2max on day 27 vs day 84–85; NMR on day 23 vs day 85–86; training starts on day 35.
  - 4w: VO2max on day 62 vs day 91.
- **Change variables.** The script asserts these identities to 1e-6:
  - `calculated.variables.vo2_max_change` = `vo2_max_2 − vo2_max_1`;
  - `pct_body_fat_change` = `nmr_fat_2 − nmr_fat_1`, a change in percentage points;
  - `pct_body_lean_change` = `nmr_lean_2 − nmr_lean_1`.
- **1w and 2w animals have no post-training VO2max or NMR measurement** (0 of 30 in each group). They cannot be part of this analysis.

## 2. Did training work? Per-group change: median [IQR], n [measured]

| group | n | VO2max change | % fat change (pts) | lean change (pts) | NMR weight change (g) |
|---|---|---|---|---|---|
| control | 24 | −6.03 [−8.15, −1.03] | +2.20 [0.80, 3.42] | −1.30 [−2.90, −0.70] | +19.3 [14.2, 25.7] |
| 4w* | 36 | +1.96 [0.35, 3.90] | −1.20 [−2.23, −0.40] | +1.85 [1.10, 3.00] | +3.0 [−5.3, 8.4] |
| 8w | 27 | +12.22 [6.24, 14.70] | −2.20 [−4.90, 0.05] | +2.50 [0.70, 3.55] | −1.2 [−15.6, 8.4] |
| 1w, 2w | 0/30 each | — | — | — | — |

- **8w vs control, within each sex (Mann-Whitney).** This is the date-matched contrast.

  | | females (14 vs 12) | males (13 vs 12) |
  |---|---|---|
  | VO2max, median difference | +19.3, p = 1.3e-4 | +17.1, p = 4e-5 |
  | fat change, median difference | −3.35 pts, p = 1.0e-4 | −5.9 pts, p = 2.4e-5 |

  Baseline VO2max and baseline fat do not differ between 8w and control (p ≥ 0.05).
- **Training worked physiologically.**
- **Fat change differs sharply by sex within 8w:** F +0.05 [−1.0, 0.9] vs M −4.9 [−5.4, −4.5]. Because of this, the sex-adjusted correlations are the primary ones.
- *4w vs control: the 4w animals are from a different arrival cohort, were tested and sacrificed in Nov–Dec 2018 against Sep–Oct for controls, and have a shorter pre-to-post interval (29 d vs 57 d). The 4w row is not a clean dose step.

## 3. Molecular score and pre-specified test

**Score.** For each tissue × omic:
- Model: logreg_l2 with C = 0.1 fixed (no tuning), on median impute → top-2000 variance → z-score, plus a sex indicator, as in Task B.
- Trained on control vs 8w, using the Task B animals (those with all three omics, via `build_blocks` from `07_fusion_vs_baselines.py`).
- Folds: `grouped_kfold(sex_group, n_splits = min(4, n // 3) = 4)`.
- Score = the **out-of-fold logit**, averaged over 10 fold seeds.
- Secondary scores, fitted on training-fold animals only, so an animal is never in its own reference centroid:
  - projection onto the control→8w centroid axis;
  - Euclidean distance from the control centroid.

**OOF AUROC of the primary score (control vs 8w, n = 17–20).**

| tissue | TRNSCRPT | PROT | METAB |
|---|---|---|---|
| CORTEX | **0.38** | 0.89 | 0.74 |
| HEART | 1.00 | 1.00 | 0.89 |
| KIDNEY | 0.86 | 0.95 | 0.98 |
| LIVER | 0.96 | 0.97 | 0.96 |
| LUNG | 0.77 | 0.84 | 0.97 |
| SKM-GN | 0.90 | 1.00 | 0.93 |
| WAT-SC | 0.79 | 0.91 | 0.90 |

- CORTEX TRNSCRPT does not separate the groups, so a within-group score there is meaningless.
- **Score stability:** the fold-seed sd of an animal's score is a median 0.37× the between-animal sd among 8w animals (max 1.55×). The within-group ranking is therefore partly fold noise [measured].

**Primary family, declared in the script header before any correlation was run.**
- Animals: 8w only.
- Measure: sex-stratified rank correlation (ranks within sex; permutation within sex, 20,000 draws).
- Pairs: primary score vs {VO2max change, fat change}, across 7 tissues × 3 omics.
- **42 tests.** BH and Bonferroni are applied over those 42.
- Expected sign under a training response: + for VO2max, − for fat.

**Power** (`power_reference.csv`) [measured by simulation]:
- n = 10 plain Spearman: critical |rho| = 0.64 (two-sided 0.05). About 0.80 is needed for 80% power (Fisher-z approximation).
- The stratified test with 5 F + 5 M has critical |rho| = **0.70**; with 4 F + 5 M (SKM-GN) it is 0.74.

## 4. Results

**Primary: within 8w, sex-stratified rho (n = 10; SKM-GN n = 9)** [measured]

| tissue | omic | rho VO2max | p | rho fat | p |
|---|---|---|---|---|---|
| CORTEX | TRNSCRPT | +0.30 | 0.47 | −0.50 | 0.19 |
| CORTEX | PROT | −0.25 | 0.55 | +0.65 | 0.08 |
| CORTEX | METAB | 0.00 | 1.00 | +0.20 | 0.64 |
| HEART | TRNSCRPT | 0.00 | 1.00 | +0.20 | 0.65 |
| HEART | PROT | −0.30 | 0.46 | +0.65 | 0.08 |
| HEART | METAB | +0.05 | 0.94 | +0.15 | 0.73 |
| KIDNEY | TRNSCRPT | −0.10 | 0.84 | −0.10 | 0.84 |
| KIDNEY | PROT | 0.00 | 1.00 | −0.30 | 0.46 |
| KIDNEY | METAB | +0.30 | 0.46 | −0.50 | 0.19 |
| LIVER | TRNSCRPT | +0.15 | 0.74 | −0.45 | 0.25 |
| LIVER | PROT | +0.10 | 0.85 | +0.55 | 0.14 |
| LIVER | METAB | +0.25 | 0.54 | −0.40 | 0.31 |
| LUNG | TRNSCRPT | 0.00 | 1.00 | −0.30 | 0.47 |
| LUNG | PROT | +0.55 | 0.14 | −0.60 | 0.11 |
| LUNG | METAB | 0.00 | 1.00 | −0.40 | 0.31 |
| SKM-GN | TRNSCRPT | −0.34 | 0.41 | −0.18 | 0.66 |
| SKM-GN | PROT | −0.19 | 0.64 | +0.54 | 0.16 |
| SKM-GN | METAB | −0.45 | 0.27 | −0.01 | 0.96 |
| WAT-SC | TRNSCRPT | +0.10 | 0.84 | −0.20 | 0.64 |
| WAT-SC | PROT | −0.40 | 0.31 | +0.10 | 0.85 |
| WAT-SC | METAB | +0.45 | 0.25 | −0.30 | 0.46 |

**Primary family summary.**
- **0 of 42 have p < 0.05.** The smallest p is 0.076: CORTEX PROT and HEART PROT fat. Both run in the *opposite* direction to a training response (+0.65).
- Min BH q = 0.94. Bonferroni gives 1 for all.
- 23 of 42 are in the expected direction (binomial p = 0.64). The mean rho in the expected direction is +0.03.
- No |rho| reaches the 0.70 critical value.

**Secondary results.** Counts are of p < 0.05 within each family of 126 tests: 7 tissues × 3 omics × 6 variables (four changes plus the two baselines). About 6.3 are expected by chance. See `score_physiology_correlations.csv`.
- **Pooled control + 8w (uninformative by design):** 55–63 of 126 at p < 0.05. For the primary pairs it is 34 of 42, with median rho in the expected direction 0.62 (n = 17–20).
  - This only restates that the score separates the groups and training changes physiology.
  - **Any group separator would do this, batch included, so it is not evidence for training over batch** [interpreted].
- **Within 8w, unadjusted Spearman:** 4 of 126. The unadjusted within-8w rho leans in the expected direction (mean +0.18). Sex-stratified it is +0.02, so the lean is the sex difference in fat change, not individual response [interpreted].
- **Negative control, within controls** (same model, same score): 9 of 126 (sex-stratified). For the primary pairs it is 3 of 42, all *negative* rho:
  - KIDNEY METAB VO2max: −0.75, p = 0.033;
  - LIVER PROT fat: −0.73, p = 0.037;
  - SKM-GN PROT VO2max: −0.79, n = 8, p = 0.043.

  This is the chance rate, the same as the 8w arm.
- **Within 4w** (control-vs-4w model; the cohort caveat above applies): 1 of 42 primary pairs at p < 0.05. WAT-SC METAB fat is −0.78 (n = 10, p = 0.022), with q not significant. 24 of 42 are in the expected direction.
- Centroid projection and distance scores show the same null pattern (2–13 of 126 within 8w).

Figures:
- `score_vs_physiology_8w.png`: 8w animals, score vs change, 7 tissues × 3 omics × 2 variables.
- `within_group_rho_forest.png`: primary rho for 8w vs controls, with the 0.70 critical band.

## 5. Verdict

- **No tissue shows a within-group link between score and physiology that n = 10 can distinguish from zero** [measured]. This holds for all 42 primary tests: none reaches the 0.70 critical value, and none is significant after correction. The direction is at chance level and matches the negative control.
- **This is absence of evidence, not evidence of absence** [interpreted]:
  - At n = 10, only true |rho| ≈ 0.8 would be detected reliably.
  - Fold noise in the OOF score is a sizeable fraction of the between-animal spread.
- **It therefore cannot adjudicate training vs batch** [interpreted]. A real but moderate individual-response link (|rho| of about 0.3–0.5) would look exactly like this.
- **The pooled correlations are large, but group membership explains them** [interpreted]. They add nothing to the batch question.

**What would give more power (not run):**
- Use all animals with TRNSCRPT only (about 13–14 8w animals per tissue).
- Pool tissues per animal, for example the mean within-sex score rank across tissues, which gives one test with n ≈ 27 8w animals. That has about 80% power for rho of about 0.5.
