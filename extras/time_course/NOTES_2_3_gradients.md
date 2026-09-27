# Parts 2–3: molecular and covariate duration gradients (control vs 1w/2w/4w/8w)

## 0. Pre-specification (written before any fixed-arm, PC-separation or covariate number for 1w/2w/4w was looked at)

What had already been seen when this was written: `taskB_duration_summary.csv` (best-arm / best-single / best-per-omic AUROC
per tissue × duration), `batch_covariate_auroc.csv` (8w only), and the RNA plate layout (Lib_barcode_well by group, see §4).
Not seen: per-arm AUROC / log-loss for any fixed arm, any unsupervised separation, any 1w/2w/4w covariate classifier.

**Fixed arm, all tissues:** `single:TRNSCRPT/logreg_l2` (tuned C, prefilter 2000, sex covariate, same folds as Task B).
Reasons, none of them data-driven: (i) same assay and same libraries as the Part-3 RNA covariates, so the molecular-vs-covariate
comparison is on one set of vials; (ii) logistic regression gives probabilities, so held-out log-loss is meaningful (centroid
log-loss is not); (iii) one arm for all tissues avoids per-tissue selection (best-of-9 is optimistically biased).

**Decisive metric (Part 2):** unsupervised TRNSCRPT separation `D` on the Task B animals of each comparison:
log2-CPM (Task B filter) → centre within sex → top 2000 genes by variance (label-free) → z-score → PCA (5 PCs, fit on the
comparison's animals, labels unused) → `D = ||centroid_Nw − centroid_control|| / sqrt(pooled within-group mean squared
distance to own centroid)`. Null: 1000 permutations of group within sex (PCA fixed, it is label-free). Effect
`e = D_obs − mean(D_null)`, noise `s = sd(D_null)`. PROT and METAB get the same statistic (reported, not decisive).
Secondary: held-out log-loss and AUROC of the fixed arm (fold mean ± sd).

**Decisive covariate metric (Part 3):** logistic regression (as in 07_batch_check, C = 0.1) on `rna_depth_qc` = the union of
07_batch_check's `trnscrpt_depth` and `trnscrpt_qc` sets (reads, splices, RIN, duplication, mapping, composition
fractions, 5′–3′ bias), same fold construction, within-sex permutation null. The `trnscrpt_library` set is reported twice:
as defined in 07_batch_check (includes plate row / column) and without plate position, because plate position was
already known to encode sacrifice chronology (§4).

**Profile rule** (applied to `e` for Part 2 and to covariate AUROC for Part 3; weeks = 1, 2, 4, 8):
rho = Spearman(weeks, metric) over the 4 points; noise = sqrt(s_1w² + s_8w²) for `e` (null sd) and mean(fold sd at 1w, 8w) for AUROC.
1. INCREASING: rho ≥ 0.8 and metric(8w) − metric(1w) > noise.
2. DECREASING: 8w is the lowest of the four and metric(1w) − metric(8w) > noise.
3. FLAT: max − min over the 4 points ≤ noise.
4. otherwise NON-MONOTONE.
Reading: molecular INCREASING with covariate FLAT/DECREASING → training; both moving together → batch/cohort;
every 1w/2w/4w point is cohort- and season-confounded with the controls (sacrificed 4–6 months apart), 8w is not.

Deviation from §0 recorded after the fact: the first profile pass classified HEART and SKM-GN TRNSCRPT as NON-MONOTONE
because `scipy.stats.spearmanr` returned 0.7999999… for a true ρ = 0.8; the rule was applied with ρ rounded to 9 decimals
(`classify_profile`), which is the rule as written. Both verdicts hinge on ρ being exactly 0.8 (one 2w/4w swap), so treat them as fragile.

## 1. Scripts, outputs, checks

- `investigations/time_course/tc2_molecular_gradient.py` → `part2_matrix.csv` (tissue × duration, all metrics), `part2_profiles.csv`,
  `part2_pc_separation.csv` (all 3 omics), `part2_fixed_arm_per_fold.csv` (from taskB_by_duration.csv), `part2_fixed_arm_null.csv`
  (200 within-sex permutations per tissue × duration, tuned exactly as Task B; ≈ 20 s each, ≈ 10 min in total), `part2_gradient.png`.
  Check [measured]: the recomputed fixed arm reproduces taskB_by_duration.csv for all 28 comparisons (max |ΔAUROC| 1e-16, |Δlog-loss| 3e-14).
- `tc3_covariate_gradient.py` → `part3_covariate_auroc.csv`, `part3_profiles.csv`, `part3_vs_part2.csv`, `part3_gradient.png`,
  `part3_repro_check.csv`. Logistic null: 200 within-sex permutations for depth, qc and rna_depth_qc (not library sets, cost).
  Check [measured]: on 07_batch_check's own animals the 8w covariate AUROCs reproduce batch_covariate_auroc.csv in 41 of 42 rows exactly;
  KIDNEY trnscrpt_library RF differs by 0.016 (0.745 vs 0.760), not explained.
- `tc3b_depth_puzzle.py` → `depth_constants.csv`, `depth_layout.csv`, `depth_within_sex_spearman.csv`, `depth_within_sex_summary.csv`,
  `depth_cross_tissue.csv`, `depth_decomposition.csv`, `depth_composition_by_duration.csv`, `depth_layout.png`.
- Phenotype columns come from the pipeline's `data/raw/pheno.csv` (dotted names such as `vo2.max.test.vo2_max_1` and `nmr.testing.nmr_fat_2`), not from the c4.0 package.
- Part 2 uses the Task B animals, which have all three omics: 17–20 per comparison. Part 3 uses the same animals and the same fold construction.
- Permutation p has a floor of 1/201 ≈ 0.005 for 200 permutations and 1/1001 for 1000.

## 2. Tissue × duration matrix (* = cohort- and season-confounded contrast; 8w is the only date-matched one)

Columns: best-single AUROC ± fold sd (published, best of 9 arms, optimistic); fixed arm TRNSCRPT/logreg_l2 AUROC ± fold sd
(permutation p); fixed-arm held-out log-loss ± fold sd (null median, p); **TRNSCRPT PC separation e = D − null mean** (z, p);
PROT e (p); METAB e (p); **covariate-only rna_depth_qc logistic AUROC ± fold sd (p)**.

| tissue | dur | n ctrl+trained | best single | fixed AUROC | fixed log-loss | **TRN e** | PROT e | METAB e | **RNA covariates AUROC** |
|---|---|---|---|---|---|---|---|---|---|
| CORTEX | 1w* | 9+10 | 0.96 ± 0.08 | 0.67 ± 0.24 (p 0.194) | 0.72 ± 0.20 (null med 0.79, p 0.269) | 0.12 (z 1.1, p 0.149) | 0.18 (0.079) | 0.30 (0.034) | 0.65 ± 0.34 (p 0.244) |
| CORTEX | 2w* | 9+10 | 0.96 ± 0.08 | 0.96 ± 0.08 (p 0.010) | 0.59 ± 0.10 (null med 0.78, p 0.020) | 0.29 (z 2.5, p 0.016) | 0.24 (0.035) | 0.24 (0.063) | 0.83 ± 0.14 (p 0.050) |
| CORTEX | 4w* | 9+10 | 0.83 ± 0.19 | 0.50 ± 0.36 (p 0.493) | 0.78 ± 0.21 (null med 0.78, p 0.473) | 0.07 (z 0.5, p 0.266) | 0.25 (0.027) | 0.38 (0.023) | 0.62 ± 0.28 (p 0.294) |
| CORTEX | 8w | 9+10 | 1.00 ± 0.00 | 0.58 ± 0.10 (p 0.398) | 0.71 ± 0.05 (null med 0.80, p 0.229) | -0.09 (z -0.7, p 0.769) | 0.23 (0.018) | 0.07 (0.296) | 0.52 ± 0.24 (p 0.408) |
| HEART | 1w* | 10+10 | 1.00 ± 0.00 | 1.00 ± 0.00 (p 0.005) | 0.52 ± 0.06 (null med 0.82, p 0.005) | 0.37 (z 2.9, p 0.003) | 0.20 (0.071) | 0.61 (0.001) | 0.75 ± 0.22 (p 0.090) |
| HEART | 2w* | 10+10 | 1.00 ± 0.00 | 0.92 ± 0.17 (p 0.005) | 0.40 ± 0.17 (null med 0.83, p 0.005) | 0.53 (z 4.0, p 0.003) | 0.29 (0.015) | 0.74 (0.001) | 0.92 ± 0.10 (p 0.020) |
| HEART | 4w* | 10+10 | 1.00 ± 0.00 | 0.92 ± 0.17 (p 0.005) | 0.41 ± 0.19 (null med 0.83, p 0.005) | 0.52 (z 4.1, p 0.001) | 0.62 (0.001) | 0.72 (0.001) | 0.62 ± 0.28 (p 0.229) |
| HEART | 8w | 10+10 | 1.00 ± 0.00 | 1.00 ± 0.00 (p 0.005) | 0.30 ± 0.09 (null med 0.83, p 0.005) | 0.58 (z 4.4, p 0.001) | 0.65 (0.001) | 0.60 (0.002) | 0.96 ± 0.08 (p 0.005) |
| KIDNEY | 1w* | 10+10 | 1.00 ± 0.00 | 1.00 ± 0.00 (p 0.005) | 0.53 ± 0.08 (null med 0.89, p 0.010) | 0.39 (z 3.2, p 0.002) | -0.15 (0.872) | 0.49 (0.005) | 0.83 ± 0.24 (p 0.015) |
| KIDNEY | 2w* | 10+10 | 1.00 ± 0.00 | 1.00 ± 0.00 (p 0.005) | 0.53 ± 0.14 (null med 0.86, p 0.005) | 0.14 (z 2.1, p 0.033) | 0.17 (0.096) | 0.31 (0.045) | 0.67 ± 0.24 (p 0.154) |
| KIDNEY | 4w* | 10+10 | 0.83 ± 0.14 | 0.67 ± 0.27 (p 0.109) | 0.64 ± 0.19 (null med 0.84, p 0.035) | 0.22 (z 1.7, p 0.050) | 0.10 (0.259) | 0.19 (0.110) | 0.75 ± 0.10 (p 0.065) |
| KIDNEY | 8w | 10+10 | 1.00 ± 0.00 | 0.96 ± 0.08 (p 0.005) | 0.51 ± 0.14 (null med 0.83, p 0.005) | 0.34 (z 2.8, p 0.002) | 0.44 (0.003) | 0.45 (0.008) | 0.88 ± 0.25 (p 0.010) |
| LIVER | 1w* | 10+9 | 1.00 ± 0.00 | 0.92 ± 0.10 (p 0.010) | 0.42 ± 0.07 (null med 0.81, p 0.005) | 0.31 (z 2.6, p 0.004) | -0.00 (0.493) | 0.41 (0.002) | 0.81 ± 0.24 (p 0.050) |
| LIVER | 2w* | 10+10 | 1.00 ± 0.00 | 0.92 ± 0.10 (p 0.010) | 0.49 ± 0.11 (null med 0.83, p 0.010) | 0.31 (z 2.5, p 0.008) | 0.44 (0.007) | 0.66 (0.001) | 0.88 ± 0.16 (p 0.010) |
| LIVER | 4w* | 10+10 | 0.96 ± 0.08 | 0.75 ± 0.32 (p 0.055) | 0.60 ± 0.24 (null med 0.86, p 0.020) | 0.18 (z 1.6, p 0.070) | 0.06 (0.305) | 0.44 (0.002) | 0.83 ± 0.33 (p 0.020) |
| LIVER | 8w | 10+10 | 1.00 ± 0.00 | 1.00 ± 0.00 (p 0.005) | 0.29 ± 0.18 (null med 0.84, p 0.005) | 0.44 (z 4.1, p 0.001) | 0.51 (0.002) | 0.59 (0.003) | 0.59 ± 0.41 (p 0.239) |
| LUNG | 1w* | 10+10 | 0.96 ± 0.08 | 0.92 ± 0.17 (p 0.005) | 0.38 ± 0.16 (null med 0.86, p 0.005) | 0.55 (z 3.8, p 0.004) | 0.23 (0.047) | 0.19 (0.130) | 0.71 ± 0.08 (p 0.149) |
| LUNG | 2w* | 10+10 | 1.00 ± 0.00 | 0.96 ± 0.08 (p 0.010) | 0.47 ± 0.15 (null med 0.85, p 0.010) | 0.45 (z 2.9, p 0.008) | 0.35 (0.019) | 0.35 (0.031) | 0.88 ± 0.25 (p 0.025) |
| LUNG | 4w* | 10+10 | 1.00 ± 0.00 | 1.00 ± 0.00 (p 0.005) | 0.34 ± 0.05 (null med 0.85, p 0.005) | 0.47 (z 3.4, p 0.002) | 0.14 (0.131) | 0.19 (0.120) | 0.62 ± 0.32 (p 0.219) |
| LUNG | 8w | 10+10 | 1.00 ± 0.00 | 0.83 ± 0.24 (p 0.025) | 0.63 ± 0.25 (null med 0.87, p 0.040) | 0.30 (z 2.0, p 0.037) | 0.27 (0.039) | 0.40 (0.018) | 0.46 ± 0.08 (p 0.547) |
| SKM-GN | 1w* | 8+10 | 1.00 ± 0.00 | 1.00 ± 0.00 (p 0.010) | 0.27 ± 0.11 (null med 0.85, p 0.005) | 0.42 (z 3.2, p 0.002) | 0.28 (0.035) | 0.29 (0.032) | 0.77 ± 0.16 (p 0.080) |
| SKM-GN | 2w* | 8+10 | 1.00 ± 0.00 | 1.00 ± 0.00 (p 0.010) | 0.32 ± 0.10 (null med 0.80, p 0.005) | 0.57 (z 4.2, p 0.001) | 0.35 (0.016) | 0.59 (0.001) | 0.88 ± 0.14 (p 0.045) |
| SKM-GN | 4w* | 8+10 | 1.00 ± 0.00 | 1.00 ± 0.00 (p 0.005) | 0.40 ± 0.08 (null med 0.83, p 0.005) | 0.48 (z 3.5, p 0.001) | 0.59 (0.001) | 0.43 (0.013) | 1.00 ± 0.00 (p 0.010) |
| SKM-GN | 8w | 8+9 | 1.00 ± 0.00 | 0.92 ± 0.17 (p 0.025) | 0.37 ± 0.26 (null med 0.83, p 0.015) | 0.70 (z 4.8, p 0.002) | 0.69 (0.001) | 0.55 (0.006) | 0.83 ± 0.19 (p 0.060) |
| WAT-SC | 1w* | 10+10 | 0.83 ± 0.24 | 0.83 ± 0.24 (p 0.040) | 0.62 ± 0.37 (null med 1.09, p 0.030) | 0.26 (z 1.8, p 0.050) | 0.19 (0.085) | 0.23 (0.063) | 0.92 ± 0.17 (p 0.005) |
| WAT-SC | 2w* | 10+10 | 1.00 ± 0.00 | 0.75 ± 0.29 (p 0.100) | 0.83 ± 0.51 (null med 1.06, p 0.174) | 0.30 (z 2.3, p 0.017) | 0.37 (0.009) | 0.57 (0.001) | 0.50 ± 0.24 (p 0.333) |
| WAT-SC | 4w* | 10+10 | 0.96 ± 0.08 | 0.67 ± 0.36 (p 0.174) | 1.02 ± 0.94 (null med 0.99, p 0.537) | -0.11 (z -1.5, p 0.961) | 0.40 (0.002) | 0.25 (0.046) | 0.25 ± 0.32 (p 0.821) |
| WAT-SC | 8w | 10+10 | 0.96 ± 0.08 | 0.79 ± 0.25 (p 0.045) | 1.04 ± 1.09 (null med 0.94, p 0.627) | 0.24 (z 2.1, p 0.021) | 0.59 (0.001) | 0.71 (0.001) | 0.69 ± 0.14 (p 0.100) |
## 3. Profiles (pre-specified rule, §0)

| tissue | **TRN e** (decisive) | PROT e | METAB e | fixed −log-loss | fixed AUROC | **rna_depth_qc** (decisive) | depth | qc | library (with plate position) | library, no position |
|---|---|---|---|---|---|---|---|---|---|---|
| CORTEX | **DECREASING** (ρ -0.8) | FLAT | DECREASING | NON-MONOTONE | NON-MONOTONE | **NON-MONOTONE** | FLAT | NON-MONOTONE | FLAT | NON-MONOTONE |
| HEART | **INCREASING** (ρ 0.8) | INCREASING | FLAT | INCREASING | NON-MONOTONE | **NON-MONOTONE** | INCREASING | FLAT | DECREASING | FLAT |
| KIDNEY | **NON-MONOTONE** (ρ -0.2) | INCREASING | NON-MONOTONE | NON-MONOTONE | NON-MONOTONE | **FLAT** | INCREASING | FLAT | DECREASING | NON-MONOTONE |
| LIVER | **NON-MONOTONE** (ρ 0.2) | INCREASING | NON-MONOTONE | NON-MONOTONE | NON-MONOTONE | **FLAT** | FLAT | DECREASING | DECREASING | DECREASING |
| LUNG | **DECREASING** (ρ -0.8) | NON-MONOTONE | FLAT | DECREASING | FLAT | **DECREASING** | INCREASING | DECREASING | DECREASING | NON-MONOTONE |
| SKM-GN | **INCREASING** (ρ 0.8) | INCREASING | NON-MONOTONE | FLAT | FLAT | **NON-MONOTONE** | INCREASING | FLAT | DECREASING | NON-MONOTONE |
| WAT-SC | **NON-MONOTONE** (ρ -0.6) | INCREASING | INCREASING | FLAT | FLAT | **NON-MONOTONE** | NON-MONOTONE | NON-MONOTONE | DECREASING | FLAT |
## 4. Depth puzzle and plate layout

**Constant within tissue [measured, `depth_constants.csv`].** In all 7 tissues, each of these has a single value across all 50 study vials: RNA plate, RNA extraction date, library batch, kit, robot, flowcell, flowcell run, sequencing date, sequencing batch and machine. The **lane** is constant too: every library of a tissue was run on the same lane pair ("1;2" for CORTEX and LUNG, "3;4" for the others). The one exception is KIDNEY `Lib_prep_date`: one 1w vial is dated 08/30/2019 and the other 49 are dated 08/26/2019. Within control+8w, everything is constant. So lane and flowcell cannot explain anything within a tissue.

**Plate layout = sacrifice chronology [measured, `depth_layout.csv`, `depth_layout.png`].** In every tissue, the barcode well order (column-major) follows sacrifice date and time: Spearman ρ = 0.996, n = 50. The animal→relative-well map is identical across tissues for 96 % of animals.
- Control and 8w share plate columns 1–3: 3/5, 5/3 and 2/2 animals (control/8w) per column; χ² p = 0.61, n = 20.
- 4w occupies columns 3–4, 1w columns 4–6 and 2w columns 5–7.

Consequences:
- **Within control+8w [measured]:** there is no layout mechanism. The two groups are interleaved in the same columns, and within each sex the well position simply follows the sacrifice order.
- **For 1w/2w/4w vs control [measured], with the cohort caveat in the same breath:** plate position alone separates the groups. `trnscrpt_library` with row/column reaches AUROC 0.92–1.00 at 1w and 2w in every tissue, and falls to 0.15–0.75 at 8w. That is a concrete library-processing axis (plate column) that is perfectly aliased with sacrifice cohort, so edge or position effects cannot be told apart from cohort, season or duration.
- Without plate position, the library covariates show no such pattern (`trnscrpt_library_noposition`).

**What does depth separate in control vs 8w?**
- HEART is the only tissue where depth separates the groups [measured, `depth_decomposition.csv`]:
  - 8w libraries have about 10 % fewer reads: median 33.9 M vs 37.4 M, pooled AUROC 0.17, p = 0.014, n = 10+10.
  - The direction is the same in both sexes: within-sex Spearman ρ(reads, 8w) = −0.59 in each sex, n = 10, p = 0.071.
  - Spliced reads per mapped read are unchanged (0.566 vs 0.567), so this is read count, not composition.
  - The "depth" set's stronger AUROC in HEART comes from the splice-junction counts, which scale with reads (0.93).
- Within control+8w and within sex, reads were correlated with 20 variables in 14 tissue × sex cells, 280 tests, n = 10 each [measured, `depth_within_sex_summary.csv`]:
  - Variables: sacrifice day, sacrifice order, time of death, well index/row/column, RNA extraction concentration, RIN, 260/280, 260/230, library RNA/DNA concentration, molarity, fragment size, terminal body weight, NMR weight and fat % pre/post, VO2max pre/post, time to freeze, and group.
  - No variable has more than 2 of 14 cells at p < 0.05 (about 0.7 expected by chance). No variable has a consistent sign across cells; the lowest sign-test p is RIN, 10/13 positive, p = 0.09.
  - The largest single cells are consistent with multiplicity: KIDNEY males, well row ρ = −0.86 (p = 0.001); HEART females, library molarity ρ = −0.84; HEART males, sacrifice order ρ = 0.77 and pre-training VO2max ρ = −0.78.
- Depth is not a property of the animal or of the index [measured, `depth_cross_tissue.csv`]:
  - The per-animal depth rank does not carry over between tissues: median cross-tissue ρ = 0.02 over 21 pairs, n = 50.
  - This holds even for HEART–LIVER–LUNG, whose libraries use the same index sequences for the same animals (96–100 % identical index 1; ρ = 0.08–0.21).
- **Answer: no explanation found** for the lower HEART depth at 8w. It is not lane, flowcell, plate column, index, sacrifice order, RNA quality or physiology within these data [measured]. It appears only in the one date-matched contrast (HEART reads AUROC 0.56, 0.57, 0.42 and 0.17 at 1w/2w/4w/8w, the first three cohort-confounded) [measured]. A pooling or normalization step specific to that sub-plate cannot be excluded [interpreted].
- **SKM-GN "QC" separation is composition, not depth [measured]:**
  - In 8w vs control, pct_chrM is higher (AUROC 0.93; 20.0 vs 16.3 %), pct_mRNA higher (0.90) and pct_intronic lower (0.10), all p ≤ 0.003, n = 10+10. The direction is the same in each sex (5+5; p = 0.008 in males, 0.056 in females).
  - Spliced reads per mapped read fall mechanically with chrM (ρ = −0.94).
  - The chrM shift grows with duration: AUROC 0.60, 0.69, 0.93 and 0.93 at 1w/2w/4w/8w. The 1w, 2w and 4w values are cohort-confounded with the controls.
  - [interpreted] A higher mitochondrial read fraction in trained muscle is what mitochondrial biogenesis produces, so here the "QC covariate" is plausibly a readout of training, not a batch.

## 5. Per-tissue verdicts

The decisive metrics were fixed in §0: TRNSCRPT e for the molecular profile and rna_depth_qc logistic AUROC for the covariate profile. Every 1w/2w/4w statement below compares animals sacrificed 2–6 months after the controls, from different arrival cohorts and seasons.

- **HEART.**
  - The molecular separation rises with duration [measured]. TRNSCRPT e is 0.37 at 1w, 0.53 at 2w, 0.52 at 4w and 0.58 at 8w (z 2.9 to 4.4, all p ≤ 0.003, n = 20 per contrast; INCREASING by the rule, but only just: ρ = 0.8 and Δ = 0.21 against noise 0.19). The 1w–4w points are cohort- and season-confounded.
  - PROT rises from 0.20 to 0.65 (INCREASING). Held-out log-loss of the fixed arm falls from 0.52 to 0.30 (INCREASING), with 1w/2w/4w again cohort-confounded.
  - The RNA covariates do not follow this gradient [measured]: 0.75, 0.92, 0.63 and 0.96 (NON-MONOTONE), cohort-confounded except at 8w. At 8w they do separate the groups (0.96, p = 0.005), through read depth (§4, unexplained).
  - **Reading [interpreted]:** a training-consistent gradient in RNA and protein, with the cohort confound on 1w/2w/4w. The 8w RNA contrast cannot be cleanly attributed to training because of the unexplained depth difference. PROT is not subject to that covariate and shows the clearest monotone rise.
- **SKM-GN.**
  - TRNSCRPT e is 0.42, 0.57, 0.48 and 0.70 (INCREASING, ρ = 0.8, Δ = 0.28 against noise 0.20, n = 17–18) [measured]. PROT rises from 0.28 to 0.69 (INCREASING). Both carry the cohort confound at 1w/2w/4w.
  - The RNA covariates are high throughout (0.77, 0.88, 1.00, 0.83; NON-MONOTONE; cohort-confounded except at 8w, where p = 0.06). They are driven by the chrM/mRNA composition shift, which itself grows with duration (§4) [measured].
  - **Reading [interpreted]:** training, with the caveat that the "covariate" is plausibly the same biology (mitochondrial RNA fraction). The fixed-arm log-loss is FLAT (0.27 to 0.37), so the classifier metric does not show the gradient, while the unsupervised one does.
- **KIDNEY.**
  - TRNSCRPT e is 0.39, 0.14, 0.22 and 0.34 (NON-MONOTONE). 1w is the largest separation, and 1w is cohort-confounded [measured]. METAB is also NON-MONOTONE (0.49, 0.31, 0.19, 0.45).
  - PROT alone rises, from −0.15 to 0.44 (INCREASING). 1w/2w/4w are cohort-confounded.
  - The RNA covariates are FLAT and high (0.83, 0.67, 0.75, 0.88; 8w p = 0.010), cohort-confounded except at 8w.
  - **Reading [interpreted]:** no RNA duration gradient. The RNA and metabolite separation at 1w, which is at least as large as at 8w and cohort-confounded, points to cohort, season or batch contributing to kidney RNA/METAB separation. PROT is the only kidney layer consistent with a training gradient.
- **LIVER.**
  - TRNSCRPT is NON-MONOTONE (0.31, 0.31, 0.18, 0.44). PROT is INCREASING (0.00 to 0.51). 1w/2w/4w are cohort-confounded.
  - The RNA covariates are FLAT (0.81, 0.88, 0.83, 0.59) but lowest at 8w (p = 0.24), while molecular separation is highest at 8w [measured].
  - **Reading [interpreted]:** the 8w separation is the cleanest (covariates at null); the earlier cohort-confounded contrasts mix duration with cohort.
- **LUNG.**
  - TRNSCRPT is DECREASING (0.55, 0.45, 0.47, 0.30), and fixed-arm −log-loss is also DECREASING. The RNA covariates are DECREASING too (0.71, 0.88, 0.62, 0.46). 1w/2w/4w are cohort-confounded [measured].
  - **Reading [interpreted]:** molecular and covariate separations move together, which is the batch/cohort pattern. The largest lung RNA separations are in the cohort-confounded contrasts. At 8w a smaller separation remains (e = 0.30, p = 0.037; fixed AUROC 0.83, p = 0.025) with covariates at null.
- **CORTEX.**
  - TRNSCRPT is DECREASING (0.12, 0.29, 0.07, −0.09). At 8w there is no RNA separation (fixed AUROC 0.58, p = 0.40). 1w/2w/4w are cohort-confounded [measured].
  - PROT is FLAT and modest (0.18–0.25, p 0.02–0.08).
  - **Reading [interpreted]:** cortex RNA shows no training effect in the matched contrast. What RNA separation there is appears only in cohort-confounded contrasts, i.e. batch/cohort.
- **WAT-SC.**
  - TRNSCRPT is NON-MONOTONE (0.26, 0.30, −0.11, 0.24). PROT is INCREASING (0.19 to 0.59) and METAB INCREASING (0.23 to 0.71). 1w/2w/4w are cohort-confounded.
  - The RNA covariates separate most at 1w (0.92, p = 0.005, cohort-confounded) and not at 8w (0.69, p = 0.10) [measured].
  - **Reading [interpreted]:** the protein and metabolite gradients are training-consistent. RNA shows none.

**Across tissues [interpreted]:**
- The published best-single-omic AUROC is 1.00 in 19 of 28 cells, and the best-of-all-arms AUROC in 21 of 28 [measured], so neither can order durations.
- The unsupervised measure and PROT e do order them, and PROT is INCREASING in 5 of 7 tissues [measured]. PROT has no RNA-library covariate, but its own TMT plex/channel layout by duration was not checked here.
- An increasing profile over 1w→8w remains a between-cohort sequence, not a trajectory: 1w/2w/4w each carry their own cohort, season and plate-column offset.

## 6. Not done / limits

- The noise used in the rule is the permutation-null sd of D, not a bootstrap sd of the observed D. Treat the INCREASING calls for HEART and SKM-GN as marginal.
- Permutation null for the library covariate sets: not run, for cost. The library-with-position AUROC ≈ 1 at 1w/2w needs no null.
- The TMT plex × duration layout for PROT was not examined. It is the analogue of §4 for the PROT gradient and should be checked before calling PROT INCREASING "training".
