# MoTrPAC evaluation pipeline (formerly "the starter kit")

Leakage-safe, calibrated evaluation of the MoTrPAC 6-month rat endurance-training data: animal-grouped
splits, tuned simple baselines, compact tissue panels, conformal prediction sets and certificates, shift
tests, BodyMap/GTEx transfer, fusion and RNA–protein discordance. Pre-hackathon results are complete:
`results/SUMMARY.md` (numbers) and `../../docs/findings/` (write-up). Rules: `CLAUDE.md`, `docs/EVALUATION_RULES.md`.

## Quick start

```bash
conda activate motrpac-py                      # see ../../docs/setup/ for the environment
cd code/pipeline
PYTHONPATH=src python -m pytest -q             # 18 leakage / I/O / conformal tests
make inventory QUICK=1 PY=python               # any phase: make <phase> [QUICK=1]
```

`data/raw/` already holds the CSV export of `MotrpacRatTraining6moData` 2.0.0 (made once with
`Rscript R/export_motrpac.R data/raw`, which needs the `motrpac-r` env); `data/external/` holds GTEx v8
and the rat BodyMap. `make synthetic && make all QUICK=1` smoke-tests on synthetic data (writes to
`data/raw/_SYNTHETIC`; delete it afterwards).

## Layout

```
CLAUDE.md               Project rules Claude Code reads automatically (evaluation rules, conventions)
README.md               This file
Makefile                make <phase>; QUICK=1 for small/fast runs
requirements.txt        Python deps (core deps only; optional extras listed inside)
R/                      One-time export from MotrpacRatTraining6moData → CSV
src/motrpac/            Library: io, splits, models, conformal, discordance, plots, report
scripts/                Numbered phases 00–10, one script each
docs/                   DATA_GUIDE, EVALUATION_RULES, EXTERNAL_VALIDATION, GTEX_TRANSFER
tests/                  Leakage and I/O tests (pytest)
data/raw/               R-package CSV export (1.9 GB)
data/external/          GTEx v8 + rat BodyMap (5.8 GB)
results/                Outputs (git-ignored)
```

## How the phases fit together

| Phase | Script | Question it answers |
|---|---|---|
| 01 | `R/export_motrpac.R` | Get the data out of R into flat CSVs, one file per assay × tissue |
| 02 | `scripts/02_inventory.py` | What is actually there: tissues × assays × samples, animal overlap across assays |
| 03 | `scripts/03_eda.py` | Does tissue dominate variance? Sex? Time point? Any batch structure? Missingness? |
| 04 | `scripts/04_fingerprint_baselines.py` | Tissue classification with animal-grouped CV; tuned simple baselines; confusable pairs |
| 05 | `scripts/05_compact_panels.py` | How small can a panel be? Stability of selected features across folds |
| 06 | `scripts/06_conformal_certify.py` | Prediction sets with coverage guarantee; certified panel size (LTT-style) |
| 07 | `scripts/07_fusion_vs_baselines.py` | Does multi-omic fusion beat the best single-omic tuned baseline? |
| 08 | `scripts/08_shift_tests.py` | Does the fingerprint/guarantee survive a held-out sex, time point, or training state? |
| 09 | `scripts/09_discordance.py` | Transcript–protein discordance: how much, where, and is it predictable? |
| 10 | `scripts/10_make_report.py` | Collect everything into `results/REPORT.md` |

## Lessons already baked in (from the synthetic smoke test)

- Univariate F-test feature selection is the wrong tool for multiclass panels: the top-k F-scores
  are k markers of the same easy tissue, so a "20-gene panel" left most tissues unrepresented
  (~0.3 accuracy). The default selector is class-aware round-robin (`RoundRobinSelector`);
  `make panels` runs both so the failure is documented in the report.
- The certified panel size is larger than the point-estimate size (k = 15 vs 8 on synthetic data
  at α = δ = 0.05) — that gap is the honest-reporting message, not a bug.
- LAC prediction sets can be empty at large α when the classifier is confident; APS never are but
  are larger. Report both.
- Sex-specific tissues (OVARY, TESTES) make the held-out-sex experiment a label-shift problem;
  `08_shift_tests.py` reports unseen classes separately.

## What is verified vs. what to confirm

Verified against the MoTrPAC documentation (Sept 2026):
- Rat data are public via the `MotrpacRatTraining6moData` R package (v2.0.0), installed from GitHub.
- Object naming (`<ASSAY>_<TISSUE>_NORM_DATA`, `..._RAW_COUNTS`, `..._DA`, `PHENO`, `FEATURE_TO_GENE`,
  `RAT_TO_HUMAN_GENE`) and the sample-level layout (4 leading columns, then one column per viallabel).
- Study design: 6-month-old rats, 1/2/4/8 weeks of training vs. sedentary controls, both sexes,
  ~3–6 animals per sex per time point per assay, 18 solid tissues + blood + plasma.

Since verified (see `../../docs/data/`): the package equals portal release c1.0; human acute-exercise
results are public at summary level (`MotrpacHumanPreSuspensionAnalysis` 0.2.4); cross-tissue proteomics is
not comparable (per-tissue TMT references; missingness identifies tissue). Still to confirm with organizers:
which data release and slice the hackathon uses (`../../docs/findings/QUESTIONS_FOR_ORGANIZERS.md`,
`../../docs/scoping/RECOMMENDATIONS.md` §7).

## License

MIT — the hackathon requires open-source outputs.
