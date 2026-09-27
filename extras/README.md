# extras — additional analyses, not part of the submission

These scripts and notebooks were written during the project, but the submission site does not use their outputs.
They are **not covered by the tests or CI** and are kept for completeness and reuse. They import the library from
`src/tfp` (`pip install -e .` or `PYTHONPATH=src`) and read `data/` and `results/` from the repository root.

| path | what it does | needs |
|---|---|---|
| `scripts/09_discordance.py` | transcript–protein discordance of the training response per tissue (three definitions); writes `results/09_discordance/` | `data/raw`, `tfp_extras/discordance.py` |
| `scripts/10_make_report.py` | a one-page `results/SUMMARY.md` from every phase's CSVs | `results/` |
| `scripts/12_bodymap_prepare.py` | the GEO GSE53960 archive as a fallback expression matrix; the submission uses the Bioconductor package `bodymapRat` (`R/export_bodymap.R`) | the GEO supplementary archive |
| `scripts/quantile_fix_diff.py` | cell-by-cell diff of two result trees (used for the 2026-09-25 conformal-quantile correction) | `--snapshot DIR` (a pre-fix results tree) |
| `tfp_extras/discordance.py` | the discordance helpers formerly in `src/tfp` | — |
| `time_course/tc2_molecular_gradient.py` | training-duration gradients of the within-tissue molecular separation (pre-specified; `NOTES_2_3_gradients.md`) | `data/raw`, `results/07_fusion`, joblib; ≈ 10 min |
| `time_course/tc3_covariate_gradient.py`, `tc3b_depth_puzzle.py` | the same gradients on processing covariates only; the sequencing-depth difference in heart | `data/raw`, `results/07_fusion` |
| `time_course/tc4_physiology.py`, `tc4b_pooled_across_tissues.py` | physiology (VO2max, body composition) as an anchor for the training response; within-group correlations (`NOTES_4_physiology.md`) | `data/raw` |
| `notebooks/` | `01_replication.ipynb`, `02_transfer.ipynb`: narrative notebooks built by `_build/build.py` from `_build/sections/`, executed on 2026-09-26 against the results of that date, with `expected_values.csv` (287 self-check keys). **Frozen:** the library pasted into them predates later changes to `src/tfp`, and three cells read files outside the repository. Read them; do not expect them to re-execute. | nbformat, jupyter |

Why these stay out of the site: the duration gradients compare animals sacrificed months apart (cohort, season and
plate position travel with training duration), the depth puzzle is one cell of the covariate figure the site does
show, the within-group physiology correlations are an underpowered null (n = 10 per group), and the discordance
analysis is about the exercise response rather than tissue identity. The design dates, the physiology group tests and
the fingerprint-by-duration analysis from the same investigation were promoted into the submission path
(`scripts/15_time_course_design.py`, `scripts/15_fingerprint_by_duration.py`) and appear on the site's Exercise page.

Promoting a script: `git mv extras/<path> scripts/NN_name.py`, add a Makefile target and, if the site reads its
outputs, provenance entries in `scripts/30_export_site_data.py`; the frozen snapshot then picks the new files up.
