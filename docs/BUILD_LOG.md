# Build log — submission site, explorer and repo

Times are local (Spark, UTC-5 host clock as reported by `date`). Every entry: done / next / blockers.

## 2026-09-26 17:30–18:20 — orientation and plan
**Done.** Read CLAUDE.md, EVALUATION_RULES, DATA_GUIDE, all of `src/motrpac`, scripts 04/05/06/12/13, the
result tables of phases 02–15, `QUANTILE_FIX_CHANGES.md`, `expected_values.csv` (published = pre-fix,
reference = post-fix). Inventoried `results/`: phases 02–09, 12–14 (+14_cpm), 15_time_course present;
16–26 absent; no pre-registration file anywhere in the workspace. Confirmed the toolchain: motrpac-py
(pandas 3.0.6, sklearn 1.9.1), Node 20, system Chrome + Playwright (npm, scratch install works).
Probed the RNA batch metadata: 17 plates, 17 library batches, 4 flowcells; only Ovary–Testes shares a
level of every processing variable → "1 of 171" reproduces, and it is the sex contrast.
Initialised git (`main` = pipeline as received), branch `submission-build`. Wrote the spec (verbatim
brief) and the plan under `docs/superpowers/`.
**Next.** Task 1 (phase-16 identifiability recompute), Tasks 2–3 (`--save-scores` in 06/12/13, regen runs
in the background), Task 4 (export script + provenance test).
**Blockers.** None. Note: the in-distribution coverage headline is 0.908 (pooled) / 0.916 (one vial per
animal) post-fix, not the 0.96 in the brief; the site will use the post-fix values.
