# Build log — submission site, explorer and repo

Times are local (Spark, UTC-5 host clock as reported by `date`). Every entry: done / next / blockers.

## 2026-09-26 17:30–18:20 — orientation and plan
**Done.** Read CLAUDE.md, EVALUATION_RULES, DATA_GUIDE, all of `src/tfp`, scripts 04/05/06/12/13, the
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

## 2026-09-26 18:20–19:20 — Phase A: recomputes and regenerations
**Done.** `tfp.batch` + `scripts/16_identifiability.py` (results/16_identifiability/): RNA-seq 1 of 171 tissue pairs
estimable (OVARY|TESTES); METHYL/ATAC 0 of 28; TMT layers 0 (plex = tissue × label); immunoassay plates hold 2–4 tissues
(16 of 136 pairs). QC-only baseline reproduced the notebook exactly (0.8732 / 0.9488 / 0.9755). `--save-scores` added
to phases 06, 12, 13 (+ `transfer.save_transfer_scores`, per-draw recalibration thresholds); regeneration runs under
results/31_site_regen/ reproduce the published tables cell-for-cell (7 tests). Export script `scripts/30_export_site_data.py`
written; JS conformal port + fixture test written.
**Next.** Run the export, check anchors, reconciliation doc + banners (Task 5), then the design system (Task 6).
**Blockers.** None.

## 2026-09-26 19:20–2026-09-27 00:10 — Phase B and C: design system and all nine pages
**Done.** Export run: 235 provenance entries, 84 copied tables, all 55 sanity anchors reconcile (the spec's
"heart → skeletal muscle 0.907" is the SKM-GN + SKM-VL fraction, exported as such); `docs/NUMBERS_RECONCILIATION.md`
(102 of 182 comparable numbers moved with the fix) and banners on the three stale documents. Design system
(theme.css, charts.js template, site.js components, vendored Plotly 2.35.2), the JS conformal port (2,713
assertions, 540 infinite-threshold cases), and pages: Home, Explore (tissue card with client-side sets incl.
BodyMap thymus/uterus and recalibrated thresholds, gene explorer, panel builder, calculator), Transfer,
Fingerprint, Identifiability, Beyond, Methods, Limitations, About. Render pass: 36 renders (9 pages × 2 sizes ×
2 themes), zero console errors, no horizontal overflow. Fixed along the way: a protein-id symbol mapping in the
gene export (Pgk2 was NP_001012130.1), select/caption overflow on phones, the duplicate symIdx in transfer.js.
**Next.** Task 13 (link check, offline mode, sizes, load time, a11y), Task 14 (rename src/motrpac → src/tfp,
README, LICENSE, CITATION, compliance doc, workflows, Makefile targets, abstract, summary figure), final review.
**Blockers.** None. Optional 6b (bridge variance) still open.

## 2026-09-27 00:10–00:40 — verification, repo polish, bridge measurement
**Done.** Full verification pass: 36 renders clean (normal) + 36 clean with the CDN blocked (vendored Plotly), 46 Python
tests, JS fixture test, link check clean, load 241 ms, no missing alt text, site 7.7 MB. Rename src/motrpac → src/tfp
(57 files; the executed section notebook reports the inline library identical to src/tfp). README, abstract (199
words), compliance doc, CITATION.cff, LICENSE data note, GitHub Pages + tests workflows, Makefile targets, summary
figure (300 dpi). Optional 6b done: `scripts/16_identifiability.py --bridge` measures batch on the reference pools
from the portal count files: for the two gastrocnemius-derived pools on 6 plates at both sites, Σ V_batch / Σ V_tissue
= 0.017 / 0.017 over all genes (0.010 on the panel genes the pool expresses); wired into the Identifiability and Home
pages with the definition and the muscle-pool caveat. Whole-branch review dispatched to a fresh-context reviewer.
**Next.** Review findings, final export, tag `hackathon-submission-v1`, final report.
**Blockers.** None.
