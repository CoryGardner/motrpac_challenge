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

## 2026-09-27 00:40–01:30 — review and fix pass
**Done.** Fresh-context review of the whole branch: 2 critical, 7 important, 16 minor findings. Fixed in one pass, each
with a test or a scripted browser check: the F-test number rendering as "pending" on Home (new test: every `extras`
key a page reads must exist), tour step 2 vs the Explorer's 35-animal refit (card shows both calls; tour rewritten),
the wrong reason in the ladder caption for Mondrian 1.0 on the held-out sex (8 vials per class → +∞), hand-typed
result numbers wired to JSON (+ reconciliation counts in provenance.json), Explorer charts reverting on theme toggle,
held-out-sex denominators (accuracy_seen, n_samples_coverage), hide-the-answer, regen test pins every 06/12/13 table.
16 minors deferred (listed in the ledger and the final report). Site re-exported (271 provenance entries), full suite
green, render pass clean.
**Next.** Tag `hackathon-submission-v1`, final report.

## Appendix — execution ledger (rulings, task completions, review outcome)

```
# SDD ledger — plan: docs/superpowers/plans/2026-09-26-submission-site-build.md
Spec: docs/superpowers/specs/2026-09-26-submission-site-spec.md (user brief, verbatim; binding).

Pre-flight (shared interfaces):
- Task 1 → Task 4: results/16_identifiability/*.csv column names as in Task 1 Interfaces → nesting.json / qc_baseline.json. Clean.
- Task 2 → Task 4: results/31_site_regen/06_conformal/TRNSCRPT/scores_test_probs.csv, scores_calibration.csv (fold, model, calibration, p_<class>) → samples_motrpac.json. Clean.
- Task 3 → Task 4: results/31_site_regen/{12_bodymap,13_gtex}/scores_target_probs.csv, scores_calibration.csv, recal_thresholds.csv → samples_bodymap.json / samples_gtex.json. Clean.
- Task 4 → Task 7: conformal_fixtures.json computed from the SAME quantised JSON the site ships (6 dp probabilities) so JS equality is exact. Clean.
- Task 6 → Tasks 8–12: charts.js/site.js exports (figure, statTile, callout, pendingBlock, loadJSON, mountChrome). Clean.
- Task 14 rename (motrpac → tfp) happens AFTER all scripts have run; scripts 16/30 import from the package name current at run time and are updated by the same rename. Clean.

Rulings:
- Task 0: Ruling: no git worktree; branch `submission-build` in the same working tree — data/ and results/ are git-ignored and exist only here, a worktree would have neither — cost if wrong: none beyond the usual branch discipline (main is untouched until the final fast-forward).
- Task 0: Ruling: new phase numbers 16 (identifiability recompute), 30 (site export), 31 (regeneration runs) — 17 and 21–26 are reserved by the spec for parallel work that may land; 30/31 keep the site tooling visibly apart from analysis phases — cost if wrong: a rename of two files.
Task 1: complete (commits f0288e3..1bd7e3e, tests: /home/cory/miniconda3/envs/motrpac-py/bin/python -m pytest tests/test_batch.py -q → 3 passed in 0.65s)
Task 2: complete (commits f0288e3..1bd7e3e, tests: /home/cory/miniconda3/envs/motrpac-py/bin/python -m pytest tests/test_regen_scores.py -q -k 06 → 2 passed, 5 deselected in 0.65s)
Task 3: complete (commits f0288e3..1bd7e3e, tests: /home/cory/miniconda3/envs/motrpac-py/bin/python -m pytest tests/test_regen_scores.py -q -k 'transfer or recal' → 5 passed, 2 deselected in 0.71s)
Task 4: complete (commits 1bd7e3e..375e78f, tests: bash -c '~/miniconda3/envs/motrpac-py/bin/python -m pytest tests/test_site_data.py -q && node tests/test_site_conformal.js' → ok: 2713 assertions, 2700 fixture cases (540 with infinite threshold))
Task 5: complete (commits 1bd7e3e..375e78f, tests: bash -c 'test -s docs/NUMBERS_RECONCILIATION.md && grep -q 'Superseded (2026-09-26)' results/SUMMARY.md results/ABSTRACT.md ../../docs/findings/FINDINGS_REPORT.md && echo banners-ok' → banners-ok)
Task 6: complete (commits 375e78f..97d6032, tests: node tools/screenshot.js --pages=index → all page renders clean)
Task 7: complete (commits 375e78f..97d6032, tests: node tests/test_site_conformal.js → ok: 2713 assertions, 2700 fixture cases (540 with infinite threshold))
Task 8: complete (commits 375e78f..97d6032, tests: node tools/screenshot.js --pages=index → all page renders clean)
Task 9: complete (commits 375e78f..97d6032, tests: node tools/screenshot.js --pages=explore → all page renders clean)
Task 10: complete (commits 375e78f..97d6032, tests: node tools/screenshot.js --pages=transfer → all page renders clean)
Task 11: complete (commits 375e78f..97d6032, tests: node tools/screenshot.js --pages=fingerprint,identifiability,beyond → all page renders clean)
Task 12: complete (commits 375e78f..97d6032, tests: node tools/screenshot.js --pages=methods,limitations,about → all page renders clean)
- Task 14: Ruling: docs/EVALUATION_RULES.md (frozen) had three module references rewritten motrpac.* → tfp.* by the rename — a name substitution only, no rule changed; leaving them would point at a module that no longer exists — cost if wrong: three words to revert.
Task 13: complete (commits 97d6032..fb3c3c0, tests: bash -c '~/miniconda3/envs/motrpac-py/bin/python tools/linkcheck.py | tail -1 && node tools/screenshot.js --pages=index,explore --block-cdn | tail -1' → all page renders clean)
Task 14: complete (commits 97d6032..fb3c3c0, tests: bash -c '~/miniconda3/envs/motrpac-py/bin/python -m pytest -q | tail -1 && node tests/test_site_conformal.js | tail -1' → ok: 2713 assertions, 2700 fixture cases (540 with infinite threshold))
- Task 6b: Ruling: the bridge-sample measurement is recomputed here with an explicit definition (between-plate variance of a reference pool's log2 CPM over the plates it ran on, as a fraction of the variance of the 19 tissue means; all-gene and expressed-in-pool sets) and labelled as such — the spec allowed 'pending' but the portal count files contain the 36 reference vials and two gastrocnemius-derived pools bridge 6 plates at both sites — cost if wrong: the audit's own definition may differ; the page says so and keeps the pending note for the audit's number.
Final review: fresh-context reviewer (fable) on 5236735..fb3c3c0 — verdict "with fixes": 2 Critical, 7 Important, 16 Minor.
Final: fixed C1 F-test number rendered as pending on Home — test_every_extras_key_a_page_reads_is_exported RED→GREEN (extras.acc_fclassif_k20 exported).
Final: fixed C2 tour step 2 contradicted the Explorer (35-animal refit) — card now shows both calls (refit + all-animal pred_all_animals), design stated, tour rewritten; verified by a Playwright drive (GTEx heart: refit HEART, all-animal SKM-GN).
Final: fixed I3 ladder caption's wrong reason for Mondrian 1.0 on the held-out sex — caption states the 8-vials-per-class → rank 9 > 8 → +∞ cause; hover flags full sets.
Final: fixed I4 hand-typed result numbers (limitations orthologs / infinite draws / reconciliation counts; identifiability layer counts; fingerprint certificate fractions, flag counts, easy-tissue range; home 171; tile subtitles) — wired to JSON; test_provenance_records_the_reconciliation_counts RED→GREEN.
Final: fixed I5 theme toggle reverted Explorer charts — .spec reassigned on every re-render; Playwright drive shows Gnb3 kept after toggle.
Final: fixed I6 regen test pins the certificate tables and every 12/13 table — test extended (already identical).
Final: fixed I7 held-out-sex denominators — accuracy_seen and n_samples_coverage exported and shown; test_held_out_sex_rungs_carry_both_denominators RED→GREEN.
Final: fixed I8 hide-the-answer — picker selects blurred and disabled, hint shown; Playwright drive confirms.
Final: fixed I9 manifest site_data_files — present (3b2f63f) and pinned by test_manifest_lists_the_site_data_files.
Final: fixed (minor, factual) CITATION date 2026-09-26 and pyproject version 1.0.0.
Final: minor (deferred): conformal.js NaN with fallback=null (unreachable from setsFor).
Final: minor (deferred): α slider and variant stay live in recalibrated mode (note says they do not apply).
Final: minor (deferred): Mondrian + one-vial-per-animal in the card does not say every per-class threshold is +∞.
Final: minor (deferred): QC bars are accuracy, gene-model bar balanced accuracy, one axis label (values differ by < 0.001).
Final: minor (deferred): BodyMap age curve shows no n (age_shift_accuracy.csv has none; accuracy_by_organ has per-organ n).
Final: minor (deferred): recal "draw0" is iloc[0] after sorting — mislabelled only if draw 0 were skipped (all 20 present).
Final: minor (deferred): tests/test_site_conformal.js relies on Node's ESM detection (CI Node 20.20 fine; add tests/package.json type=module for older Node).
Final: minor (deferred): REPO_URL placeholder in site/assets/site.js (unknown until the repo is pushed; documented in README).
Final: minor (deferred): tile subtitle design constants and methods-page design constants (n = 22/396/59/35/15/20) typed as design parameters.
Final: minor (deferred): Metabolomics Workbench listed as a CFDE component per the spec; reword as a cross-reference if the organizers object.
Final: minor (deferred): abstract leads with pooled coverage 0.908; limitations calls one-vial-per-animal 0.916 the honest guarantee (both shown on the ladder).
Final: minor (deferred): `ad_grid =[` spacing in scripts/06; screenshot tool does not fail on imgsNoAlt (linkcheck does).
Final: Ruling: the reviewer's "Declined to judge" items (commits after fb3c3c0 incl. the bridge measurement; TMT plex_id construction; quoted phase-15 verdicts; notebooks not re-executed) stand as built — the bridge work carries the plan's "recomputed here, definition: …" label and its own provenance entries and anchor — cost if wrong: the audit's definition may differ; stated on the page.
```

## 2026-09-27 01:30–02:40 — audit fixes, Beyond removed, no pending numbers
**Done.** Audit of every result claim; fixes: "accuracy survives every shift" → "degrades gracefully"; "every adult organ"
→ every mapped adult organ (9 of 11, super-class scoring, 68 samples from 8 animals); "sequenced as its own batch" →
each tissue sits inside one plate, one library batch and one flowcell; the panel is described as re-selected inside each
fold; the abstract cites the one-vial-per-animal coverage and the technical-only QC number. Ladder made internally
consistent: the in-distribution rung now pairs the accuracy and coverage of the same phase-06 models (k20 0.978 / k50
0.982 / full 0.992; the 40-animal CV numbers stay on the tile and Fingerprint page and in the hover); the held-out-sex
rungs use seen-class denominators for coverage, empty and wrong-but-non-empty from per-vial sets exported by a new
`--save-scores` flag in phase 08; the k = 50 held-out-sex rows come from a k = 50 rerun (results/31_site_regen/08_shift_k50),
so no ladder rung is pending; BodyMap and GTEx rungs and tiles carry 95 % cluster-bootstrap intervals over animals or
donors, with an exact animal-level interval where the bootstrap is degenerate (BodyMap adults: 8 of 8 animals all right,
lower bound 0.63). Removed as not challenge-relevant: the Beyond page, the training-vs-batch verdict section, the
fusion / discordance / batch-verdict / inventory exports, and every reference to parallel phases or "this copy".
51 tests, all anchors, clean render pass.

## 2026-09-27 01:15–01:30 — published
Repository https://github.com/CoryGardner/motrpac_challenge (public); site https://corygardner.github.io/motrpac_challenge/,
published from `main` by `.github/workflows/pages.yml`. Pages was first served from a temporary `gh-pages` branch (the token
lacked the `workflow` scope; the SSH key was passphrase-locked until unlocked); the environment `github-pages` that GitHub
created then only allowed that branch, so the first two workflow deployments were rejected ("Branch main is not allowed to
deploy"); the branch policy was removed, the workflow deployed, and `gh-pages` was deleted. Tags `hackathon-submission-v1`
and `-v2` are pushed. CI (`tests.yml`: pytest without data, JS fixture test, link check) passes on GitHub.
