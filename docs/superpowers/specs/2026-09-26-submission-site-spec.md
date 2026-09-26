# Submission build spec — site, explorer and repo (Molecular Tissue Fingerprints)

This is the brief the build was made from, saved verbatim on 2026-09-26 so the plan and the code can be
read against it. Repo root for every relative path below is `code/pipeline/` (this repository).

---

# Build the submission: site, explorer, and repo — Molecular Tissue Fingerprints

You are building the final hackathon submission from a cold start, working independently for as long
as it takes. There are no check-ins. Log progress to `docs/BUILD_LOG.md` every ~30 minutes of work
(what's done, what's next, blockers). If you are blocked on one item for more than 30 minutes, record
it, move on, and return at the end if time remains. Stop only when the final report is written.

## The challenge, verbatim

> **Molecular Tissue Fingerprints — Can molecular signatures identify a tissue reliably?**
> Use one or more omics layers to predict tissue identity and determine whether a compact,
> interpretable signature is sufficient. Data: individual-sample rat endurance-training data across
> a selected set of tissues, with GTEx or another tissue resource for external context.

Rules of the competition: at least one MoTrPAC component, at least one CFDE component, and the
deliverable is an open-source repo. The audience is judges with mixed backgrounds (biologists,
bioinformaticians, computer scientists) who will spend five to ten minutes on it.

## Our answer, which the site must deliver in the first thirty seconds

**Yes — a 20-gene panel identifies 19 rat tissues at 0.976 balanced accuracy and names every adult
organ correctly in another laboratory's rats. But "reliably" has three parts, and only the first
survives on its own.**

1. **Accuracy.** A compact signature is sufficient: k = 20 → 0.976, k = 50 → 0.993, full model 0.995.
   A ten-gene stable core of textbook markers. The selector, not the classifier, was the hard part
   (round-robin 0.976 vs F-test 0.399 at the same k).
2. **The guarantee.** A conformal prediction set promises the true tissue is in the set 90% of the
   time. Under shift the classifier keeps naming tissues correctly while the guarantee collapses —
   coverage 0.96 in-distribution → 0.618 across labs → 0.364 across species — and it collapses by
   *abstaining* (empty sets), not by confident mistakes. Three target animals repair it within species
   (0.943 coverage at ~1 tissue per set); across species they do not (11.7 of 19 tissues per set).
3. **Identifiability.** Within a single multi-tissue study the tissue axis cannot be separated from
   processing: every tissue was extracted, library-prepped and sequenced as its own batch, across
   every omic layer (1 estimable tissue pair out of 171 in RNA-seq, and it is the sex contrast). QC
   covariates alone classify tissue at 0.975. So within-study accuracy is not evidence that the
   signature is biology. The evidence is external: the panel transfers to a lab where none of these
   batches exist. And where batch could actually be measured (one bridging standard), it was ~1.7%
   of the variance that separates tissues.

Then two things beyond identity: the failures are explained gene by gene (a juvenile testis has no
spermatids yet, so Pgk2 is absent), and the fingerprint's real value is as an *instrument* — the
axis you subtract to see everything else.

## 0. Orient — read before writing anything

- `CLAUDE.md`, `docs/EVALUATION_RULES.md` (FROZEN), `docs/DATA_GUIDE.md`, `docs/DATA_INVENTORY.md`
- `src/motrpac/` — every module. `splits.py` is FROZEN (add functions only).
- `notebooks/02_transfer.ipynb` and `notebooks/expected_values.csv` — the post-fix verified numbers.
- `results/` — inventory every phase directory that exists. Phases 04–14 are the core; 15–20 may hold
  training-axis work; 21–26 may hold the identifiability audit and the exercise decomposition. Some
  may be absent in this copy; that is expected.
- `results/QUANTILE_FIX_CHANGES.md` and `backup/pipeline_history/results_pre_quantile_fix_2026-09-25/`
  — the conformal quantile was fixed on 2026-09-25. `results/` is post-fix and is the truth. Any
  number in `FINDINGS_REPORT.md`, `SUMMARY.md` or `ABSTRACT.md` that disagrees with `results/` is
  stale; `results/` wins, always.

Hard rules:
- **No number on the site is typed by hand.** Every number is exported from `results/` (or
  recomputed from data with the existing pipeline) by `scripts/NN_export_site_data.py` into
  `site/data/*.json`, with a provenance manifest. A number you cannot source renders as "pending"
  with the reason. Never fabricate, never estimate, never round a stale value into a fresh one.
- The animal is the unit; every accuracy shown carries its per-fold spread and the n behind it.
- Data stays out of git: derived JSON for the site is fine (panel-gene values and per-sample class
  scores from public data); full matrices are not.
- Do not modify existing `results/` directories. New scripts take the next free phase numbers.
- Do not present anything as final that depends on work still running elsewhere — see §3.6.

## 1. Phase A — Establish the truth (do this first; nothing else is valid until it is done)

**A1. Inventory.** Write `site/data/manifest.json` listing every results phase present, its files,
the git hash of the repo, and the generation timestamp. Every later export cites this manifest.

**A2. Numbers reconciliation.** Write `docs/NUMBERS_RECONCILIATION.md`: every headline number the
site will display, its source file and row, its post-fix value, and (where it moved) its pre-fix
value. Recompute the in-distribution conformal coverage — the 0.962 quoted in older documents is
PRE-fix (the fix lowers the threshold by one rank at 22 calibration points) — and use the post-fix
value everywhere. Where the three stale documents disagree with `results/`, list the disagreement;
do not edit those documents, but add a banner at the top of each stating that `results/` and this
reconciliation supersede them.

**A3. Per-sample exports for interactivity.** The explorer needs, for every sample, its class
scores and the calibration scores that produce prediction sets. Check what phases 04/06/12/13 already
persisted. Where per-sample outputs are missing, regenerate them by adding a `--save-scores` flag to
the existing scripts (do not duplicate model logic), then **verify that the regenerated aggregates
reproduce `expected_values.csv` to the recorded tolerance** before using them. Export:
- MoTrPAC TRNSCRPT, 5-fold animal-grouped CV: out-of-fold class probabilities (19 classes) for every
  study vial, for k = 20, k = 50 and the full model; plus, per fold, the LAC calibration scores under
  the 18 fit / 22 calibration / 10 test design, both pooled-vial and one-vial-per-animal, per class
  (needed for marginal, Mondrian and floored sets). Include `pid`, tissue, sex, group.
- Rat BodyMap (phase 12 design): class probabilities for all 316 samples, organ, age, assumed animal
  id, and the MoTrPAC calibration scores per class from the 15 calibration animals.
- GTEx (phase 13 design): class probabilities for all 2,485 samples, tissue, donor, and the rat
  calibration scores.
- Recalibration thresholds: if per-draw thresholds for n = 3 and n = 5 target individuals are
  recoverable with the pipeline's seed, export the first draw's thresholds per (target, model, n);
  otherwise export the mean coverage / set size from `recalibration.csv` and mark the explorer's
  recalibration control as summary-only.
- Panel-gene expression per sample, quantised to 2 decimals: log2 CPM in MoTrPAC (899 vials) and
  BodyMap (316), log2 TPM in GTEx (2,485, through 1:1 orthologs), for the union of: the k = 20 and
  k = 50 panels, the 10-gene stable core, all 51 ever-selected genes, the developmental markers
  (Pgk2, Prm1, Tnp1, Tnp2, Acrv1, Mybph, Klf1, Alas2, Hbq1b, Hbb, Fcrl5), and every gene in
  `panel_gene_check.csv`. Cap ~100 genes. Include each gene's annotation: marker tissue, bootstrap
  frequency, effect size, `pct_mrna` correlation, training-regulated flag, fails-in-GTEx flag.
- Aggregate tables as JSON: panel curve (both selectors), stable core, in-distribution confusion,
  BodyMap age accuracy, BodyMap conformal and recalibration, GTEx per-tissue accuracy and confusion
  (k20/k50/full), GTEx conformal and recalibration, representation results incl. the CPM test,
  certification results and the animals-needed table, QC-only baselines (technical and composition
  SEPARATE), batch nesting per variable, fusion and batch-verdict tables, discordance summary.

Keep every JSON under 3 MB and the whole `site/` under 15 MB. Quantise floats; never ship a full
gene × sample matrix.

**Sanity anchors** — these are post-fix self-check values. `results/` is truth; if an export disagrees
with one of these, investigate before proceeding, do not paper over it:
MoTrPAC full model 0.995 ± 0.005, k20 0.976 ± 0.008, F-test k20 0.399 · BodyMap adults k20 1.000,
all ages 0.9104, by age 0.765 / 0.970 / 1.000 / 0.908, native 0.9925 · BodyMap adult marginal
coverage k20 0.618 (empty 0.382), floored 0.691, recal n3 k20 0.943 coverage at 0.998 set size ·
GTEx k20 0.654, k50 0.781, full 0.855, native 0.979, heart k20 0.033 (→ skeletal muscle 0.907),
ovary 0.000 (→ HEART 0.953), marginal coverage k20 0.364 (empty 0.616), full 0.062 (empty 0.938),
recal n3 k20 0.954 coverage at 11.70 set size with 45% infinite draws · QC-only technical
0.873 ± 0.027, composition 0.949 ± 0.011, all 0.975 ± 0.020; 17 plates, 17 library batches,
4 flowcells · shared genes 21,040 of 21,193; orthologs 14,609 (14,569 in GTEx) · discordance
1,948 / 1,899 / 4, AUROC 0.780 with flag / 0.742 without · fusion 0 of 7 beat single omic, 7 of 7
beat null.

## 2. Phase B — Site architecture and design system

**Architecture.** A static site in `site/`, viewable with `python -m http.server -d site 8000` and
deployable to GitHub Pages unchanged: relative paths only, `.nojekyll`, no server routing, no build
step required to view. Plain HTML + CSS + vanilla JS (ES modules are fine). Charts with Plotly.js,
pinned version, loaded from CDN with a vendored fallback in `site/vendor/` so the demo works offline.
Multi-page with one shared nav; every page has the same header, footer (repo link, results git hash,
generation date) and a "Source" caption under every figure naming the `results/` file it came from.

**Design system** — define once in `site/assets/theme.css` as CSS custom properties, and once in
`site/assets/charts.js` as a Plotly template every chart uses. Do not deviate per page.

Typeface: `system-ui, -apple-system, "Segoe UI", sans-serif` everywhere, including hero numbers.
Prose max width 72ch, line-height 1.55. Generous whitespace; no boxes-within-boxes.

Surfaces and ink (light / dark — support both via `prefers-color-scheme` and a toggle):
page `#f9f9f7` / `#0d0d0d`; chart surface `#fcfcfb` / `#1a1a19`; primary ink `#0b0b0b` / `#ffffff`;
secondary ink `#52514e` / `#c3c2b7`; muted `#898781`; gridline `#e1e0d9` / `#2c2c2a`;
axis `#c3c2b7` / `#383835`.

Categorical palette, fixed order, never cycled (light / dark):
1 blue `#2a78d6`/`#3987e5` · 2 orange `#eb6834`/`#d95926` · 3 aqua `#1baf7a`/`#199e70` ·
4 yellow `#eda100`/`#c98500` · 5 magenta `#e87ba4`/`#d55181` · 6 green `#008300`/`#008300` ·
7 violet `#4a3aa7`/`#9085e9` · 8 red `#e34948`/`#e66767`.
Sequential (one hue, blue): `#cde2fb #9ec5f4 #6da7ec #3987e5 #256abf #184f95 #0d366b`.
Diverging: blue ↔ red with neutral midpoint `#f0efec` / `#383835`.
Status (icon + label always, never colour alone): good `#0ca30c`, warning `#fab219`, serious
`#ec835a`, critical `#d03b3b`.

**Nineteen tissues cannot be told apart by hue.** Tissue identity is carried by labels and tooltips,
never by colour alone. Where colour is used on tissue-level marks, it encodes ORGAN SYSTEM, eight
groups in the palette order above: 1 brain (CORTEX, HIPPOC, HYPOTH) · 2 muscle (SKM-GN, SKM-VL,
HEART) · 3 adipose (WAT-SC, BAT) · 4 gut (COLON, SMLINT) · 5 gonad (OVARY, TESTES) · 6 circulation
& immune (BLOOD, SPLEEN, VENACV) · 7 visceral (LIVER, KIDNEY, LUNG) · 8 endocrine (ADRNL). Per-tissue
bar charts use one hue with the axis carrying identity. Confusion matrices are sequential heatmaps.

Chart rules, enforced by the template: one y-axis per chart, never dual; no rainbows; thin marks,
2px lines, ≥ 8px markers; hairline gridlines only; legend present for ≥ 2 series and direct labels
when ≤ 4; hover tooltips on every mark; text in ink tokens, never series colour; a data table view
toggle on every chart. Every figure is a block with: **title = the claim**, subtitle = what is
plotted, the chart, and a caption with "Source:" and "What it does not show:" — the house style.

Components to build once: nav, stat tile (hero number, label, source), figure block, callout
(note / caveat / pending), sample-card, table with sticky header, theme toggle.

## 3. Phase C — Pages, in this priority order

Build in this order so the essential pages exist even if time runs short: Home → Explore → Transfer
→ Fingerprint → Identifiability → Beyond → Methods → Limitations → About.

**3.1 `index.html` — Home.** H1 is the challenge question: "Can a molecular signature identify a
tissue reliably?" Subtitle is our one-sentence answer. Then four stat tiles: 0.976 (20-gene balanced
accuracy, 19 tissues), 1.000 (adult BodyMap organs, another lab), 0.618 (coverage there, the
guarantee that did not travel), 1 of 171 (estimable tissue pairs within study). Then the **transfer
ladder** — the site's signature figure: paired accuracy and coverage bars per shift (in-distribution,
held-out sex, different lab, different species k20, different species full), dashed line at 0.90,
interactive controls for model and conformal variant. Then three short sections — Accuracy, The
guarantee, Identifiability — each one paragraph, one chart, one link. End with the "Five-minute
tour": an ordered list of five links with one line each on what to click.

**3.2 `explore.html` — The Explorer.** This is the interactive centrepiece; spend real time on it.

(a) **The tissue card.** A sample picker: source (MoTrPAC held-out vial / BodyMap by age / GTEx by
tissue) plus a "random sample" button, and a "hide the answer" toggle for demos. Controls: α slider
0.05–0.30, model (k20 / k50 / full), conformal variant (marginal / Mondrian / floored), calibration
(source-calibrated / recalibrated on target animals when thresholds are available). Output: the
prediction set as chips; a badge — **Confident** (one tissue), **Ambiguous** (several), **Abstains**
(empty) — with a one-sentence plain-English explanation; the top-5 class probabilities as a bar; and
the panel-gene strip: this sample's within-dataset z for each panel gene against the tissue
reference profiles. Prediction sets are computed **client-side** from the exported scores, so the
α slider is exact: port `conformal_quantile` to JS faithfully, including the rank rule
⌈(n+1)(1−α)⌉ and the +∞ case when that rank exceeds n, and the floored-Mondrian rule (per-class
threshold never below the marginal one). Make sure BodyMap thymus and uterus samples are in the
picker — the fingerprint abstaining on an organ it has never seen is the best moment on the site.

(b) **The gene explorer.** Gene picker over the ~100 exported genes. Three linked panels: MoTrPAC
log2 CPM by tissue (strip + median, sex as marker shape), BodyMap by organ × age, GTEx by tissue.
An annotation card: marker tissue, bootstrap frequency, effect size, risk flags, fails-in-GTEx.
Preload Pgk2 (collapses in the 2-week testis) and Gnb3 (the heart marker that points at ovary in
human) as the two featured examples.

(c) **The panel builder.** k slider over the grid {1,2,3,5,8,10,15,20,30,50,100}: genes selected at
that k, CV balanced accuracy ± SD, and which tissues have a marker yet — a 19-tissue coverage strip
that fills as k grows. Shows why the jump from 10 to 15 happens.

(d) **The animals-needed calculator.** α and δ sliders → n ≥ ln δ / ln(1 − α) with the four
reference points (59 / 29 / 22 / 11). One sentence: the guarantee is paid for in animals, not genes.

**3.3 `transfer.html` — Does it transfer?** The ladder with full controls; the empty-set stack
(covered / wrong-but-non-empty / empty per shift); the recalibration-cost scatter (coverage vs mean
set size, one point per target × model × n, target region shaded — BodyMap points inside, GTEx
three-donor points at ~11.7); the BodyMap age curve beside small multiples of the developmental
markers (log2 CPM vs age, panel marker distinguished from confirmatory genes); GTEx per-tissue
accuracy and the row-normalised confusion heatmap with a k20/k50/full toggle; the representation
result (z-score vs rank vs pairs, the CPM test, the ovary exception); the panel-survival table.
State plainly, in a callout, that "three donors restore coverage" is FALSE across species and why
(45% of draws have too few samples for a finite threshold).

**3.4 `fingerprint.html` — The signature.** Panel curve, round robin vs F-test, log x, fold SD
bands; the ten-gene stable core as an annotated table; in-distribution per-tissue accuracy and
confusion; the hard tissues (BAT, the two muscles, VENACV) and why (weak single markers, 1.26–2.44
log2 CPM above the next tissue; LIVER's instability is redundancy, six near-equivalent candidates).

**3.5 `identifiability.html` — What is the fingerprint identifying?** Frame it with the organizers'
own question: a set of analytes, or a set of analytes given it came from a particular centre? The
nesting table per layer (levels per tissue, tissues per level, Cramér's V) as a heatmap; the
estimable-pairs verdict per layer; QC-only classification with technical and composition shown
SEPARATELY and the caveat that composition features (`pct_chrM`, `pct_globin`, chrX/chrY) are read
biologically by MoTrPAC itself; the bridge-sample measurement if available; the two failure modes
(signal present but confounded — RNA, metabolomics; signal removed by quantification — TMT
proteomics, where values are ratios to a per-tissue reference); and the resolution: external
transfer. Do NOT present ATAC as an exception — its extraction dates cross tissues but each flowcell
holds one tissue, which closes it. If `results/21_identifiability/` is present with the
estimable-pairs framing, use it. If absent, recompute the RNA and proteomics nesting from metadata
(minutes) and show the other layers as "from the parallel audit, pending merge" with the facts
stated. If a batch-verdict-v2 (phase 15) exists, use its verdicts; otherwise show the 3-of-7
unresolvable table with a note that the QC set behind it contains biologically-read covariates and
a revision is pending.

**3.6 `beyond.html` — The fingerprint as an instrument.** This page depends on work running in
parallel. Build it as a function of what exists:
- If a training-transfer phase (likely 17) exists: before displaying its leave-one-tissue-out result,
  read its split function and confirm it blocks on BOTH `pid` and tissue and has a passing unit
  test; show the per-tissue AUROC chart against its null, badged "preliminary" unless the
  permutation null is present. If blocking cannot be confirmed, do not show the number.
- If decomposition phases (22–26) exist: forest-plot picker for the top tissue-agnostic genes (lead
  with the heat-shock genes), the I² distribution, the cross-tissue response correlation map, the
  leave-one-tissue-out AUROC with the QC-covariate control arm, the RNA-vs-protein shared-effect
  scatter.
- Whatever is absent renders as a clean "in progress" section quoting the pre-registration file if
  it exists, so the page is honest rather than empty. `make site` must regenerate this page when
  the results land.

**3.7 `methods.html`.** An SVG pipeline diagram (median imputation → variance prefilter → z-score →
round-robin selector → classifier, all inside the fold); animal-grouped splits; the selector and why
the F-test fails on multiclass; conformal sets with the quantile rule written out; the certificate;
the transfer protocol; the frozen evaluation rules (link to the file); the two-notebook self-check
scheme; software versions.

**3.8 `limitations.html`.** Honest and specific: pre-split gene filter (known deviation, unmeasured);
assumed BodyMap animal ids; super-class scoring in transfer; single-sex tissues and sex confounded
with cohort; only 8w matched on collection date; one release (rn6 c1.0; rn7 changes the
training-regulated sets); small n; the quantile fix history and which documents are stale; the
untuned QC baseline.

**3.9 `about.html`.** Team (placeholder block for names — leave clearly marked for editing), the
challenge statement, data sources with links and licences (MoTrPAC c1.0 / R package, rat BodyMap
GSE53960, GTEx v8 open-access expression), competition compliance (MoTrPAC is the dataset; CFDE
components: GTEx, and Metabolomics Workbench — a CFDE Data Coordinating Center holding MoTrPAC
metabolomics as PR001020), how to reproduce, how to cite, licence.

## 4. Phase D — Verification (not optional; a broken page is worse than no page)

- **JS conformal port test.** Export Python-computed prediction sets for a fixed sample subset at
  α ∈ {0.05, 0.10, 0.20}, all three variants, and assert the JS implementation reproduces them
  exactly, including infinite thresholds. Run under Node in `tests/test_site_conformal.js`.
- **Provenance test.** `tests/test_site_data.py`: every value in `site/data/*.json` matches its
  cited `results/` file within tolerance; every stat tile on the home page has a manifest entry.
- **Render pass.** Screenshot every page at 1440×900 and 390×844 with headless Chromium (Playwright
  if installed, else `chromium --headless --screenshot`) into `site/_screenshots/` (git-ignored), in
  both themes. **Open and look at every screenshot.** Fix label collisions, overflow, unreadable
  ticks, empty charts. Zero console errors on every page.
- Link check (internal links resolve; external are well-formed), total size, largest JSON, load
  time under 2 s locally, keyboard focus on all controls, alt text on all images.

## 5. Phase E — Repo polish and deploy-readiness

- Rename `src/motrpac/` (it collides with MoTrPAC's real R package names) — e.g. `src/tfp/`. Update
  imports, the notebook builder, the Makefile; run the test suite.
- Root `README.md` rewritten as the submission README: title, the one-paragraph answer, a screenshot
  of the home page, how to view the site (one command), a key-results table generated from
  `results/`, repo map, how to reproduce (`make` targets), data access (what is and is not in the
  repo, and where to get the 128 GB), licence, citation, compliance.
- `LICENSE` (MIT for code; note data licences), `CITATION.cff`, `docs/COMPETITION_COMPLIANCE.md`.
- `.github/workflows/pages.yml` publishing `site/` to GitHub Pages; `.github/workflows/tests.yml`
  running the tests that do not need the raw data.
- `make site` = export data → run site tests → (optional) screenshots. Document it.
- `docs/ABSTRACT_SUBMISSION.md`: a 200-word abstract whose numbers are filled from `results/`.
- `figures/summary_figure.png` at 300 dpi: a three-panel composite (transfer ladder · stable core ·
  identifiability heatmap) for slides, using the same design system.
- Commit in logical chunks with clear messages; tag `hackathon-submission-v1` at the end.

## 6. Time budget (guidance, not a contract)

A ≤ 3 h · B ≤ 1.5 h · C ≤ 6 h with Home, Explore and Transfer first · D ≤ 2 h · E ≤ 2 h. If you must
cut, cut from the bottom of §3 and never from §1 or §4.

## 7. Final report

1. How to view the site locally, and what deploying to GitHub Pages requires.
2. Page-by-page status: done / partial / pending, with reasons.
3. Verification results: tests passed, screenshots reviewed, sizes, load times.
4. The numbers reconciliation summary — every headline number that changed from the stale documents.
5. Anything that could not be done validly, and why.
6. The five-minute judge tour, as it appears on the home page.
7. The abstract.

Do not tune anything to make a number look better. If a regenerated aggregate disagrees with the
self-check, that is a finding to report, not a discrepancy to hide.
