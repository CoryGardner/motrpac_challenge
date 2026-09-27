   PREFLIGHT (≤ 10 min, before Phase 0; do not wait for me afterwards):
   1. `date -u` → write `START <timestamp>` as the first line of docs/MULTIOMIC_LOG.md.
   2. Confirm branch = multiomic-overnight; `ls data/raw data/external` resolves to the existing MoTrPAC tree; `df -h .` shows ≥ 300 GB free where downloads will go.
   3. Reachability, one line each (HTTP code from `curl -sSIL --max-time 20`): PRIDE (FTP and HTTPS), Cell/Elsevier supplementary host, EMBO/MSB, UniProt REST, Ensembl REST, Metabolomics Workbench REST, MetaboLights.
   4. Whether MoTrPAC RII (reporter-ion intensity) proteomics files exist locally; if not, whether their download URL answers.
   Write all of it under "Preflight" in docs/MULTIOMIC_LOG.md, then start Phase 0. An unreachable host is noted now and handled by the drop policy later — it is not a reason to stop.

# Overnight autonomous run: can the proteome and metabolome carry a tissue fingerprint that transfers?

You are working unattended for 10+ hours from a cold start in this repo. Nobody will answer
questions. You may use as much time, compute and disk as you need (a few hundred GB is fine) and
download any public data you can reach without a login. Iterate, decide, and keep going. The
morning deliverable is a report a person can read in ten minutes.

## Ground rules

- `git checkout -b multiomic-overnight`. Never commit to main; never modify `site/`, `results_frozen/`,
  `README.md`, `docs/EVALUATION_RULES.md` or `src/tfp/splits.py` (you may ADD functions to a new
  module). All new outputs go to `results_multiomic/<phase>/` and `data/external_multiomic/`.
  The submission must be exactly as you found it if this branch is deleted.
- Read first: `CLAUDE.md`, `docs/EVALUATION_RULES.md`, `docs/DATA_GUIDE.md`, `docs/DATA_INVENTORY.md`
  (§4 proteomics and metabolomics, esp. the quant-id RII files), `src/tfp/` (io, models, conformal,
  transfer, batch), `scripts/04_prot_diagnostic.py`, `scripts/16_identifiability.py`,
  `results_frozen/16_identifiability/`, `results_frozen/07_fusion/`.
- Same rigour as the submission: the animal (or donor) is the unit of every split; every mean with
  its spread and n; permutation nulls where a null is needed; pre-register predictions before
  running each phase; never type a number into markdown — README files read their CSVs.
- **Reporting contract.** Maintain `docs/MULTIOMIC_REPORT.md` with, at the top, a dated
  "Significant findings so far" list (one line each, with the number, its n and its source file),
  then one section per phase (question · data · design · result · what it does not show). Rewrite
  the top list after EVERY phase. Keep a timestamped `docs/MULTIOMIC_LOG.md` of decisions,
  downloads (URL, size, sha256, date), failures and time spent. If the run dies, these two files
  are the output.
- Time-box every download and every phase (boxes below). When blocked > 30 min on one thing, log
  it, move on, return at the end if time remains. A dataset that cannot be matched (IDs, units,
  tissues) is a recorded negative result, not a reason to force it.
- Never attempt raw mass-spectrometry reprocessing (no .raw/.mzML searches). Processed tables only.
- Do not present a within-study number as evidence of biology; the audit already showed why.

## The question, and what "solved" means

The track allows one or more omic layers; we used one. The proteome fails as distributed (TMT
ratios to a per-tissue reference erase the cross-tissue axis) and the metabolome has no external
check. "Solved" means, in order of strength:
1. a protein or metabolite fingerprint exists in MoTrPAC on a scale where cross-tissue comparison is
   defined (the RII rescue);
2. it transfers to an independently processed dataset, with accuracy AND conformal coverage
   reported, and recalibration on a few target individuals measured;
3. the RNA fingerprint and the protein/metabolite fingerprint agree — same tissues easy, same pairs
   confused, RNA-selected markers behave the same at the protein level;
4. at least one external dataset has a design in which tissue is NOT nested in batch, so the
   fingerprint is identifiable there — the positive counterexample our audit lacks;
5. a two-layer fingerprint is more robust under shift than either layer alone (fusion judged by
   transfer, not by within-study accuracy, which is a ceiling task).

## Phase 0 — setup and pre-registration (≤ 30 min)

Create the branch, dirs, report skeleton and `docs/PREREGISTRATION_MULTIOMIC.md` stating, before any
data is touched: (a) RII-normalised MoTrPAC proteomics recovers tissue as the dominant axis (tissue
R² on PC1 > 0.5, from 0.001 on ratios); (b) a protein panel selected on RII classifies the 7
proteomics tissues at ≥ 0.95 balanced accuracy under animal-grouped CV; (c) proteins of the RNA
panel genes that are quantified show tissue specificity at the protein level with the same marker
tissue in ≥ 70% of cases; (d) an RII-selected protein panel transfers to a human proteome atlas at
accuracy well above chance and at coverage below 0.90 with MoTrPAC calibration (the RNA pattern);
(e) a metabolite panel selected on MoTrPAC HILIC+ transfers to an external rodent tissue metabolome
above chance; (f) at least one external atlas has tissue crossed with its batch variable.

## Phase 1 — the RII rescue (local data only, ≤ 2.5 h; guaranteed to produce a result)

`results_multiomic/01_rii/`. From the portal quant-id proteomics folders (prot-pr; also ph and, in
c2.0, ac if time), read the RII tables for all 7 tissues. Per plex: drop reference channels, filter
proteins by minimum peptides, normalise each channel to its column total (or median), log2, then
combine plexes within a tissue and stack across tissues. Report, on the stacked matrix: tissue R² of
PC1–3 (against 0.001/0.004 on ratios), the missingness-only baseline, the per-tissue-mean-removed
classifier (against 0.193), and a RoundRobinSelector protein panel curve under the existing
animal-grouped folds with the pipeline unchanged. State plainly that plex is nested in tissue here
too (the audit), so this is within-study evidence only — its job is to show the signal exists.
Also compute the cross-tissue RNA–protein correlation per gene on the same animals (this was
impossible on ratios): Spearman of tissue means across the 7 tissues, distribution over genes, and
the value for each RNA panel gene whose protein is quantified. That is the "does the fingerprint
transfer across layers" number Cory asked for.

## Phase 2 — data discovery (≤ 2 h total; log every attempt)

Try, in this order, and stop each attempt at 20 min or 20 GB. Prefer processed supplementary tables
and repository REST APIs; PRIDE/ProteomeXchange FTP and Metabolomics Workbench REST need no login.
Record for each: URL, what it contains, samples × features, tissue list, units, the batch/plex
metadata available, and whether it can be matched to our features.

Proteome atlases (human, cross-species through the existing 1:1 rat–human ortholog table;
UniProt → gene symbol/Ensembl via the UniProt ID-mapping REST if needed):
- Jiang et al. 2020, Cell 183:269, "A Quantitative Proteome Map of the Human Body" — 32 tissues,
  14 GTEx donors, TMT, matched GTEx RNA-seq; PRIDE PXD016999; processed abundances in the
  supplementary tables. The matched RNA makes it the key dataset for Phase 5.
- Wang et al. 2019, Mol Syst Biol 15:e8503, "A deep proteome and transcriptome abundance atlas of
  29 healthy human tissues" — label-free iBAQ with matched RNA-seq; PRIDE PXD010154; EV tables are
  open access; also ProteomicsDB.
- A mouse multi-tissue proteome (e.g. Geiger et al. 2013 MCP, 28 tissues; or any newer atlas you
  find) — rat–mouse matching by gene symbol, with the caveat recorded.
- Search PRIDE and ProteomeXchange for any RAT multi-tissue proteome with processed tables; if one
  exists it is the BodyMap-equivalent and jumps to the front.
Metabolome atlases (match by RefMet name first, then HMDB/KEGG/InChIKey; MoTrPAC's
metadata-metabolites files carry RefMet names):
- Sato et al. 2022, Cell Metab 34:329, "Atlas of exercise metabolism reveals time-dependent
  signatures of metabolic homeostasis" — mouse, several tissues, sedentary and exercised: the best
  fit, because it tests both the tissue axis and exercise invariance. Supplementary tables; check
  Metabolomics Workbench.
- Metabolomics Workbench REST (`/rest/study/...`): search for mouse or rat studies with ≥ 4 tissues
  and named metabolites; rank by tissue overlap with MoTrPAC's 19.
- MetaboLights likewise.
If nothing external is reachable, say so at the top of the report and spend the time on Phases 1,
5-local and 7.

## Phase 3 — protein fingerprint transfer (≤ 2.5 h; needs a Phase 2 proteome)

`results_multiomic/03_prot_transfer/`. Reuse `src/tfp/transfer.py` unchanged: panels selected on
MoTrPAC RII (all animals), features matched to the target, z-scored within dataset, super-classes
where a target tissue maps to several of ours, every target split grouped on donor. Report per
target tissue: accuracy, main wrong call, conformal coverage with MoTrPAC calibration, empty-set
fraction, and recalibration on 3 and 5 donors with set size. Then the ladder for proteins beside
the RNA ladder on the same axis. Also the reverse direction where a target has enough samples:
select on the atlas, score MoTrPAC RII. Note that human atlases are adults of both sexes and
post-mortem; say what that mixes in.

## Phase 4 — metabolite fingerprint transfer (≤ 2.5 h; needs a Phase 2 metabolome)

`results_multiomic/04_metab_transfer/`. Panel on MoTrPAC HILIC+ (19 tissues) and on the deep
platforms (9 tissues), matched to the target by RefMet name; same transfer protocol; report matched
feature counts first and stop the leg if the overlap is < 30 named metabolites. If the target is
Sato 2022: also fit on sedentary mice only and test on exercised mice (the invariance test), and
compare to our RNA result (0.961 / 0.903).

## Phase 5 — fusion judged by transfer, on matched RNA + protein (≤ 2 h)

`results_multiomic/05_fusion_transfer/`. Where an atlas has RNA and protein on the same samples
(Jiang 2020; Wang 2019): score the RNA panel, the protein panel, and a late-fusion of the two (mean
of class probabilities, and a stacked logistic regression fit on MoTrPAC only) on the same human
samples. The question is not "is fusion more accurate in MoTrPAC" (it is a ceiling task; 0 of 7
there) but "is a two-layer fingerprint more robust under a species shift than either layer" —
accuracy AND coverage AND empty-set fraction. Also: within the atlas, is the confusion structure the
same across layers (same pairs confused)? Report the correlation of the two confusion matrices'
off-diagonals. If no matched atlas was obtained, run the local version: RNA-vs-RII confusion
structure on MoTrPAC's 7 shared tissues, animal-grouped.

## Phase 6 — identifiability of the external designs (≤ 1 h; metadata only)

`results_multiomic/06_external_identifiability/`. For every dataset obtained, apply the audit's
framing from `scripts/16_identifiability.py` to its sample metadata: for each batch variable
available (TMT plex, MS run/batch, acquisition date, plate), levels per tissue, tissues per level,
Cramér's V, and estimable tissue pairs. A design in which tissues were spread across plexes or runs
is the positive counterexample: report it prominently, with the numbers, and state what it means
for our audit's conclusion ("nesting is a choice of design, not a property of the assay").

## Phase 7 — synthesis (≤ 1 h, and repeat at the end regardless of how far you got)

- Rewrite "Significant findings so far" at the top of `docs/MULTIOMIC_REPORT.md`: numbered, one
  line each, strongest first, each with n, spread or null, and the CSV it comes from.
- Add a section "What this changes in the submission, if anything": which sentences on the site or
  in the README would change, quoted, with the replacement — but do NOT edit them.
- Add "For the talk": at most three sentences a presenter could say tomorrow, only for findings
  that survived their null and their caveats.
- Add "Not done / not possible", with reasons, and the download log summary.
- Draft `site/_drafts/multiomic.html` content as Markdown only (no site changes): the figures you
  would show, with their source CSVs.
- Commit everything on the branch in logical chunks; do not merge.

## Order and drop policy

0 → 1 → 2 → (3, 4, 5 in whatever order the data arrived, most complete dataset first) → 6 → 7.
If time is short, Phase 1, Phase 6 and Phase 7 are the ones that must exist. A clean negative
("no reachable external proteome could be matched; here is why") is a valid Phase 2 outcome and
must be reported as one, not papered over.

Definition of a significant finding, for the top list: an external transfer accuracy with its
coverage and n; an external design with crossed batches; the RII tissue-R² recovery; the
cross-tissue RNA–protein correlation of the panel genes; a fusion-under-shift comparison. Anything
that is only within-study accuracy is context, not a finding.
