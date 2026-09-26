# EXTERNAL_VALIDATION — the rat BodyMap as an independent test of the tissue fingerprint

Purpose: test the MoTrPAC-trained fingerprint (panel + classifier + conformal calibration) on rat
tissue data produced by a different lab, at different ages, with a different library protocol,
**without** the orthology loss of the human leg (`docs/GTEX_TRANSFER.md`). Same species, so a
failure here is a failure of the fingerprint or of its calibration, not of the mapping.

## Dataset: Rat BodyMap, GEO GSE53960

Yu Y. et al., *A rat RNA-Seq transcriptomic BodyMap across 11 organs and 4 developmental stages*,
Nature Communications 5:3230 (2014). Fischer 344 rats, **11 organs** (adrenal gland, brain, heart,
kidney, liver, lung, skeletal muscle, spleen, thymus, testes, uterus), **4 ages** (2, 6, 21, 104
weeks), **both sexes**, 4 biological replicates per organ × sex × age (320 biological samples;
sex-specific organs only in one sex). rRNA-depleted total-RNA libraries, Illumina HiSeq 2000.
MoTrPAC animals are 6-month-old (≈ 26-week) F344 rats with poly-A mRNA libraries (globin-depleted
for blood), so age, library chemistry and lab are all shifted at once — that is the point.

## Two ways to get it (done 2026-09-17: `bodymapRat` via `R/export_bodymap.R`, hub file 15.8 MB; the GEO archive is also on disk under `data/external/bodymap/`)

| Source | What you get | Size | Needs |
|---|---|---|---|
| Bioconductor ExperimentHub package **`bodymapRat`** (v1.28.0, CC BY 4.0) | `SummarizedExperiment` of STAR gene counts, **32,637 Ensembl genes (ENSRNOG) × 652 RNA-seq runs** (technical replicates / lanes included); `colData` has organ, sex (F/M), stage (weeks), `techRep`, `rnaRIN`, instrument, flowcell/lane, SRA/GEO ids | counts matrix on the order of tens of MB (ExperimentHub cache; exact size shown on first download) | R 4.4 (present), `BiocManager` (present), `BiocManager::install("bodymapRat")`, network |
| GEO supplementary **`GSE53960_RAW.tar`** | 320 per-sample files of **AceView gene symbols with the authors' expression values** (not counts, not Ensembl; `scripts/12_bodymap_prepare.py` builds a matrix; only 13,853 symbols match MoTrPAC) — a fallback, not the primary source | **84 MB** | `curl`; already extracted |

`bodymapRat` is the source used: STAR counts, Ensembl IDs (21,040 of the 21,193 MoTrPAC genes present, 99.3%), 652 runs summed to 316 biological samples (4 of GEO's 320 are absent: aged kidney F/M, 6-week male spleen). The animal id is the replicate index within organ × stage × sex (`sex_stage_replicate`), an assumption — no animal id exists in GEO or SRA. Export it
with R to `data/external/bodymap_counts.csv` (genes × runs, `feature_ID` = ENSRNOG) and
`data/external/bodymap_meta.csv` (run, organ, sex, stage_weeks, animal/biological-replicate id,
techRep, rnaRIN). Collapse technical replicates (sum counts per biological sample) before use.
Check the Ensembl release: MoTrPAC counts use the Rnor_6.0 annotation of the package; ID overlap
with the panel genes must be reported, and the (few) panel genes missing from BodyMap listed.

## Organ mapping to MoTrPAC tissue codes

| BodyMap organ | MoTrPAC code(s) | How to score |
|---|---|---|
| Adrenal | ADRNL | exact |
| Brain (whole) | CORTEX, HIPPOC, HYPOTH | no exact match: count any brain-region prediction as correct ("brain" super-class), report the region distribution |
| Heart | HEART | exact |
| Kidney | KIDNEY | exact |
| Liver | LIVER | exact |
| Lung | LUNG | exact |
| Muscle (skeletal; site to confirm from the SRA metadata) | SKM-GN, SKM-VL | "skeletal muscle" super-class; report the GN/VL split |
| Spleen | SPLEEN | exact |
| Testes | TESTES | exact (male only) |
| Thymus | — | **out-of-distribution**: no MoTrPAC class; the right outcome is an empty or very large prediction set, never a confident label |
| Uterus | — | **out-of-distribution** (OVARY is not uterus); same treatment as thymus |

MoTrPAC tissues with no BodyMap counterpart: BAT, WAT-SC, BLOOD, COLON, SMLINT, VENACV, OVARY,
PLASMA (no transcripts anyway). They stay in the classifier as possible labels.

## Protocol (`scripts/12_bodymap_validate.py`; results in `results/12_bodymap/` and the REPORT section "12 · External validation"; the variance prefilter is applied to raw log values before z-scoring — an earlier run had it after, which made the 5,000-gene prefilter arbitrary; `scripts/14_transfer_representations.py` compares z-score, rank and pair representations and the transferability-aware selector)

1. Fit the panel and classifier on **all** MoTrPAC animals (BodyMap is fully external, so no
   MoTrPAC animal is held out for accuracy). Keep a MoTrPAC animal-grouped calibration split for
   the "rat-calibrated" prediction sets.
2. Normalize BodyMap exactly like MoTrPAC counts: log2 CPM on the total library size
   (`motrpac.io.log_cpm`), restricted to the panel genes. Run twice: (a) no re-standardization
   (same species, same units — the honest first attempt), (b) per-gene z-score within each
   dataset, and report both.
3. Report per mapped organ: accuracy, the confusion matrix against the 19 MoTrPAC classes, and,
   for thymus and uterus, the set-size distribution (they should not be covered by singletons).
4. Coverage: (a) with MoTrPAC calibration, (b) after re-calibrating on a BodyMap calibration split
   grouped by animal, on the remaining BodyMap animals. The gap is the headline.
5. Stratify by age (2 / 6 / 21 / 104 weeks) and sex: does the guarantee hold for juveniles and
   aged rats, i.e. outside the 6-month window the panel was selected in?

## What the run showed about the juveniles and the untested markers (2026-09-18; `results/12_bodymap/juvenile_marker_check.csv`, `panel_gene_check.csv`)

- At 2 weeks the 20-gene panel calls testes SKM-VL (0 of 4 correct) and spleen BLOOD (0 of 8), against 4 of 4 and
  8 of 8 at 21 weeks. The panel's only testis marker, Pgk2, is a spermatid-specific glycolytic isozyme that is not
  expressed before puberty (0.2 log2 CPM in 2-week testes vs 9.7 in adults), while the SKM-VL marker Mybph is higher
  in the juvenile testis (4.6 vs 2.7), so the juvenile testis carries muscle evidence and no testis evidence. The
  juvenile spleen, still erythropoietic, expresses more of the blood marker Hbq1b (3.3 vs 1.7) and less of the spleen
  marker Fcrl5 (4.6 vs 5.8), so it is called blood. Both are developmental-stage effects on single markers: the panel
  is an adult fingerprint, and a one-marker-per-tissue panel has no redundancy to absorb them.
- Seven of the 20 panel markers have no BodyMap organ — Akr1c3 (OVARY), Trim29 (WAT-SC), Hbq1b (BLOOD), Otop1 (BAT),
  Ccl25 (SMLINT), Abo3 (COLON), Gdf10 (VENACV) — so ovary, adipose, blood, BAT, small intestine, colon and vena cava
  are **untested within species**; their only external check is the human transfer (`docs/GTEX_TRANSFER.md`), where
  adipose and the gonads are the tissues that fail.
- Conformal variants on the target: `results/12_bodymap/coverage_by_organ.csv` has the per-organ coverage of the
  marginal, Mondrian and marginal-floor Mondrian sets (each class quantile floored at the marginal one) under MoTrPAC
  calibration; the comparison with the in-distribution and held-out-sex numbers is in the "Conditional coverage"
  block of `results/SUMMARY.md`.
- Accuracy convention: sample-weighted over the mapped samples (every sample counts once), as for GTEx.

## Confounds to state up front

- Library chemistry: rRNA-depleted total RNA (BodyMap) vs poly-A mRNA (MoTrPAC) changes the
  intronic/non-coding read share; genes with long 3' UTR or high intronic signal shift most.
  Panel genes flagged as `pct_mrna`-correlated in `results/05_panels/TRNSCRPT/candidate_panel_annotated.csv`
  are the ones to watch.
- Age: 2-week rats are juveniles; 104-week rats are aged. Genes that are training-regulated in
  MoTrPAC (flagged in the same file) may also be age- or activity-sensitive.
- Whole brain vs cortex/hippocampus/hypothalamus; unknown muscle site; different Ensembl release;
  technical replicates must be collapsed before any split.
