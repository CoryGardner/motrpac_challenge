# DATA_GUIDE — the MoTrPAC rat endurance-training data

Sources: `MotrpacRatTraining6moData` R package v2.0.0 (https://motrpac.github.io/MotrpacRatTraining6moData/;
installed from GitHub commit `f831a4f`, tag v2.1.0, whose DESCRIPTION also reads 2.0.0 — `R/install_deps.R` pins it),
its vignette, and the companion `MotrpacRatTraining6mo` analysis package. The Nature 2024 paper (MoTrPAC Study
Group, *Temporal dynamics of the multi-omic response to endurance exercise training*, Nature 629:174–183)
explains the normalization and differential-analysis choices in its supplementary methods.

## 1. Study design

- 6-month-old rats (Fischer 344), **both sexes**.
- Progressive treadmill endurance training for **1, 2, 4, or 8 weeks**; tissues collected
  **48 h after the last bout**. Sex-matched **sedentary controls**.
- **Whole blood, plasma, and 18 solid tissues**; most assays done in a subset of tissues.
- **3–6 animals per sex per time point per assay.** So per tissue you have roughly
  2 sexes × 5 groups × ~5 = **~50–60 animals**. This is the single most important number for
  planning evaluation: everything is small-n, and folds must be counted in animals, not samples.

## 2. Identifiers (this is where leakage happens)

| ID | Meaning | Notes |
|---|---|---|
| `pid` | animal ("participant") | 8-digit. **Group splits on this.** |
| `bid` | collection event — **one per animal in this study** (147 ↔ 147 pids), not per tissue | 5-digit; `bid == viallabel[:5]` |
| `viallabel` | one vial (one assay's aliquot of a biospecimen) | 11-digit; column names in sample-level tables |
| `labelid` | specimen label | rarely needed |

Different assays on the same tissue use *different vials* of the *same biospecimen*, so to
combine transcriptomics and proteomics from the same animal+tissue you join on `bid` (or `pid`
when tissues differ). Never join on `viallabel`.

## 3. Tissue and assay codes

Tissues (`TISSUE_ABBREV`): ADRNL adrenal · BAT brown adipose · BLOOD whole blood · COLON ·
CORTEX cerebral cortex · HEART · HIPPOC hippocampus · HYPOTH hypothalamus · KIDNEY · LIVER ·
LUNG · OVARY · PLASMA · SKM-GN gastrocnemius · SKM-VL vastus lateralis · SMLINT small intestine ·
SPLEEN · TESTES · VENACV vena cava · WAT-SC subcutaneous white adipose.

Assays (`ASSAY_ABBREV`): TRNSCRPT RNA-seq · PROT global proteomics · PHOSPHO · ACETYL · UBIQ ·
METAB metabolomics/lipidomics · IMMUNO multiplex immunoassays · ATAC · METHYL (RRBS).
ATAC and METHYL data are **not** in the R package (they are on the MoTrPAC Data Hub); only their sample
metadata (`ATAC_META`, `METHYL_META` → `data/raw/meta/ATAC.csv`, `METHYL.csv`) is exported. Not used by this
pipeline except their sample metadata (phase 16, the identifiability audit).

Hyphenated codes (`SKM-GN`, `SKM-VL`, `WAT-SC`) become `SKMGN`, `SKMVL`, `WATSC` in R object
names and in our file names. Values inside the `tissue` column keep the hyphen.

Coverage, from `make inventory` on package v2.0.0 (2026-09-17): TRNSCRPT in 19 tissues (18 solid +
BLOOD; 50 samples = 50 animals per tissue, OVARY 24, TESTES 25; 13.8k–17.4k genes after the
package's expression filter, 32,883 in raw counts); PROT and PHOSPHO in 7 (CORTEX, HEART, KIDNEY,
LIVER, LUNG, SKM-GN, WAT-SC; 57–60 samples, i.e. 6 per sex × group, more animals than
transcriptomics); ACETYL and UBIQ in HEART, LIVER; METAB in 19 (PLASMA instead of BLOOD;
45–52 samples, VENACV 32); IMMUNO in 17 (30 samples). The **fusion set** is the tissues with all three
of TRNSCRPT + PROT + METAB: the 7 PROT tissues, with 47–50 animals shared; `02_inventory.py` prints them.
Unmatched samples: 0.

## 4. Data objects and our export layout

| R object pattern | Content | Exported to |
|---|---|---|
| `PHENO` | 5,955 samples × 510 phenotype/training/terminal variables, rows keyed by `viallabel` | `data/raw/pheno.csv` |
| `<ASSAY>_<TISSUE>_NORM_DATA` | normalized sample-level data; cols `feature`, `feature_ID`, `tissue`, `assay`, then viallabels | `data/raw/norm/<ASSAY>__<TISSUE>.csv` |
| `TRNSCRPT_<TISSUE>_RAW_COUNTS` | RNA-seq raw counts, same layout | `data/raw/counts/TRNSCRPT__<TISSUE>.csv` |
| `<ASSAY>_<TISSUE>_DA` | differential analysis per feature × sex × comparison_group (DESeq2 / limma) | `data/raw/da/<ASSAY>__<TISSUE>.csv` |
| `<ASSAY>_META` | assay-level sample metadata (batch, plex, QC) | `data/raw/meta/<ASSAY>.csv` |
| `METAB_NORM_DATA_NESTED[[platform]][[tissue]]` | metabolomics: 13 platforms, each leaf features × viallabels (feature ID in rownames); **not** per-tissue objects | `data/raw/norm/METAB__<TISSUE>.csv` (flattened, see below) + `meta/METAB_VIALS.csv`, `meta/METAB_FEATURES.csv` |
| `IMMUNO_NORM_DATA_NESTED[[panel]][[tissue]]` | immunoassays: 6 panels, each leaf samples × analytes with a `viallabel` column | `data/raw/norm/IMMUNO__<TISSUE>.csv` (flattened the same way) + `meta/IMMUNO_VIALS.csv`, `meta/IMMUNO_FEATURES.csv` |
| `FEATURE_TO_GENE` | feature_ID → gene (4M rows; export is subset to features present) | `data/raw/feature_to_gene.csv` |
| `RAT_TO_HUMAN_GENE` | rat ↔ human orthologs | `data/raw/rat_to_human_gene.csv` |
| `METAB_FEATURE_ID_MAP` | metabolite feature annotation | `data/raw/metab_feature_id_map.csv` |
| `TRAINING_REGULATED_FEATURES` | features regulated by training at 5% FDR | `data/raw/training_regulated_features.csv` |
| `OUTLIERS` | samples excluded from DA | `data/raw/outliers.csv` |
| `TISSUE_ABBREV`, `ASSAY_ABBREV`, `GROUP_COLORS` | code tables | `data/raw/codes/*.csv` |

`R/export_motrpac.R` discovers objects by regex, so it survives minor naming differences
(e.g., `HEART_PROT_DA` vs `PROT_HEART_DA`), and writes `data/raw/manifest.json` listing what it found.

In `*_NORM_DATA`, the `feature` column is non-NA **only for training-regulated features**; use
`feature_ID` as the key.

**Flattened METAB / IMMUNO tables.** Each platform (metabolomics) or panel (immunoassay) received its
own vial of a biospecimen, so viallabels differ across platforms for the same `bid`. The exporter
averages vials within a platform, joins platforms on `bid`, and writes **one column per biospecimen,
named by one representative viallabel** (one present in PHENO when possible) so `bid == viallabel[:5]`
and the PHENO lookup keep working; `meta/<ASSAY>_VIALS.csv` maps every original vial to its platform,
bid and exported column. Features measured on several platforms appear as **duplicate `feature_ID`
rows** (platforms in sorted name order; `tfp.io` suffixes the repeats with `__dupN`);
`meta/<ASSAY>_FEATURES.csv` says which platform each row came from. Tissues outside the nine core
metabolomics tissues have only the `metab-u-hilicpos` platform, so METAB feature counts differ a lot
by tissue and the cross-tissue inner join is small. `METAB_<TISSUE>_DA` uses `feature_ID_da`, which
differs from the sample-data ID for meta-regression-merged metabolites (see `metab_feature_id_map.csv`).
There is no `METAB_META` object; `TRAINING_REGULATED_NORM_DATA` and `*_NORM_DATA_05FDR` are
deliberately not exported (manifest `skipped` lists every object with a reason).

## 5. Key `PHENO` columns

- `pid`, `bid`, `viallabel`
- `sex` (`male`/`female`; raw form `registration___sex` is `"1"` female, `"2"` male)
- `group` (`control`, `1w`, `2w`, `4w`, `8w`) and `study_group_timepoint`
- `specimen_processing___sampletypedescription` (tissue name, long form)
- physiology worth having around: `vo2.max.test.vo2_max`, `nmr.testing.nmr_fat`,
  `terminal.weight.bw`, `calculated.variables.pct_body_fat_change`, `time_to_freeze`

`tfp.io.load_pheno()` standardizes these into `pid, bid, sex, group, tissue_long` and keeps
the rest.

In package v2.0.0 as installed, PHENO is **6,156 vials × 509 variables from 147 animals**, `pid`/`bid`
are integers, `viallabel` is a character column, `sex`/`group`/`tissue` are already in the clean form,
and the raw-style names use dots (`registration.sex`, `specimen.processing.sampletypedescription`),
not `___`. `load_pheno()` accepts both spellings. `bid == viallabel[:5]` holds for every row.
**Export bug:** `data/raw/pheno.csv` has `labelid = 0` in all rows (integer64 lost on export); use `viallabel`/`bid`/`pid`.
**UBIQ:** `UBIQ_*_NORM_DATA`/`_DA` are the portal's *protein-corrected* files (PHOSPHO/ACETYL are not).

## 6. Gotchas that will cost time if not known

1. **Cross-tissue comparability of proteomics.** TMT proteomics is normalized per tissue,
   typically as log-ratios to a per-tissue reference channel, so absolute levels are not
   comparable across tissues the way log-CPM transcripts are. A cross-tissue proteomic
   "fingerprint" trained on `PROT_*_NORM_DATA` may be learning the normalization, not biology.
   Checked in phase 04 (`04_prot_diagnostic.py`): proteomics is treated as within-tissue only.
   For the cross-tissue classification task, transcriptomics from **raw counts re-normalized
   uniformly** (log2 CPM per sample, done by `tfp.io.load_counts`) is the safest layer.
   Metabolomics has a similar issue (platform-specific, per-tissue runs).
2. **Feature-ID spaces differ across tissues for some assays.** Transcripts share Ensembl gene
   IDs across tissues (good). Proteomics feature IDs may be protein accessions, sometimes with
   isoform suffixes; metabolite IDs are shared for named metabolites but not for unknowns.
   `stack_tissues()` reports intersection vs. union sizes — look at them.
3. **Tissue identity from transcriptomics is trivially separable with many features.** The
   interesting outputs are (a) the smallest panel that still works, (b) which tissue pairs are
   confusable (SKM-GN vs SKM-VL; CORTEX vs HIPPOC vs HYPOTH; WAT-SC vs BAT), (c) whether a
   panel selected on one sex/time point/training state holds on the other, and (d) whether the
   panel's genes have human orthologs and transfer.
4. **Samples from one animal are not independent** — one rat contributes up to ~18 tissue
   samples. Group every split by `pid`. This also means conformal calibration at the sample
   level is only approximately exchangeable; `tfp.conformal` has a `one_per_group` option
   that subsamples one sample per animal for the calibration set.
5. **Outliers.** `OUTLIERS` lists samples the consortium excluded from differential analysis.
   Default: keep them for classification but flag them; report both if results differ.
6. **DA tables are per sex.** Comparisons are `1w/2w/4w/8w vs control` within sex. Discordance
   analyses on DA statistics must be done within sex × time point.
7. **Sizes.** `FEATURE_TO_GENE` is ~4M rows; never load it whole in Python — the export subsets it.
   Transcript matrices are ~30k features × ~60 samples per tissue; stacked across 18 tissues
   that's ~1,000 samples — fine in memory.

## 7. Human data (for the human-transfer leg)

MoTrPAC's human protocol collects blood, vastus lateralis muscle, and subcutaneous adipose
before and after acute endurance or resistance exercise (tissue papers appeared as preprints in
March 2026; Data Hub: https://motrpac-data.org/). This pipeline's human leg uses **GTEx** for
tissue-identity transfer (see `docs/GTEX_TRANSFER.md`) and `RAT_TO_HUMAN_GENE` for orthology;
the MoTrPAC human tissues are the natural next test once their processed data are in hand.

## 8. Citing

Publications using these data must acknowledge MoTrPAC and cite the data package version:
"Data used in the preparation of this article were obtained from the Molecular Transducers of
Physical Activity Consortium (MoTrPAC) MotrpacRatTraining6moData R package [version]."
The manifest records the version.

## 9. Portal inputs for phase 16 (`MOTRPAC_PORTAL`)

Phase 16 (`scripts/16_identifiability.py`) has two optional inputs from the MoTrPAC Data Hub
(https://motrpac-data.org/), release **c1.0**, that are not part of the R-package export (20 files, ≈ 165 MB):

- the RNA-seq library-QC table, all tissues:
  `rat-training-06/c1.0/transcriptomics/qa-qc/motrpac_pass1b-06_transcript-rna-seq_qa-qc-metrics.csv`;
- the 19 per-tissue RSEM gene-count tables, read by `--bridge` (batch measured directly on the bridging
  reference pools, whose vials are in the portal count files):
  `rat-training-06/c1.0/transcriptomics/t<code>-<tissue>/transcript-rna-seq/motrpac_pass1b-06_<tissue>_transcript-rna-seq_rsem-genes-count.txt`
  for `t30-blood-rna`, `t52-hippocampus`, `t53-cortex`, `t54-hypothalamus`, `t55-gastrocnemius`, `t56-vastus-lateralis`,
  `t58-heart`, `t59-kidney`, `t60-adrenal`, `t61-colon`, `t62-spleen`, `t63-testes`, `t64-ovaries`, `t66-lung`,
  `t67-small-intestine`, `t68-liver`, `t69-brown-adipose`, `t70-white-adipose`, `t99-vena-cava`
  (`PORTAL_TISSUE_DIRS` in the script).

Download them into one directory that holds `rat-training-06/` and point `MOTRPAC_PORTAL` at it
(`export MOTRPAC_PORTAL=<dir>`, or `make identifiability MOTRPAC_PORTAL=<dir>`). `make portal-check` lists the
20 paths and says which are present. Without them phase 16 still runs: the QC-only baseline then uses the same QC
columns from the package's `TRNSCRPT_META` export (`data/raw/meta/TRNSCRPT.csv`), and the bridge block is reported
as unavailable.
