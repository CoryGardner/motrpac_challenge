# Third-party data and components

The code, the site and the derived tables written by this repository's scripts are MIT-licensed (`LICENSE`). The
repository also **contains third-party data**, in small, subset or transformed form; each item keeps its own terms.
The full raw datasets are not in the repository.

## Rat BodyMap (CC BY 4.0)

- **What:**
  - `tests/fixtures/bodymap_counts_subset.csv.gz`: the raw counts of 4 rat BodyMap samples, every gene (about 32,600 rows).
  - `site/data/expr_bodymap.json`: log2 CPM of at most 100 genes for the 316 BodyMap samples.
  - `site/data/samples_bodymap.json`: the per-sample predictions of the models for the same samples.
  - The BodyMap result tables under `results_frozen/12_bodymap/` and `results_frozen/31_site_regen/12_bodymap/`.
- **Source:** Yu et al., *Nature Communications* 5:3230 (2014); GEO GSE53960; obtained through the Bioconductor
  ExperimentHub package `bodymapRat` 1.28.0.
- **Licence:** CC BY 4.0, https://creativecommons.org/licenses/by/4.0/.
- **Changes:**
  - technical runs summed per sample;
  - subset to a few samples (fixture) or to at most 100 genes (site);
  - counts converted to log2 CPM on the total library;
  - models scored on them.

## GTEx v8 (GTEx data-use policy)

- **What:**
  - `site/data/expr_gtex.json`: log2(TPM + 1) of at most 100 genes, 1:1 orthologs, on a donor-capped subset.
  - `site/data/samples_gtex.json`: the per-sample predictions for that subset.
  - The GTEx result tables under `results_frozen/13_gtex/` and `results_frozen/31_site_regen/13_gtex/`.
- **Source:** GTEx v8 open-access gene expression and sample attributes (release 2017-06-05, RNASeQCv1.1.9);
  GTEx Consortium, *Science* 369:1318–1330 (2020).
- **Terms:** the GTEx data-use policy, https://gtexportal.org/home/license.
- **Changes:**
  - subset to 17 tissues, at most 150 donors per tissue, and at most 100 genes;
  - mapped to rat 1:1 orthologs;
  - log2(TPM + 1).

## MoTrPAC rat endurance training, release c1.0 (MoTrPAC data-use terms)

- **What:**
  - `site/data/expr_motrpac.json`: log2 CPM of at most 100 genes for the 899 RNA-seq vials.
  - `site/data/samples_motrpac.json`: per-vial predictions.
  - The other MoTrPAC-derived tables in `site/data/` and `results_frozen/`: per-fold and per-sample scores, QC summaries, and
    per-tissue means of the panel genes.
- **Source:** the R package `MotrpacRatTraining6moData` 2.0.0, GitHub commit
  `f831a4fe421ec11687640452484a8247137aa74a`. The package is MIT-licensed; the data carry the consortium's data-use terms,
  https://motrpac-data.org/. MoTrPAC Study Group, *Nature* 629:174–183 (2024).
- **Acknowledgement** (from `docs/DATA_GUIDE.md` §8): "Data used in the preparation of this article were obtained from
  the Molecular Transducers of Physical Activity Consortium (MoTrPAC) MotrpacRatTraining6moData R package [version]."
  `site/data/manifest.json` records the version.
- **Changes:** exported to CSV, counts converted to log2 CPM, subset to the genes the site shows.

## Multiomic follow-up (`results_multiomic/`)

The follow-up commits only derived tables. These are:

- per-protein and per-metabolite summaries;
- variance partitions and PCA scores;
- panel genes;
- per-sample transfer predictions and accuracies;
- recalibration and design tables.

No source matrix is committed. The sources are listed below; the download URLs and sha256 are in
`docs/REPRODUCIBILITY.md` and `results_multiomic/02_discovery/download_log.csv`.

- **MoTrPAC proteomics reporter-ion files** (Data Hub quant-id folders, `prot-pr`, releases c1.0 and c2.0): the
  consortium's data-use terms.
  - Derived: `results_multiomic/01_rii/` (per-vial PCA scores `pca_scores_*.csv`, variance partitions, panel and
    stability tables, per-tissue summaries) and `rii_meta.csv` (vial metadata: tissue, sex, plex, channel).
- **Jiang et al. 2020**, *Cell* 183:269, doi:10.1016/j.cell.2020.08.036, the GTEx tissue proteome (PRIDE PXD016999),
  supplementary tables: the journal's terms (Elsevier supplementary material).
  - Derived: `results_multiomic/03_prot_transfer/`, `rawppm/`, `05_fusion_transfer/`, `08_verification/`: per-sample
    predictions and per-tissue accuracies on the 44 or 94 mapped samples.
- **Wang et al. 2019**, *Mol Syst Biol* 15:e8503, doi:10.15252/msb.20188503 (PRIDE PXD010154): CC BY 4.0.
  - Derived: `results_multiomic/03_prot_transfer/wang2019/`.
- **Geiger et al. 2013**, *Mol Cell Proteomics* 12:1709, doi:10.1074/mcp.M112.024919, supplementary table: the
  journal's terms.
  - Derived: `results_multiomic/03_prot_transfer/geiger2013/`.
- **Sato et al. 2022**, *Cell Metab* 34:329, doi:10.1016/j.cmet.2021.12.016, supplementary tables: the journal's terms.
  - Derived: `results_multiomic/04_metab_transfer/hilic_sato/`, `deep_sato/`, and `08_verification/sato_invariance_design.csv`.
- **Metabolomics Workbench ST003188**, "A metabolic atlas of mouse aging" (Mullen Lab, University of Southern
  California; project PR001984, doi:10.21228/M88J0W): CC BY 4.0.
  - Derived: `results_multiomic/04_metab_transfer/deep_mw/`, `hilic_mw/`.

## Software

- **Plotly.js 2.35.2** (`site/vendor/plotly-cartesian-2.35.2.min.js`): MIT, Plotly, Inc.; the licence text is
  `site/vendor/LICENSE.plotly.txt`.
