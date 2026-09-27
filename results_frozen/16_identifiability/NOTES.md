# Phase 16 — identifiability recompute

Metadata only; no expression value is used. Study vials only (viallabel starts with 9).

**Estimable pair.** Two tissues are estimable when they share a level of every processing variable listed in `variables_used` (estimable_pairs.csv), so that a contrast between them exists inside one batch. Everything else in nesting_<ASSAY>.csv is descriptive: levels per tissue, tissues per level, shared levels, Cramér's V between tissue and the variable, and the number of tissue pairs sharing at least one level.

**QC-only baseline.** Copied from notebooks/_build/sections/s14_qc_baseline.py: multinomial logistic regression (C = 1) on the consortium's per-library QC numbers, no gene, median imputation and standardisation fit inside the fold, on the phase-04 animal-grouped folds. `technical` and `composition` are kept separate because the composition fractions (mitochondrial, globin, rRNA, intronic reads, chrX/chrY) are read biologically by MoTrPAC itself and are not pure processing.

**TMT layers.** The plex labels S1–S6 repeat in every tissue; a plex is 11 channels (10 samples + the tissue's reference pool), so the physical plex is (tissue, label) and is nested in tissue by construction. The channel is shared across everything and carries no tissue information; in HEART and LIVER it encodes sex (phase 03).

**Immunoassays.** The Luminex plates are the one layer whose plates hold several tissues (the plate names list the tissue codes), so tissue pairs on a shared plate are estimable there.

Layers: TRNSCRPT (recomputed), METHYL (recomputed), ATAC (recomputed), PROT (recomputed), PHOSPHO (recomputed), ACETYL (recomputed), UBIQ (recomputed), IMMUNO (recomputed), METAB (unavailable).

**Bridging standards (--bridge).** V_batch = variance of a pool's log2 CPM across the plates it was run on (one vial per plate); V_tissue = variance of the 19 tissue means of study-vial log2 CPM; ratio per gene, median and ratio of sums over gene sets 36 reference vials, 6 pools run on more than one plate (bridge_variance.csv, bridge_variance_per_gene.csv).
