# Phase 6 — identifiability of the external designs (metadata only)

Built by `scripts/multiomic/06_external_identifiability.py` on 2026-09-27 07:49 UTC; numbers from the CSVs here and from `results_frozen/16_identifiability/`.

## Design comparison — `design_comparison.csv`

| dataset | batch_variable | n_samples | n_tissues | n_levels | max_tissues_per_level | cramers_v | n_pairs_estimable | n_pairs_total | source | frac_pairs_estimable |
|---|---|---|---|---|---|---|---|---|---|---|
| MoTrPAC TRNSCRPT | RNA_extr_plate_ID | 899 | 19 | 17 | 2 | 1.000 | 1 | 171 | results_frozen/16_identifiability | 0.006 |
| MoTrPAC METHYL | DNA_extr_plate_ID | 400 | 8 | 8 | 1 | 1.000 | 0 | 28 | results_frozen/16_identifiability | 0.000 |
| MoTrPAC ATAC | Nuclei_extr_date | 500 | 8 | 25 | 2 | 0.918 | 0 | 28 | results_frozen/16_identifiability | 0.000 |
| MoTrPAC PROT | plex_id | 420 | 7 | 42 | 1 | 1.000 | 0 | 21 | results_frozen/16_identifiability | 0.000 |
| MoTrPAC PHOSPHO | plex_id | 420 | 7 | 42 | 1 | 1.000 | 0 | 21 | results_frozen/16_identifiability | 0.000 |
| MoTrPAC ACETYL | plex_id | 108 | 2 | 12 | 1 | 1.000 | 0 | 1 | results_frozen/16_identifiability | 0.000 |
| MoTrPAC UBIQ | plex_id | 120 | 2 | 12 | 1 | 1.000 | 0 | 1 | results_frozen/16_identifiability | 0.000 |
| MoTrPAC IMMUNO | plate_id | 1375 | 17 | 20 | 4 | 0.612 | 16 | 136 | results_frozen/16_identifiability | 0.118 |
| Jiang2020 | tmt_run | 448 | 33 | 56 | 8 | 0.324 | 424 | 528 | results_multiomic/06_external_identifiability | 0.803 |
| Wang2019 | ms_experiment | 50 | 31 | 47 | 1 | 1.000 | 0 | 465 | results_multiomic/06_external_identifiability | 0.000 |
| Sato2022 | round | 191 | 8 | 6 | 7 | 0.447 | 21 | 28 | results_multiomic/06_external_identifiability | 0.750 |
| MW_ST003188 | batch | 840 | 12 | 12 | 1 | 1.000 | 0 | 66 | results_multiomic/06_external_identifiability | 0.000 |

## Per variable — `nesting_all.csv`

| dataset | variable | n_levels | n_tissues | median_levels_per_tissue | max_levels_per_tissue | max_tissues_per_level | n_levels_shared | tissues_in_one_level | cramers_v | n_pairs_sharing_level | n_pairs_total |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Jiang2020 | tmt_run | 56 | 33 | 14.000 | 28 | 8 | 56 | 0 | 0.324 | 424 | 528 |
| Jiang2020 | tmt_tag | 8 | 33 | 6.000 | 8 | 24 | 8 | 1 | 0.263 | 503 | 528 |
| Jiang2020 | donor | 14 | 32 | 7.000 | 11 | 18 | 14 | 1 | 0.310 | 431 | 496 |
| Wang2019 | ms_experiment | 47 | 31 | 1.000 | 9 | 1 | 0 | 24 | 1.000 | 0 | 465 |
| Sato2022 | round | 6 | 8 | 3.000 | 3 | 7 | 3 | 0 | 0.447 | 21 | 28 |
| Sato2022 | run_day | 4 | 2 | 4.000 | 4 | 2 | 4 | 0 | 0.000 | 1 | 1 |
| Sato2022 | animal | 24 | 8 | 24.000 | 24 | 8 | 24 | 0 | 0.027 | 28 | 28 |
| MW_ST003188 | batch | 12 | 12 | 1.000 | 1 | 1 | 0 | 12 | 1.000 | 0 | 66 |
| MW_ST003188 | mouse | 70 | 12 | 70.000 | 70 | 12 | 70 | 0 | 0.000 | 66 | 66 |

## Estimable pairs — `estimable_pairs.csv`

| dataset | n_samples | n_tissues | n_pairs_total | n_pairs_estimable | variables_used | frac_pairs_estimable | note |
|---|---|---|---|---|---|---|---|
| Jiang2020 | 448 | 33 | 528 | 424 | tmt_run | 0.803 | TMT run = plex (10 channels + 1 reference); the batch variable of a TMT design |
| Wang2019 | 50 | 31 | 465 | 0 | ms_experiment | 0.000 | one MS experiment (label-free run) per tissue sample |
| Sato2022 | 191 | 8 | 28 | 21 | round | 0.750 | Metabolon ROUND (per tissue table) and RUN DAY where given; ROUND labels repeat across tissue tables the way MoTrPAC's S1–S6 do, so a shared label is not evidence of a shared run |
| MW_ST003188 | 840 | 12 | 66 | 0 | batch | 0.000 | Batch from the mwtab SUBJECT_SAMPLE_FACTORS |

Not assessable: Geiger2013: one pooled sample per tissue, no per-sample run/batch variable in the supplementary tables (`unavailable.csv`).

## Reading

- **Jiang 2020 is the positive counterexample.** Its TMT runs each hold up to 8 different tissues (Cramér's V between run and tissue 0.324, 56 runs, 424 of 528 tissue pairs share a run); a donor also spans up to 11 runs. Tissue is therefore separable from plex inside that study — the same TMT chemistry as MoTrPAC, a different allocation of samples to plexes.
- MoTrPAC proteomics (plex_id, V = 1.000, 0 of 21 pairs), MoTrPAC RNA-seq (1 of 171 pairs), Wang 2019 (one run per tissue), MW ST003188 (one batch per organ) are the nested designs.
- Nesting is a choice of design, not a property of the assay: TMT can hold ten tissues in one plex, and when it does the audit's confound disappears.
