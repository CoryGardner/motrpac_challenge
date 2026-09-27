# Phase 2 — data discovery

Built by `scripts/multiomic/02_discovery.py` on 2026-09-27 07:24 UTC. Every attempt is in `attempts.csv`; every download (URL, bytes, sha256, time) in `download_log.csv` (19 files, 237 MB); usable datasets in `datasets.csv`; tissue maps in `tissue_maps.csv`.

## Attempts (prompt order)

| order | dataset | source | status |
|---|---|---|---|
| 1 | Jiang 2020, Cell 183:269 — Quantitative Proteome Map of the Human Body | PRIDE PXD016999 (FTP/HTTPS) | no processed table |
| 1 | Jiang 2020 | Elsevier CDN supplementary tables (ars.els-cdn.com, PII S0092867420310783, mmc1–mmc8) | obtained |
| 1 | Jiang 2020 | GTEx portal eGTEx 'Proteomics' fileset (gtexportal.org API) | alternative source |
| 2 | Wang 2019, Mol Syst Biol 15:e8503 — 29 healthy human tissues | PRIDE PXD010154 (2019/07) | over the box; ranged read of proteinGroups.txt possible |
| 2 | Wang 2019 | embopress.org → link.springer.com article page | blocked |
| 2 | Wang 2019 | Europe PMC REST supplementaryFiles bundle (PMC6379049) | see below |
| 3 | Geiger 2013, MCP 12:1709 — 28 mouse tissues (SILAC mouse) | Elsevier CDN (PII S1535947620310860, mmc1.zip 9.8 MB) | obtained |
| 3 | newer mouse multi-tissue proteome | PRIDE search API (keyword 'tissue atlas', Mus musculus) | not pursued |
| 4 | rat multi-tissue proteome | PRIDE search API (Rattus norvegicus; keywords 'tissue', 'atlas organs'), ProteomeCentral PROXI (species 10116) | negative — no rat BodyMap-equivalent proteome found |
| 5 | Sato 2022, Cell Metab 34:329 — Atlas of exercise metabolism (mouse) | Metabolomics Workbench (title search over all 4,589 studies) | absent |
| 5 | Sato 2022 | Elsevier CDN (PII S1550413121006355, mmc2–mmc8) | obtained |
| 6 | Metabolomics Workbench mouse/rat multi-tissue studies | REST /rest/study/study_id/ST/summary (all studies) filtered to Mus musculus / Rattus norvegicus with tissue words in the title | ST003188 obtained (datatable, mwtab); others recorded |
| 6 | MoTrPAC's own metabolomics on Metabolomics Workbench (ST0026xx) | REST | not applicable |
| 7 | MetaboLights | https://www.ebi.ac.uk/metabolights/ws/studies | not searched |
| 8 | MoTrPAC Data Hub (GCS bucket) | storage.googleapis.com/motrpac-data-hub | not needed |

## Datasets obtained

| dataset | layer | samples | features | individuals | tissues | units | batch_metadata | matched_to_motrpac | usable |
|---|---|---|---|---|---|---|---|---|---|
| Jiang 2020 human proteome map | protein (TMT, human) | 201 | 12627 | 14 | 32 | log2 relative abundance to the run's pooled reference (cleaned); raw reporter intensities also given | TMT run (56), TMT tag (10), donor | 3097 of 3570 RII genes by 1:1 ortholog | yes — Phases 3, 5, 6 |
| Wang 2019 human tissue atlas (29 tissues) | protein (label-free iBAQ, human) | 29 | 13640 | one donor per tissue (Table EV1 A) | 29 | iBAQ | MS run per tissue (one sample per tissue) | 3123 of 3570 RII genes by 1:1 ortholog | yes — Phase 3 secondary target (one sample per tissue), Phase 5 (matched RNA), Phase 6 |
| Geiger 2013 mouse tissue proteome (28 tissues) | protein (SILAC H/L ratios to one SILAC-mouse standard; mouse) | 28 | 7349 | pooled mice per tissue | 28 | normalised H/L SILAC ratio (common heavy standard) | none per sample (one run per tissue) | 2800 RII genes by gene symbol (rat↔mouse, caveat: symbol match, not orthology table) | yes — Phase 3 secondary target (covers KIDNEY and WAT-SC, n = 1 per tissue) |
| Sato 2022 atlas of exercise metabolism (mouse) | metabolite (Metabolon HD4 untargeted; mouse) | 191 | 547–804 named per tissue (1159 union) | 24 | 8 | peak area (raw) / median-scaled | ROUND, RUN DAY (per sample), platform | RefMet-matched names: 100 of 1159 queried, 7 in MoTrPAC (any platform), 7 in HILIC+; crude name overlap 191 | yes — Phase 4 (target + exercise-invariance test), Phase 6 |
| Metabolomics Workbench ST003188 — A metabolic atlas of mouse aging (Mullen lab, USC, 2025) | metabolite (targeted RP-negative triple-quad; mouse) | 840 | 190 | 70 | 12 | AU (peak area) | Batch (12 levels — one per organ, nested) | 109 RefMet names in MoTrPAC (any platform), 46 in HILIC+ | yes — Phase 4 (second target), Phase 6 (nested design) |

## Tissue maps

| dataset | external_tissue | n | motrpac_tissue |
|---|---|---|---|
| Jiang 2020 | Stomach | 11 | OOD |
| Jiang 2020 | Thyroid | 11 | OOD |
| Jiang 2020 | Colon - Sigmoid | 11 | OOD |
| Jiang 2020 | Heart - Atrial Appendage | 11 | HEART |
| Jiang 2020 | Muscle - Skeletal | 11 | SKM-GN (RNA: SKM-GN/SKM-VL) |
| Jiang 2020 | Esophagus - Gastroesophageal Junction | 10 | OOD |
| Jiang 2020 | Artery - Aorta | 10 | VENACV (RNA only, imperfect) |
| Jiang 2020 | Pancreas | 8 | OOD |
| Jiang 2020 | Spleen | 8 | SPLEEN (RNA only) |
| Jiang 2020 | Skin - Not Sun Exposed (Suprapubic) | 8 | OOD |
| Jiang 2020 | Lung | 8 | LUNG |
| Jiang 2020 | Esophagus - Muscularis | 8 | OOD |
| Jiang 2020 | Adrenal Gland | 7 | ADRNL (RNA only) |
| Jiang 2020 | Vagina | 7 | OOD |
| Jiang 2020 | Skin - Sun Exposed (Lower leg) | 7 | OOD |
| Jiang 2020 | Heart - Left Ventricle | 7 | HEART |
| Jiang 2020 | Breast - Mammary Tissue | 7 | OOD |
| Jiang 2020 | Uterus | 6 | OOD |
| Jiang 2020 | Nerve - Tibial | 6 | OOD |
| Jiang 2020 | Liver | 5 | LIVER |
| Jiang 2020 | Ovary | 4 | OVARY (RNA only) |
| Jiang 2020 | Esophagus - Mucosa | 4 | OOD |
| Jiang 2020 | Colon - Transverse | 4 | COLON (RNA only) |
| Jiang 2020 | Artery - Tibial | 4 | OOD |
| Jiang 2020 | Small Intestine - Terminal Ileum | 3 | SMLINT (RNA only) |
| Jiang 2020 | Testis | 3 | TESTES (RNA only) |
| Jiang 2020 | Brain - Cerebellum | 3 | OOD |
| Jiang 2020 | Brain - Cortex | 2 | CORTEX |
| Jiang 2020 | Pituitary | 2 | OOD |
| Jiang 2020 | Artery - Coronary | 2 | OOD |
| Jiang 2020 | Prostate | 2 | OOD |
| Jiang 2020 | Minor Salivary Gland | 1 | OOD |
| Geiger 2013 | Adrenal gland | 1 | ADRNL (RNA only) |
| Geiger 2013 | Brain cortex | 1 | CORTEX |
| Geiger 2013 | Brain medulla | 1 | OOD |
| Geiger 2013 | Brown fat | 1 | BAT (RNA only) |
| Geiger 2013 | Cerebellum | 1 | OOD |
| Geiger 2013 | Colon | 1 | COLON (RNA only) |
| Geiger 2013 | Diaphragm | 1 | OOD |
| Geiger 2013 | Duodenum | 1 | OOD |
| Geiger 2013 | Embryonic tissue | 1 | OOD |
| Geiger 2013 | Eye | 1 | OOD |
| Geiger 2013 | Heart | 1 | HEART |
| Geiger 2013 | Ileum | 1 | OOD |
| Geiger 2013 | Jejunum | 1 | OOD |
| Geiger 2013 | Kidney cortex | 1 | KIDNEY |
| Geiger 2013 | Kidney medulla | 1 | OOD |
| Geiger 2013 | Liver | 1 | LIVER |
| Geiger 2013 | Lung | 1 | LUNG |
| Geiger 2013 | Midbrain | 1 | OOD |
| Geiger 2013 | Muscle | 1 | SKM-GN |
| Geiger 2013 | Olfactory bulb | 1 | OOD |
| Geiger 2013 | Ovary | 1 | OVARY (RNA only) |
| Geiger 2013 | Pancreas | 1 | OOD |
| Geiger 2013 | Salivary gland | 1 | OOD |
| Geiger 2013 | Spleeen | 1 | SPLEEN (RNA only) |
| Geiger 2013 | Stomach | 1 | OOD |
| Geiger 2013 | Thymus | 1 | OOD |
| Geiger 2013 | Uterus | 1 | OOD |
| Geiger 2013 | White fat | 1 | WAT-SC |
| Sato 2022 | BAT | 24 | BAT |
| Sato 2022 | eWAT | 24 | WAT-SC (visceral vs subcutaneous: imperfect) |
| Sato 2022 | HEART | 24 | HEART |
| Sato 2022 | HYPOTHALAMUS | 24 | HYPOTH |
| Sato 2022 | iWAT | 23 | WAT-SC |
| Sato 2022 | LIVER | 24 | LIVER |
| Sato 2022 | MUSCLE | 24 | SKM-GN/SKM-VL |
| Sato 2022 | SERUM | 24 | PLASMA (serum vs plasma) |
| MW ST003188 | Bladder | 70 | OOD |
| MW ST003188 | Brain | 70 | CORTEX/HIPPOC/HYPOTH |
| MW ST003188 | Heart | 70 | HEART |
| MW ST003188 | Kidney | 70 | KIDNEY |
| MW ST003188 | Liver | 70 | LIVER |
| MW ST003188 | Lung | 70 | LUNG |
| MW ST003188 | Muscle (Quad) | 70 | SKM-GN/SKM-VL |
| MW ST003188 | Pancreas | 70 | OOD |
| MW ST003188 | Plasma | 70 | PLASMA |
| MW ST003188 | Spleen | 70 | SPLEEN |
| MW ST003188 | Thymus | 70 | OOD |
| MW ST003188 | Tongue | 70 | OOD |

## Negative results

- No rat multi-tissue proteome with processed tables was found on PRIDE or ProteomeXchange (keyword searches); the BodyMap-equivalent proteome does not exist in the public repositories searched.
- Sato 2022 is not on Metabolomics Workbench; its supplementary Metabolon tables were used instead.
- Wang 2019: the processed MaxQuant bundle on PRIDE is 28.2 GB, over the 20 GB box; the Europe PMC supplementary bundle is the route (see `datasets.csv` for its state).
- Kidney and adipose are absent from Jiang 2020, so two of the seven proteomics tissues have no human protein target; Geiger 2013 (mouse, one pooled sample per tissue) covers both.
