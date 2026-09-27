START 2026-09-27T06:45:07Z

# MULTIOMIC_LOG

Timestamped log of decisions, downloads (URL, size, sha256, date), failures and time spent for the overnight multiomic run on branch multiomic-overnight. Prompt: prompts/overnight.md.

## Preflight

- 2026-09-27T06:45Z branch `multiomic-overnight` confirmed (`git branch --show-current`); working tree clean apart from the untracked `data` symlink and `prompts/`.
- `ls data/raw data/external` resolves: `data` → `/home/cory/projects/MoTrPAC_back_to_transcriptome/code/pipeline/data` (raw/, external/ with bodymap and gtex subsets). `df -h .` → 1.3 TB free on `/` (≥ 300 GB ok). Downloads go to `data/external_multiomic/` (same filesystem).
- Python: `/home/cory/miniconda3/envs/motrpac-py/bin/python` (3.12.14; numpy 2.5.3, pandas 3.0.6, scikit-learn 1.9.1, scipy 1.18.1 = requirements.txt pins). `tfp` is not installed in that env, so every script runs with `PYTHONPATH=src`. 20 cores, 121 GB RAM.
- Reachability (`curl -sSIL --max-time 20`, HTTP code; HEAD unless noted):
  - PRIDE FTP `ftp://ftp.pride.ebi.ac.uk/pride/data/archive/` → 250 (FTP ok); HTTPS `https://ftp.pride.ebi.ac.uk/pride/data/archive/` → 200; PRIDE REST `/pride/ws/archive/v2/projects/PXD016999` → 200.
  - Cell/Elsevier: `https://www.cell.com/` → 403 (HEAD and GET, bot wall); `https://ars.els-cdn.com/` → 404 on the root but **200 on GET** for `content/image/1-s2.0-S0092867420310783-mmc2.xlsx` (Jiang 2020 supplement) → supplementary tables are downloadable.
  - EMBO/MSB `https://www.embopress.org/` → 200 (article pages redirect to link.springer.com).
  - UniProt REST `https://rest.uniprot.org/uniprotkb/P69905.json` → 200.
  - Ensembl REST `https://rest.ensembl.org/info/ping` → 200.
  - Metabolomics Workbench REST → 406 on HEAD, **200 on GET** (`/rest/study/study_id/ST000001/summary`); the all-study summary (`/rest/study/study_id/ST/summary`, 4,589 studies, 2.7 MB) downloaded fine.
  - MetaboLights `https://www.ebi.ac.uk/metabolights/ws/studies` → 200.
  - GTEx portal API `https://gtexportal.org/api/v2/dataset/fileList` → 200 (lists an eGTEx "Proteomics" fileset).
  - MoTrPAC Data Hub GCS bucket `storage.googleapis.com/motrpac-data-hub` → 403/401 anonymous (not needed, see next line).
- MoTrPAC RII proteomics files: **present locally**, no download needed. `/home/cory/projects/MoTrPAC_back_to_transcriptome/data/quant-id/rat-training-06/{c1.0,c2.0}/proteomics-untargeted/<tissue>/prot-{pr,ph,ac,ub}/*rii-results*.txt` — c1.0 prot-pr for all 7 tissues (cortex 251,907 / gastrocnemius 106,937 / heart 150,980 / kidney 186,016 / lung 265,008 / liver 170,492 / white adipose 248,639 peptide rows; 42–80 MB each), prot-ph for all 7, prot-ac and prot-ub for heart and liver; c2.0 has prot-pr, prot-ph, prot-ac for all 7 (v2.0). Peptide-level reporter-ion intensities with `protein_id, redundant_ids, is_contaminant, peptide_score, sequence, gene_symbol, entrez_id` plus one column per viallabel and `Ref_S1..Ref_S6`; `*vial-metadata.txt` gives `tmt_plex` (S1–S6) and `tmt11_channel` per vial. This tree is read-only; nothing is copied into the repo.
- Preflight finished 2026-09-27T06:53Z (≈ 8 min). Starting Phase 0.

## Phase 0 — setup and pre-registration

- 2026-09-27T06:53Z start. Created `results_multiomic/{00_setup,01_rii,02_discovery,03_prot_transfer,04_metab_transfer,05_fusion_transfer,06_external_identifiability,07_synthesis}`, `scripts/multiomic/` (new scripts live here; the numbered pipeline scripts are untouched), `data/external_multiomic/` (resolves to `/home/cory/projects/MoTrPAC_back_to_transcriptome/code/pipeline/data/external_multiomic`, outside the repo, like the existing `data/` tree).
- The `data` symlink is excluded from git via `.git/info/exclude` (not `.gitignore`, which is a tracked submission file).
- Wrote `docs/PREREGISTRATION_MULTIOMIC.md` (predictions (a)–(f) with pass rules, fixed analysis choices, drop policy) and the `docs/MULTIOMIC_REPORT.md` skeleton.
- Decision: new external downloads started in the background during Phase 0/1 so they run while Phase 1 computes on local data; each is logged under Phase 2 with URL, size, sha256 and time (`data/external_multiomic/_logs/downloads.tsv`, fetch scripts beside it). Started: Jiang 2020 mmc1–mmc8 (Elsevier CDN), Sato 2022 mmc2–mmc8, Metabolomics Workbench ST003188 (datatable, data, mwtab), Europe PMC supplementary bundle for Wang 2019 (PMC6379049).
- Discovery notes so far (details in Phase 2): PRIDE PXD016999 (Jiang 2020) holds only raw/.msf files plus `ModifiedSampleInfo_v2.xlsx` (40 KB, fetched) — the processed abundances are the Cell supplementary tables (mmc3 123 MB = protein level, mmc4 29 MB = RNA) and the same tables on the GTEx portal eGTEx "Proteomics" fileset. PRIDE PXD010154 (Wang 2019) holds MaxQuant txt bundles: `30healthy_human_tissues_fullproteome_Ensembl_txt.zip` is 28.2 GB (over the 20 GB box) but the server accepts HTTP range requests, so `proteinGroups.txt` will be read out of the zip by ranged reads. Metabolomics Workbench ST003188 "A metabolic atlas of mouse aging" (2025-11-18): 903 samples, 12 organs (Plasma, Bladder, Brain, Heart, Kidney, Liver, Lung, Muscle (Quad), Pancreas, Spleen, Thymus, Tongue; 73 each), 190 named metabolites with RefMet names, triple-quad RP-negative, units AU, both sexes, ages 1–24 months. Sato 2022 is not on Metabolomics Workbench (title search); its supplements are on the Elsevier CDN (mmc2.zip 27 MB, mmc3.xlsx 20 MB, mmc4–8 small).
