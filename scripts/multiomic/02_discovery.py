#!/usr/bin/env python
"""Multiomic Phase 2 — data discovery: what external proteome and metabolome atlases could be reached without a login,
what they contain, and whether they can be matched to MoTrPAC's features.

Reads the download log (data/external_multiomic/_logs/downloads.tsv), the converted tables under data/external_multiomic/
and MoTrPAC's own annotation, and writes results_multiomic/02_discovery/: attempts.csv (every attempt, including the
failures), datasets.csv (one row per usable dataset: URL, contents, samples × features, tissues, units, batch metadata,
matched features), tissue_maps.csv, README.md, REPORT_SECTION.md, STATUS.json. Nothing here is a finding; the phase
records what exists and what does not.
"""
from __future__ import annotations

import json
import re
import time
from pathlib import Path

import numpy as np
import pandas as pd

from tfp import config as C, report
from tfp.transfer import one_to_one_orthologs

ROOT = C.ROOT
EXT = ROOT / "data" / "external_multiomic"
OUT = ROOT / "results_multiomic" / "02_discovery"
RII = ROOT / "results_multiomic" / "01_rii"


def norm(s):
    return re.sub(r"[^a-z0-9]", "", str(s).lower())


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    dl = pd.read_csv(EXT / "_logs" / "downloads.tsv", sep="\t", dtype=str)
    dl["bytes"] = pd.to_numeric(dl["bytes"], errors="coerce")
    dl.to_csv(OUT / "download_log.csv", index=False)

    # ---- attempts: every host/dataset tried, in the prompt's order, with the outcome ------------------------------------
    attempts = [
        {"order": 1, "dataset": "Jiang 2020, Cell 183:269 — Quantitative Proteome Map of the Human Body", "source": "PRIDE PXD016999 (FTP/HTTPS)", "outcome": "raw .raw/.msf only + ModifiedSampleInfo_v2.xlsx (40 KB); no processed table on PRIDE",
         "status": "no processed table"},
        {"order": 1, "dataset": "Jiang 2020", "source": "Elsevier CDN supplementary tables (ars.els-cdn.com, PII S0092867420310783, mmc1–mmc8)", "outcome": "all 8 files downloaded (169 MB in 4 s); mmc3 = protein raw/normalised/relative abundances + tissue medians; mmc2 = TMT run × tag × donor × tissue design; mmc4 = matched RNA log TPM",
         "status": "obtained"},
        {"order": 1, "dataset": "Jiang 2020", "source": "GTEx portal eGTEx 'Proteomics' fileset (gtexportal.org API)", "outcome": "same tables listed (Table_S1 protein 116.9 MB, Table_S2 RNA 28 MB); not downloaded (Elsevier copies used)", "status": "alternative source"},
        {"order": 2, "dataset": "Wang 2019, Mol Syst Biol 15:e8503 — 29 healthy human tissues", "source": "PRIDE PXD010154 (2019/07)", "outcome": "MaxQuant txt bundles: 30healthy_human_tissues_fullproteome_Ensembl_txt.zip = 28.2 GB (> 20 GB box), 9-tissue zip 8.6 GB; server accepts range requests",
         "status": "over the box; ranged read of proteinGroups.txt possible"},
        {"order": 2, "dataset": "Wang 2019", "source": "embopress.org → link.springer.com article page", "outcome": "bot wall (3 KB page, no supplement links)", "status": "blocked"},
        {"order": 2, "dataset": "Wang 2019", "source": "Europe PMC REST supplementaryFiles bundle (PMC6379049)", "outcome": "see datasets.csv / log — bundle download of the EV tables", "status": "see below"},
        {"order": 3, "dataset": "Geiger 2013, MCP 12:1709 — 28 mouse tissues (SILAC mouse)", "source": "Elsevier CDN (PII S1535947620310860, mmc1.zip 9.8 MB)", "outcome": "downloaded; Table S2 = label-free intensities of 7,315 proteins × 28 tissues (one pooled sample per tissue), gene names",
         "status": "obtained"},
        {"order": 3, "dataset": "newer mouse multi-tissue proteome", "source": "PRIDE search API (keyword 'tissue atlas', Mus musculus)", "outcome": "3 hits; PXD082651 (2026, non-canonical peptides, lifespan multi-tissue) not pursued (non-canonical peptide focus)", "status": "not pursued"},
        {"order": 4, "dataset": "rat multi-tissue proteome", "source": "PRIDE search API (Rattus norvegicus; keywords 'tissue', 'atlas organs'), ProteomeCentral PROXI (species 10116)", "outcome": "100 rat projects with 'tissue', none a multi-tissue atlas; 'atlas organs' 0 hits; PROXI 0 datasets", "status": "negative — no rat BodyMap-equivalent proteome found"},
        {"order": 5, "dataset": "Sato 2022, Cell Metab 34:329 — Atlas of exercise metabolism (mouse)", "source": "Metabolomics Workbench (title search over all 4,589 studies)", "outcome": "not deposited there", "status": "absent"},
        {"order": 5, "dataset": "Sato 2022", "source": "Elsevier CDN (PII S1550413121006355, mmc2–mmc8)", "outcome": "downloaded (53 MB); mmc3.xlsx = Metabolon HD4 tables per tissue (BAT, eWAT, HEART, HYPOTHALAMUS, iWAT, LIVER, MUSCLE, SERUM), raw peak areas, animal id, treatment, time after exercise, round/run day",
         "status": "obtained"},
        {"order": 6, "dataset": "Metabolomics Workbench mouse/rat multi-tissue studies", "source": "REST /rest/study/study_id/ST/summary (all studies) filtered to Mus musculus / Rattus norvegicus with tissue words in the title",
         "outcome": "645 candidates; ST003188 'A metabolic atlas of mouse aging' (903 samples, 12 organs, 190 RefMet-named metabolites, 70 mice × 12 organs, both sexes, 5 ages) is the largest multi-organ study with named metabolites; also ST000121 (anesthesia, C57BL/6J tissues, 245 samples, 2015), ST004034 (multi-tissue lipidome, 275), ST003197 (UQ/RQ panel, 104), ST000314 (NMR 4 tissues, 92)",
         "status": "ST003188 obtained (datatable, mwtab); others recorded"},
        {"order": 6, "dataset": "MoTrPAC's own metabolomics on Metabolomics Workbench (ST0026xx)", "source": "REST", "outcome": "same data as the local portal files; not an external check", "status": "not applicable"},
        {"order": 7, "dataset": "MetaboLights", "source": "https://www.ebi.ac.uk/metabolights/ws/studies", "outcome": "reachable (200); not searched further within the 2 h box because two matchable rodent metabolomes were already in hand", "status": "not searched"},
        {"order": 8, "dataset": "MoTrPAC Data Hub (GCS bucket)", "source": "storage.googleapis.com/motrpac-data-hub", "outcome": "403/401 anonymous; not needed — the quant-id tree is local", "status": "not needed"},
    ]
    att = pd.DataFrame(attempts)

    # ---- datasets obtained: shapes, tissues, units, batch metadata, matchability -------------------------------------
    rows, tissue_rows = [], []
    feats = pd.read_csv(RII / "rii_features.csv", index_col=0, dtype=str)
    inner = feats[feats["in_inner"] == "True"]
    orth = one_to_one_orthologs()
    orth_inner = inner.dropna(subset=["ensembl_gene"]).merge(orth, left_on="ensembl_gene", right_on="RAT_ENSEMBL_ID")
    fm = pd.read_csv(C.RAW_DIR / "metab_feature_id_map.csv", dtype=str)
    mo_ref_all = set(fm["metabolite_refmet"].dropna().map(norm)); mo_ref_hil = set(fm.loc[fm["dataset"] == "metab-u-hilicpos", "metabolite_refmet"].dropna().map(norm))
    mo_name = set(fm["metabolite_name"].dropna().map(norm))
    # Jiang 2020
    J = EXT / "jiang2020"
    if (J / "protein_relative.parquet").exists():
        rel = pd.read_parquet(J / "protein_relative.parquet"); cm = pd.read_csv(J / "protein_relative_columns.csv", dtype=str).set_index("column")
        des = pd.read_csv(J / "experimental_design.csv", dtype=str)
        samp = [c for c in rel.columns if c.startswith("GTEX")]
        tis = cm["row2"].reindex(samp)
        jg = set(rel["gene.id"].astype(str))
        m = orth_inner[orth_inner["HUMAN_ORTHOLOG_ENSEMBL_ID"].isin(jg)]
        rna = pd.read_parquet(J / "rna_log_tpm.parquet", columns=["gene.id"]); rc = pd.read_csv(J / "rna_log_tpm_columns.csv", dtype=str)
        rsamp = [c for c in rc["column"] if str(c).startswith("GTEX")]
        key = lambda s: "-".join(s.split("-")[:3])
        matched_rna = len({key(s) for s in samp} & {key(s) for s in rsamp})
        rows.append({"dataset": "Jiang 2020 human proteome map", "layer": "protein (TMT, human)", "url": "https://ars.els-cdn.com/content/image/1-s2.0-S0092867420310783-mmc3.xlsx (+ mmc2 design, mmc4 RNA)",
                     "contents": "Table S2: raw reporter intensities per TMT channel (56 runs × 10 channels incl. one pooled reference per run), normalised abundances, cleaned relative abundances (log2 to run reference); Table S1: run × tag × donor × tissue; Table S3: matched RNA log TPM",
                     "samples": len(samp), "features": int(rel.shape[0]), "individuals": int(des["Individual ID"].dropna().nunique()), "tissues": int(tis.nunique()), "tissue_list": ";".join(sorted(tis.dropna().unique())),
                     "units": "log2 relative abundance to the run's pooled reference (cleaned); raw reporter intensities also given", "batch_metadata": "TMT run (56), TMT tag (10), donor", "matched_to_motrpac": f"{m['RAT_ENSEMBL_ID'].nunique()} of {inner['ensembl_gene'].nunique()} RII genes by 1:1 ortholog",
                     "matched_rna_protein_samples": matched_rna, "usable": "yes — Phases 3, 5, 6"})
        for t, n in tis.value_counts().items():
            tissue_rows.append({"dataset": "Jiang 2020", "external_tissue": t, "n": int(n), "motrpac_tissue": {"Muscle - Skeletal": "SKM-GN (RNA: SKM-GN/SKM-VL)", "Heart - Atrial Appendage": "HEART", "Heart - Left Ventricle": "HEART", "Lung": "LUNG", "Liver": "LIVER", "Brain - Cortex": "CORTEX",
                                                                                                                "Adrenal Gland": "ADRNL (RNA only)", "Colon - Transverse": "COLON (RNA only)", "Small Intestine - Terminal Ileum": "SMLINT (RNA only)", "Spleen": "SPLEEN (RNA only)", "Testis": "TESTES (RNA only)", "Ovary": "OVARY (RNA only)", "Artery - Aorta": "VENACV (RNA only, imperfect)"}.get(t, "OOD")})
    # Wang 2019
    W = EXT / "wang2019"
    wz = W / "PMC6379049_supplementaryFiles.zip"
    wang_status = "not obtained"
    if wz.exists():
        import zipfile
        try:
            z = zipfile.ZipFile(wz); names = [i.filename for i in z.infolist()]
            wang_status = f"Europe PMC bundle readable: {len(names)} files ({'; '.join(names[:8])})"
        except Exception as e:
            wang_status = f"Europe PMC bundle not readable ({wz.stat().st_size / 1e6:.0f} MB, {type(e).__name__}); PRIDE zip 28.2 GB over the box"
    wrow = {"dataset": "Wang 2019 human tissue atlas (29 tissues)", "layer": "protein (label-free iBAQ, human)", "url": "https://ftp.pride.ebi.ac.uk/pride/data/archive/2019/07/PXD010154/ ; https://www.ebi.ac.uk/europepmc/webservices/rest/PMC6379049/supplementaryFiles",
            "contents": "MaxQuant txt bundles on PRIDE (28.2 GB / 8.6 GB); EV tables via Europe PMC", "samples": np.nan, "features": np.nan, "individuals": np.nan, "tissues": 29, "tissue_list": "", "units": "iBAQ", "batch_metadata": "MS run per tissue (one sample per tissue)",
            "matched_to_motrpac": "", "matched_rna_protein_samples": np.nan, "usable": wang_status}
    wp = W / "wang2019_tissue_ibaq.csv"
    if wp.exists():
        wv = pd.read_csv(wp, index_col=0)
        n_w = orth_inner[orth_inner["HUMAN_ORTHOLOG_ENSEMBL_ID"].isin(set(wv.index.astype(str)))]["RAT_ENSEMBL_ID"].nunique()
        wrow.update({"samples": wv.shape[1], "features": wv.shape[0], "tissue_list": ";".join(wv.columns), "individuals": "one donor per tissue (Table EV1 A)",
                     "matched_to_motrpac": f"{n_w} of {inner['ensembl_gene'].nunique()} RII genes by 1:1 ortholog", "matched_rna_protein_samples": int(len(set(wv.columns) & set(pd.read_csv(W / 'wang2019_tissue_rna_fpkm.csv', index_col=0, nrows=1).columns))),
                     "contents": "Table EV1 C: gene-level label-free intensities per tissue; EV1 A: donor sex, age, MS experiment id per tissue; EV2 B: RNA FPKM per tissue (matched)",
                     "usable": "yes — Phase 3 secondary target (one sample per tissue), Phase 5 (matched RNA), Phase 6"})
    rows.append(wrow)
    # Geiger 2013
    Gd = EXT / "geiger2013"
    if (Gd / "tableS1_ratio_HL_normalized.csv").exists():
        gv = pd.read_csv(Gd / "tableS1_ratio_HL_normalized.csv", index_col=0); ga = pd.read_csv(Gd / "tableS1_annotation.csv", index_col=0, dtype=str)
        gv = gv.drop(columns=[c for c in gv.columns if str(c).strip() == "" or str(c).startswith("Unnamed")], errors="ignore")
        g1 = ga["Gene names"].astype(str).str.split(";").str[0].str.strip().str.upper()
        n_match = len(set(g1) & set(inner["gene_symbol"].dropna().str.upper()))
        rows.append({"dataset": "Geiger 2013 mouse tissue proteome (28 tissues)", "layer": "protein (SILAC H/L ratios to one SILAC-mouse standard; mouse)", "url": "https://ars.els-cdn.com/content/image/1-s2.0-S1535947620310860-mmc1.zip",
                     "contents": "Table S1: MaxQuant normalised H/L ratios per tissue against the common SILAC-mouse standard (one pooled sample per tissue), gene names; Table S2: top-100 intensities", "samples": gv.shape[1], "features": gv.shape[0], "individuals": "pooled mice per tissue", "tissues": gv.shape[1], "tissue_list": ";".join(gv.columns),
                     "units": "normalised H/L SILAC ratio (common heavy standard)", "batch_metadata": "none per sample (one run per tissue)", "matched_to_motrpac": f"{n_match} RII genes by gene symbol (rat↔mouse, caveat: symbol match, not orthology table)",
                     "matched_rna_protein_samples": np.nan, "usable": "yes — Phase 3 secondary target (covers KIDNEY and WAT-SC, n = 1 per tissue)"})
        for t in gv.columns:
            tissue_rows.append({"dataset": "Geiger 2013", "external_tissue": t, "n": 1, "motrpac_tissue": {"Brain cortex": "CORTEX", "Heart": "HEART", "Kidney cortex": "KIDNEY", "Liver": "LIVER", "Lung": "LUNG", "Muscle": "SKM-GN", "White fat": "WAT-SC", "Brown fat": "BAT (RNA only)", "Adrenal gland": "ADRNL (RNA only)", "Spleeen": "SPLEEN (RNA only)", "Ovary": "OVARY (RNA only)", "Colon": "COLON (RNA only)"}.get(t, "OOD")})
    # Sato 2022
    Sd = EXT / "sato2022"
    if (Sd / "sheets_summary.csv").exists():
        ss = pd.read_csv(Sd / "sheets_summary.csv")
        names = set(pd.read_csv(Sd / "named_biochemicals_union.csv", dtype=str)["name"].dropna())
        rm_parts = sorted(Sd.glob("_refmet_shard*.csv")) + ([Sd / "refmet_match.csv"] if (Sd / "refmet_match.csv").exists() else [])
        rm = pd.concat([pd.read_csv(p, dtype=str) for p in rm_parts], ignore_index=True).drop_duplicates("query") if rm_parts else pd.DataFrame(columns=["query", "refmet_name"])
        sato_ref = set(rm["refmet_name"].dropna().map(norm))
        crude = {norm(n) for n in names}
        rows.append({"dataset": "Sato 2022 atlas of exercise metabolism (mouse)", "layer": "metabolite (Metabolon HD4 untargeted; mouse)", "url": "https://ars.els-cdn.com/content/image/1-s2.0-S1550413121006355-mmc3.xlsx",
                     "contents": "per-tissue Metabolon tables: raw peak areas (OrigScale) and normalised, sample metadata (animal id, treatment sedentary/exercise, time after exercise, round, run day), chemical annotation (KEGG, HMDB, PubChem, CAS)",
                     "samples": int(ss["n_samples"].sum()), "features": f"{int(ss['n_named'].min())}–{int(ss['n_named'].max())} named per tissue ({len(names)} union)", "individuals": int(ss["n_animals"].max()), "tissues": int(len(ss)), "tissue_list": ";".join(ss["tissue_sheet"]),
                     "units": "peak area (raw) / median-scaled", "batch_metadata": "ROUND, RUN DAY (per sample), platform", "matched_to_motrpac": f"RefMet-matched names: {len(rm)} of {len(names)} queried, {len(sato_ref & mo_ref_all)} in MoTrPAC (any platform), {len(sato_ref & mo_ref_hil)} in HILIC+; crude name overlap {len(crude & (mo_ref_all | mo_name))}",
                     "matched_rna_protein_samples": np.nan, "usable": "yes — Phase 4 (target + exercise-invariance test), Phase 6"})
        for _, r in ss.iterrows():
            tissue_rows.append({"dataset": "Sato 2022", "external_tissue": r["tissue_sheet"], "n": int(r["n_samples"]), "motrpac_tissue": {"BAT": "BAT", "eWAT": "WAT-SC (visceral vs subcutaneous: imperfect)", "iWAT": "WAT-SC", "HEART": "HEART", "HYPOTHALAMUS": "HYPOTH", "LIVER": "LIVER", "MUSCLE": "SKM-GN/SKM-VL", "SERUM": "PLASMA (serum vs plasma)"}.get(r["tissue_sheet"], "OOD")})
    # Metabolomics Workbench ST003188
    Md = EXT / "mw_ST003188"
    if (Md / "samples.csv").exists():
        sm = pd.read_csv(Md / "samples.csv", dtype=str); mets = pd.read_csv(Md / "metabolites.csv", dtype=str)
        study = sm[sm["mouse"].notna() & ~sm["source"].isin(["Pool", "Blank"])]
        mw_ref = set(mets["refmet_name"].dropna().map(norm))
        rows.append({"dataset": "Metabolomics Workbench ST003188 — A metabolic atlas of mouse aging (Mullen lab, USC, 2025)", "layer": "metabolite (targeted RP-negative triple-quad; mouse)", "url": "https://www.metabolomicsworkbench.org/rest/study/analysis_id/AN005236/datatable",
                     "contents": "peak areas of 190 RefMet-named metabolites; factors organ, sex, age; mwtab gives Batch per sample", "samples": int(len(study)), "features": int(mets.shape[0]), "individuals": int(study["mouse"].nunique()), "tissues": int(study["organ"].nunique()), "tissue_list": ";".join(sorted(study["organ"].unique())),
                     "units": "AU (peak area)", "batch_metadata": "Batch (12 levels — one per organ, nested)", "matched_to_motrpac": f"{len(mw_ref & mo_ref_all)} RefMet names in MoTrPAC (any platform), {len(mw_ref & mo_ref_hil)} in HILIC+",
                     "matched_rna_protein_samples": np.nan, "usable": "yes — Phase 4 (second target), Phase 6 (nested design)"})
        for t, n in study["organ"].value_counts().items():
            tissue_rows.append({"dataset": "MW ST003188", "external_tissue": t, "n": int(n), "motrpac_tissue": {"Plasma": "PLASMA", "Brain": "CORTEX/HIPPOC/HYPOTH", "Heart": "HEART", "Kidney": "KIDNEY", "Liver": "LIVER", "Lung": "LUNG", "Muscle (Quad)": "SKM-GN/SKM-VL", "Spleen": "SPLEEN"}.get(t, "OOD")})
    ds = pd.DataFrame(rows)
    ds.to_csv(OUT / "datasets.csv", index=False)
    att.to_csv(OUT / "attempts.csv", index=False)
    tm = pd.DataFrame(tissue_rows)
    tm.to_csv(OUT / "tissue_maps.csv", index=False)
    total_mb = dl["bytes"].sum() / 1e6
    lines = [f"# Phase 2 — data discovery", "", f"Built by `scripts/multiomic/02_discovery.py` on {time.strftime('%Y-%m-%d %H:%M UTC', time.gmtime())}. Every attempt is in `attempts.csv`; every download (URL, bytes, sha256, time) in `download_log.csv` "
             f"({len(dl)} files, {total_mb:.0f} MB); usable datasets in `datasets.csv`; tissue maps in `tissue_maps.csv`.", "",
             "## Attempts (prompt order)", "", report.df_to_md(att[["order", "dataset", "source", "status"]]), "",
             "## Datasets obtained", "", report.df_to_md(ds[["dataset", "layer", "samples", "features", "individuals", "tissues", "units", "batch_metadata", "matched_to_motrpac", "usable"]]), "",
             "## Tissue maps", "", report.df_to_md(tm, max_rows=120), "",
             "## Negative results", "",
             "- No rat multi-tissue proteome with processed tables was found on PRIDE or ProteomeXchange (keyword searches); the BodyMap-equivalent proteome does not exist in the public repositories searched.",
             "- Sato 2022 is not on Metabolomics Workbench; its supplementary Metabolon tables were used instead.",
             "- Wang 2019: the processed MaxQuant bundle on PRIDE is 28.2 GB, over the 20 GB box; the Europe PMC supplementary bundle is the route (see `datasets.csv` for its state).",
             "- Kidney and adipose are absent from Jiang 2020, so two of the seven proteomics tissues have no human protein target; Geiger 2013 (mouse, one pooled sample per tissue) covers both.", ""]
    (OUT / "README.md").write_text("\n".join(lines))
    sec = ["- question · which external proteome and metabolome atlases are reachable without a login, and can they be matched to MoTrPAC's features?",
           f"- data · {len(dl)} files, {total_mb:.0f} MB downloaded (`results_multiomic/02_discovery/download_log.csv`: URL, bytes, sha256, time); {len(att)} attempts logged (`attempts.csv`).",
           "- design · prompt order: Jiang 2020 → Wang 2019 → mouse atlas → rat atlas search → Sato 2022 → Metabolomics Workbench → MetaboLights; 20 min / 20 GB per attempt; processed tables only.",
           "- result · " + "; ".join(f"**{r.dataset.split(' (')[0].split(' — ')[0]}** ({r.layer}; {r.samples} samples, {r.tissues} tissues; {r.matched_to_motrpac})" for r in ds.itertuples()) + ". "
           "Negative: no rat multi-tissue proteome exists in PRIDE/ProteomeXchange; Sato 2022 is not on Metabolomics Workbench; Wang 2019's PRIDE bundle (28 GB) is over the box (`datasets.csv`, `attempts.csv`).",
           "- what it does not show · nothing about transfer; the tissue maps (`tissue_maps.csv`) carry imperfect matches (atrial appendage → HEART, serum → PLASMA, epididymal fat → WAT-SC) that the later phases inherit.",
           "- files · `results_multiomic/02_discovery/README.md`."]
    (OUT / "REPORT_SECTION.md").write_text("\n".join(sec) + "\n")
    (OUT / "STATUS.json").write_text(json.dumps({"status": "DONE", "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}))
    print(ds[["dataset", "samples", "tissues", "matched_to_motrpac", "usable"]].to_string(index=False))
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
