#!/usr/bin/env python
"""Phase 16 — identifiability: is the tissue axis separable from processing inside this study?

Three recomputations from metadata alone (no expression values):
  a) Batch nesting per omic layer — for every processing variable in the export (extraction plate,
     library batch, flowcell, lane, dates, site; TMT plex and channel; immunoassay plate), how many
     levels a tissue spans, how many tissues a level holds, Cramér's V with tissue, and how many of
     the tissue pairs share at least one level (`tfp.batch.nesting_table`).
  b) Estimable tissue pairs per layer — pairs that share a level of EVERY processing variable of the
     layer, so that a within-batch contrast exists (`tfp.batch.estimable_pairs`).
  c) The QC-only baseline of notebook 02 §13 (`notebooks/_build/sections/s14_qc_baseline.py`), copied
     verbatim in its feature lists and model: a multinomial logistic regression on the consortium's
     RNA-seq library QC numbers, no gene, on the phase-04 animal-grouped folds, technical and
     composition feature sets kept separate.

Layers without a processing variable in the export (metabolomics: platforms only) are recorded as
`unavailable` — never guessed.

Outputs (results/16_identifiability/): nesting_<ASSAY>.csv, estimable_pairs.csv, batch_counts.csv,
qc_only_per_fold.csv, qc_only_summary.csv, layers.json, NOTES.md
"""
from __future__ import annotations

import json
import re
import warnings

import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, balanced_accuracy_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from tfp import batch, cli, config as C, io
from tfp.splits import assert_no_group_leak, grouped_kfold

# free-text tissue names used in the metadata tables → pipeline tissue codes
_TISSUE_ALIASES = {
    "heart": "HEART", "liver": "LIVER", "whiteadipose": "WAT-SC", "gastrocnemius": "SKM-GN", "brownadipose": "BAT",
    "kidney": "KIDNEY", "hippocampus": "HIPPOC", "lung": "LUNG", "cortex": "CORTEX", "hypothalmus": "HYPOTH",
    "hypothalamus": "HYPOTH", "vastuslateralis": "SKM-VL", "spleen": "SPLEEN", "adrenal": "ADRNL",
    "paxgenerna": "BLOOD", "blood": "BLOOD", "colon": "COLON", "smallintestine": "SMLINT", "aorta": "VENACV",
    "venacava": "VENACV", "testes": "TESTES", "ovaries": "OVARY", "ovary": "OVARY", "plasma": "PLASMA",
}


def tissue_code(text: str) -> str | None:
    key = re.sub(r"[^a-z]", "", str(text).lower().replace("rat-", "").replace("powder", ""))
    return _TISSUE_ALIASES.get(key)


# processing variables per layer, and the subset that defines an estimable pair (all must be shared)
LAYERS = {
    "TRNSCRPT": {"file": "TRNSCRPT.csv", "tissue": "Tissue",
                 "variables": ["GET_site", "RNA_extr_plate_ID", "RNA_extr_date", "Lib_prep_date", "Lib_batch_ID",
                               "Seq_date", "Seq_flowcell_ID", "Seq_flowcell_lane", "Seq_batch"],
                 "estimable": ["RNA_extr_plate_ID", "Lib_batch_ID", "Seq_flowcell_ID"]},
    "METHYL": {"file": "METHYL.csv", "tissue": "Tissue",
               "variables": ["GET_site", "DNA_extr_plate_ID", "DNA_extr_date", "Lib_prep_date", "Lib_batch_ID",
                             "Seq_date", "Seq_flowcell_ID", "Seq_flowcell_lane", "Seq_batch"],
               "estimable": ["DNA_extr_plate_ID", "Lib_batch_ID", "Seq_flowcell_ID"]},
    "ATAC": {"file": "ATAC.csv", "tissue": "general.description",
             "variables": ["GET_site", "Sample_batch", "Nuclei_extr_date", "Tagmentation_date", "PCR_date", "Seq_date",
                           "Seq_flowcell_ID", "Seq_flowcell_lane", "Seq_batch"],
             "estimable": ["Nuclei_extr_date", "Tagmentation_date", "PCR_date", "Seq_flowcell_ID"]},
    # TMT: the plex labels S1–S6 repeat in every tissue, but a plex is 11 channels = 10 samples + that tissue's
    # reference pool, so a physical plex never holds two tissues; the plex identity is (tissue, label).
    "PROT": {"file": "PROT.csv", "tissue": "tissue", "variables": ["plex_id", "tmt11_channel"], "estimable": ["plex_id"]},
    "PHOSPHO": {"file": "PHOSPHO.csv", "tissue": "tissue", "variables": ["plex_id", "tmt11_channel"], "estimable": ["plex_id"]},
    "ACETYL": {"file": "ACETYL.csv", "tissue": "tissue", "variables": ["plex_id", "tmt11_channel"], "estimable": ["plex_id"]},
    "UBIQ": {"file": "UBIQ.csv", "tissue": "tissue", "variables": ["plex_id", "tmt11_channel"], "estimable": ["plex_id"]},
    "IMMUNO": {"file": "IMMUNO.csv", "tissue": "tissue", "variables": ["plate_id", "panel_name"], "estimable": ["plate_id"]},
}
UNAVAILABLE = {"METAB": "the export carries the platform (dataset) per vial but no run, plate or batch variable"}

# ---- QC-only baseline: feature lists copied verbatim from notebooks/_build/sections/s14_qc_baseline.py ----
QC_TECH = ["RIN", "r_260_280", "r_260_230", "pct_adapter_detected", "pct_trimmed", "pct_GC", "pct_dup_sequence",
           "pct_umi_dup", "pct_multimapped", "median_5_3_bias", "reads_log10", "avg_input_read_length"]
QC_COMP = ["pct_rRNA", "pct_globin", "pct_chrM", "pct_chrX", "pct_chrY", "pct_mrna", "pct_coding", "pct_utr",
           "pct_intronic", "pct_intergenic"]
QC_SETS = {"technical": QC_TECH, "composition": QC_COMP, "all": QC_TECH + QC_COMP}
QC_TABLE = C.ROOT.parents[1] / "data/quant-id/rat-training-06/c1.0/transcriptomics/qa-qc/motrpac_pass1b-06_transcript-rna-seq_qa-qc-metrics.csv"


def qc_model():
    return Pipeline([("imp", SimpleImputer(strategy="median")), ("sc", StandardScaler()),
                     ("clf", LogisticRegression(C=1.0, max_iter=10000))])


def load_layer(name: str, spec: dict, study_only: bool = True) -> pd.DataFrame | None:
    p = C.META_DIR / spec["file"]
    if not p.exists():
        return None
    m = pd.read_csv(p, dtype=str, low_memory=False)
    if "viallabel" in m.columns and study_only:
        m = m[m["viallabel"].astype(str).str.startswith("9")]
    if name == "ATAC":
        m["tissue_code"] = m[spec["tissue"]].str.split("_").str[0].map(tissue_code)
    else:
        m["tissue_code"] = m[spec["tissue"]].map(lambda t: tissue_code(t) or (t if t in C.TISSUES else None))
    unmapped = m.loc[m["tissue_code"].isna(), spec["tissue"]].unique().tolist()
    if unmapped:
        raise ValueError(f"{name}: unmapped tissue names {unmapped}")
    if "tmt_plex" in m.columns:
        m["plex_id"] = m["tissue_code"] + ":" + m["tmt_plex"].astype(str)
    return m


def qc_only_baseline(out) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """Notebook 02 §13 (s14) reproduced: same matrix (phase-04 stacked counts), same folds, same model."""
    pheno = io.load_pheno()
    om = io.stack_tissues("TRNSCRPT", source=cli.resolve_source("TRNSCRPT", "auto"), join="inner", pheno=pheno, verbose=False)
    folds = list(grouped_kfold(om.meta, "tissue", 5, C.SEED))
    if QC_TABLE.exists():
        q = pd.read_csv(QC_TABLE, dtype=str)
        q = q[q["vial_label"].str.startswith("9")].drop_duplicates("vial_label").set_index("vial_label")
        source = str(QC_TABLE.relative_to(C.ROOT.parents[1]))
    else:  # same columns, exported from the R package metadata
        q = pd.read_csv(C.META_DIR / "TRNSCRPT.csv", dtype=str, low_memory=False).drop_duplicates("viallabel").set_index("viallabel")
        source = "data/raw/meta/TRNSCRPT.csv (consortium QC table not found)"
    idx = om.meta.index.astype(str)
    missing = int((~idx.isin(q.index)).sum())
    if missing:
        raise ValueError(f"{missing} study vials have no QC row")
    qc = q.loc[idx].copy()
    qc["reads_log10"] = np.log10(qc["reads"].astype(float))
    y = om.meta["tissue"].astype(str).to_numpy()
    rows = []
    for name, cols in QC_SETS.items():
        X = qc[cols].apply(pd.to_numeric, errors="coerce").to_numpy(float)
        for f, (tr, te) in enumerate(folds):
            assert_no_group_leak(om.meta, tr, te)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                m = qc_model().fit(X[tr], y[tr])
            p = m.predict(X[te])
            rows.append({"features": name, "n_features": len(cols), "fold": f, "n_test_vials": len(te),
                         "n_test_animals": int(om.meta.iloc[te]["pid"].nunique()),
                         "accuracy": accuracy_score(y[te], p), "balanced_accuracy": balanced_accuracy_score(y[te], p)})
    pf = pd.DataFrame(rows)
    summ = pf.groupby("features", sort=False).agg(n_features=("n_features", "first"), acc_mean=("accuracy", "mean"),
                                                   acc_sd=("accuracy", "std"), bal_acc_mean=("balanced_accuracy", "mean"),
                                                   bal_acc_sd=("balanced_accuracy", "std"), n_folds=("fold", "nunique"),
                                                   n_test_animals_mean=("n_test_animals", "mean")).reset_index()
    info = {"qc_source": source, "n_vials": int(om.n_samples), "n_animals": int(om.n_animals),
            "n_tissues": int(om.meta["tissue"].nunique()), "chance": 1 / om.meta["tissue"].nunique(),
            "model": "LogisticRegression(C=1.0, multinomial) after median imputation and standardisation, fit inside each fold",
            "folds": "grouped_kfold(meta, 'tissue', 5, SEED): the phase-04 folds"}
    return pf, summ, info


PORTAL_RNA = C.ROOT.parents[1] / "data/quant-id/rat-training-06/c1.0/transcriptomics"


def bridge_variance(out) -> tuple[pd.DataFrame, pd.DataFrame, dict] | None:
    """Batch measured directly on the bridging standards.

    The consortium sequenced reference-standard RNA pools (vial labels starting with 8; `Sample_category` = ref in
    the QC table) alongside the study vials. A pool is identified by its bid (label[:5]) and type (label[5:7], 99 or
    88); a pool that was run on several extraction plates is the same RNA measured in several batches. For each such
    pool and gene, the between-plate variance of its log2 CPM (n − 1 denominator over its vials, one per plate) is
    the batch variance V_batch. The tissue-separating variance V_tissue is the variance of the 19 tissue means of
    study-vial log2 CPM (the pipeline's stacked matrix). The ratio V_batch / V_tissue per gene, summarised over gene
    sets as the median and as the ratio of sums, is what "batch as a fraction of the variance that separates tissues"
    means here. Study vials and reference vials use the same unit: log2(CPM + 1) on the total library over all genes.
    Reads the portal per-tissue RSEM count files; returns None when they are absent."""
    files = sorted(PORTAL_RNA.glob("t*/transcript-rna-seq/*rsem-genes-count.txt"))
    if not files:
        return None
    cols, frames = [], []
    for f in files:
        hdr = open(f).readline().rstrip("\n").split("\t")
        refs = [c for c in hdr[1:] if c.startswith("8")]
        if not refs:
            continue
        df = pd.read_csv(f, sep="\t", usecols=["gene_id"] + refs, index_col="gene_id")
        frames.append(df)
        cols += refs
    counts = pd.concat(frames, axis=1)
    lcpm = np.log2(counts.div(counts.sum(axis=0), axis=1) * 1e6 + 1.0)   # genes × ref vials
    qc = pd.read_csv(QC_TABLE, dtype=str) if QC_TABLE.exists() else None
    plate = qc.set_index("vial_label")["RNA_extr_plate_ID"].to_dict() if qc is not None else {}
    site = qc.set_index("vial_label")["GET_site"].to_dict() if qc is not None else {}
    pool_tissue = qc.set_index("vial_label")["Tissue"].to_dict() if qc is not None else {}
    # study vials: the pipeline's stacked log2 CPM matrix (21,193 genes) and its tissue means
    pheno = io.load_pheno()
    om = io.stack_tissues("TRNSCRPT", source=cli.resolve_source("TRNSCRPT", "auto"), join="inner", pheno=pheno, verbose=False)
    genes = [g for g in om.X.columns if g in lcpm.index]
    X = om.X[genes]
    tissue_means = X.groupby(om.meta["tissue"].to_numpy()).mean()
    v_tissue = tissue_means.var(axis=0, ddof=1)                          # variance of the 19 tissue means
    v_within = X.groupby(om.meta["tissue"].to_numpy()).var(ddof=1).mean(axis=0)  # mean within-tissue variance
    sets = {"all_genes": genes}
    d05 = C.RESULTS_DIR / "05_panels" / "TRNSCRPT"
    d12 = C.RESULTS_DIR / "12_bodymap"
    if (d12 / "panel_gene_check.csv").exists():
        sets["panel_k20"] = [g for g in pd.read_csv(d12 / "panel_gene_check.csv")["feature_ID"] if g in set(genes)]
    if (d05 / "candidate_panel_annotated.csv").exists():
        sets["stable_core"] = [g for g in pd.read_csv(d05 / "candidate_panel_annotated.csv")["feature_ID"] if g in set(genes)]
    if (d05 / "stability_k20_annotated.csv").exists():
        sets["stability_51"] = [g for g in pd.read_csv(d05 / "stability_k20_annotated.csv")["feature_ID"] if g in set(genes)]
    pools = {}
    for v in lcpm.columns:
        pools.setdefault((v[:5], v[5:7]), []).append(v)
    rows, per_gene = [], []
    for (bid, ptype), vials in sorted(pools.items()):
        plates = sorted({plate.get(v, "?") for v in vials})
        if len(vials) < 2 or len(plates) < 2:
            continue
        vb = lcpm.loc[genes, vials].var(axis=1, ddof=1)
        # a pool can only measure batch on the genes it expresses: restrict to mean log2 CPM ≥ 1 in the pool
        expressed = set(lcpm.loc[genes, vials].mean(axis=1).loc[lambda s: s >= 1.0].index)
        pool_sets = dict(sets)
        pool_sets["all_genes_expressed_in_pool"] = [g for g in genes if g in expressed]
        if "panel_k20" in sets:
            pool_sets["panel_k20_expressed_in_pool"] = [g for g in sets["panel_k20"] if g in expressed]
        for name, gs in pool_sets.items():
            if not gs:
                continue
            r = (vb[gs] / v_tissue[gs]).replace([np.inf, -np.inf], np.nan)
            rows.append({"pool_bid": bid, "pool_type": ptype, "pool_tissue": pool_tissue.get(vials[0], "?"), "n_vials": len(vials), "n_plates": len(plates),
                         "sites": ";".join(sorted({site.get(v, "?") for v in vials})), "gene_set": name, "n_genes": len(gs),
                         "median_ratio_batch_over_tissue": float(np.nanmedian(r)), "sum_ratio_batch_over_tissue": float(vb[gs].sum() / v_tissue[gs].sum()),
                         "median_v_batch": float(vb[gs].median()), "median_v_tissue": float(v_tissue[gs].median()), "median_v_within_tissue": float(v_within[gs].median()),
                         "median_ratio_batch_over_within": float(np.nanmedian((vb[gs] / v_within[gs]).replace([np.inf, -np.inf], np.nan)))})
        for g in sets.get("panel_k20", []):
            per_gene.append({"pool_bid": bid, "pool_type": ptype, "feature_ID": g, "v_batch": float(vb[g]), "v_tissue": float(v_tissue[g]), "v_within_tissue": float(v_within[g]),
                             "ratio_batch_over_tissue": float(vb[g] / v_tissue[g]) if v_tissue[g] > 0 else np.nan, "n_plates": len(plates),
                             "mean_log2cpm_in_pool": float(lcpm.loc[g, vials].mean()), "expressed_in_pool": g in expressed})
    summary, pg = pd.DataFrame(rows), pd.DataFrame(per_gene)
    if len(pg):
        sym = io.map_to_gene_symbols(pg["feature_ID"].unique())
        pg["gene_symbol"] = pg["feature_ID"].map(sym)
    summary.to_csv(out / "bridge_variance.csv", index=False)
    pg.to_csv(out / "bridge_variance_per_gene.csv", index=False)
    info = {"n_ref_vials": int(lcpm.shape[1]), "n_pools": len(pools), "n_bridging_pools": int(summary[["pool_bid", "pool_type"]].drop_duplicates().shape[0]) if len(summary) else 0,
            "n_genes_shared": len(genes), "source": str(PORTAL_RNA.relative_to(C.ROOT.parents[1])) + "/t*/transcript-rna-seq/*rsem-genes-count.txt (reference vials) + the pipeline's stacked study matrix",
            "definition": "V_batch = variance of a pool's log2 CPM across the plates it was run on (one vial per plate); V_tissue = variance of the 19 tissue means of study-vial log2 CPM; ratio per gene, median and ratio of sums over gene sets"}
    return summary, pg, info


def main() -> None:
    ap = cli.common_parser("Identifiability: batch nesting, estimable pairs, QC-only baseline")
    ap.add_argument("--skip-qc", action="store_true", help="skip the QC-only baseline (metadata tables only)")
    ap.add_argument("--bridge", action="store_true", help="also measure batch directly on the bridging reference pools (needs the portal count files)")
    args = ap.parse_args()
    cli.banner("16_identifiability", args)
    out = cli.outdir("16_identifiability", args.out)

    layers, est_rows, count_rows = {}, [], []
    for name, spec in LAYERS.items():
        m = load_layer(name, spec)
        if m is None:
            layers[name] = {"status": "unavailable", "reason": f"data/raw/meta/{spec['file']} not present"}
            continue
        variables = [v for v in spec["variables"] if v in m.columns]
        tab = batch.nesting_table(m, "tissue_code", variables)
        tab.insert(0, "assay", name)
        tab.to_csv(out / f"nesting_{name}.csv", index=False)
        est_vars = [v for v in spec["estimable"] if v in m.columns]
        pairs, total = batch.estimable_pairs(m, "tissue_code", est_vars)
        est_rows.append({"assay": name, "n_tissues": int(m["tissue_code"].nunique()), "n_samples": int(len(m)),
                         "n_pairs_total": total, "n_pairs_estimable": len(pairs),
                         "estimable_pairs": ";".join(f"{a}|{b}" for a, b in pairs), "variables_used": ";".join(est_vars)})
        for v in variables:
            count_rows.append({"assay": name, "variable": v, "n_levels": int(m[v].dropna().nunique())})
        layers[name] = {"status": "recomputed", "n_samples": int(len(m)), "n_tissues": int(m["tissue_code"].nunique()),
                        "variables": variables, "estimable_variables": est_vars, "n_pairs_total": total,
                        "n_pairs_estimable": len(pairs), "estimable_pairs": pairs, "source": f"data/raw/meta/{spec['file']}"}
        print(f"  {name}: {len(m)} samples, {m['tissue_code'].nunique()} tissues; estimable pairs {len(pairs)}/{total} {pairs}")
    for name, reason in UNAVAILABLE.items():
        layers[name] = {"status": "unavailable", "reason": reason}
    pd.DataFrame(est_rows).to_csv(out / "estimable_pairs.csv", index=False)
    pd.DataFrame(count_rows).to_csv(out / "batch_counts.csv", index=False)

    qc_info = None
    if not args.skip_qc:
        pf, summ, qc_info = qc_only_baseline(out)
        pf.to_csv(out / "qc_only_per_fold.csv", index=False)
        summ.to_csv(out / "qc_only_summary.csv", index=False)
        print(summ.to_string(index=False, float_format="%.4f"))
    bridge_info = None
    if args.bridge:
        b = bridge_variance(out)
        if b is None:
            bridge_info = {"status": "unavailable", "reason": "portal RNA-seq count files (data/quant-id/.../transcript-rna-seq) not present"}
            print("  bridge: portal count files not found; skipped")
        else:
            summary, pg, bridge_info = b
            bridge_info["status"] = "recomputed"
            print(summary[summary["gene_set"].isin(["panel_k20", "all_genes"])].to_string(index=False, float_format="%.4f"))
    (out / "layers.json").write_text(json.dumps({"layers": layers, "qc_only": qc_info, "bridge": bridge_info}, indent=1))

    notes = [
        "# Phase 16 — identifiability recompute\n",
        "Metadata only; no expression value is used. Study vials only (viallabel starts with 9).\n",
        "**Estimable pair.** Two tissues are estimable when they share a level of every processing variable listed in "
        "`variables_used` (estimable_pairs.csv), so that a contrast between them exists inside one batch. "
        "Everything else in nesting_<ASSAY>.csv is descriptive: levels per tissue, tissues per level, shared levels, "
        "Cramér's V between tissue and the variable, and the number of tissue pairs sharing at least one level.\n",
        "**QC-only baseline.** Copied from notebooks/_build/sections/s14_qc_baseline.py: multinomial logistic regression "
        "(C = 1) on the consortium's per-library QC numbers, no gene, median imputation and standardisation fit inside "
        "the fold, on the phase-04 animal-grouped folds. `technical` and `composition` are kept separate because the "
        "composition fractions (mitochondrial, globin, rRNA, intronic reads, chrX/chrY) are read biologically by MoTrPAC "
        "itself and are not pure processing.\n",
        "**TMT layers.** The plex labels S1–S6 repeat in every tissue; a plex is 11 channels (10 samples + the tissue's "
        "reference pool), so the physical plex is (tissue, label) and is nested in tissue by construction. The channel "
        "is shared across everything and carries no tissue information; in HEART and LIVER it encodes sex (phase 03).\n",
        "**Immunoassays.** The Luminex plates are the one layer whose plates hold several tissues (the plate names list "
        "the tissue codes), so tissue pairs on a shared plate are estimable there.\n",
        "Layers: " + ", ".join(f"{k} ({v['status']})" for k, v in layers.items()) + ".\n",
        ("**Bridging standards (--bridge).** " + bridge_info["definition"] + f" {bridge_info['n_ref_vials']} reference vials, "
         f"{bridge_info['n_bridging_pools']} pools run on more than one plate (bridge_variance.csv, bridge_variance_per_gene.csv).\n")
        if bridge_info and bridge_info.get("status") == "recomputed" else "",
    ]
    (out / "NOTES.md").write_text("\n".join(notes))
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
