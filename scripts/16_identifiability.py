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


def main() -> None:
    ap = cli.common_parser("Identifiability: batch nesting, estimable pairs, QC-only baseline")
    ap.add_argument("--skip-qc", action="store_true", help="skip the QC-only baseline (metadata tables only)")
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
    (out / "layers.json").write_text(json.dumps({"layers": layers, "qc_only": qc_info}, indent=1))

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
    ]
    (out / "NOTES.md").write_text("\n".join(notes))
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
