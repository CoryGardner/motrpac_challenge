#!/usr/bin/env python
"""Phase 05c — training response of the panel genes in their marker tissues (consortium DA tables).

Question: does endurance training move the tissue-marker genes enough to matter for a tissue fingerprint?
Per gene, the largest training log2 fold change in its marker tissue (consortium DESeq2, sex × duration vs the
sex-matched sedentary controls, `data/raw/da/TRNSCRPT__<TISSUE>.csv`) is set against the tissue effect the panel
rests on: the marker-tissue mean log2 CPM minus the highest mean among the other tissues (`effect_size` and
`next_highest_tissue` of the phase-05 annotation, computed on the stacked counts matrix).

Gene sets. `k20` = the 20-gene panel of the all-animal fit (results/12_bodymap/panel_survival.csv, row
"k=20 (all-animal fit)"; symbols mapped to feature IDs as scripts/30_export_site_data.py does, with the annotated
tables taking precedence over feature_to_gene). `core` = results/05_panels/TRNSCRPT/candidate_panel_annotated.csv
(the 10-gene stable core). `stability_51` = results/05_panels/TRNSCRPT/stability_k20_annotated.csv (every gene
selected in at least one of the 50 k = 20 fits). Marker tissue and tissue effect come from the annotated files.
Membership in the consortium 5 % FDR training-regulated set (data/raw/training_regulated_features.csv, assay
TRNSCRPT) is looked up in the marker tissue and counted over all tissues; per tissue it is checked against
`selection_fdr < 0.05` in the DA table.

Caveats. The DA n is small (the TRNSCRPT sample metadata holds 4–5 vials per tissue × sex × group; the number is
recomputed and written to the background table), so a null log2FC is weak evidence of no response. The tissue
effect is a difference of tissue means, not a per-animal estimate, and is not on the same footing as a DESeq2
coefficient; the ratio of the two is descriptive. Nothing here changes any panel.

Outputs (results/05_panels/TRNSCRPT/): panel_training_response.csv (one row per gene), panel_training_background.csv
(one row per marker tissue), panel_training_summary.csv (one row per gene set), panel_training_NOTES.md.
Under 30 s; `--quick` changes nothing.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from tfp import cli, config as C, io, report

DA_COLS = ["feature_ID", "sex", "comparison_group", "logFC", "logFC_se", "adj_p_value", "selection_fdr"]
CONTRASTS = [(s, g) for s in ("female", "male") for g in ("1w", "2w", "4w", "8w")]


def read_da(tissue: str) -> pd.DataFrame:
    """The columns of `io.load_da("TRNSCRPT", tissue)` this phase needs (same file; usecols keeps it fast)."""
    p = C.DA_DIR / f"TRNSCRPT__{C.tissue_token(tissue)}.csv"
    da = pd.read_csv(p, usecols=DA_COLS, dtype={"feature_ID": str}, low_memory=False)
    for c in ("sex", "comparison_group"):
        da[c] = da[c].astype(str).str.lower().str.strip()
    for c in ("logFC", "logFC_se", "adj_p_value", "selection_fdr"):
        da[c] = pd.to_numeric(da[c], errors="coerce")
    return da


def k20_ids(panels_dir: Path, results: Path, stab: pd.DataFrame) -> list[str]:
    """The site's k = 20 panel: symbols of the all-animal fit → feature IDs (scripts/30_export_site_data.py, gene_set)."""
    surv = pd.read_csv(results / "12_bodymap" / "panel_survival.csv").set_index("panel")
    symbols = surv.loc["k=20 (all-animal fit)", "genes"].split(";")
    f2g = io.load_feature_to_gene()
    f2g = f2g.dropna(subset=["feature_ID"])
    f2g = f2g[f2g["feature_ID"].astype(str).str.startswith("ENSRNOG")]
    symcol = "gene_symbol" if "gene_symbol" in f2g.columns else "symbol"
    sym_to_id: dict[str, str] = {}
    for fid, s in zip(f2g["feature_ID"], f2g[symcol]):          # first mapping wins
        if isinstance(s, str) and s not in sym_to_id:
            sym_to_id[s] = fid
    gc12_path = results / "12_bodymap" / "panel_gene_check.csv"
    if gc12_path.exists():                                        # ids named in the results tables win
        gc12 = pd.read_csv(gc12_path, dtype={"feature_ID": str})
        sym_to_id.update({s: f for f, s in zip(gc12["feature_ID"], gc12["gene_symbol"]) if isinstance(s, str)})
    sym_to_id.update({s: f for f, s in zip(stab["feature_ID"], stab["gene_symbol"]) if isinstance(s, str)})
    ids = [sym_to_id.get(s) for s in symbols]
    assert all(ids) and len(set(ids)) == 20, ("k20 panel symbols did not all map to distinct feature IDs", symbols, ids)
    site = C.ROOT / "site" / "data" / "genes.json"
    if site.exists():                                             # informational cross-check with the exported site data
        exported = set(json.loads(site.read_text()).get("sets", {}).get("k20", []))
        print(f"  k20 ids equal the site export's sets.k20: {exported == set(ids)}" if exported else "  (site export has no k20 set)")
    return ids


def n_per_sex_group(tissue_codes: list[str]) -> pd.DataFrame:
    """Study vials per tissue × sex × group in the exported RNA-seq tables (the DESeq2 n before outlier removal):
    only the header row of each counts (else norm) table is read, joined to PHENO for sex and group."""
    ph = io.load_pheno()[["sex", "group"]]
    rows = []
    for t in tissue_codes:
        p = C.COUNTS_DIR / f"TRNSCRPT__{C.tissue_token(t)}.csv"
        if not p.exists():
            p = C.NORM_DIR / f"TRNSCRPT__{C.tissue_token(t)}.csv"
        if not p.exists():
            rows.append({"tissue": t, "n_per_sex_group_min": np.nan, "n_per_sex_group_max": np.nan, "n_vials": np.nan})
            continue
        vials = [c for c in pd.read_csv(p, nrows=0).columns if c not in C.LEADING_COLS and str(c).startswith("9")]
        m = ph.reindex(vials).dropna()
        ct = m.groupby(["sex", "group"]).size()
        rows.append({"tissue": t, "n_per_sex_group_min": int(ct.min()), "n_per_sex_group_max": int(ct.max()), "n_vials": len(m)})
    return pd.DataFrame(rows)


def main() -> None:
    ap = cli.common_parser("Training response of the panel genes in their marker tissues")
    args = ap.parse_args()
    cli.banner("05_panel_training_response", args)
    out = cli.outdir("05_panels/TRNSCRPT", args.out)
    results = C.RESULTS_DIR
    stab = pd.read_csv(results / "05_panels" / "TRNSCRPT" / "stability_k20_annotated.csv", dtype={"feature_ID": str})
    core = pd.read_csv(results / "05_panels" / "TRNSCRPT" / "candidate_panel_annotated.csv", dtype={"feature_ID": str})
    k20 = k20_ids(out, results, stab)
    sets = {"k20": list(k20), "core": core["feature_ID"].tolist(), "stability_51": stab["feature_ID"].tolist()}
    ann = pd.concat([stab, core[~core["feature_ID"].isin(stab["feature_ID"])]], ignore_index=True).set_index("feature_ID")
    genes = list(dict.fromkeys(sets["core"] + sets["k20"] + sets["stability_51"]))
    missing = [g for g in genes if g not in ann.index]
    if missing:
        gc12 = pd.read_csv(results / "12_bodymap" / "panel_gene_check.csv", dtype={"feature_ID": str}).set_index("feature_ID")
        print(f"  {len(missing)} genes are not in the annotated panel files; marker tissue from panel_gene_check.csv, no tissue effect: {missing}")
        extra = gc12.loc[[g for g in missing if g in gc12.index], ["gene_symbol", "marker_tissue"]]
        ann = pd.concat([ann, extra])
    print(f"  {len(genes)} genes: k20 {len(sets['k20'])}, core {len(sets['core'])}, stability_51 {len(sets['stability_51'])}")

    reg = pd.read_csv(C.RAW_DIR / "training_regulated_features.csv", dtype=str, usecols=["feature_ID", "assay", "tissue"])
    reg = reg[reg["assay"] == "TRNSCRPT"][["feature_ID", "tissue"]].drop_duplicates()
    reg_tissues = reg.groupby("feature_ID")["tissue"].agg(lambda s: sorted(set(s))).to_dict()
    reg_by_tissue = reg.groupby("tissue")["feature_ID"].agg(set).to_dict()

    tissues = sorted(set(ann.loc[genes, "marker_tissue"].dropna()))
    per_gene_max = {}
    background = []
    da_rows = {}
    for t in tissues:
        da = read_da(t)
        tested = set(da["feature_ID"])
        abs_max = da.assign(a=da["logFC"].abs()).groupby("feature_ID")["a"].max()
        per_gene_max[t] = abs_max
        regulated = reg_by_tissue.get(t, set())
        sel = set(da.loc[da["selection_fdr"] < 0.05, "feature_ID"])
        background.append({"tissue": t, "n_genes_tested": len(tested), "n_regulated_5pct": len(regulated & tested),
                           "frac_regulated_5pct": len(regulated & tested) / len(tested),
                           "regulated_set_equals_selection_fdr_lt_0.05": regulated == sel,
                           "median_abs_logfc_regulated": float(abs_max.reindex(sorted(regulated & tested)).median()) if regulated & tested else np.nan,
                           "median_abs_logfc_all_genes": float(abs_max.median()),
                           "n_contrasts": int(da.groupby(["sex", "comparison_group"]).ngroups)})
        da_rows[t] = da[da["feature_ID"].isin(genes)].set_index(["feature_ID", "sex", "comparison_group"])
    bg = pd.DataFrame(background).merge(n_per_sex_group(tissues), on="tissue", how="left")
    bg.to_csv(out / "panel_training_background.csv", index=False)

    rows = []
    for g in genes:
        a = ann.loc[g]
        t = a["marker_tissue"]
        d = da_rows.get(t)
        sub = d.loc[g] if (d is not None and g in d.index.get_level_values(0)) else None
        def val(sex, grp, col):
            try:
                return float(sub.loc[(sex, grp), col]) if sub is not None else np.nan
            except KeyError:
                return np.nan
        max_abs = float(sub["logFC"].abs().max()) if sub is not None else np.nan
        which = sub["logFC"].abs().idxmax() if sub is not None else None
        eff = float(a["effect_size"]) if "effect_size" in a.index and pd.notna(a["effect_size"]) else np.nan
        reg_t = reg_tissues.get(g, [])
        rows.append({"feature_ID": g, "symbol": a["gene_symbol"], "in_k20": g in set(sets["k20"]), "in_core": g in set(sets["core"]),
                     "in_stability_51": g in set(sets["stability_51"]), "marker_tissue": t,
                     "next_highest_tissue": a.get("next_highest_tissue", np.nan), "tissue_effect_log2cpm": eff,
                     "max_abs_logfc_marker": max_abs, "contrast_of_max_abs_logfc": f"{which[0]}_{which[1]}" if which is not None else None,
                     "logfc_8w_female": val("female", "8w", "logFC"), "logfc_8w_female_se": val("female", "8w", "logFC_se"),
                     "logfc_8w_male": val("male", "8w", "logFC"), "logfc_8w_male_se": val("male", "8w", "logFC_se"),
                     "min_adj_p_marker": float(sub["adj_p_value"].min()) if sub is not None else np.nan,
                     "min_selection_fdr_marker": float(sub["selection_fdr"].min()) if sub is not None else np.nan,
                     "tested_in_marker_da": sub is not None,
                     "regulated_marker_5pct": t in reg_t, "n_tissues_regulated_anywhere": len(reg_t), "regulated_tissues": ";".join(reg_t),
                     "ratio_tissue_over_training": eff / max_abs if (np.isfinite(eff) and max_abs and np.isfinite(max_abs)) else np.nan})
    resp = pd.DataFrame(rows).sort_values(["marker_tissue", "symbol"]).reset_index(drop=True)
    resp.to_csv(out / "panel_training_response.csv", index=False)

    summ = []
    for name, ids in sets.items():
        s = resp[resp["feature_ID"].isin(ids)]
        summ.append({"gene_set": name, "n_genes": len(s), "n_regulated_marker_5pct": int(s["regulated_marker_5pct"].sum()),
                     "regulated_marker_genes": ";".join(s.loc[s["regulated_marker_5pct"], "symbol"]),
                     "n_regulated_anywhere": int((s["n_tissues_regulated_anywhere"] > 0).sum()),
                     "max_abs_logfc_marker": float(s["max_abs_logfc_marker"].max()),
                     "min_max_abs_logfc_marker": float(s["max_abs_logfc_marker"].min()),
                     "min_tissue_effect": float(s["tissue_effect_log2cpm"].min()), "max_tissue_effect": float(s["tissue_effect_log2cpm"].max()),
                     "median_ratio_tissue_over_training": float(s["ratio_tissue_over_training"].median()),
                     "min_ratio_tissue_over_training": float(s["ratio_tissue_over_training"].min())})
    summary = pd.DataFrame(summ)
    summary.to_csv(out / "panel_training_summary.csv", index=False)

    show = ["symbol", "marker_tissue", "in_k20", "in_core", "tissue_effect_log2cpm", "max_abs_logfc_marker", "contrast_of_max_abs_logfc",
            "logfc_8w_female", "logfc_8w_male", "min_adj_p_marker", "regulated_marker_5pct", "n_tissues_regulated_anywhere", "ratio_tissue_over_training"]
    with pd.option_context("display.width", 250, "display.max_rows", 200):
        print(resp[show].to_string(index=False))
        print(bg.to_string(index=False))
        print(summary.to_string(index=False))
    n_lo, n_hi = int(bg["n_per_sex_group_min"].min()), int(bg["n_per_sex_group_max"].max())
    notes = "\n\n".join([
        "# Phase 05c: training response of the panel genes in their marker tissues",
        __doc__.split("\n", 1)[1].strip(),
        f"## Run\n\nDESeq2 n per sex × group in these tables (TRNSCRPT sample metadata, marker tissues): {n_lo}–{n_hi} vials. "
        "Every marker tissue's consortium 5 % FDR set equals `selection_fdr < 0.05` in its DA table: "
        f"{bool(bg['regulated_set_equals_selection_fdr_lt_0.05'].all())}. `max_abs_logfc_marker` is the largest |log2FC| over the "
        f"{int(bg['n_contrasts'].max())} sex × duration contrasts; `ratio_tissue_over_training` = `tissue_effect_log2cpm` / `max_abs_logfc_marker`.",
        "### Per gene set [measured]\n\n" + report.df_to_md(summary, floatfmt=".3f"),
        "### Per gene [measured]\n\n" + report.df_to_md(resp[show], floatfmt=".3f", max_rows=100),
        "### Per marker tissue [measured]\n\n" + report.df_to_md(bg, floatfmt=".3f"),
    ])
    (out / "panel_training_NOTES.md").write_text(notes + "\n")
    report.add_section("05c · Training response of the panel genes (consortium DA)", notes.split("\n", 2)[2],
                       params={"n_genes": len(genes), "tissues": len(tissues)})
    print(f"wrote {out / 'panel_training_response.csv'}, panel_training_background.csv, panel_training_summary.csv, panel_training_NOTES.md")


if __name__ == "__main__":
    main()
