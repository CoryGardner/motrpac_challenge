#!/usr/bin/env python
"""Phase 03b — the EDA read-out: headline numbers pulled from results/03_eda*/ CSVs and the
per-tissue tables, with the interpretation that goes to the hackathon. Run after 03_eda.py
(default run in results/03_eda, inner-join run in results/03_eda_inner, core-9 metabolomics run in
results/03_eda_metab_core9). Writes results/03_eda/readout_*.csv and appends one REPORT.md section.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd

from tfp import cli, config as C, io, report

_spec = importlib.util.spec_from_file_location("eda", Path(__file__).with_name("03_eda.py"))
eda = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(eda)

R = C.RESULTS_DIR
ASSAYS = ["TRNSCRPT", "PROT", "METAB"]


def read(path: Path) -> pd.DataFrame | None:
    return pd.read_csv(path) if path.exists() else None


def main() -> None:
    ap = cli.common_parser("EDA read-out")
    args = ap.parse_args()
    cli.banner("03_eda_readout", args)
    out = cli.outdir("03_eda", args.out)
    pheno = io.load_pheno()
    body = []

    # ---- 1. feature joins across tissues -------------------------------------------------
    join_rows, spec_rows = [], []
    for a in ASSAYS:
        source = cli.resolve_source(a, "auto")
        files = io.list_sample_files(source)
        files = files[files["assay"] == a]
        sets = {}
        for _, r in files.iterrows():
            om = io.load_counts(r["tissue"], pheno, assay=a) if source == "counts" else io.load_norm(a, r["tissue"], pheno)
            sets[r["tissue"]] = set(om.X.columns)
        union = set().union(*sets.values())
        inner = set.intersection(*sets.values())
        outer_n = len(read(R / "03_eda" / f"feature_stats_{a}.csv"))
        inner_run = read(R / "03_eda_inner" / f"feature_stats_{a}.csv")
        join_rows.append({"assay": a, "source": source, "tissues": len(sets), "union": len(union),
                          "outer_min_present_0.8": outer_n, "inner": len(inner),
                          "inner_run": len(inner_run) if inner_run is not None else np.nan,
                          "outer_features_absent_in_some_tissue": outer_n - len(inner)})
        if a == "TRNSCRPT":
            n_t = pd.Series({f: sum(f in s for s in sets.values()) for f in union})
            spec_rows = [{"expressed_in_n_tissues": "1–3", "genes": int((n_t <= 3).sum())},
                         {"expressed_in_n_tissues": "4–18", "genes": int(((n_t > 3) & (n_t < len(sets))).sum())},
                         {"expressed_in_n_tissues": f"all {len(sets)}", "genes": int((n_t == len(sets)).sum())}]
    joins = pd.DataFrame(join_rows)
    joins.to_csv(out / "readout_joins.csv", index=False)
    spec = pd.DataFrame(spec_rows)
    core9 = read(R / "03_eda_metab_core9" / "feature_stats_METAB.csv")
    fp9 = read(R / "03_eda_metab_core9" / "feature_platforms_METAB.csv")
    fp19 = read(R / "03_eda" / "feature_platforms_METAB.csv")
    body.append("### 1. Which features survive a cross-tissue join\n\n" + report.df_to_md(joins) + "\n\n"
                "`outer_min_present_0.8` is what `03_eda.py` stacks by default; `inner` is what phases 04–08 use. "
                "Every feature counted in the last column is absent from at least one tissue and would be "
                "**median-imputed** there under the outer join, which turns missingness into a tissue label — a leak "
                "of the class into the features, not biology. That is the reason to keep the inner join for every "
                "cross-tissue classifier.\n\n"
                f"METAB collapses because only `metab-u-hilicpos` was run on all 19 tissues: the outer join keeps "
                f"{joins.loc[joins.assay == 'METAB', 'outer_min_present_0.8'].item()} features, of which "
                f"{int(fp19.loc[fp19.platform == 'metab-u-hilicpos', 'n_features'].sum()) if fp19 is not None else '?'} "
                f"are hilicpos-only and the rest hilicpos plus another platform; the inner join keeps "
                f"{joins.loc[joins.assay == 'METAB', 'inner'].item()}. Restricting to the 9 tissues that have 10–13 "
                f"platforms (SKM-GN, HEART, KIDNEY, LIVER, LUNG, BAT, WAT-SC, HIPPOC, PLASMA) an inner join keeps "
                f"**{len(core9) if core9 is not None else '?'}** features from "
                f"{fp9['platform'].str.split('|').explode().nunique() if fp9 is not None else '?'} platforms. "
                "Use: inner join over all 19 tissues (hilicpos, identical provenance everywhere) for the tissue-identity "
                "task; inner join over the core tissues for fusion (the 7 fusion tissues are all core). The stacked PCA "
                "is the same story under either join (tissue R² on PC1–3 differs by < 0.02), so nothing in this EDA "
                "depends on the choice; the classifiers do.\n\n"
                "TRNSCRPT has the mirror-image problem: `load_counts` filters genes per tissue, so the inner join drops "
                "tissue-specific genes — exactly the genes a compact panel would pick:\n\n" + report.df_to_md(spec) + "\n\n"
                "Decision for T4/T5: stack raw counts with an outer join where an absent gene is **0 counts, not the "
                "median**, and filter on the stacked matrix, so tissue-specific genes stay available and no imputation "
                "artefact is created.")

    # ---- 2. what the leading PCs are ------------------------------------------------------
    vrows = []
    for a in ASSAYS:
        v = read(R / "03_eda" / f"variance_partition_{a}.csv")
        vi = read(R / "03_eda_inner" / f"variance_partition_{a}.csv")
        for i in range(3):
            row = {"assay": a, "PC": v.loc[i, "PC"], "explained": v.loc[i, "explained"]}
            for c in ("R2_tissue", "R2_sex", "R2_group", "R2_pid"):
                row[c] = v.loc[i, c] if c in v.columns else np.nan
            row["R2_tissue_inner_join"] = vi.loc[i, "R2_tissue"] if vi is not None else np.nan
            vrows.append(row)
    vt = pd.DataFrame(vrows)
    vt.to_csv(out / "readout_variance_pc1-3.csv", index=False)
    p = vt[vt.assay == "PROT"].iloc[0]
    body.append("### 2. What the leading PCs are, per assay (stacked tissues, outer join)\n\n"
                + report.df_to_md(vt, floatfmt=".3f") + "\n\n"
                f"Transcripts and metabolites are tissue on PC1–3 (R² ≥ 0.97 and ≥ 0.97). Proteomics is not: tissue "
                f"explains {p['R2_tissue']:.3f} of PC1 and {vt[vt.assay == 'PROT'].iloc[1]['R2_tissue']:.3f} of PC2, "
                f"while sex explains {p['R2_sex']:.2f} and the animal {p['R2_pid']:.2f}. This is the normalization: "
                "`PROT_*_NORM_DATA` are log2 ratios to a per-tissue reference pool, median-centred per sample, so every "
                "tissue sits at zero and only within-tissue biology is left. A cross-tissue proteomic fingerprint built "
                "on these tables would learn residual covariance and missingness, not protein abundance — the data "
                "answer to Q15. Metabolomics keeps tissue identity despite per-sample centring because the *profile* "
                "differs by tissue.")

    # ---- 3. time course within tissue ----------------------------------------------------
    wrows, top = [], []
    for a in ASSAYS:
        w = read(R / "03_eda" / f"within_tissue_pca_{a}.csv")
        w = w.assign(assay=a, ratio=w["max_R2_group_PC1-5"] / w["null95_group"])
        wrows.append({"assay": a, "tissues": len(w), "group_visible": int(w["group_visible"].sum()),
                      "sex_visible": int((w["max_R2_sex_PC1-5"] > w["null95_sex"]).sum()),
                      "median_max_R2_group": w["max_R2_group_PC1-5"].median(),
                      "median_null95_group": w["null95_group"].median(),
                      "not_visible": ", ".join(w.loc[~w["group_visible"], "tissue"])})
        top.append(w.sort_values("max_R2_group_PC1-5", ascending=False).head(5)[
            ["assay", "tissue", "n", "n_features", "max_R2_group_PC1-5", "null95_group", "PC_with_max_group_R2", "max_R2_sex_PC1-5"]])
    ws = pd.DataFrame(wrows)
    ws.to_csv(out / "readout_within_tissue_summary.csv", index=False)
    topdf = pd.concat(top)
    body.append("### 3. Is the time course visible within a tissue?\n\n" + report.df_to_md(ws, floatfmt=".2f") + "\n\n"
                "Strongest tissues per assay (max R² of training group over within-tissue PCs 1–5, own full feature "
                "table; `null95` = 95th percentile of the same maximum with labels shuffled):\n\n"
                + report.df_to_md(topdf, floatfmt=".2f") + "\n\n"
                "Read: the time course is visible, but never on the axes that carry most variance — in transcripts it "
                "sits on PC3 in HEART and on PC5 elsewhere; sex is the dominant within-tissue axis in LIVER, KIDNEY, "
                "WAT-SC, HEART, ADRNL (R² 0.8–0.98). Metabolomics shows the training signal more broadly than "
                "transcripts. A `group_visible` at 1.0–1.3× the null (e.g. OVARY, TESTES, COLON, CORTEX, SKM-VL) is "
                "at the edge of chance and should not be read as a result. See the within-tissue PC1/PC2 figures "
                "coloured by group in `results/03_eda/within_pca_*_by_group.png`.")

    # ---- 4. batch --------------------------------------------------------------------------
    wt = read(R / "03_eda" / "within_tissue_pca_TRNSCRPT.csv")
    qc = wt[["tissue", "PC1_explained", "max_R2_sex_PC1-5", "max_R2_group_PC1-5",
             "max_R2_pct_mrna_PC1-5", "max_R2_pct_mrna|design_PC1-5", "R2_pct_mrna~sex", "R2_pct_mrna~group",
             "max_R2_pct_uniquely_mapped_PC1-5", "max_R2_pct_uniquely_mapped|design_PC1-5",
             "max_R2_pct_rRNA_PC1-5", "max_R2_pct_rRNA|design_PC1-5", "R2_pct_rRNA~sex"]]
    qc_flag = qc[(qc["max_R2_pct_mrna|design_PC1-5"] >= 0.3) | (qc["max_R2_pct_uniquely_mapped|design_PC1-5"] >= 0.3)
                 | (qc["max_R2_pct_rRNA|design_PC1-5"] >= 0.3)].sort_values("max_R2_pct_mrna|design_PC1-5", ascending=False)
    bp_t = read(R / "03_eda" / "batch_partition_TRNSCRPT.csv").head(3)
    n_const = {c: int((wt[f"n_{c}"] == 1).sum()) for c in ("GET_site", "RNA_extr_plate_ID", "Lib_batch_ID", "Seq_flowcell_ID")}

    wp = read(R / "03_eda" / "within_tissue_pca_PROT.csv")
    ch = wp[["tissue", "n", "max_R2_group_PC1-5", "max_R2_tmt11_channel_PC1-5", "max_R2_tmt11_channel|design_PC1-5",
             "R2_tmt11_channel~sex", "R2_tmt11_channel~group", "max_R2_tmt_plex_PC1-5", "max_R2_tmt_plex|design_PC1-5"]]
    def locate(assay: str, cov: str, source: str) -> pd.DataFrame:
        """Per tissue: the within-tissue PC on which `cov` has the largest design-residualized R², that
        PC's variance share, and the sex / group R² on the same PC (to show it is not design)."""
        files = io.list_sample_files(source)
        rows = []
        for t in sorted(files.loc[files["assay"] == assay, "tissue"]):
            om = io.load_counts(t, pheno, assay=assay) if source == "counts" else io.load_norm(assay, t, pheno)
            keep = (om.X.notna().mean(axis=0) >= 0.8).to_numpy()
            om = om.select_features(list(om.X.columns[keep]))
            if om.n_samples < 8 or om.n_features < 5:
                continue
            St, evt = eda.pca_of(om, n_components=5, seed=args.seed)
            bc = eda.batch_covariates(assay, om.X.index)
            if cov not in bc.columns or bc[cov].nunique() < 2:
                continue
            cell = (om.meta["sex"].astype(str) + "/" + om.meta["group"].astype(str)).to_numpy()
            r_res, r_raw = eda.r2_of(eda.residualize(St, cell), bc[cov]), eda.r2_of(St, bc[cov])
            j = int(np.nanargmax(r_res))
            r_sex = eda.r2_categorical(St, om.meta["sex"])[j] if om.meta["sex"].nunique() > 1 else 0.0
            r_grp = eda.r2_categorical(St, om.meta["group"])[j] if om.meta["group"].nunique() > 1 else 0.0
            rows.append({"tissue": t, "PC": f"PC{j+1}", "PC_variance_share": evt[j], f"R2_{cov}": r_raw[j],
                         f"R2_{cov}|design": r_res[j], "R2_sex_same_PC": r_sex, "R2_group_same_PC": r_grp})
        return pd.DataFrame(rows)

    loc = locate("PROT", "tmt11_channel", "norm")
    loc.to_csv(out / "readout_prot_channel_location.csv", index=False)
    loc_t = locate("TRNSCRPT", "pct_mrna", cli.resolve_source("TRNSCRPT", "auto"))
    loc_t.to_csv(out / "readout_trnscrpt_pct_mrna_location.csv", index=False)
    loc_t_flag = loc_t[loc_t["R2_pct_mrna|design"] >= 0.3].sort_values("R2_pct_mrna|design", ascending=False)
    bp_p = read(R / "03_eda" / "batch_partition_PROT.csv").head(3)
    bp_m = read(R / "03_eda" / "batch_partition_METAB.csv").head(3)

    body.append("### 4. Batch structure\n\n"
                "**TRNSCRPT.** Sequencing site, RNA-extraction plate, library batch and flowcell are constant within "
                f"every tissue ({', '.join(f'{c}: {n}/{len(wt)} tissues' for c, n in n_const.items())}): each tissue was "
                "processed as one unit, so batch and tissue identity are perfectly confounded on the stacked matrix "
                "(the nested `tissue|plate` factor adds 0.000 over tissue on PC1–3) and cannot be separated by any "
                "analysis of these data. Within a tissue there is no batch variation to test, but the library QC "
                "metrics are informative — the leading within-tissue PC is a **library-composition axis** in several "
                "tissues, and it survives removal of the sex × group design means:\n\n"
                + report.df_to_md(qc_flag, floatfmt=".2f") + "\n\n"
                "Where that axis sits (tissues whose design-residualized `pct_mrna` R² on some PC is ≥ 0.3; the sex and "
                "group R² on the *same* PC show it is not the design):\n\n" + report.df_to_md(loc_t_flag, floatfmt=".2f")
                + "\n\nRead: in ADRNL, BAT, VENACV and SMLINT a leading within-tissue PC is mRNA fraction / mapping "
                "rate, not biology; in BAT, VENACV and SKM-GN that metric is also associated with training group "
                "(`~group` 0.35–0.37), so within-tissue 'training' signals there must be checked against `pct_mrna`. "
                "In HEART, LIVER and WAT-SC the same metrics track sex (`~sex` 0.5–0.8) and vanish after the design "
                "is removed — sex-linked RNA composition, not batch. Stacked PC1–3 with the batch columns:\n\n"
                + report.df_to_md(bp_t, floatfmt=".2f") + "\n\n"
                "**PROT.** TMT plex explains 0.00 of every PC, within and across tissues: it was regressed out by "
                "`limma::removeBatchEffect` and any biology aligned with plex went with it. TMT **channel** is a "
                "different story: after removing the sex × group means it still explains up to "
                f"{ch['max_R2_tmt11_channel|design_PC1-5'].max():.2f} of a within-tissue PC "
                f"({ch.loc[ch['max_R2_tmt11_channel|design_PC1-5'].idxmax(), 'tissue']}), more than training group does "
                "in KIDNEY, LUNG, CORTEX and SKM-GN. In HEART and LIVER channel *is* sex (`~sex` = 1.0: the plex layout "
                "put the sexes in fixed channels), elsewhere channel is independent of sex and group:\n\n"
                + report.df_to_md(ch, floatfmt=".2f") + "\n\n"
                "Where the channel effect sits (PC with the largest design-residualized channel R², and that PC's share "
                "of within-tissue variance):\n\n" + report.df_to_md(loc, floatfmt=".2f") + "\n\n"
                "Read plainly: on the stacked proteomics PCs, tissue explains 0.001 / 0.004 / 0.21 of PC1–3 and the "
                "channel-within-tissue factor 0.05 / 0.68 / 0.40 — channel explains more of PC2 and PC3 than tissue "
                "does, although part of that is sex riding on channel in HEART and LIVER. For within-tissue "
                "proteomics tasks (fusion task B, discordance) `tmt11_channel` must be a covariate or a blocking "
                "factor; for cross-tissue tasks the normalization issue (section 2) dominates everything else.\n\n"
                + report.df_to_md(bp_p, floatfmt=".2f") + "\n\n"
                "**METAB.** The package has no sample-level run or batch metadata for metabolomics (no `METAB_META`), "
                "so the only covariate available is how many platforms measured a sample, which is a property of the "
                "tissue (constant within 13 of 17 tissues) and adds nothing over tissue on the stacked PCs "
                "(`tissue|n_platforms` = tissue). A metabolomics batch check would need run-order or site metadata, "
                "which the package does not carry.\n\n" + report.df_to_md(bp_m, floatfmt=".2f"))

    # ---- 5. single-sex tissues ---------------------------------------------------------------
    ss = []
    for a in ASSAYS:
        t = read(R / "03_eda" / f"tissue_by_sex_{a}.csv")
        for _, r in t.iterrows():
            if (r.get("female", 0) == 0) != (r.get("male", 0) == 0):
                ss.append({"assay": a, "tissue": r["tissue"], "female": int(r.get("female", 0)), "male": int(r.get("male", 0))})
    ssd = pd.DataFrame(ss)
    body.append("### 5. Single-sex tissues and the sex-shift experiment\n\n" + report.df_to_md(ssd) + "\n\n"
                "Confirmed in the stacked metadata: OVARY is female-only and TESTES male-only in TRNSCRPT and METAB; "
                "PROT has neither. For T7 (`08_shift_tests.py`, held-out sex): the source never contains the target's "
                "gonad class, so under train-male/test-female OVARY is an unseen class (24 of 449 target samples) and "
                "under train-female/test-male TESTES is (25 of 450). No prediction set can cover them; they are "
                "reported as unseen and excluded from 'seen' coverage, and the fix made this session guarantees they "
                "never enter calibration either. The honest framing of the sex shift is therefore label shift plus an "
                "open-set problem for 2 of 19 classes, not pure covariate shift.")

    report.add_section("03 · EDA read-out (real data)", "\n\n".join(body),
                       params={"inputs": "results/03_eda, results/03_eda_inner, results/03_eda_metab_core9"})
    print(f"wrote {out}/readout_*.csv and the REPORT.md section")


if __name__ == "__main__":
    main()
