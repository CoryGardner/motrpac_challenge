#!/usr/bin/env python
"""Time course, part 3b: what does RNA sequencing depth separate in control vs 8w?

1. Verify that RNA plate, extraction date, library batch / date, flowcell, lane, sequencing date / batch / machine are constant
   within each of the 7 fusion tissues (all 50 study vials and the control+8w subset).
2. Plate layout: Lib_barcode_well (row A–H, column) by group — were durations / control vs 8w placed in different columns?
   Column-major well index vs sacrifice order.
3. Within control+8w and WITHIN SEX (sex is blocked by sacrifice date there): Spearman rho (n reported) of reads and reads_raw
   vs sacrifice date, sacrifice order (date + time of death), time of death, well index / row / column, RNA extraction
   concentration, RIN, 260/280, 260/230, library RNA / DNA concentration, molarity, fragment size, and physiology (terminal body
   weight, NMR weight / fat % pre and post, VO2max pre and post, time to freeze), plus group (8w vs control, as rho with a 0/1 code).
4. Is depth a property of the animal or the index? Cross-tissue per-animal depth correlation, and whether an animal carries
   the same index sequence in different tissues.
Uses the pipeline's data/raw/pheno.csv (dotted column names, e.g. vo2.max.test.vo2_max_1) and data/raw/meta/TRNSCRPT.csv
(= the c1.0 RNA QC table, study vials only).
Run from code/pipeline:  PYTHONPATH=src python investigations/time_course/tc3b_depth_puzzle.py
Outputs: results/15_time_course/2_3_gradients/{depth_constants.csv, depth_layout.csv, depth_within_sex_spearman.csv,
depth_within_sex_summary.csv, depth_cross_tissue.csv, depth_decomposition.csv, depth_composition_by_duration.csv, depth_layout.png}
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

from motrpac import config as C

import os
import sys
sys.path.insert(0, os.path.dirname(__file__))
from tc2_molecular_gradient import OUT, TISSUES  # noqa: E402

CONST = ["RNA_extr_plate_ID", "RNA_extr_date", "Lib_batch_ID", "Lib_prep_date", "Lib_kit_id", "Lib_robot", "Seq_flowcell_ID",
         "Seq_flowcell_run", "Seq_flowcell_lane", "Seq_date", "Seq_batch", "Seq_machine_ID"]
LIB_NUM = ["RNA_extr_conc", "RIN", "r_260_280", "r_260_230", "Lib_RNA_conc", "Lib_DNA_conc", "Lib_molarity", "Lib_frag_size"]
PHYS = {"terminal.weight.bw": "body_weight_terminal", "nmr.testing.nmr_weight_2": "nmr_weight_post", "nmr.testing.nmr_fat_2": "nmr_fat_pct_post",
        "nmr.testing.nmr_fat_1": "nmr_fat_pct_pre", "vo2.max.test.vo2_max_1": "vo2max_pre", "vo2.max.test.vo2_max_2": "vo2max_post",
        "time_to_freeze": "time_to_freeze"}
DEPTH = ["reads", "reads_raw"]


def hours(s):
    td = pd.to_timedelta(s.astype(str).where(s.notna()), errors="coerce")
    return td.dt.total_seconds() / 3600


def main():
    ph = pd.read_csv(C.RAW_DIR / "pheno.csv", dtype=str, low_memory=False)
    tm = pd.read_csv(C.META_DIR / "TRNSCRPT.csv", dtype=str, low_memory=False).drop_duplicates("viallabel")
    tm = tm[tm["viallabel"].str.startswith("9")]
    keep = ["viallabel", "pid", "sex", "group", "tissue", "key.d_sacrifice", "specimen.collection.t_death"] + list(PHYS)
    m = tm.merge(ph[keep].drop_duplicates("viallabel"), on="viallabel", how="inner")
    m = m[m["tissue"].isin(TISSUES)].copy()
    for c in DEPTH + LIB_NUM:
        m[c] = pd.to_numeric(m[c], errors="coerce")
    for c, new in PHYS.items():
        m[new] = pd.to_numeric(m[c], errors="coerce")
    m["sac_day"] = (pd.to_datetime(m["key.d_sacrifice"], format="%d%b%Y") - pd.Timestamp("2018-01-01")).dt.days
    m["t_death_h"] = hours(m["specimen.collection.t_death"])
    w = m["Lib_barcode_well"].str.upper().str.extract(r"^([A-H])(\d{1,2})$")
    m["well_row"] = w[0].map({c: i + 1 for i, c in enumerate("ABCDEFGH")})
    m["well_col"] = w[1].astype(float)
    m["well_col_rel"] = m["well_col"] - m.groupby("tissue")["well_col"].transform("min") + 1
    m["well_index"] = (m["well_col_rel"] - 1) * 8 + m["well_row"]

    # 1. constants
    crow = []
    for t in TISSUES:
        for subset, d in (("all_groups", m[m["tissue"] == t]), ("control+8w", m[(m["tissue"] == t) & m["group"].isin(["control", "8w"])])):
            crow.append({"tissue": t, "subset": subset, "n_vials": len(d), **{c: int(d[c].nunique(dropna=False)) for c in CONST},
                         "lane_value": ";".join(sorted(d["Seq_flowcell_lane"].dropna().unique()))})
    const = pd.DataFrame(crow)
    const.to_csv(OUT / "depth_constants.csv", index=False)
    print(const.to_string(index=False))

    # 2. layout
    lrow = []
    for t in TISSUES:
        d = m[m["tissue"] == t]
        for g, dg in d.groupby("group"):
            lrow.append({"tissue": t, "group": g, "n": len(dg), "cols_rel": ",".join(str(int(x)) for x in sorted(dg["well_col_rel"].unique())),
                         "well_index_min": int(dg["well_index"].min()), "well_index_max": int(dg["well_index"].max()),
                         "sac_first": dg["key.d_sacrifice"].iloc[dg["sac_day"].argmin()], "sac_last": dg["key.d_sacrifice"].iloc[dg["sac_day"].argmax()]})
        # order of wells vs sacrifice chronology, all 50 vials
        r, p = stats.spearmanr(d["well_index"], d["sac_day"] * 24 + d["t_death_h"].fillna(12))
        lrow.append({"tissue": t, "group": "ALL (rho well_index vs sacrifice datetime)", "n": len(d), "cols_rel": f"rho={r:.3f}, p={p:.2g}"})
        c8 = d[d["group"].isin(["control", "8w"])]
        ct = pd.crosstab(c8["well_col_rel"], c8["group"])
        lrow.append({"tissue": t, "group": "control+8w: col × group (Fisher/chi2 p)", "n": len(c8),
                     "cols_rel": "; ".join(f"col{int(i)}: " + ",".join(f"{g}={v}" for g, v in row.items()) for i, row in ct.iterrows())
                     + f" | p={stats.chi2_contingency(ct)[1]:.2g}"})
    lay = pd.DataFrame(lrow)
    lay.to_csv(OUT / "depth_layout.csv", index=False)
    print(lay.to_string(index=False))

    # 3. within-sex correlations in control+8w
    c8 = m[m["group"].isin(["control", "8w"])].copy()
    c8["sac_dt"] = c8["sac_day"] * 24 + c8["t_death_h"]
    c8["sac_order"] = c8.groupby(["tissue", "sex"])["sac_dt"].rank()
    c8["is_8w"] = (c8["group"] == "8w").astype(float)
    xvars = (["is_8w", "sac_day", "sac_order", "t_death_h", "well_index", "well_row", "well_col_rel"] + LIB_NUM + list(PHYS.values()))
    srow = []
    for t in TISSUES:
        for sx in ("male", "female"):
            d = c8[(c8["tissue"] == t) & (c8["sex"] == sx)]
            for dv in DEPTH:
                for xv in xvars:
                    ok = d[[dv, xv]].dropna()
                    if len(ok) < 5 or ok[xv].nunique() < 2:
                        srow.append({"tissue": t, "sex": sx, "depth": dv, "variable": xv, "n": len(ok), "rho": np.nan, "p": np.nan})
                        continue
                    r, p = stats.spearmanr(ok[dv], ok[xv])
                    srow.append({"tissue": t, "sex": sx, "depth": dv, "variable": xv, "n": len(ok), "rho": r, "p": p})
    sp = pd.DataFrame(srow)
    sp.to_csv(OUT / "depth_within_sex_spearman.csv", index=False)
    # summary across the 14 tissue × sex cells: median rho, # p<0.05 (expected 0.7 by chance), sign consistency
    summ = sp.dropna(subset=["rho"]).groupby(["depth", "variable"]).agg(
        n_cells=("rho", "size"), median_n=("n", "median"), median_rho=("rho", "median"), n_pos=("rho", lambda v: int((v > 0).sum())),
        n_p05=("p", lambda v: int((v < 0.05).sum()))).reset_index()
    # sign test across cells
    summ["sign_test_p"] = [stats.binomtest(int(r.n_pos), int(r.n_cells)).pvalue for r in summ.itertuples()]
    summ = summ.sort_values(["depth", "n_p05", "sign_test_p"], ascending=[True, False, True])
    summ.to_csv(OUT / "depth_within_sex_summary.csv", index=False)
    print(summ[summ["depth"] == "reads"].round(3).to_string(index=False))
    sig = sp[(sp["p"] < 0.05) & (sp["depth"] == "reads")].sort_values("p")
    print("reads, cells with p<0.05:\n", sig.round(3).to_string(index=False))

    # 4. cross-tissue: per-animal depth consistency and index identity
    wide = m.pivot_table(index="pid", columns="tissue", values="reads")
    z = wide.rank(pct=True)
    xrow = []
    for i, a in enumerate(TISSUES):
        for b in TISSUES[i + 1:]:
            ok = z[[a, b]].dropna()
            r, p = stats.spearmanr(ok[a], ok[b])
            same_idx = m[m["tissue"].isin([a, b])].pivot_table(index="pid", columns="tissue", values="Lib_index_1", aggfunc="first").dropna()
            xrow.append({"tissue_a": a, "tissue_b": b, "n_animals": len(ok), "rho_reads": r, "p": p,
                         "frac_same_index1": float((same_idx[a] == same_idx[b]).mean()) if len(same_idx) else np.nan,
                         "same_plate_cols": bool(m.loc[m.tissue == a, "well_col"].min() == m.loc[m.tissue == b, "well_col"].min())})
    xt = pd.DataFrame(xrow)
    # well-index effect: does the same relative well give similar depth across tissues? (depth rank per tissue, averaged per well)
    m["reads_rank"] = m.groupby("tissue")["reads"].rank(pct=True)
    wr = m.pivot_table(index="well_index", columns="tissue", values="reads_rank")
    # animal effect vs well effect (they coincide if layout identical); report whether the animal→well map is identical across tissues
    same_map = m.pivot_table(index="pid", columns="tissue", values="well_index").nunique(axis=1).le(1).mean()
    xt["frac_animals_same_rel_well_all_tissues"] = same_map
    xt.to_csv(OUT / "depth_cross_tissue.csv", index=False)
    print(xt.round(3).to_string(index=False))
    print(f"median cross-tissue rho of per-animal depth rank: {xt['rho_reads'].median():.3f} over {len(xt)} pairs; "
          f"fraction of animals at the same relative well in every tissue: {same_map:.2f}")

    # 5. decomposition: is the "depth" separation read count or spliced-read composition? (control vs 8w, pooled and within sex)
    for c in ("uniquely_mapped", "num_splices", "pct_chrM", "pct_mrna", "pct_intronic", "pct_uniquely_mapped", "pct_rRNA", "pct_globin"):
        c8[c] = pd.to_numeric(c8[c], errors="coerce")
    c8["splices_per_mapped"] = c8["num_splices"] / c8["uniquely_mapped"]
    drow = []
    for t in TISSUES:
        for sx in ("pooled", "male", "female"):
            d = c8[c8["tissue"] == t] if sx == "pooled" else c8[(c8["tissue"] == t) & (c8["sex"] == sx)]
            y = d["group"] == "8w"
            for v in ("reads_raw", "reads", "uniquely_mapped", "num_splices", "splices_per_mapped", "pct_chrM", "pct_mrna", "pct_intronic",
                      "pct_uniquely_mapped", "pct_rRNA", "pct_globin", "RIN"):
                a, b = d.loc[y, v].dropna(), d.loc[~y, v].dropna()
                if len(a) < 2 or len(b) < 2:
                    continue
                u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
                drow.append({"tissue": t, "sex": sx, "variable": v, "n_8w": len(a), "n_control": len(b),
                             "auroc_8w_higher": u / (len(a) * len(b)), "p": p, "median_8w": a.median(), "median_control": b.median()})
        d = c8[c8["tissue"] == t]
        for v in ("pct_chrM", "pct_mrna", "pct_intronic"):
            ok = d[["splices_per_mapped", v]].dropna()
            r, p = stats.spearmanr(ok["splices_per_mapped"], ok[v])
            drow.append({"tissue": t, "sex": "pooled", "variable": f"rho(splices_per_mapped, {v})", "n_8w": int((d["group"] == "8w").sum()),
                         "n_control": int((d["group"] == "control").sum()), "auroc_8w_higher": r, "p": p})
    dec = pd.DataFrame(drow)
    dec.to_csv(OUT / "depth_decomposition.csv", index=False)
    print(dec[dec["tissue"].isin(["HEART", "KIDNEY", "SKM-GN"])].round(3).to_string(index=False))

    # 6. the same composition / depth variables for every duration vs control (pooled sexes; 1w/2w/4w carry the cohort confound)
    allm = m.copy()
    for c in ("uniquely_mapped", "num_splices", "pct_chrM", "pct_mrna", "pct_intronic"):
        allm[c] = pd.to_numeric(allm[c], errors="coerce")
    allm["splices_per_mapped"] = allm["num_splices"] / allm["uniquely_mapped"]
    drow = []
    for t in TISSUES:
        for dur in ("1w", "2w", "4w", "8w"):
            d = allm[(allm["tissue"] == t) & allm["group"].isin(["control", dur])]
            y = d["group"] == dur
            for v in ("reads", "splices_per_mapped", "pct_chrM", "pct_mrna", "pct_intronic", "RIN"):
                a, b = d.loc[y, v].dropna(), d.loc[~y, v].dropna()
                if len(a) < 2 or len(b) < 2:
                    continue
                u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
                drow.append({"tissue": t, "duration": dur, "variable": v, "n_trained": len(a), "n_control": len(b),
                             "auroc_trained_higher": u / (len(a) * len(b)), "p": p})
    bd = pd.DataFrame(drow)
    bd.to_csv(OUT / "depth_composition_by_duration.csv", index=False)
    print(bd.pivot_table(index=["tissue", "variable"], columns="duration", values="auroc_trained_higher").round(2).to_string())

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, 2, figsize=(11, 4))
    col = {"control": "#444444", "1w": "#e69f00", "2w": "#56b4e9", "4w": "#009e73", "8w": "#d55e00"}
    d = m[m["tissue"] == "HEART"]
    for g, dg in d.groupby("group"):
        axs[0].scatter(dg["well_index"], dg["sac_day"], c=col[g], label=g, s=25)
    axs[0].set_xlabel("relative plate well index (column-major), HEART RNA plate")
    axs[0].set_ylabel("sacrifice day (days since 2018-01-01)")
    axs[0].legend(fontsize=7, frameon=False)
    axs[0].set_title("Library plate order follows sacrifice chronology (same map in all 7 tissues)", fontsize=8)
    for t in ("HEART", "SKM-GN", "KIDNEY"):
        d = c8[c8["tissue"] == t]
        for sx, mk in (("male", "o"), ("female", "^")):
            dd = d[d["sex"] == sx]
            axs[1].scatter(dd["sac_order"] + (0.15 if t == "SKM-GN" else 0.3 if t == "KIDNEY" else 0), dd["reads"] / 1e6, marker=mk,
                           c=[col[g] for g in dd["group"]], s=22, alpha=0.8)
    axs[1].set_xlabel("sacrifice order within sex (control+8w), HEART / SKM-GN / KIDNEY offset")
    axs[1].set_ylabel("reads (millions)")
    axs[1].set_title("Depth vs sacrifice order within sex (o male, ^ female; grey control, red 8w)", fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "depth_layout.png", dpi=150)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
