#!/usr/bin/env python
"""Write a small synthetic dataset in the exact export layout (docs/DATA_GUIDE.md §4).

Purpose: smoke-test every phase before the real export exists, and give Claude Code something
concrete to develop against. The biology is fake but the *structure* is faithful: animals with
pid/bid/viallabel, both sexes, 5 groups, tissue-specific markers, a confusable tissue pair,
RNA–protein pairs with a discordant subset, DA tables per sex × time point.

Usage: python scripts/make_synthetic_data.py [--out data/raw] [--seed 1]
A marker file `_SYNTHETIC` is written so later phases can flag synthetic results.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from motrpac import config as C

TRN_TISSUES = ["SKM-GN", "SKM-VL", "HEART", "LIVER", "KIDNEY", "WAT-SC", "BAT", "CORTEX"]
PROT_TISSUES = ["SKM-GN", "HEART", "LIVER", "KIDNEY", "WAT-SC"]
METAB_TISSUES = ["SKM-GN", "HEART", "LIVER", "KIDNEY", "WAT-SC", "CORTEX"]


def bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(ranked, 1.0)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(C.RAW_DIR))
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--per-cell", type=int, default=6, help="animals per sex × group")
    ap.add_argument("--n-genes", type=int, default=1500)
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)
    out = Path(args.out)
    for d in ("norm", "counts", "da", "meta", "codes"):
        (out / d).mkdir(parents=True, exist_ok=True)

    # ---- animals ------------------------------------------------------------------------
    animals = []
    pid_counter = 10_000_000
    for sex in ("female", "male"):
        for grp in C.GROUP_ORDER:
            for _ in range(args.per_cell):
                pid_counter += 1
                animals.append({"pid": str(pid_counter), "sex": sex, "group": grp})
    animals = pd.DataFrame(animals)
    n_animals = len(animals)
    animal_effect = rng.normal(0, 0.3, size=(n_animals, args.n_genes))

    # ---- gene model ---------------------------------------------------------------------
    G = args.n_genes
    gene_ids = [f"ENSRNOG{i:011d}" for i in range(1, G + 1)]
    gene_syms = [f"GENE{i:04d}" for i in range(1, G + 1)]
    base = rng.normal(5, 2, size=G)
    tissue_effect = {t: np.zeros(G) for t in TRN_TISSUES}
    marker_blocks = {}
    block = max(4, min(60, (G - 50) // len(TRN_TISSUES)))   # markers per tissue, scales with --n-genes
    ptr = 0
    for t in TRN_TISSUES:
        if t == "SKM-VL":
            # confusable pair: shares SKM-GN markers, differs in only a few genes
            tissue_effect[t] = tissue_effect["SKM-GN"].copy()
            nb = max(2, block // 8)
            tissue_effect[t][ptr:ptr + nb] += 2.5
            marker_blocks[t] = list(range(ptr, ptr + nb))
            ptr += nb
        else:
            tissue_effect[t][ptr:ptr + block] += 4.0
            marker_blocks[t] = list(range(ptr, ptr + block))
            ptr += block
    sex_genes = list(range(ptr, ptr + 20)); ptr += 20
    train_genes = list(range(ptr, ptr + 30)); ptr += 30
    time_factor = {"control": 0.0, "1w": 0.3, "2w": 0.5, "4w": 0.8, "8w": 1.0}

    # ---- protein / metabolite maps -------------------------------------------------------
    n_prot = min(600, int(0.4 * G))
    prot_gene_idx = rng.choice(G, size=n_prot, replace=False)
    prot_ids = [f"NP_{i:06d}" for i in range(1, n_prot + 1)]
    discordant = np.zeros(n_prot, dtype=bool)
    discordant[rng.choice(n_prot, size=int(0.3 * n_prot), replace=False)] = True   # 30% of proteins ignore RNA
    metab_ids = [f"metab_{i:03d}" for i in range(1, 201)]
    metab_tissue_effect = {t: rng.normal(0, 1.0, size=200) for t in METAB_TISSUES}

    # ---- generate samples --------------------------------------------------------------
    pheno_rows, bid_counter = [], 10_000
    norm = {("TRNSCRPT", t): {} for t in TRN_TISSUES}
    counts = {t: {} for t in TRN_TISSUES}
    prot = {t: {} for t in PROT_TISSUES}
    metab = {t: {} for t in METAB_TISSUES}
    all_tissues = sorted(set(TRN_TISSUES) | set(PROT_TISSUES) | set(METAB_TISSUES))
    for ai, a in animals.iterrows():
        for t in all_tissues:
            bid_counter += 1
            bid = str(bid_counter)
            if t in TRN_TISSUES:
                x = base + tissue_effect[t] + animal_effect[ai] + rng.normal(0, 0.5, size=G)
                if a["sex"] == "male":
                    x[sex_genes] += 1.5
                x[train_genes] += 0.8 * time_factor[a["group"]]
                if t in ("SKM-GN", "SKM-VL"):
                    x[train_genes] += 0.6 * time_factor[a["group"]]  # muscle responds more
                x = np.clip(x, 0, None)
                vl = f"{bid}{rng.integers(100000, 999999)}"
                norm[("TRNSCRPT", t)][vl] = x
                lib = rng.integers(15_000_000, 30_000_000)
                counts[t][vl] = rng.poisson((2 ** x - 1) / 1e6 * lib)
                pheno_rows.append({"viallabel": vl, "bid": bid, "pid": a["pid"], "sex": a["sex"],
                                   "group": a["group"], "tissue": t, "assay": "TRNSCRPT"})
                rna_for_prot = x
            else:
                rna_for_prot = None
            if t in PROT_TISSUES:
                r = rna_for_prot[prot_gene_idx] if rna_for_prot is not None else base[prot_gene_idx]
                p = 0.6 * (r - base[prot_gene_idx]) + rng.normal(0, 0.6, size=n_prot)
                p[discordant] = rng.normal(0, 0.6, size=discordant.sum())
                vl = f"{bid}{rng.integers(100000, 999999)}"
                prot[t][vl] = p
                pheno_rows.append({"viallabel": vl, "bid": bid, "pid": a["pid"], "sex": a["sex"],
                                   "group": a["group"], "tissue": t, "assay": "PROT"})
            if t in METAB_TISSUES:
                m = metab_tissue_effect[t] + rng.normal(0, 0.7, size=200)
                if a["sex"] == "male":
                    m[:10] += 0.8
                vl = f"{bid}{rng.integers(100000, 999999)}"
                metab[t][vl] = m
                pheno_rows.append({"viallabel": vl, "bid": bid, "pid": a["pid"], "sex": a["sex"],
                                   "group": a["group"], "tissue": t, "assay": "METAB"})

    # ---- write sample-level tables -----------------------------------------------------
    def write_table(d: dict, ids: list[str], tissue: str, assay: str, folder: str, regulated: set | None = None):
        df = pd.DataFrame(d, index=ids)
        lead = pd.DataFrame({"feature": [f"{assay};{tissue};{i}" if regulated and i in regulated else np.nan for i in ids],
                             "feature_ID": ids, "tissue": tissue, "assay": assay}, index=ids)
        pd.concat([lead, df], axis=1).to_csv(out / folder / f"{assay}__{C.tissue_token(tissue)}.csv", index=False)

    reg_genes = {gene_ids[i] for i in train_genes}
    for t in TRN_TISSUES:
        write_table(norm[("TRNSCRPT", t)], gene_ids, t, "TRNSCRPT", "norm", reg_genes)
        write_table(counts[t], gene_ids, t, "TRNSCRPT", "counts")
    for t in PROT_TISSUES:
        write_table(prot[t], prot_ids, t, "PROT", "norm")
    for t in METAB_TISSUES:
        write_table(metab[t], metab_ids, t, "METAB", "norm")

    # ---- PHENO ---------------------------------------------------------------------------
    ph = pd.DataFrame(pheno_rows)
    ph["study_group_timepoint"] = np.where(ph["group"] == "control", "control - 8w", "training - " + ph["group"])
    ph["specimen_processing___sampletypedescription"] = ph["tissue"]
    ph["registration___sex"] = ph["sex"].map({"female": "1", "male": "2"})
    ph["vo2.max.test.vo2_max"] = rng.normal(70, 8, size=len(ph)).round(1)
    ph["terminal.weight.bw"] = rng.normal(300, 40, size=len(ph)).round(1)
    ph.drop(columns=["tissue", "assay"]).to_csv(out / "pheno.csv", index=False)

    # ---- DA tables (t-tests within sex, group vs control) --------------------------------
    def da_table(d: dict, ids: list[str], tissue: str, assay: str, statname: str):
        df = pd.DataFrame(d, index=ids)
        meta = ph.set_index("viallabel").loc[df.columns]
        rows = []
        for sex in ("female", "male"):
            ctrl = df.loc[:, ((meta["sex"] == sex) & (meta["group"] == "control")).values]
            for grp in ("1w", "2w", "4w", "8w"):
                trt = df.loc[:, ((meta["sex"] == sex) & (meta["group"] == grp)).values]
                tt = stats.ttest_ind(trt.to_numpy(), ctrl.to_numpy(), axis=1, equal_var=False)
                p = np.nan_to_num(tt.pvalue, nan=1.0)
                rows.append(pd.DataFrame({"feature_ID": ids, "tissue": tissue, "assay": assay, "sex": sex,
                                          "comparison_group": grp,
                                          "logFC": (trt.mean(axis=1) - ctrl.mean(axis=1)).values,
                                          statname: np.nan_to_num(tt.statistic), "p_value": p,
                                          "adj_p_value": bh(p)}))
        pd.concat(rows).to_csv(out / "da" / f"{assay}__{C.tissue_token(tissue)}.csv", index=False)

    for t in TRN_TISSUES:
        da_table(norm[("TRNSCRPT", t)], gene_ids, t, "TRNSCRPT", "zscore")
    for t in PROT_TISSUES:
        da_table(prot[t], prot_ids, t, "PROT", "tscore")

    # ---- annotation tables ----------------------------------------------------------------
    f2g = pd.concat([
        pd.DataFrame({"feature_ID": gene_ids, "gene_symbol": gene_syms, "ensembl_gene": gene_ids, "assay": "TRNSCRPT"}),
        pd.DataFrame({"feature_ID": prot_ids, "gene_symbol": [gene_syms[i] for i in prot_gene_idx],
                      "ensembl_gene": [gene_ids[i] for i in prot_gene_idx], "assay": "PROT"}),
    ])
    f2g.to_csv(out / "feature_to_gene.csv", index=False)
    pd.DataFrame({"RAT_ENSEMBL_ID": gene_ids, "RAT_SYMBOL": gene_syms,
                  "HUMAN_ORTHOLOG_SYMBOL": [s if i % 10 else np.nan for i, s in enumerate(gene_syms)]}
                 ).to_csv(out / "rat_to_human_gene.csv", index=False)
    pd.DataFrame({"feature_ID": [gene_ids[i] for i in train_genes], "assay": "TRNSCRPT"}
                 ).to_csv(out / "training_regulated_features.csv", index=False)
    pd.DataFrame(columns=["viallabel", "tissue", "assay", "reason"]).to_csv(out / "outliers.csv", index=False)
    pd.DataFrame({"name": C.TISSUES, "value": C.TISSUES}).to_csv(out / "codes" / "tissue_abbrev.csv", index=False)
    pd.DataFrame({"name": C.ASSAYS, "value": C.ASSAYS}).to_csv(out / "codes" / "assay_abbrev.csv", index=False)
    pd.DataFrame({"name": ["TRNSCRPT", "PROT", "METAB"], "value": ["#1f77b4", "#ff7f0e", "#2ca02c"]}
                 ).to_csv(out / "codes" / "assay_colors.csv", index=False)
    manifest = {"package": "SYNTHETIC", "version": f"seed{args.seed}", "note": "fake data for smoke tests",
                "n_animals": int(n_animals), "tissues_trnscrpt": TRN_TISSUES, "tissues_prot": PROT_TISSUES,
                "tissues_metab": METAB_TISSUES, "discordant_protein_fraction": 0.3,
                "confusable_pair": ["SKM-GN", "SKM-VL"]}
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2))
    (out / "_SYNTHETIC").write_text("This directory holds SYNTHETIC data. Delete it before running the real export.\n")
    print(f"synthetic data written to {out}: {n_animals} animals, {len(ph)} vials")


if __name__ == "__main__":
    main()
