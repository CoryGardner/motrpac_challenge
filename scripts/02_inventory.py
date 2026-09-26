#!/usr/bin/env python
"""Phase 02 — what is actually on disk: tissues × assays, samples, animals, overlap, fusion set."""
from __future__ import annotations

import json

import numpy as np
import pandas as pd

from motrpac import cli, config as C, io, plots, report


def main() -> None:
    ap = cli.common_parser("Inventory of exported MoTrPAC tables")
    args = ap.parse_args()
    cli.banner("02_inventory", args)
    out = cli.outdir("02_inventory", args.out)
    pheno = io.load_pheno()
    manifest = io.load_manifest()

    rows = []
    for kind in ("norm", "counts"):
        files = io.list_sample_files(kind)
        for _, r in files.iterrows():
            feat, vals = io.read_sample_table(r["path"])
            meta = io.sample_meta_for(vals.columns, pheno, r["tissue"], r["assay"])
            rows.append({"kind": kind, "assay": r["assay"], "tissue": r["tissue"], "n_features": vals.shape[0],
                         "n_samples": vals.shape[1], "n_animals": int(meta["pid"].nunique()),
                         "n_unmatched": int(meta["pid"].isna().sum()),
                         "n_female": int((meta["sex"] == "female").sum()), "n_male": int((meta["sex"] == "male").sum()),
                         "groups": ",".join(f"{g}:{n}" for g, n in meta["group"].value_counts().sort_index().items()),
                         "frac_nan": float(vals.isna().to_numpy().mean()),
                         "n_regulated_flagged": int(feat["feature"].notna().sum()) if "feature" in feat else 0})
    inv = pd.DataFrame(rows).sort_values(["kind", "assay", "tissue"]).reset_index(drop=True)
    inv.to_csv(out / "inventory.csv", index=False)
    print(inv.to_string())

    # tissue × assay presence (norm) and the fusion tissue set
    norm = inv[inv["kind"] == "norm"]
    present = norm.pivot_table(index="tissue", columns="assay", values="n_samples", aggfunc="sum", fill_value=0)
    present.to_csv(out / "tissue_by_assay.csv")
    fig = plots.missingness_heatmap(present, "samples per tissue × assay (normalized data)", out / "tissue_by_assay.png")
    core = [a for a in ("TRNSCRPT", "PROT", "METAB") if a in present.columns]
    fusion_tissues = list(present.index[(present[core] > 0).all(axis=1)]) if core else []

    # animal-level overlap across assays within each fusion tissue
    overlap_rows = []
    for t in fusion_tissues:
        mats = {a: io.load_norm(a, t, pheno) for a in core}
        bids = {a: set(m.meta["bid"]) for a, m in mats.items()}
        common = set.intersection(*bids.values())
        overlap_rows.append({"tissue": t, **{f"n_{a}": len(b) for a, b in bids.items()}, "n_common_bid": len(common)})
    overlap = pd.DataFrame(overlap_rows)
    overlap.to_csv(out / "fusion_overlap.csv", index=False)

    # DA table inventory
    da_rows = []
    for p in sorted((C.DA_DIR).glob("*.csv")) if C.DA_DIR.exists() else []:
        assay, tok = p.stem.split("__")
        da = pd.read_csv(p, nrows=5)
        n = sum(1 for _ in open(p)) - 1
        da_rows.append({"assay": assay, "tissue": C.tissue_from_token(tok), "n_rows": n, "columns": ",".join(da.columns[:12])})
    da_inv = pd.DataFrame(da_rows)
    da_inv.to_csv(out / "da_inventory.csv", index=False)

    # animals overall
    animals = pheno.drop_duplicates("pid")
    design = animals.groupby(["sex", "group"]).size().unstack(fill_value=0)
    design = design.reindex(columns=[g for g in C.GROUP_ORDER if g in design.columns])
    design.to_csv(out / "animals_by_sex_group.csv")

    body = [
        f"Data package: `{manifest.get('package', '?')}` version `{manifest.get('version', '?')}`"
        + (" — **SYNTHETIC**" if cli.is_synthetic() else ""),
        f"Animals in PHENO: **{animals.shape[0]}**; vials: {pheno.shape[0]}.",
        "Animals by sex × group:\n\n" + report.df_to_md(design.reset_index()),
        "Sample-level tables:\n\n" + report.df_to_md(inv.drop(columns=["groups"])),
        report.figure_md(fig, "Samples per tissue × assay"),
        f"Tissues with all of {', '.join(core)}: **{', '.join(fusion_tissues) or 'none'}**",
        ("Animal overlap within fusion tissues (join on bid):\n\n" + report.df_to_md(overlap)) if len(overlap) else "",
        ("DA tables:\n\n" + report.df_to_md(da_inv)) if len(da_inv) else "No DA tables found.",
        "Check against the paper: ~18 tissues TRNSCRPT, ~7 PROT, most tissues METAB; 3–6 animals per sex × group per assay.",
    ]
    report.add_section("02 · Inventory", "\n\n".join(b for b in body if b), params={"raw_dir": str(C.RAW_DIR)})
    (out / "summary.json").write_text(json.dumps({"fusion_tissues": fusion_tissues, "n_animals": int(animals.shape[0]),
                                                  "tables": int(len(inv))}, indent=2))
    print(f"fusion tissues: {fusion_tissues}")


if __name__ == "__main__":
    main()
