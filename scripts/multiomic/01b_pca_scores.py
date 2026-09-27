#!/usr/bin/env python
"""Multiomic Phase 1b — per-vial PCA scores of the two proteomics matrices behind
results_multiomic/01_rii/variance_partition.csv (reporter-ion log2 ppm) and variance_partition_ratio.csv (the
distributed ratios), for the "one picture" chart on the Multiomic page. No new analysis: the matrices are rebuilt with
the loading code of scripts/multiomic/01_rii_rescue.py (same defaults: prot-pr, release c1.0, >= 2 peptides, channel
total; proteins quantified in every tissue; the distributed ratios via io.stack_tissues(join="inner")) and projected with
its pca_scores() and seed.

Gate: the tissue R² of PC1 and PC2 recomputed from the written scores must equal the published values to 4 decimals;
otherwise the script exits non-zero and writes nothing.
Writes results_multiomic/01_rii/pca_scores_rii.csv and pca_scores_ratio.csv
(viallabel, tissue, sex, plex_id, PC1, PC2, explained_PC1, explained_PC2).
Usage: python scripts/multiomic/01b_pca_scores.py
"""
from __future__ import annotations

import importlib.util
import sys

import numpy as np
import pandas as pd

from tfp import config as C, io, rii

ROOT = C.ROOT
OUT = ROOT / "results_multiomic" / "01_rii"


def _rescue():
    spec = importlib.util.spec_from_file_location("rii_rescue", ROOT / "scripts" / "multiomic" / "01_rii_rescue.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main() -> int:
    R = _rescue()
    seed = C.SEED
    pheno = io.load_pheno()
    # the reporter-ion matrix (01_rii_rescue.py section 1, defaults)
    om_outer, _ = rii.stack_rii(rii.PROT_TISSUES, "prot-pr", "c1.0", pheno, join="outer", min_peptides=2, norm="total")
    common = [c for c in om_outer.X.columns if om_outer.X[c].notna().groupby(om_outer.meta["tissue"].to_numpy()).any().all()]
    X_rii, meta_rii = om_outer.X[common], om_outer.meta
    # the distributed ratio matrix (01_rii_rescue.py section 2)
    om_ratio = io.stack_tissues("PROT", source="norm", join="inner", pheno=pheno, verbose=False)
    meta_r = om_ratio.meta.copy()
    pm = pd.read_csv(C.META_DIR / "PROT.csv", dtype=str).drop_duplicates("viallabel").set_index("viallabel")
    meta_r["plex_id"] = meta_r["tissue"] + ":" + pm["tmt_plex"].reindex(meta_r.index).astype(str).to_numpy()

    frames, problems = {}, []
    for tag, X, meta, pub in (("rii", X_rii, meta_rii, "variance_partition.csv"), ("ratio", om_ratio.X, meta_r, "variance_partition_ratio.csv")):
        S, evr = R.pca_scores(X.to_numpy(dtype=float), seed=seed)
        vp = pd.read_csv(OUT / pub)
        df = pd.DataFrame({"viallabel": X.index.astype(str), "tissue": meta["tissue"].astype(str).to_numpy(), "sex": meta["sex"].astype(str).to_numpy(),
                           "plex_id": meta["plex_id"].astype(str).to_numpy(), "PC1": S[:, 0], "PC2": S[:, 1], "explained_PC1": evr[0], "explained_PC2": evr[1]})
        r2 = R.r2_categorical(df[["PC1", "PC2"]].to_numpy(), df["tissue"])
        r2s = R.r2_categorical(df[["PC1", "PC2"]].to_numpy(), df["sex"])
        for j, pc in enumerate(("PC1", "PC2")):
            row = vp[vp["PC"] == pc].iloc[0]
            for name, got, want in ((f"{tag} {pc} tissue R²", r2[j], row["R2_tissue"]), (f"{tag} {pc} sex R²", r2s[j], row["R2_sex"]), (f"{tag} {pc} explained", evr[j], row["explained"])):
                ok = round(float(got), 4) == round(float(want), 4)
                print(f"  {name:28s} recomputed {got:.6f}  published {want:.6f}  {'ok' if ok else 'DIFFERENT'}")
                if not ok:
                    problems.append((name, float(got), float(want)))
        print(f"  {tag}: {len(df)} vials × {X.shape[1]} proteins")
        frames[tag] = df
    if problems:
        print("GATE FAILED — not writing the scores:", problems)
        return 1
    for tag, df in frames.items():
        df.to_csv(OUT / f"pca_scores_{tag}.csv", index=False)
    print(f"gate passed; wrote {OUT / 'pca_scores_rii.csv'} and {OUT / 'pca_scores_ratio.csv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
