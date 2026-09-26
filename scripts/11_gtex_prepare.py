#!/usr/bin/env python
"""Phase 11 — prepare a GTEx subset for the human-transfer leg (docs/GTEX_TRANSFER.md).

Inputs in data/external/gtex/ (download from https://gtexportal.org/home/downloads/adult-gtex):
  --tpm    gene TPM GCT file (gzipped ok), e.g. GTEx_Analysis_*_RNASeQCv*_gene_tpm.gct.gz
  --reads  gene read-count GCT file (optional), e.g. GTEx_Analysis_*_RNASeQCv*_gene_reads.gct.gz
  --attrs  sample attributes, e.g. GTEx_Analysis_*_SampleAttributesDS.txt
Outputs (data/external/): gtex_tpm_subset.csv (samples × genes, log2(TPM+1)), gtex_meta.csv
(SAMPID, donor, SMTSD, rat_tissue), gtex_gene_symbols.csv, and with --reads gtex_cpm_subset.csv
(samples × genes, log2(CPM+1) with the library size = the sample's read total over every gene in
the GCT — the same unit as the MoTrPAC and BodyMap count matrices, motrpac.io.log_cpm). Only the
tissues that map to the MoTrPAC rat set are kept, ≤ --max-per-tissue donors per tissue, and the
CPM matrix uses exactly the samples of the TPM subset.

An output that already exists is not rebuilt unless --force is given (the inputs are ~1–2 GB
each and the parse takes minutes). File names change between GTEx releases; pass them explicitly.
"""
from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np
import pandas as pd

from motrpac import config as C

GTEX_TO_RAT = {
    "Whole Blood": "BLOOD", "Muscle - Skeletal": "SKM-GN", "Adipose - Subcutaneous": "WAT-SC",
    "Heart - Left Ventricle": "HEART", "Liver": "LIVER", "Kidney - Cortex": "KIDNEY", "Lung": "LUNG",
    "Brain - Cortex": "CORTEX", "Brain - Hippocampus": "HIPPOC", "Brain - Hypothalamus": "HYPOTH",
    "Colon - Transverse": "COLON", "Small Intestine - Terminal Ileum": "SMLINT", "Spleen": "SPLEEN",
    "Adrenal Gland": "ADRNL", "Ovary": "OVARY", "Testis": "TESTES", "Artery - Aorta": "VENACV",
}


def read_gct_columns(path: str, wanted: set[str]) -> pd.DataFrame:
    """GCT: 2 header lines, then columns Name, Description, <samples>; keep the wanted samples."""
    with pd.read_csv(path, sep="\t", skiprows=2, chunksize=2000, dtype={"Name": str, "Description": str}) as reader:
        parts = []
        for chunk in reader:
            cols = ["Name", "Description"] + [c for c in chunk.columns if c in wanted]
            parts.append(chunk[cols])
    df = pd.concat(parts)
    df["Name"] = df["Name"].str.split(".").str[0]   # strip Ensembl version
    return df.drop_duplicates("Name")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tpm", required=True)
    ap.add_argument("--reads", default=None, help="gene read-count GCT; writes gtex_cpm_subset.csv")
    ap.add_argument("--attrs", required=True)
    ap.add_argument("--out", default=str(C.EXTERNAL_DIR))
    ap.add_argument("--max-per-tissue", type=int, default=150, help="subsample donors per tissue to keep it light")
    ap.add_argument("--force", action="store_true", help="rebuild outputs that already exist")
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    tpm_out, meta_out, cpm_out = out / "gtex_tpm_subset.csv", out / "gtex_meta.csv", out / "gtex_cpm_subset.csv"

    if tpm_out.exists() and meta_out.exists() and not args.force:
        meta = pd.read_csv(meta_out, dtype=str)
        print(f"{tpm_out.name} and {meta_out.name} exist ({len(meta)} samples); not rebuilt (use --force)")
    else:
        attrs = pd.read_csv(args.attrs, sep="\t", dtype=str, usecols=lambda c: c in ("SAMPID", "SMTSD", "SMAFRZE"))
        attrs = attrs[attrs["SMTSD"].isin(GTEX_TO_RAT)]
        if "SMAFRZE" in attrs.columns:
            attrs = attrs[attrs["SMAFRZE"] == "RNASEQ"]
        attrs["donor"] = attrs["SAMPID"].str.split("-").str[:2].str.join("-")
        attrs["rat_tissue"] = attrs["SMTSD"].map(GTEX_TO_RAT)
        rng = np.random.default_rng(C.SEED)
        keep = []
        for t, d in attrs.groupby("SMTSD"):
            donors = d["donor"].unique()
            if len(donors) > args.max_per_tissue:
                donors = rng.choice(donors, size=args.max_per_tissue, replace=False)
            keep.append(d[d["donor"].isin(donors)])
        meta = pd.concat(keep)
        t0 = time.time()
        tpm = read_gct_columns(args.tpm, set(meta["SAMPID"]))
        tpm[["Name", "Description"]].set_index("Name")["Description"].to_csv(out / "gtex_gene_symbols.csv", header=["symbol"])
        X = np.log2(tpm.set_index("Name").drop(columns=["Description"]).T.astype(float) + 1.0)
        X.index.name = "SAMPID"
        X.to_csv(tpm_out)
        meta = meta.set_index("SAMPID").loc[X.index]
        meta.to_csv(meta_out)
        meta = meta.reset_index()
        print(f"wrote {X.shape[0]} samples × {X.shape[1]} genes (log2 TPM+1) in {time.time() - t0:.0f}s; tissues: {meta['SMTSD'].nunique()}")

    if args.reads:
        if cpm_out.exists() and not args.force:
            print(f"{cpm_out.name} exists; not rebuilt (use --force)")
            return
        t0 = time.time()
        wanted = set(meta["SAMPID"])
        reads = read_gct_columns(args.reads, wanted)
        counts = reads.set_index("Name").drop(columns=["Description"]).T.astype(np.float64)   # samples × all GCT genes
        counts = counts.loc[[s for s in meta["SAMPID"] if s in counts.index]]
        lib = counts.sum(axis=1)
        cpm = np.log2(counts.div(lib.replace(0, np.nan), axis=0) * 1e6 + 1.0)
        cpm.index.name = "SAMPID"
        cpm.to_csv(cpm_out)
        print(f"wrote {cpm.shape[0]} samples × {cpm.shape[1]} genes (log2 CPM+1; library size = read total over all "
              f"{counts.shape[1]} GCT genes, median {lib.median() / 1e6:.1f} M reads) in {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
