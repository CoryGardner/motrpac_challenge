"""RII (reporter-ion intensity) proteomics from the MoTrPAC portal quant-id folders — the scale on which a
cross-tissue comparison is defined.

The R-package `PROT_*_NORM_DATA` tables are log2 ratios of every channel to its plex's reference pool, median-centred
per sample; the reference pool is tissue-specific, so the tissue axis is removed by construction (phase 03/04: tissue
R² of PC1 = 0.0008). The portal's `*_rii-results.txt` files hold the raw reporter-ion intensity of every peptide in
every channel, including the reference channels (`Ref_S1..Ref_S6`). Rolled up to proteins and normalised to the
channel's total protein signal (parts per million), the values are comparable across plexes and across tissues in the
sense "fraction of this sample's quantified protein mass" — still a within-study measurement (a plex is one tissue),
but one on which tissue means can differ.

Added for the overnight multiomic run (branch multiomic-overnight); nothing here is imported by the frozen pipeline.
"""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd

from . import config as C, io

QUANT_ID = Path(os.environ.get("MOTRPAC_QUANT_ID",
                               "/home/cory/projects/MoTrPAC_back_to_transcriptome/data/quant-id/rat-training-06"))

# tissue code → portal tissue directory (the 7 proteomics tissues)
TISSUE_DIRS = {"CORTEX": "t53-cortex", "SKM-GN": "t55-gastrocnemius", "HEART": "t58-heart", "KIDNEY": "t59-kidney",
               "LUNG": "t66-lung", "LIVER": "t68-liver", "WAT-SC": "t70-white-adipose"}
PROT_TISSUES = list(TISSUE_DIRS)
LEAD_COLS = ["protein_id", "redundant_ids", "is_contaminant", "peptide_score", "sequence", "gene_symbol", "entrez_id",
             "organism_name"]


def rii_paths(tissue: str, assay: str = "prot-pr", release: str = "c1.0", root: Path | None = None) -> tuple[Path, Path]:
    """(rii-results file, vial-metadata file) for one tissue × PTM assay × portal release."""
    root = Path(root) if root else QUANT_ID
    d = root / release / "proteomics-untargeted" / TISSUE_DIRS[tissue] / assay
    suffix = "" if release == "c1.0" else "_v2.0"
    stem = f"motrpac_pass1b-06_{TISSUE_DIRS[tissue]}_{assay}"
    return d / f"{stem}_rii-results{suffix}.txt", d / f"{stem}_vial-metadata{suffix}.txt"


def load_rii_peptides(tissue: str, assay: str = "prot-pr", release: str = "c1.0",
                      root: Path | None = None) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Peptide-level RII table (lead columns + one column per vial + Ref_S*) and the vial metadata
    (vial_label, tmt11_channel, tmt_plex) for one tissue."""
    p_rii, p_vm = rii_paths(tissue, assay, release, root)
    pep = pd.read_csv(p_rii, sep="\t", low_memory=False)
    lead = [c for c in LEAD_COLS if c in pep.columns]
    pep = pep.rename(columns={c: str(c) for c in pep.columns})
    val_cols = [c for c in pep.columns if c not in lead and c not in ("ptm_id", "confident_site", "flanking_sequence",
                                                                     "ptm_score", "ptm_peptide", "sequence_id")]
    pep[val_cols] = pep[val_cols].apply(pd.to_numeric, errors="coerce")
    vm = pd.read_csv(p_vm, sep="\t", dtype=str)
    vm["vial_label"] = vm["vial_label"].str.strip()
    return pep, vm


def rii_protein_matrix(tissue: str, assay: str = "prot-pr", release: str = "c1.0", min_peptides: int = 2,
                       norm: str = "total", drop_contaminants: bool = True, root: Path | None = None,
                       min_peptides_shared: bool = False) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Peptides → proteins → per-channel normalised log2 values for one tissue.

    Per plex (tissue × S1–S6): the plex's 10 study channels are taken (the reference channel `Ref_Sk` is dropped);
    peptide rows flagged as contaminants are dropped; peptide intensities are summed per `protein_id` (NaN-skipping;
    a protein with no quantified peptide in a channel is NaN there); proteins with fewer than `min_peptides` peptides
    quantified anywhere in the plex are dropped; every channel is divided by its column total over the kept proteins
    (norm='total', × 1e6 → ppm) or by its median over the proteins quantified in all 10 channels (norm='median'),
    then log2. Plexes are combined within the tissue by an outer join on protein_id.

    Returns (log2 matrix proteins × vials, per-plex summary, protein annotation from the RII lead columns).
    """
    pep, vm = load_rii_peptides(tissue, assay, release, root)
    if drop_contaminants and "is_contaminant" in pep.columns:
        is_c = pep["is_contaminant"].astype(str).str.upper().isin(["TRUE", "1", "T"])
        pep = pep[~is_c]
    ann = (pep[[c for c in ("protein_id", "gene_symbol", "entrez_id", "redundant_ids") if c in pep.columns]]
           .drop_duplicates("protein_id").set_index("protein_id"))
    blocks, rows = [], []
    for plex, g in vm.groupby("tmt_plex"):
        # study vials only: reference-pool channels are listed in the vial metadata as Ref_Sk (or labels starting
        # with 8) and are dropped here — they are the per-tissue reference the ratio files divide by
        vials = [v for v in g["vial_label"] if v in pep.columns and str(v)[:1] == "9"]
        if not vials:
            continue
        sub = pep[["protein_id"] + vials]
        quantified = sub[vials].notna().any(axis=1)
        n_pep = sub.loc[quantified].groupby("protein_id").size()
        summed = sub.groupby("protein_id")[vials].sum(min_count=1)
        keep = n_pep[n_pep >= min_peptides].index
        summed = summed.loc[summed.index.isin(keep)]
        if norm == "total":
            tot = summed.sum(axis=0, skipna=True)
            rel = summed.div(tot, axis=1) * 1e6
        elif norm == "median":
            full = summed.dropna(axis=0)
            med = full.median(axis=0)
            rel = summed.div(med, axis=1)
        else:
            raise ValueError(f"unknown norm {norm!r}")
        log2 = np.log2(rel.where(rel > 0))
        blocks.append(log2)
        rows.append({"tissue": tissue, "tmt_plex": plex, "n_vials": len(vials), "n_peptide_rows": int(len(sub)),
                     "n_proteins_quantified": int(len(n_pep)), "n_proteins_kept": int(len(keep)),
                     "n_proteins_dropped_min_peptides": int(len(n_pep) - len(keep)),
                     "median_channel_total": float(summed.sum(axis=0, skipna=True).median()),
                     "frac_nan_kept": float(summed.isna().mean().mean())})
    mat = pd.concat(blocks, axis=1, join="outer")
    mat.index.name = "protein_id"
    return mat, pd.DataFrame(rows), ann.reindex(mat.index)


def stack_rii(tissues: list[str] | None = None, assay: str = "prot-pr", release: str = "c1.0", pheno: pd.DataFrame | None = None,
              join: str = "inner", min_peptides: int = 2, norm: str = "total", root: Path | None = None,
              verbose: bool = True) -> tuple[io.OmicsMatrix, pd.DataFrame]:
    """RII matrices of several tissues stacked into one OmicsMatrix (samples × proteins, log2 ppm), with sample metadata
    from PHENO plus `tmt_plex`, `tmt11_channel` and `plex_id` (tissue:plex — the physical plex). join='inner' keeps
    proteins present (in at least one vial) in every tissue; 'outer' keeps the union. Returns (matrix, per-plex summary)."""
    tissues = tissues or PROT_TISSUES
    pheno = io.load_pheno() if pheno is None else pheno
    mats, metas, anns, plex_rows = [], [], [], []
    for t in tissues:
        mat, plexes, ann = rii_protein_matrix(t, assay, release, min_peptides, norm, root=root)
        vm = pd.read_csv(rii_paths(t, assay, release, root)[1], sep="\t", dtype=str)
        vm["vial_label"] = vm["vial_label"].str.strip()
        vm = vm.set_index("vial_label")
        meta = io.sample_meta_for(mat.columns, pheno, t, "PROT-RII")
        meta["tmt_plex"] = vm["tmt_plex"].reindex(meta.index).to_numpy()
        meta["tmt11_channel"] = vm["tmt11_channel"].reindex(meta.index).to_numpy()
        meta["plex_id"] = t + ":" + meta["tmt_plex"].astype(str)
        keep = ~meta["pid"].isna()
        X = mat.T.loc[keep.to_numpy()]
        X.index.name = "viallabel"
        mats.append(X)
        metas.append(meta.loc[keep])
        anns.append(ann)
        plex_rows.append(plexes)
        if verbose:
            print(f"  RII {t}: {X.shape[0]} vials, {X.shape[1]} proteins ({plexes['n_proteins_kept'].min()}–{plexes['n_proteins_kept'].max()} per plex), "
                  f"NaN fraction {X.isna().mean().mean():.3f}")
    if join == "inner":
        common = set(mats[0].columns)
        for m in mats[1:]:
            common &= set(m.columns)
        common = [c for c in mats[0].columns if c in common]
        X = pd.concat([m[common] for m in mats], axis=0)
    else:
        X = pd.concat(mats, axis=0, join="outer")
    meta = pd.concat(metas, axis=0)
    ann = pd.concat(anns, axis=0)
    ann = ann[~ann.index.duplicated(keep="first")].reindex(X.columns)
    ann["assay"] = f"RII-{assay}-{release}"
    union = len(set().union(*[set(m.columns) for m in mats]))
    om = io.OmicsMatrix(X, meta, ann, "PROT-RII", name=f"RII-{assay}-{release}-{join}")
    om.notes.append(f"stack RII {assay} {release}: {len(mats)} tissues, {X.shape[0]} vials, {X.shape[1]} proteins after join={join} "
                    f"(union {union}); min_peptides={min_peptides}, norm={norm}, log2; NaN fraction {X.isna().mean().mean():.3f}")
    if verbose:
        print("  " + om.notes[-1])
    return om, pd.concat(plex_rows, ignore_index=True)


def protein_to_ensembl(protein_ids, f2g: pd.DataFrame | None = None) -> pd.DataFrame:
    """protein_id → (gene_symbol, ensembl_gene) from data/raw/feature_to_gene.csv; every mapping row is kept
    (a protein can map to more than one Ensembl gene; callers pick the one present in the RNA matrix)."""
    f2g = io.load_feature_to_gene() if f2g is None else f2g
    ids = pd.Index(pd.Series(protein_ids).astype(str))
    m = f2g[f2g["feature_ID"].isin(ids)][["feature_ID", "gene_symbol", "ensembl_gene", "entrez_gene"]].drop_duplicates()
    return m.rename(columns={"feature_ID": "protein_id"})
