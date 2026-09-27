"""Loaders for the exported MoTrPAC tables (see docs/DATA_GUIDE.md for the layout).

Everything returns an `OmicsMatrix`: samples × features with sample metadata that always
carries `pid` (animal) so downstream code can group splits correctly.
"""
from __future__ import annotations

import json
import re
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

from . import config as C


# ---------------------------------------------------------------------------------------
# Container
# ---------------------------------------------------------------------------------------
@dataclass
class OmicsMatrix:
    X: pd.DataFrame            # samples (index=viallabel) × features (columns=feature_ID)
    meta: pd.DataFrame         # index=viallabel; columns include pid, bid, sex, group, tissue, assay
    features: pd.DataFrame     # index=feature_ID; feature (regulated flag), assay, gene_symbol (if mapped)
    assay: str
    name: str = ""
    notes: list[str] = field(default_factory=list)

    @property
    def n_samples(self) -> int:
        return self.X.shape[0]

    @property
    def n_features(self) -> int:
        return self.X.shape[1]

    @property
    def n_animals(self) -> int:
        return int(self.meta["pid"].nunique())

    def groups(self) -> np.ndarray:
        return self.meta["pid"].to_numpy()

    def y(self, col: str) -> np.ndarray:
        return self.meta[col].to_numpy()

    def subset(self, mask) -> "OmicsMatrix":
        mask = np.asarray(mask, dtype=bool)
        return OmicsMatrix(self.X.loc[mask], self.meta.loc[mask], self.features, self.assay,
                           self.name, list(self.notes))

    def select_features(self, feature_ids: Iterable[str]) -> "OmicsMatrix":
        ids = [f for f in feature_ids if f in self.X.columns]
        return OmicsMatrix(self.X[ids], self.meta, self.features.loc[ids], self.assay,
                           self.name, list(self.notes))

    def __repr__(self) -> str:
        return (f"OmicsMatrix({self.name or self.assay}: {self.n_samples} samples, "
                f"{self.n_features} features, {self.n_animals} animals)")


# ---------------------------------------------------------------------------------------
# Low-level readers
# ---------------------------------------------------------------------------------------
_FILE_RE = re.compile(r"^([A-Z]+)__([A-Z]+)\.csv$")


def list_sample_files(kind: str = "norm", raw_dir: Path | None = None) -> pd.DataFrame:
    """Scan data/raw/<kind>/ for <ASSAY>__<TISSUE>.csv files."""
    raw_dir = Path(raw_dir) if raw_dir else C.RAW_DIR
    d = raw_dir / kind
    rows = []
    if d.exists():
        for p in sorted(d.glob("*.csv")):
            m = _FILE_RE.match(p.name)
            if not m:
                continue
            assay, tok = m.group(1), m.group(2)
            rows.append({"assay": assay, "tissue_token": tok, "tissue": C.tissue_from_token(tok),
                         "path": str(p), "kind": kind})
    return pd.DataFrame(rows, columns=["assay", "tissue_token", "tissue", "path", "kind"])


def _dedupe_index(ids: pd.Series) -> pd.Index:
    """Make feature IDs unique (append __dupN) without dropping anything."""
    ids = ids.astype(str)
    if ids.is_unique:
        return pd.Index(ids)
    counts = {}
    out = []
    for v in ids:
        n = counts.get(v, 0)
        out.append(v if n == 0 else f"{v}__dup{n}")
        counts[v] = n + 1
    return pd.Index(out)


def read_sample_table(path: str | Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Read a features × samples CSV. Returns (feature_table indexed by feature_ID, values)."""
    df = pd.read_csv(path, dtype={c: str for c in C.LEADING_COLS}, low_memory=False)
    lead = [c for c in C.LEADING_COLS if c in df.columns]
    if "feature_ID" not in lead:
        raise ValueError(f"{path}: no feature_ID column; columns start with {list(df.columns[:6])}")
    feat = df[lead].copy()
    vals = df.drop(columns=lead)
    vals.columns = vals.columns.astype(str)
    vals = vals.apply(pd.to_numeric, errors="coerce")
    idx = _dedupe_index(feat["feature_ID"])
    feat.index = idx
    vals.index = idx
    return feat, vals


# ---------------------------------------------------------------------------------------
# Phenotype table
# ---------------------------------------------------------------------------------------
def _first_present(df: pd.DataFrame, candidates: list[str]) -> str | None:
    for c in candidates:
        if c in df.columns:
            return c
    return None


def load_pheno(raw_dir: Path | None = None) -> pd.DataFrame:
    """PHENO with standardized columns: viallabel (index), pid, bid, sex, group, tissue_long.

    Column names in the package use either '.' or '___' separators depending on version;
    both are handled. Sex is 'male'/'female'; group is control/1w/2w/4w/8w.
    """
    raw_dir = Path(raw_dir) if raw_dir else C.RAW_DIR
    path = raw_dir / "pheno.csv"
    ph = pd.read_csv(path, dtype=str, low_memory=False)
    if "viallabel" not in ph.columns:
        raise ValueError("pheno.csv has no viallabel column")
    ph["viallabel"] = ph["viallabel"].str.strip()
    ph["bid"] = ph["bid"].str.strip() if "bid" in ph.columns else ph["viallabel"].str[:5]
    if "pid" not in ph.columns:
        raise ValueError("pheno.csv has no pid column")
    ph["pid"] = ph["pid"].str.strip()

    # sex
    sex_col = _first_present(ph, ["sex", "registration___sex", "registration.sex"])
    if sex_col is None:
        raise ValueError("no sex column in pheno")
    sex = ph[sex_col].str.strip().str.lower()
    sex = sex.replace({"1": "female", "2": "male", "f": "female", "m": "male"})
    ph["sex"] = sex

    # group / time point
    grp_col = _first_present(ph, ["group", "study_group_timepoint", "sacrificetime"])
    if grp_col is None:
        raise ValueError("no group column in pheno")
    g = ph[grp_col].astype(str).str.strip().str.lower()
    if grp_col != "group":
        # 'control - 1w' -> 'control', 'training - 8w' -> '8w'
        def _parse(v: str) -> str:
            if "control" in v:
                return "control"
            m = re.search(r"(\d)w", v)
            return f"{m.group(1)}w" if m else v
        g = g.map(_parse)
    ph["group"] = g

    tl = _first_present(ph, ["specimen_processing___sampletypedescription",
                             "specimen.processing.sampletypedescription", "tissue"])
    ph["tissue_long"] = ph[tl] if tl else np.nan

    ph = ph.drop_duplicates("viallabel").set_index("viallabel")
    return ph


# ---------------------------------------------------------------------------------------
# Sample metadata assembly
# ---------------------------------------------------------------------------------------
META_COLS = ["pid", "bid", "sex", "group", "tissue", "assay"]


def sample_meta_for(viallabels: Iterable[str], pheno: pd.DataFrame, tissue: str, assay: str) -> pd.DataFrame:
    """Join sample columns (viallabels) to PHENO; fall back to bid-level join if needed."""
    vl = pd.Index([str(v) for v in viallabels], name="viallabel")
    meta = pd.DataFrame(index=vl)
    hit = pheno.reindex(vl)
    meta["pid"] = hit["pid"].values
    meta["bid"] = hit["bid"].values
    meta["sex"] = hit["sex"].values
    meta["group"] = hit["group"].values
    missing = meta["pid"].isna()
    if missing.any():
        # bid is the first 5 characters of the viallabel; PHENO has one row per vial, so
        # any bid should resolve to an animal even if this exact vial is absent.
        by_bid = pheno.drop_duplicates("bid").set_index("bid")
        bids = vl[missing].str[:5]
        fb = by_bid.reindex(bids)
        meta.loc[missing, "pid"] = fb["pid"].values
        meta.loc[missing, "bid"] = bids.values
        meta.loc[missing, "sex"] = fb["sex"].values
        meta.loc[missing, "group"] = fb["group"].values
    still = meta["pid"].isna().sum()
    if still:
        warnings.warn(f"{tissue}/{assay}: {still} samples have no PHENO match and will be dropped")
    meta["tissue"] = tissue
    meta["assay"] = assay
    return meta


# ---------------------------------------------------------------------------------------
# Matrix loaders
# ---------------------------------------------------------------------------------------
def load_norm(assay: str, tissue: str, pheno: pd.DataFrame | None = None,
              raw_dir: Path | None = None) -> OmicsMatrix:
    """Normalized sample-level data for one assay × tissue."""
    raw_dir = Path(raw_dir) if raw_dir else C.RAW_DIR
    pheno = load_pheno(raw_dir) if pheno is None else pheno
    path = raw_dir / "norm" / f"{assay}__{C.tissue_token(tissue)}.csv"
    feat, vals = read_sample_table(path)
    meta = sample_meta_for(vals.columns, pheno, C.tissue_from_token(C.tissue_token(tissue)), assay)
    keep = ~meta["pid"].isna()
    X = vals.T.loc[keep.values]
    X.index.name = "viallabel"
    features = feat.drop(columns=["feature_ID"], errors="ignore")
    features["assay"] = assay
    return OmicsMatrix(X, meta.loc[keep], features, assay, name=f"{assay}/{tissue}")


def _read_counts_table(tissue: str, pheno: pd.DataFrame, assay: str = "TRNSCRPT",
                       raw_dir: Path | None = None) -> OmicsMatrix:
    """Raw counts for one tissue as an OmicsMatrix (samples × genes, NaN → 0); unfiltered, not logged."""
    raw_dir = Path(raw_dir) if raw_dir else C.RAW_DIR
    path = raw_dir / "counts" / f"{assay}__{C.tissue_token(tissue)}.csv"
    feat, counts = read_sample_table(path)
    counts = counts.fillna(0.0)
    meta = sample_meta_for(counts.columns, pheno, C.tissue_from_token(C.tissue_token(tissue)), assay)
    keep = ~meta["pid"].isna()
    X = counts.T.loc[keep.values]
    X.index.name = "viallabel"
    features = feat.drop(columns=["feature_ID"], errors="ignore")
    features["assay"] = assay
    return OmicsMatrix(X, meta.loc[keep], features, assay, name=f"{assay}-rawcounts/{tissue}")


def log_cpm(counts: pd.DataFrame, log: bool = True) -> pd.DataFrame:
    """samples × genes counts → CPM (library size = row sum over the genes given) → log2(CPM + 1)."""
    lib = counts.sum(axis=1).replace(0, np.nan)
    cpm = counts.div(lib, axis=0) * 1e6
    return np.log2(cpm + 1.0) if log else cpm


def load_counts(tissue: str, pheno: pd.DataFrame | None = None, assay: str = "TRNSCRPT",
                min_count: int = 10, min_frac: float = 0.2, log: bool = True,
                raw_dir: Path | None = None, filter_genes: bool = True) -> OmicsMatrix:
    """Raw RNA-seq counts → filtered log2 CPM for ONE tissue (uniform across tissues, unlike *_NORM_DATA).

    filter_genes=True keeps genes with >= min_count counts in >= min_frac of this tissue's samples
    (library size is the sum over the kept genes). For cross-tissue matrices use
    stack_tissues(source="counts"), whose default `counts_filter="stacked"` keeps every gene that
    passes this rule in at least one tissue, so tissue-restricted genes are not lost.
    """
    raw_dir = Path(raw_dir) if raw_dir else C.RAW_DIR
    pheno = load_pheno(raw_dir) if pheno is None else pheno
    raw = _read_counts_table(tissue, pheno, assay, raw_dir)
    counts = raw.X  # samples × genes
    if filter_genes:
        expressed = (counts >= min_count).mean(axis=0) >= min_frac
    else:
        expressed = pd.Series(True, index=counts.columns)
    counts = counts.loc[:, expressed.to_numpy()]
    vals = log_cpm(counts, log=log)
    om = OmicsMatrix(vals, raw.meta, raw.features.loc[vals.columns], assay, name=f"{assay}-counts/{tissue}")
    om.notes.append(f"counts filtered: {int(expressed.sum())}/{len(expressed)} features with "
                    f">= {min_count} counts in >= {min_frac:.0%} of samples; log2 CPM"
                    if filter_genes else f"counts unfiltered: {len(expressed)} features; log2 CPM")
    return om


def stack_tissues(assay: str, tissues: list[str] | None = None, source: str = "norm",
                  join: str = "inner", min_present: float = 0.8, pheno: pd.DataFrame | None = None,
                  raw_dir: Path | None = None, verbose: bool = True, counts_filter: str = "stacked",
                  min_count: int = 10, min_frac: float = 0.2, complete: bool = False,
                  drop_incomplete_samples: float | None = None) -> OmicsMatrix:
    """Concatenate one assay across tissues into one samples × features matrix.

    join='inner': features present in every tissue. join='outer': union, then drop features
    present in < min_present of samples. Always logs the sizes — look at them (DATA_GUIDE §6.2).

    source='counts' with counts_filter='stacked' (default): the raw counts of every tissue are
    stacked first (every counts table carries the full gene annotation, so an unexpressed gene is
    0 counts, never NaN), library sizes are taken over all genes, and a gene is kept if it has
    >= min_count counts in >= min_frac of the samples of AT LEAST ONE tissue. Tissue-restricted
    genes therefore survive and nothing is imputed; `join` is not used on this path.
    counts_filter='per_tissue' is the old behaviour: each tissue filtered on its own, then joined.
    complete=True additionally drops every feature with a missing value in any sample after the
    join (logged), for assays where the missingness pattern encodes tissue (PROT plexes).
    """
    raw_dir = Path(raw_dir) if raw_dir else C.RAW_DIR
    pheno = load_pheno(raw_dir) if pheno is None else pheno
    files = list_sample_files(source, raw_dir)
    files = files[files["assay"] == assay]
    if tissues is not None:
        files = files[files["tissue"].isin(tissues)]
    if files.empty:
        raise FileNotFoundError(f"no {source} files for assay {assay} under {raw_dir}")
    if source == "counts" and counts_filter == "stacked":
        return _stack_counts(assay, files, pheno, raw_dir, min_count, min_frac, verbose)
    if source == "counts" and counts_filter != "per_tissue":
        raise ValueError(f"counts_filter must be 'stacked' or 'per_tissue', got {counts_filter!r}")
    mats = []
    for _, r in files.iterrows():
        om = (load_counts(r["tissue"], pheno, assay=assay, raw_dir=raw_dir) if source == "counts"
              else load_norm(assay, r["tissue"], pheno, raw_dir=raw_dir))
        mats.append(om)
        if verbose:
            print(f"  loaded {om}")
            for note in om.notes:  # per-tissue filter counts (no silent feature drops)
                print(f"    {note}")
    if join == "inner":
        common = set(mats[0].X.columns)
        for om in mats[1:]:
            common &= set(om.X.columns)
        common = [c for c in mats[0].X.columns if c in common]
        X = pd.concat([om.X[common] for om in mats], axis=0)
    else:
        X = pd.concat([om.X for om in mats], axis=0, join="outer")
        present = X.notna().mean(axis=0)
        X = X.loc[:, present >= min_present]
    complete_note = ""
    meta = pd.concat([om.meta for om in mats], axis=0)
    n_empty = int(X.isna().all(axis=1).sum())
    if n_empty:
        complete_note += f"; WARNING {n_empty} samples have no value for any joined feature (imputed unless dropped)"
    if drop_incomplete_samples is not None:
        frac_nan = X.isna().mean(axis=1)
        keep_s = (frac_nan <= drop_incomplete_samples).to_numpy()
        dropped = meta.loc[~keep_s, "tissue"].value_counts().to_dict()
        X, meta = X.loc[keep_s], meta.loc[keep_s]
        complete_note += (f"; dropped {int((~keep_s).sum())} samples missing > {drop_incomplete_samples:.0%} of features "
                          f"({dropped})")
    if complete:
        n0 = X.shape[1]
        X = X.loc[:, X.notna().all(axis=0).to_numpy()]
        complete_note += f"; complete-features filter kept {X.shape[1]}/{n0} features with no missing value in any sample"
    features = pd.concat([om.features for om in mats], axis=0)
    features = features[~features.index.duplicated(keep="first")].reindex(X.columns)
    union = len(set().union(*[set(om.X.columns) for om in mats]))
    note = (f"stack {assay}/{source}: {len(mats)} tissues, {X.shape[0]} samples, "
            f"{X.shape[1]} features after join={join} (union {union}){complete_note}")
    if verbose:
        print("  " + note)
    out = OmicsMatrix(X, meta, features, assay, name=f"{assay}-{source}-stacked")
    out.notes.append(note)
    return out


def _stack_counts(assay: str, files: pd.DataFrame, pheno: pd.DataFrame, raw_dir: Path,
                  min_count: int, min_frac: float, verbose: bool) -> OmicsMatrix:
    """stack_tissues(source='counts', counts_filter='stacked'); see there."""
    mats = [_read_counts_table(r["tissue"], pheno, assay, raw_dir) for _, r in files.iterrows()]
    if verbose:
        for om in mats:
            print(f"  loaded {om} (raw counts)")
    counts = pd.concat([om.X for om in mats], axis=0, join="outer").fillna(0.0)
    masks = pd.concat([((om.X >= min_count).mean(axis=0) >= min_frac).rename(om.meta["tissue"].iloc[0])
                       for om in mats], axis=1).reindex(counts.columns).fillna(False)
    keep_any, keep_all = masks.any(axis=1), masks.all(axis=1)
    vals = log_cpm(counts, log=True).loc[:, keep_any.to_numpy()]
    meta = pd.concat([om.meta for om in mats], axis=0)
    features = pd.concat([om.features for om in mats], axis=0)
    features = features[~features.index.duplicated(keep="first")].reindex(vals.columns)
    n_any, n_all = int(keep_any.sum()), int(keep_all.sum())
    note = (f"stack {assay}/counts (stacked filter): {len(mats)} tissues, {vals.shape[0]} samples, "
            f"{n_any} genes with >= {min_count} counts in >= {min_frac:.0%} of samples in >= 1 tissue "
            f"(in every tissue: {n_all}; recovered {n_any - n_all} tissue-restricted genes that the per-tissue "
            f"inner join drops; annotation union {counts.shape[1]}); absent = 0 counts, log2 CPM on total library size")
    if verbose:
        print("  " + note)
    out = OmicsMatrix(vals, meta, features, assay, name=f"{assay}-counts-stacked")
    out.notes.append(note)
    out.notes.append("genes_expressed_in_n_tissues:" + ",".join(str(int(v)) for v in masks.sum(axis=1).loc[vals.columns]))
    return out


def align_by_animal(mats: dict[str, OmicsMatrix], key: str = "bid") -> dict[str, OmicsMatrix]:
    """Restrict several matrices to the same animals/biospecimens, in the same order.

    key='bid' for different assays on the same tissue; key='pid' across tissues.
    If an assay has more than one vial per key, values are averaged and a note is added.
    """
    keys = None
    for om in mats.values():
        k = set(om.meta[key].dropna())
        keys = k if keys is None else keys & k
    keys = sorted(keys)
    out = {}
    for name, om in mats.items():
        m = om.meta[om.meta[key].isin(keys)]
        Xk = om.X.loc[m.index].copy()
        Xk[key] = m[key].values
        dup = Xk[key].duplicated().sum()
        Xa = Xk.groupby(key).mean(numeric_only=True).reindex(keys)
        meta = m.drop_duplicates(key).set_index(key).reindex(keys)
        meta.index.name = key
        new = OmicsMatrix(Xa, meta, om.features, om.assay, name=om.name, notes=list(om.notes))
        if dup:
            new.notes.append(f"{name}: averaged {dup} duplicate vials per {key}")
        new.notes.append(f"aligned on {key}: {len(keys)} shared")
        out[name] = new
    return out


# ---------------------------------------------------------------------------------------
# Annotation tables
# ---------------------------------------------------------------------------------------
def load_manifest(raw_dir: Path | None = None) -> dict:
    raw_dir = Path(raw_dir) if raw_dir else C.RAW_DIR
    p = raw_dir / "manifest.json"
    return json.loads(p.read_text()) if p.exists() else {}


def load_feature_to_gene(raw_dir: Path | None = None) -> pd.DataFrame:
    raw_dir = Path(raw_dir) if raw_dir else C.RAW_DIR
    return pd.read_csv(raw_dir / "feature_to_gene.csv", dtype=str, low_memory=False)


def map_to_gene_symbols(feature_ids: Iterable[str], f2g: pd.DataFrame | None = None,
                        raw_dir: Path | None = None) -> pd.Series:
    """feature_ID -> gene_symbol (first mapping wins; unmapped -> NaN)."""
    f2g = load_feature_to_gene(raw_dir) if f2g is None else f2g
    sym_col = _first_present(f2g, ["gene_symbol", "symbol", "gene"])
    if sym_col is None:
        raise ValueError("feature_to_gene.csv has no gene_symbol column")
    m = f2g.dropna(subset=["feature_ID"]).drop_duplicates("feature_ID").set_index("feature_ID")[sym_col]
    ids = pd.Index([str(f).split("__dup")[0] for f in feature_ids])
    return pd.Series(m.reindex(ids).values, index=list(feature_ids), name="gene_symbol")


def load_rat_to_human(raw_dir: Path | None = None) -> pd.DataFrame:
    raw_dir = Path(raw_dir) if raw_dir else C.RAW_DIR
    return pd.read_csv(raw_dir / "rat_to_human_gene.csv", dtype=str, low_memory=False)


def load_da(assay: str, tissue: str, raw_dir: Path | None = None) -> pd.DataFrame:
    """Differential-analysis table for one assay × tissue, with a unified `stat` column
    (zscore for DESeq2 tables, tscore for limma tables, else logFC)."""
    raw_dir = Path(raw_dir) if raw_dir else C.RAW_DIR
    p = raw_dir / "da" / f"{assay}__{C.tissue_token(tissue)}.csv"
    da = pd.read_csv(p, low_memory=False, dtype={"feature_ID": str, "feature": str})
    stat = _first_present(da, ["zscore", "tscore", "logFC"])
    if stat is None:
        raise ValueError(f"{p}: no zscore/tscore/logFC column")
    da["stat"] = pd.to_numeric(da[stat], errors="coerce")
    pcol = _first_present(da, ["adj_p_value", "selection_fdr", "p_value"])
    da["padj"] = pd.to_numeric(da[pcol], errors="coerce") if pcol else np.nan
    for c in ("sex", "comparison_group"):
        if c in da.columns:
            da[c] = da[c].astype(str).str.lower().str.strip()
    return da


def load_outliers(raw_dir: Path | None = None) -> pd.DataFrame | None:
    raw_dir = Path(raw_dir) if raw_dir else C.RAW_DIR
    p = raw_dir / "outliers.csv"
    return pd.read_csv(p, dtype=str) if p.exists() else None
