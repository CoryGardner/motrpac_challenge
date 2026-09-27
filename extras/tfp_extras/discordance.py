"""Track 2 helpers: transcript–protein (or any two-layer) discordance.

Three complementary definitions, all computed within a tissue:
  1. sample-level: per-gene correlation between layers across matched animals
  2. DA-level: disagreement of differential statistics (sign/significance) per sex × time point
  3. model-level: residuals of protein ~ transcript per gene (how much protein is unexplained)
"""
from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
from scipy import stats

from tfp import config as C
from tfp.io import OmicsMatrix, align_by_animal, load_da, load_norm, map_to_gene_symbols


def gene_level(om: OmicsMatrix, symbols: pd.Series) -> pd.DataFrame:
    """Collapse features to gene symbols (mean of duplicates). Returns samples × genes."""
    sym = symbols.reindex(om.X.columns)
    keep = sym.notna()
    X = om.X.loc[:, keep.values]
    X.columns = sym[keep].values
    return X.T.groupby(level=0).mean().T


def channel_map(assay: str = "PROT", raw_dir=None) -> pd.Series | None:
    """viallabel → TMT channel from data/raw/meta/<ASSAY>.csv (None when absent)."""
    p = (raw_dir or C.RAW_DIR) / "meta" / f"{assay}.csv"
    if not p.exists():
        return None
    m = pd.read_csv(p, dtype=str).drop_duplicates("viallabel").set_index("viallabel")
    return m["tmt11_channel"] if "tmt11_channel" in m.columns else None


def residualize_channel(om: OmicsMatrix, channel: pd.Series) -> OmicsMatrix:
    """Unsupervised batch removal: subtract each feature's per-channel mean within the tissue (no
    labels used). Samples without a channel entry are centred on the overall mean."""
    ch = channel.reindex(om.X.index).fillna("unknown").to_numpy()
    X = om.X.copy()
    for c in np.unique(ch):
        m = ch == c
        X.loc[m] = X.loc[m] - X.loc[m].mean(axis=0)
    out = OmicsMatrix(X, om.meta, om.features, om.assay, om.name, list(om.notes))
    out.notes.append(f"channel residualized: {len(np.unique(ch))} channels")
    return out


def paired_layers(tissue: str, assay_a: str = "TRNSCRPT", assay_b: str = "PROT",
                  pheno=None, f2g=None, raw_dir=None, complete_b: bool = False,
                  residualize_b: pd.Series | None = None) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Gene-level matrices for two assays aligned on biospecimen (bid). Returns (A, B, meta).
    complete_b drops layer-b features with any missing value in this tissue; residualize_b (a
    viallabel → channel map) subtracts per-channel means from layer b before the gene collapse."""
    a = load_norm(assay_a, tissue, pheno, raw_dir)
    b = load_norm(assay_b, tissue, pheno, raw_dir)
    if complete_b:
        keep = b.X.notna().all(axis=0)
        n0 = b.n_features
        b = b.select_features(list(b.X.columns[keep.to_numpy()]))
        b.notes.append(f"{assay_b}/{tissue}: {b.n_features}/{n0} features with no missing value")
    if residualize_b is not None:
        b = residualize_channel(b, residualize_b)
    aligned = align_by_animal({assay_a: a, assay_b: b}, key="bid")
    a, b = aligned[assay_a], aligned[assay_b]
    ga = gene_level(a, map_to_gene_symbols(a.X.columns, f2g, raw_dir))
    gb = gene_level(b, map_to_gene_symbols(b.X.columns, f2g, raw_dir))
    genes = sorted(set(ga.columns) & set(gb.columns))
    return ga[genes], gb[genes], a.meta


def per_gene_correlation(A: pd.DataFrame, B: pd.DataFrame, method: str = "spearman",
                         min_n: int = 8) -> pd.DataFrame:
    rows = []
    for g in A.columns:
        x, y = A[g].to_numpy(dtype=float), B[g].to_numpy(dtype=float)
        ok = ~(np.isnan(x) | np.isnan(y))
        if ok.sum() < min_n:
            rows.append({"gene": g, "rho": np.nan, "p": np.nan, "n": int(ok.sum())})
            continue
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")  # constant inputs → NaN, reported as such
            if method == "spearman":
                r, p = stats.spearmanr(x[ok], y[ok])
            else:
                r, p = stats.pearsonr(x[ok], y[ok])
        rows.append({"gene": g, "rho": float(r), "p": float(p), "n": int(ok.sum())})
    return pd.DataFrame(rows)


def da_discordance(tissue: str, assay_a: str = "TRNSCRPT", assay_b: str = "PROT",
                   fdr: float = 0.1, f2g=None, raw_dir=None) -> pd.DataFrame:
    """Per gene × sex × comparison_group: statistics from both layers and discordance flags.

    discordant_sign: both layers significant (padj < fdr) with opposite signs
    discordant_sig : one layer significant, the other not (|stat| small)
    """
    da_a, da_b = load_da(assay_a, tissue, raw_dir), load_da(assay_b, tissue, raw_dir)
    for da in (da_a, da_b):
        da["gene"] = map_to_gene_symbols(da["feature_ID"], f2g, raw_dir).values
    keys = ["gene", "sex", "comparison_group"]
    agg = {"stat": "mean", "padj": "min"}
    A = da_a.dropna(subset=["gene"]).groupby(keys).agg(agg).rename(columns={"stat": "stat_a", "padj": "padj_a"})
    B = da_b.dropna(subset=["gene"]).groupby(keys).agg(agg).rename(columns={"stat": "stat_b", "padj": "padj_b"})
    m = A.join(B, how="inner").reset_index()
    sig_a, sig_b = m["padj_a"] < fdr, m["padj_b"] < fdr
    m["sig_a"], m["sig_b"] = sig_a, sig_b
    m["discordant_sign"] = sig_a & sig_b & (np.sign(m["stat_a"]) != np.sign(m["stat_b"]))
    m["discordant_sig"] = sig_a ^ sig_b
    m["concordant"] = sig_a & sig_b & (np.sign(m["stat_a"]) == np.sign(m["stat_b"]))
    return m


def ptm_site_counts(tissue: str, assay: str, f2g=None, raw_dir=None) -> pd.Series:
    """Number of quantified PTM sites per gene in this tissue (0 when the assay has no table here).
    PTM feature IDs look like NP_001030329.2_S14s: the protein accession is everything before the
    last underscore; it maps to a gene through feature_to_gene."""
    from .io import read_sample_table
    p = (raw_dir or C.RAW_DIR) / "norm" / f"{assay}__{C.tissue_token(tissue)}.csv"
    if not p.exists():
        return pd.Series(dtype=float, name=f"n_{assay.lower()}_sites")
    feat, _ = read_sample_table(p)
    acc = pd.Series([str(f).rsplit("_", 1)[0] for f in feat.index])
    genes = map_to_gene_symbols(acc.unique(), f2g, raw_dir)
    g = acc.map(genes.to_dict())
    return g.dropna().value_counts().rename(f"n_{assay.lower()}_sites")


def discordance_features(A: pd.DataFrame, B: pd.DataFrame, corr: pd.DataFrame,
                         regulated: set[str] | None = None, ptm: dict[str, pd.Series] | None = None) -> pd.DataFrame:
    """Per-gene covariates for a 'predict discordance' model: transcript mean/variance, protein
    mean/variance and abundance rank, ρ, training-regulated flag, PTM site counts."""
    df = pd.DataFrame({
        "gene": A.columns,
        "mean_a": A.mean(axis=0).values, "mean_b": B.mean(axis=0).values,
        "var_a": A.var(axis=0).values, "var_b": B.var(axis=0).values,
        "cv_a": (A.std(axis=0) / A.mean(axis=0).abs().replace(0, np.nan)).values,
        "cv_b": (B.std(axis=0) / B.mean(axis=0).abs().replace(0, np.nan)).values,
    }).merge(corr[["gene", "rho"]], on="gene", how="left")
    df["prot_abundance_rank"] = df["mean_b"].rank(pct=True)
    if regulated is not None:
        df["is_training_regulated"] = df["gene"].isin(regulated).astype(int)
    for name, counts in (ptm or {}).items():
        df[name] = df["gene"].map(counts).fillna(0).astype(float)
    return df


def predict_discordance(df: pd.DataFrame, target: str, features: list[str], n_splits: int = 5,
                        seed: int = C.SEED, importance: bool = False):
    """Gene-grouped CV AUROC for predicting a binary discordance label from gene covariates.
    importance=True also returns a per-covariate table: mean |standardized coefficient| and
    drop-one ΔAUROC for the logistic model, mean impurity importance for the forest."""
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.impute import SimpleImputer
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import roc_auc_score
    from sklearn.model_selection import GroupKFold
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import StandardScaler

    d = df.dropna(subset=[target]).copy()
    y = d[target].astype(int).to_numpy()
    X = d[features].to_numpy(dtype=float)
    g = d["gene"].astype(str).to_numpy()
    if len(np.unique(y)) < 2 or y.sum() < n_splits:
        empty = pd.DataFrame([{"model": "n/a", "fold": 0, "auroc": np.nan, "note": "too few positives"}])
        return (empty, pd.DataFrame()) if importance else empty
    models = {
        "logreg": Pipeline([("imp", SimpleImputer()), ("sc", StandardScaler()),
                            ("clf", LogisticRegression(max_iter=2000, class_weight="balanced"))]),
        "rf": Pipeline([("imp", SimpleImputer()),
                        ("clf", RandomForestClassifier(n_estimators=300, random_state=seed, class_weight="balanced", n_jobs=-1))]),
    }
    rows, imp_rows = [], []
    for fold, (tr, te) in enumerate(GroupKFold(n_splits=min(n_splits, len(np.unique(g)))).split(X, y, g)):
        if len(np.unique(y[te])) < 2:
            continue
        for name, m in models.items():
            m.fit(X[tr], y[tr])
            auc = roc_auc_score(y[te], m.predict_proba(X[te])[:, 1])
            rows.append({"model": name, "fold": fold, "auroc": auc, "n_test": len(te), "prevalence": float(y[te].mean())})
            if not importance:
                continue
            if name == "logreg":
                coef = m.named_steps["clf"].coef_.ravel()
                for j, f in enumerate(features):
                    keep = [i for i in range(len(features)) if i != j]
                    m1 = Pipeline([("imp", SimpleImputer()), ("sc", StandardScaler()),
                                   ("clf", LogisticRegression(max_iter=2000, class_weight="balanced"))]).fit(X[tr][:, keep], y[tr])
                    d_auc = auc - roc_auc_score(y[te], m1.predict_proba(X[te][:, keep])[:, 1])
                    imp_rows.append({"fold": fold, "covariate": f, "logreg_abs_coef": float(abs(coef[j])), "logreg_drop_one_dAUROC": float(d_auc)})
            else:
                for j, f in enumerate(features):
                    imp_rows.append({"fold": fold, "covariate": f, "rf_importance": float(m.named_steps["clf"].feature_importances_[j])})
    if importance:
        imp = pd.DataFrame(imp_rows)
        imp = imp.groupby("covariate").agg(logreg_abs_coef=("logreg_abs_coef", "mean"), logreg_drop_one_dAUROC=("logreg_drop_one_dAUROC", "mean"),
                                          rf_importance=("rf_importance", "mean")).reset_index() if len(imp) else imp
        return pd.DataFrame(rows), imp
    return pd.DataFrame(rows)
