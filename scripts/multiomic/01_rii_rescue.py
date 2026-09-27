#!/usr/bin/env python
"""Multiomic Phase 1 — the RII rescue: does MoTrPAC proteomics carry a tissue axis once it is put on a scale where
cross-tissue comparison is defined?

Reads the portal quant-id peptide-level reporter-ion intensities (tfp.rii) for the 7 proteomics tissues, rolls them up
to proteins per plex, normalises every channel to its total protein signal (log2 ppm), stacks the tissues, and reports:
  1. tissue R² of PC1–3 on the stacked matrix (same PCA + R² code as phase 03), against the same computation on the
     distributed ratio matrix and the frozen phase-03 numbers;
  2. the phase-04 diagnostic on all animal-grouped folds: model as fitted, missingness-indicators-only, per-tissue-mean-
     removed, against the frozen ratio-matrix diagnostic;
  3. a RoundRobinSelector protein panel curve (tfp.models.panel_curve, pipeline unchanged) with a label-permutation null
     at k = 20, confusion at k = 20/30, and the k = 20 panel selected on all animals;
  4. the cross-tissue RNA–protein Spearman correlation per gene on the same animals (tissue means over the 7 tissues),
     its distribution, a mismatched-pair null, and the value for every RNA panel gene whose protein is quantified, with
     the marker-tissue agreement that pre-registration (c) asks for.
Plex is nested in tissue here too (one plex = one tissue's samples + that tissue's reference pool), so every accuracy in
this phase is within-study evidence. Its job is to show the signal exists, not to certify it.

Outputs: results_multiomic/01_rii/ (CSVs, parquet matrices, README.md built from the CSVs).
"""
from __future__ import annotations

import argparse
import json
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.decomposition import PCA
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score
from sklearn.preprocessing import StandardScaler

from tfp import config as C, io, models, report, rii
from tfp.splits import assert_no_group_leak, grouped_kfold

ROOT = C.ROOT
OUT_ROOT = ROOT / "results_multiomic" / "01_rii"
FROZEN = C.FROZEN_DIR


# ---------------------------------------------------------------------------------------------------------------------
# helpers (phase-03 R² and PCA, copied so the numbers are comparable)
# ---------------------------------------------------------------------------------------------------------------------
def r2_categorical(scores: np.ndarray, cat) -> np.ndarray:
    cat = pd.Series(cat).astype(str).to_numpy()
    out = []
    for j in range(scores.shape[1]):
        s = scores[:, j]
        tot = np.var(s) * len(s)
        if tot == 0:
            out.append(0.0)
            continue
        between = sum(len(s[cat == c]) * (s[cat == c].mean() - s.mean()) ** 2 for c in np.unique(cat))
        out.append(float(between / tot))
    return np.array(out)


def pca_scores(X: np.ndarray, n_components: int = 10, top_var: int | None = 5000, seed: int = C.SEED):
    X = SimpleImputer(strategy="median").fit_transform(X)
    if top_var and X.shape[1] > top_var:
        v = X.var(axis=0)
        X = X[:, np.argsort(v)[::-1][:top_var]]
    X = StandardScaler().fit_transform(X)
    n = min(n_components, X.shape[0] - 1, X.shape[1])
    pca = PCA(n_components=n, random_state=seed)
    return pca.fit_transform(X), pca.explained_variance_ratio_


def variance_partition(om: io.OmicsMatrix, covs: list[str], seed: int, n_perm: int, rng) -> pd.DataFrame:
    S, evr = pca_scores(om.X.to_numpy(dtype=float), seed=seed)
    vp = pd.DataFrame({"PC": [f"PC{i + 1}" for i in range(S.shape[1])], "explained": evr})
    for c in covs:
        if c in om.meta.columns:
            vp[f"R2_{c}"] = r2_categorical(S, om.meta[c])
    # permutation null for the tissue R²: shuffle tissue labels across samples, 95th percentile per PC
    if n_perm:
        null = np.array([r2_categorical(S, rng.permutation(om.meta["tissue"].astype(str).to_numpy())) for _ in range(n_perm)])
        vp["R2_tissue_null95"] = np.quantile(null, 0.95, axis=0)
    return vp


def centred_accuracy(X, y, tr, te, classes):
    """Phase-04 diagnostic: centre every tissue on its TRAINING-fold means (label-using), then classify."""
    Xi = SimpleImputer(strategy="median").fit(X[tr]).transform(X)
    Xc = Xi.copy()
    for t in classes:
        mu = Xi[tr][y[tr] == t].mean(axis=0)
        Xc[y == t] -= mu
    sc = StandardScaler().fit(Xc[tr])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cen = LogisticRegression(C=0.1, max_iter=3000).fit(sc.transform(Xc[tr]), y[tr])
    p = cen.predict(sc.transform(Xc[te]))
    return float(np.mean(p == y[te])), float(balanced_accuracy_score(y[te], p))


def missingness_accuracy(M, y, tr, te):
    if M.shape[1] == 0 or M[tr].std(axis=0).max() == 0:
        return np.nan, np.nan
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        miss = LogisticRegression(C=0.1, max_iter=3000).fit(M[tr], y[tr])
    p = miss.predict(M[te])
    return float(np.mean(p == y[te])), float(balanced_accuracy_score(y[te], p))


def fmt(v, nd=3):
    return "" if v is None or (isinstance(v, float) and np.isnan(v)) else f"{v:.{nd}f}"


# ---------------------------------------------------------------------------------------------------------------------
def main() -> None:
    ap = argparse.ArgumentParser(description="Multiomic Phase 1 — RII rescue")
    ap.add_argument("--assay", default="prot-pr", help="portal PTM folder: prot-pr (global), prot-ph, prot-ac")
    ap.add_argument("--release", default="c1.0")
    ap.add_argument("--min-peptides", type=int, default=2)
    ap.add_argument("--norm", default="total", choices=["total", "median"])
    ap.add_argument("--n-splits", type=int, default=5)
    ap.add_argument("--prefilter", type=int, default=5000)
    ap.add_argument("--grid", default=",".join(map(str, C.PANEL_GRID)))
    ap.add_argument("--seed", type=int, default=C.SEED)
    ap.add_argument("--n-perm", type=int, default=200, help="label permutations for the PC R² null")
    ap.add_argument("--n-perm-panel", type=int, default=20, help="label permutations for the k=20 panel null (quick fits)")
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--skip-rna", action="store_true")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    t_start = time.time()
    out = Path(args.out) if args.out else (OUT_ROOT if args.assay == "prot-pr" else OUT_ROOT / args.assay.replace("prot-", ""))
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)
    grid = [int(k) for k in args.grid.split(",")]
    n_splits = 3 if args.quick else args.n_splits
    print(f"=== multiomic 01_rii_rescue === {vars(args)}")

    # ---- 1. build the stacked RII matrices ----------------------------------------------------------------------------
    pheno = io.load_pheno()
    om_outer, plexes = rii.stack_rii(rii.PROT_TISSUES, args.assay, args.release, pheno, join="outer",
                                     min_peptides=args.min_peptides, norm=args.norm)
    common = [c for c in om_outer.X.columns if om_outer.X[c].notna().groupby(om_outer.meta["tissue"].to_numpy()).any().all()]
    om = io.OmicsMatrix(om_outer.X[common], om_outer.meta, om_outer.features.loc[common], om_outer.assay,
                        name=om_outer.name.replace("outer", "inner"), notes=list(om_outer.notes))
    om.notes.append(f"inner: {len(common)} proteins quantified (>= 1 vial) in every tissue; NaN fraction {om.X.isna().mean().mean():.3f}")
    print("  " + om.notes[-1])
    plexes.to_csv(out / "plex_summary.csv", index=False)
    per_t = om_outer.meta.groupby("tissue").agg(n_vials=("pid", "size"), n_animals=("pid", "nunique"), n_plexes=("plex_id", "nunique"))
    per_t["n_proteins_any_vial"] = [int(om_outer.X.loc[om_outer.meta["tissue"] == t].notna().any(axis=0).sum()) for t in per_t.index]
    per_t["frac_nan_within_tissue_any_protein"] = [float(om_outer.X.loc[om_outer.meta["tissue"] == t, om_outer.X.loc[om_outer.meta["tissue"] == t].notna().any(axis=0)].isna().mean().mean()) for t in per_t.index]
    per_t.reset_index().to_csv(out / "tissue_summary.csv", index=False)
    join = pd.DataFrame([{"n_tissues": om.meta["tissue"].nunique(), "n_vials": om.n_samples, "n_animals": om.n_animals,
                          "n_proteins_union": om_outer.n_features, "n_proteins_inner": om.n_features,
                          "frac_nan_outer": float(om_outer.X.isna().mean().mean()), "frac_nan_inner": float(om.X.isna().mean().mean()),
                          "n_proteins_inner_complete": int(om.X.notna().all(axis=0).sum()),
                          "min_peptides": args.min_peptides, "norm": args.norm, "release": args.release, "assay": args.assay}])
    join.to_csv(out / "join_summary.csv", index=False)
    om.X.to_parquet(out / "rii_inner_log2ppm.parquet")
    om_outer.X.to_parquet(out / "rii_outer_log2ppm.parquet")
    om.meta.to_csv(out / "rii_meta.csv")
    p2g = rii.protein_to_ensembl(om_outer.X.columns)
    feats = om_outer.features.copy()
    feats["ensembl_gene"] = feats.index.map(p2g.drop_duplicates("protein_id").set_index("protein_id")["ensembl_gene"])
    feats["in_inner"] = feats.index.isin(om.X.columns)
    feats.to_csv(out / "rii_features.csv")

    # ---- 2. variance partition: RII vs ratios --------------------------------------------------------------------------
    covs = ["tissue", "sex", "group", "pid", "plex_id", "tmt11_channel"]
    vp = variance_partition(om, covs, args.seed, args.n_perm if not args.quick else 20, rng)
    vp.to_csv(out / "variance_partition.csv", index=False)
    om_ratio = io.stack_tissues("PROT", source="norm", join="inner", pheno=pheno, verbose=False)
    meta_r = om_ratio.meta.copy()
    pm = pd.read_csv(C.META_DIR / "PROT.csv", dtype=str).drop_duplicates("viallabel").set_index("viallabel")
    meta_r["tmt11_channel"] = pm["tmt11_channel"].reindex(meta_r.index).to_numpy()
    meta_r["plex_id"] = meta_r["tissue"] + ":" + pm["tmt_plex"].reindex(meta_r.index).astype(str).to_numpy()
    om_ratio = io.OmicsMatrix(om_ratio.X, meta_r, om_ratio.features, "PROT", name="PROT-ratio")
    vp_r = variance_partition(om_ratio, covs, args.seed, 0, rng)
    vp_r.to_csv(out / "variance_partition_ratio.csv", index=False)
    vp_med = None
    if args.norm == "total" and not args.quick:
        om_med, _ = rii.stack_rii(rii.PROT_TISSUES, args.assay, args.release, pheno, join="inner", min_peptides=args.min_peptides,
                                  norm="median", verbose=False)
        vp_med = variance_partition(om_med, ["tissue"], args.seed, 0, rng)
        vp_med.to_csv(out / "variance_partition_median_norm.csv", index=False)
    print(vp[["PC", "explained", "R2_tissue", "R2_plex_id", "R2_sex", "R2_pid"]].head(3).round(4).to_string(index=False))

    # ---- 3. diagnostics on every animal-grouped fold -----------------------------------------------------------------
    X = om.X.to_numpy(dtype=float)
    y = om.meta["tissue"].astype(str).to_numpy()
    g = om.groups()
    classes = sorted(np.unique(y))
    Xo = om_outer.X.loc[om.X.index].to_numpy(dtype=float)
    M_outer = np.isnan(Xo).astype(float)
    M_outer = M_outer[:, M_outer.std(axis=0) > 0]
    M_inner = np.isnan(X).astype(float)
    M_inner = M_inner[:, M_inner.std(axis=0) > 0]
    folds = list(grouped_kfold(om.meta, "tissue", n_splits, args.seed))
    rows = []
    for f, (tr, te) in enumerate(folds):
        assert_no_group_leak(om.meta, tr, te)
        t0 = time.time()
        est = models.fit_tuned("logreg_l2", X[tr], y[tr], g[tr], prefilter=args.prefilter, quick=args.quick, seed=args.seed)
        p = est.predict(X[te])
        r = {"fold": f, "n_train_animals": int(om.meta.iloc[tr]["pid"].nunique()), "n_test_animals": int(om.meta.iloc[te]["pid"].nunique()),
             "n_test": len(te), "model_as_fitted_acc": float(np.mean(p == y[te])), "model_as_fitted_bal_acc": float(balanced_accuracy_score(y[te], p))}
        r["missingness_outer_acc"], r["missingness_outer_bal_acc"] = missingness_accuracy(M_outer, y, tr, te)
        r["missingness_inner_acc"], r["missingness_inner_bal_acc"] = missingness_accuracy(M_inner, y, tr, te)
        r["per_tissue_means_removed_acc"], r["per_tissue_means_removed_bal_acc"] = centred_accuracy(X, y, tr, te, classes)
        r["chance"] = 1.0 / len(classes)
        r["seconds"] = round(time.time() - t0, 1)
        rows.append(r)
        print(f"  fold {f}: fitted {r['model_as_fitted_acc']:.3f}  miss(outer) {fmt(r['missingness_outer_acc'])}  miss(inner) {fmt(r['missingness_inner_acc'])}  "
              f"means-removed {r['per_tissue_means_removed_acc']:.3f}  ({r['seconds']}s)")
    diag = pd.DataFrame(rows)
    diag.to_csv(out / "diagnostic_accuracy.csv", index=False)
    num = [c for c in diag.columns if c.endswith("_acc") or c == "chance"]
    ds = diag[num].agg(["mean", "std"]).T.reset_index().rename(columns={"index": "quantity", "mean": "mean", "std": "sd"})
    ds["n_folds"] = len(diag)
    ds["n_test_animals_mean"] = diag["n_test_animals"].mean()
    ds["n_missingness_features_outer"] = M_outer.shape[1]
    ds["n_missingness_features_inner"] = M_inner.shape[1]
    ds.to_csv(out / "diagnostic_accuracy_summary.csv", index=False)

    # ---- 4. panel curve (pipeline unchanged) + null + all-animal k=20 panel ---------------------------------------------
    curve, selected, preds = models.panel_curve(om, label="tissue", grid=grid, kind="logreg_l2", n_splits=n_splits,
                                                prefilter=args.prefilter, quick=args.quick, seed=args.seed, verbose=False,
                                                return_predictions=True)
    curve.to_csv(out / "panel_curve.csv", index=False)
    selected.to_csv(out / "selected_by_fold.csv", index=False)
    preds.to_csv(out / "panel_predictions.csv", index=False)
    agg = curve.groupby("k").agg(bal_acc_mean=("balanced_accuracy", "mean"), bal_acc_sd=("balanced_accuracy", "std"),
                                 acc_mean=("accuracy", "mean"), macro_f1_mean=("macro_f1", "mean"), n_folds=("fold", "nunique"),
                                 n_test_animals_mean=("n_test_animals", "mean")).reset_index()
    for kk in (20, 30):
        pk = preds[preds["k"] == kk]
        if pk.empty:
            continue
        conf = pd.crosstab(pk["y_true"], pk["y_pred"]).reindex(index=classes, columns=classes, fill_value=0)
        conf.to_csv(out / f"confusion_k{kk}.csv")
        cp = pd.DataFrame([{"k": kk, "true": a, "predicted": b, "count": int(conf.loc[a, b]), "frac_of_true": float(conf.loc[a, b] / max(conf.loc[a].sum(), 1))}
                           for a in classes for b in classes if a != b and conf.loc[a, b] > 0],
                          columns=["k", "true", "predicted", "count", "frac_of_true"]).sort_values("count", ascending=False)
        cp.to_csv(out / f"confusable_pairs_k{kk}.csv", index=False)
    # label-permutation null at k = 20 (quick fits: C fixed at 0.1; the null's job is the chance level)
    k_null = 20 if 20 in grid else grid[0]
    null_vals = []
    for i in range(args.n_perm_panel if not args.quick else 3):
        yp = rng.permutation(y)
        accs = []
        for tr, te in folds:
            est = models.fit_tuned("logreg_l2", X[tr], yp[tr], g[tr], k=k_null, prefilter=args.prefilter, quick=True, seed=args.seed)
            accs.append(balanced_accuracy_score(yp[te], est.predict(X[te])))
        null_vals.append(float(np.mean(accs)))
    obs = float(agg.loc[agg["k"] == k_null, "bal_acc_mean"].iloc[0])
    null = pd.DataFrame([{"k": k_null, "observed_bal_acc_mean": obs, "null_mean": float(np.mean(null_vals)), "null_sd": float(np.std(null_vals)),
                          "null_q95": float(np.quantile(null_vals, 0.95)), "null_max": float(np.max(null_vals)), "n_perm": len(null_vals),
                          "chance": 1.0 / len(classes)}])
    null.to_csv(out / "panel_null_k20.csv", index=False)
    agg["null_q95_k20"] = np.where(agg["k"] == k_null, null["null_q95"].iloc[0], np.nan)
    agg.to_csv(out / "panel_curve_summary.csv", index=False)
    print(agg[["k", "bal_acc_mean", "bal_acc_sd"]].round(3).to_string(index=False))
    # the k = 20 panel on all animals (what a transfer phase would carry), with marker tissue and gene
    pipe = models.make_pipeline("logreg_l2", k=k_null, prefilter=args.prefilter, seed=args.seed)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pipe.fit(X, y)
    chosen = models.selected_feature_ids(pipe, np.asarray(om.X.columns))
    Z = pd.DataFrame(StandardScaler().fit_transform(SimpleImputer(strategy="median").fit_transform(X)), index=om.X.index, columns=om.X.columns)
    mu = Z[chosen].groupby(y).mean()
    panel = pd.DataFrame({"protein_id": chosen, "gene_symbol": om.features.loc[chosen, "gene_symbol"].to_numpy(),
                          "ensembl_gene": feats.loc[chosen, "ensembl_gene"].to_numpy(),
                          "marker_tissue": mu.idxmax().to_numpy(), "mean_z_in_marker": mu.max().to_numpy(),
                          "next_highest_tissue": [mu[c].drop(mu[c].idxmax()).idxmax() for c in chosen],
                          "effect_z": [float(mu[c].max() - mu[c].drop(mu[c].idxmax()).max()) for c in chosen],
                          "frac_nan": om.X[chosen].isna().mean().to_numpy()})
    panel.to_csv(out / "panel_k20_all_animals.csv", index=False)

    # ---- 4b. the same headline numbers on the complete-protein matrix (no missing value anywhere) ---------------------
    # missingness alone identifies tissue on this scale too (per-tissue searches), so every headline number is repeated on
    # the proteins quantified in EVERY vial, where no imputation happens and the NaN pattern carries nothing.
    comp_cols = [c for c in om.X.columns if om.X[c].notna().all()]
    om_c = io.OmicsMatrix(om.X[comp_cols], om.meta, om.features.loc[comp_cols], om.assay, name="RII-complete")
    vp_c = variance_partition(om_c, ["tissue", "sex", "pid", "plex_id"], args.seed, 0, rng)
    vp_c.to_csv(out / "variance_partition_complete.csv", index=False)
    Xc_ = om_c.X.to_numpy(dtype=float)
    rows_c = []
    for f, (tr, te) in enumerate(folds):
        est = models.fit_tuned("logreg_l2", Xc_[tr], y[tr], g[tr], prefilter=args.prefilter, quick=args.quick, seed=args.seed)
        p = est.predict(Xc_[te])
        rc = {"fold": f, "n_features": len(comp_cols), "n_test_animals": int(om.meta.iloc[te]["pid"].nunique()),
              "model_as_fitted_acc": float(np.mean(p == y[te])), "model_as_fitted_bal_acc": float(balanced_accuracy_score(y[te], p))}
        rc["per_tissue_means_removed_acc"], rc["per_tissue_means_removed_bal_acc"] = centred_accuracy(Xc_, y, tr, te, classes)
        rows_c.append(rc)
    diag_c = pd.DataFrame(rows_c)
    diag_c.to_csv(out / "diagnostic_accuracy_complete.csv", index=False)
    curve_c, _ = models.panel_curve(om_c, label="tissue", grid=grid, kind="logreg_l2", n_splits=n_splits, prefilter=args.prefilter,
                                    quick=args.quick, seed=args.seed, verbose=False)
    curve_c.to_csv(out / "panel_curve_complete.csv", index=False)
    agg_c = curve_c.groupby("k").agg(bal_acc_mean=("balanced_accuracy", "mean"), bal_acc_sd=("balanced_accuracy", "std"),
                                     n_folds=("fold", "nunique")).reset_index()
    agg_c["n_features"] = len(comp_cols)
    agg_c.to_csv(out / "panel_curve_complete_summary.csv", index=False)
    print(f"  complete matrix: {len(comp_cols)} proteins; PC1 tissue R2 {vp_c['R2_tissue'].iloc[0]:.3f}; "
          f"fitted acc {diag_c['model_as_fitted_acc'].mean():.3f}; k=20 bal acc {agg_c.loc[agg_c['k'] == k_null, 'bal_acc_mean'].iloc[0] if k_null in set(agg_c['k']) else float('nan'):.3f}")

    # ---- 5. RNA–protein cross-tissue correlation on the same animals ------------------------------------------------
    rp_summary = pd.DataFrame(); rp_panel = pd.DataFrame(); rp_panel_summary = pd.DataFrame()
    if not args.skip_rna:
        rna = io.stack_tissues("TRNSCRPT", tissues=rii.PROT_TISSUES, source="counts", pheno=pheno, verbose=False)
        Xp = om_outer.X
        mp = om_outer.meta
        rna_means, prot_means, n_shared = {}, {}, {}
        for t in rii.PROT_TISSUES:
            pr = set(rna.meta.loc[rna.meta["tissue"] == t, "pid"]) & set(mp.loc[mp["tissue"] == t, "pid"])
            n_shared[t] = len(pr)
            rna_means[t] = rna.X[rna.meta["tissue"].eq(t).to_numpy() & rna.meta["pid"].isin(pr).to_numpy()].mean(axis=0)
            sub = Xp[mp["tissue"].eq(t).to_numpy() & mp["pid"].isin(pr).to_numpy()]
            m = sub.mean(axis=0, skipna=True)
            m[sub.notna().sum(axis=0) < 3] = np.nan       # a tissue mean needs >= 3 quantified animals
            prot_means[t] = m
        R = pd.DataFrame(rna_means)     # genes × 7
        P = pd.DataFrame(prot_means)    # proteins × 7
        P = P.dropna(axis=0)
        pairs = p2g[p2g["protein_id"].isin(P.index) & p2g["ensembl_gene"].isin(R.index)].drop_duplicates(["protein_id", "ensembl_gene"])
        # gene level: average the log2 ppm tissue means of the proteins that map to the gene
        Pg = P.loc[pairs["protein_id"]].groupby(pairs["ensembl_gene"].to_numpy()).mean()
        n_prot = pairs.groupby("ensembl_gene")["protein_id"].agg(lambda s: ";".join(sorted(set(s))))
        Rg = R.loc[Pg.index]
        rho = np.array([spearmanr(Rg.loc[i].to_numpy(), Pg.loc[i].to_numpy()).statistic for i in Pg.index])
        pear = np.array([np.corrcoef(Rg.loc[i].to_numpy(), Pg.loc[i].to_numpy())[0, 1] for i in Pg.index])
        sym = io.map_to_gene_symbols(Pg.index)
        rp = pd.DataFrame({"ensembl_gene": Pg.index, "gene_symbol": sym.to_numpy(), "protein_ids": n_prot.reindex(Pg.index).to_numpy(),
                           "spearman": rho, "pearson": pear, "rna_marker_tissue": Rg.idxmax(axis=1).to_numpy(),
                           "prot_marker_tissue": Pg.idxmax(axis=1).to_numpy(), "rna_range_log2cpm": (Rg.max(axis=1) - Rg.min(axis=1)).to_numpy(),
                           "prot_range_log2ppm": (Pg.max(axis=1) - Pg.min(axis=1)).to_numpy()})
        rp["same_marker_tissue"] = rp["rna_marker_tissue"] == rp["prot_marker_tissue"]
        rp.sort_values("spearman", ascending=False).to_csv(out / "rna_protein_correlation.csv", index=False)
        # mismatched-pair null: RNA of gene i against the protein of a random other gene j
        n_null = 20000 if not args.quick else 2000
        ii = rng.integers(0, len(Pg), n_null)
        jj = rng.integers(0, len(Pg), n_null)
        ok = ii != jj
        Rv, Pv = Rg.to_numpy(), Pg.to_numpy()
        null_rho = np.array([spearmanr(Rv[a], Pv[b]).statistic for a, b in zip(ii[ok], jj[ok])])
        q95 = float(np.nanquantile(null_rho, 0.95))
        rng_q = pd.qcut(rp["prot_range_log2ppm"], 4, labels=["Q1 smallest protein range", "Q2", "Q3", "Q4 largest protein range"])
        by_q = rp.groupby(rng_q, observed=True)["spearman"].agg(["median", "mean", "size"]).reset_index().rename(columns={"prot_range_log2ppm": "protein_range_quartile"})
        rp_summary = pd.DataFrame([{"n_genes": len(rp), "n_proteins_with_all_7_tissue_means": int(P.shape[0]),
                                    "n_shared_animals_per_tissue": ";".join(f"{t}:{n_shared[t]}" for t in rii.PROT_TISSUES),
                                    "spearman_median": float(np.nanmedian(rho)), "spearman_q25": float(np.nanquantile(rho, 0.25)),
                                    "spearman_q75": float(np.nanquantile(rho, 0.75)), "spearman_mean": float(np.nanmean(rho)),
                                    "frac_spearman_gt_0.5": float(np.mean(rho > 0.5)), "frac_spearman_lt_0": float(np.mean(rho < 0)),
                                    "frac_same_marker_tissue": float(rp["same_marker_tissue"].mean()), "chance_same_marker": 1 / 7,
                                    "null_mismatched_median": float(np.nanmedian(null_rho)), "null_mismatched_q95": q95,
                                    "frac_genes_above_null_q95": float(np.mean(rho > q95)), "n_null_pairs": int(ok.sum()),
                                    **{f"spearman_median_{r['protein_range_quartile']}": r["median"] for _, r in by_q.iterrows()}}])
        rp_summary.to_csv(out / "rna_protein_correlation_summary.csv", index=False)
        by_q.to_csv(out / "rna_protein_correlation_by_protein_range.csv", index=False)
        # RNA panel genes
        lists = {"stable_core_T5": FROZEN / "05_panels" / "TRNSCRPT" / "candidate_panel_annotated.csv",
                 "stability_k20": FROZEN / "05_panels" / "TRNSCRPT" / "stability_k20_annotated.csv",
                 "panel_k20_all_animals": FROZEN / "12_bodymap" / "panel_gene_check.csv"}
        rows = {}
        for name, p in lists.items():
            if not p.exists():
                continue
            d = pd.read_csv(p)
            for _, r in d.iterrows():
                e = rows.setdefault(r["feature_ID"], {"ensembl_gene": r["feature_ID"], "gene_symbol": r.get("gene_symbol"), "marker_tissue_19": r.get("marker_tissue"), "lists": set()})
                e["lists"].add(name)
        pg = pd.DataFrame([{**v, "lists": ";".join(sorted(v["lists"]))} for v in rows.values()])
        pg["marker_in_prot7"] = pg["marker_tissue_19"].isin(rii.PROT_TISSUES)
        quant_any = set(p2g.loc[p2g["protein_id"].isin(om_outer.X.columns), "ensembl_gene"])
        pg["protein_quantified_any_tissue"] = pg["ensembl_gene"].isin(quant_any)
        pg["protein_ids_any"] = pg["ensembl_gene"].map(p2g[p2g["protein_id"].isin(om_outer.X.columns)].groupby("ensembl_gene")["protein_id"].agg(lambda s: ";".join(sorted(set(s)))))
        pg = pg.merge(rp[["ensembl_gene", "spearman", "pearson", "rna_marker_tissue", "prot_marker_tissue", "prot_range_log2ppm", "protein_ids"]], on="ensembl_gene", how="left")
        pg["protein_means_in_all_7"] = pg["spearman"].notna()
        # for genes whose protein has tissue means in only some tissues: marker among the tissues where it is quantified
        for i, r in pg[~pg["protein_means_in_all_7"] & pg["protein_quantified_any_tissue"]].iterrows():
            pids = [p for p in str(r["protein_ids_any"]).split(";") if p in P.index or p in Xp.columns]
            sub = Xp[pids].mean(axis=1)
            tm = sub.groupby(mp["tissue"].to_numpy()).agg(lambda s: s.mean() if s.notna().sum() >= 3 else np.nan).dropna()
            if len(tm):
                pg.loc[i, "prot_marker_tissue"] = tm.idxmax()
                pg.loc[i, "prot_tissues_quantified"] = ";".join(tm.index)
                pg.loc[i, "rna_marker_tissue"] = R.loc[r["ensembl_gene"]].idxmax() if r["ensembl_gene"] in R.index else np.nan
        pg["rna_marker_in7_equals_marker19"] = pg["rna_marker_tissue"] == pg["marker_tissue_19"]
        pg["same_marker_rna_vs_protein_in7"] = (pg["rna_marker_tissue"] == pg["prot_marker_tissue"]) & pg["prot_marker_tissue"].notna()
        pg["testable_c"] = pg["marker_in_prot7"] & pg["prot_marker_tissue"].notna()
        pg["protein_marker_equals_rna_marker19"] = pg["testable_c"] & (pg["prot_marker_tissue"] == pg["marker_tissue_19"])
        pg = pg.sort_values(["testable_c", "spearman"], ascending=[False, False])
        pg.to_csv(out / "rna_protein_panel_genes.csv", index=False)
        n_test = int(pg["testable_c"].sum())
        rp_panel_summary = pd.DataFrame([{"n_rna_panel_genes": len(pg), "n_with_protein_quantified": int(pg["protein_quantified_any_tissue"].sum()),
                                          "n_with_protein_means_in_all_7": int(pg["protein_means_in_all_7"].sum()),
                                          "n_marker_in_prot7": int(pg["marker_in_prot7"].sum()), "n_testable_c": n_test,
                                          "n_same_marker": int(pg["protein_marker_equals_rna_marker19"].sum()),
                                          "frac_same_marker": float(pg["protein_marker_equals_rna_marker19"].sum() / n_test) if n_test else np.nan,
                                          "prereg_c_threshold": 0.70,
                                          "prereg_c_pass": bool(n_test and pg["protein_marker_equals_rna_marker19"].sum() / n_test >= 0.70),
                                          "spearman_median_panel_genes": float(pg["spearman"].median()),
                                          "spearman_median_all_genes": float(np.nanmedian(rho)),
                                          "n_panel_genes_above_null_q95": int((pg["spearman"] > q95).sum())}])
        rp_panel_summary.to_csv(out / "rna_protein_panel_summary.csv", index=False)
        rp_panel = pg
        print(rp_summary.T.to_string()); print(rp_panel_summary.T.to_string())

    # ---- 6. README from the CSVs ------------------------------------------------------------------------------------
    frozen_vp = pd.read_csv(FROZEN / "03_eda" / "variance_partition_PROT.csv")
    frozen_diag = pd.read_csv(FROZEN / "04_baselines" / "PROT" / "diagnostic_accuracy.csv")
    js = join.iloc[0]
    lines = [f"# Phase 1 — the RII rescue ({args.assay}, release {args.release})", "",
             f"Built by `scripts/multiomic/01_rii_rescue.py` on {time.strftime('%Y-%m-%d %H:%M UTC', time.gmtime())}; every number below is read from a CSV in this directory.", "",
             "## Data",
             f"- Peptide-level reporter-ion intensities from the portal quant-id `{args.assay}` folders, {int(js['n_tissues'])} tissues, {int(js['n_vials'])} vials, {int(js['n_animals'])} animals (`join_summary.csv`, `tissue_summary.csv`, `plex_summary.csv`).",
             f"- Per plex: reference channel dropped, contaminants dropped, peptides summed per protein, proteins with < {args.min_peptides} quantified peptides dropped, each channel normalised to its {args.norm} (log2 ppm).",
             f"- Union {int(js['n_proteins_union'])} proteins; **{int(js['n_proteins_inner'])} quantified in every tissue** (the main matrix; NaN fraction {fmt(js['frac_nan_inner'])}, {int(js['n_proteins_inner_complete'])} with no missing value at all); outer NaN fraction {fmt(js['frac_nan_outer'])}.",
             "- **Plex is nested in tissue here too** (one plex = 10 samples of one tissue + that tissue's reference pool; `plex_id` has one tissue per level), so everything below is within-study evidence.", "",
             "## 1. Variance partition (PCA on the stacked matrix, phase-03 code) — `variance_partition.csv`, `variance_partition_ratio.csv`", ""]
    tab = pd.DataFrame({"PC": vp["PC"].head(3), "explained_RII": vp["explained"].head(3), "R2_tissue_RII": vp["R2_tissue"].head(3),
                        "R2_tissue_null95_RII": vp.get("R2_tissue_null95", pd.Series([np.nan] * 3)).head(3),
                        "R2_plex_id_RII": vp["R2_plex_id"].head(3), "R2_sex_RII": vp["R2_sex"].head(3), "R2_pid_RII": vp["R2_pid"].head(3),
                        "R2_tissue_ratio_same_code": vp_r["R2_tissue"].head(3), "R2_tissue_ratio_frozen_phase03": frozen_vp["R2_tissue"].head(3)})
    lines += [report.df_to_md(tab, floatfmt=".4f"), "",
              f"Pre-registration (a) asks for tissue R² of PC1 > 0.5: observed **{fmt(vp['R2_tissue'].iloc[0], 3)}** → {'PASS' if vp['R2_tissue'].iloc[0] > 0.5 else 'FAIL'}."
              + (f" Median-normalised sensitivity: PC1 R² {fmt(vp_med['R2_tissue'].iloc[0], 3)} (`variance_partition_median_norm.csv`)." if vp_med is not None else ""),
              "R2_plex_id equals R2_tissue up to the within-tissue plex split because plex is nested in tissue; it cannot be separated here.", "",
              "## 2. Phase-04 diagnostic on every animal-grouped fold — `diagnostic_accuracy.csv`, `diagnostic_accuracy_summary.csv`", ""]
    ds2 = ds.copy()
    ds2["frozen_ratio_fold0"] = ds2["quantity"].map({"model_as_fitted_acc": frozen_diag["model_as_fitted"].iloc[0],
                                                     "missingness_outer_acc": frozen_diag["missingness_indicators_only"].iloc[0],
                                                     "per_tissue_means_removed_acc": frozen_diag["per_tissue_means_removed"].iloc[0],
                                                     "chance": frozen_diag["chance_balanced"].iloc[0]})
    lines += [report.df_to_md(ds2, floatfmt=".3f"), "",
              f"Read: `model_as_fitted` is the pipeline's logreg_l2 on the RII matrix ({n_splits} folds, {fmt(diag['n_test_animals'].mean(), 1)} test animals per fold). "
              "`missingness_outer` classifies tissue from the NaN pattern of the union matrix alone (which proteins were quantified in which tissue — an artefact of "
              "per-tissue searches, still present on the RII scale); `missingness_inner` does the same on the every-tissue matrix (within-tissue plex gaps only). "
              "`per_tissue_means_removed` erases every protein's tissue mean using the labels: on ratios it collapsed to chance because the classifier lived on "
              "normalisation offsets; on RII the tissue means ARE the fingerprint, so the collapse is expected — the difference between the two matrices is what "
              "the means are (section 1 and 4), not whether the classifier uses them.", "",
              "## 3. Protein panel curve (RoundRobinSelector → logreg_l2, animal-grouped folds, pipeline unchanged) — `panel_curve_summary.csv`, `panel_null_k20.csv`", "",
              report.df_to_md(agg, floatfmt=".3f"), "",
              f"Label-permutation null at k = {k_null} ({int(null['n_perm'].iloc[0])} permutations, quick fits): mean {fmt(null['null_mean'].iloc[0])}, 95th percentile {fmt(null['null_q95'].iloc[0])}, "
              f"chance {fmt(null['chance'].iloc[0])}; observed {fmt(obs)}.",
              f"Pre-registration (b) asks for mean balanced accuracy ≥ 0.95 at some k ≤ 100: best {fmt(agg['bal_acc_mean'].max())} at k = {int(agg.loc[agg['bal_acc_mean'].idxmax(), 'k'])} → "
              f"{'PASS' if agg['bal_acc_mean'].max() >= 0.95 else 'FAIL'} (within-study; context, not a finding).", "",
              f"The k = {k_null} panel selected on all animals (`panel_k20_all_animals.csv`):", "", report.df_to_md(panel, floatfmt=".2f"), "",
              f"**Missingness-free check** (`variance_partition_complete.csv`, `diagnostic_accuracy_complete.csv`, `panel_curve_complete_summary.csv`): on the "
              f"{len(comp_cols)} proteins quantified in every vial (no NaN, no imputation), tissue R² of PC1 = {fmt(vp_c['R2_tissue'].iloc[0])}, PC2 = {fmt(vp_c['R2_tissue'].iloc[1])}; "
              f"logreg_l2 accuracy {fmt(diag_c['model_as_fitted_acc'].mean())} ± {fmt(diag_c['model_as_fitted_acc'].std())}, means-removed {fmt(diag_c['per_tissue_means_removed_acc'].mean())}; "
              f"panel k = {k_null} balanced accuracy {fmt(agg_c.loc[agg_c['k'] == k_null, 'bal_acc_mean'].iloc[0]) if k_null in set(agg_c['k']) else 'n/a'} ± "
              f"{fmt(agg_c.loc[agg_c['k'] == k_null, 'bal_acc_sd'].iloc[0]) if k_null in set(agg_c['k']) else 'n/a'}.", ""]
    for kk in (20, 30):
        p = out / f"confusable_pairs_k{kk}.csv"
        if p.exists():
            cp = pd.read_csv(p)
            lines += [f"Confused pairs at k = {kk} (pooled over folds; {int(cp['count'].sum())} errors):", "", report.df_to_md(cp.head(10), floatfmt=".2f") if len(cp) else "_none_", ""]
    if len(rp_summary):
        s = rp_summary.iloc[0]
        lines += ["## 4. Cross-tissue RNA–protein correlation on the same animals — `rna_protein_correlation.csv`, `rna_protein_correlation_summary.csv`", "",
                  f"Tissue means over the animals shared by RNA-seq and proteomics in each tissue ({s['n_shared_animals_per_tissue']}); Spearman across the 7 tissues per gene "
                  f"(protein means averaged over the proteins mapping to the gene; a tissue mean needs ≥ 3 quantified animals). n = {int(s['n_genes'])} genes.",
                  f"- Spearman median **{fmt(s['spearman_median'])}** (IQR {fmt(s['spearman_q25'])}–{fmt(s['spearman_q75'])}); {fmt(100 * s['frac_spearman_gt_0.5'], 1)} % of genes > 0.5; "
                  f"{fmt(100 * s['frac_spearman_lt_0'], 1)} % < 0.",
                  f"- Mismatched-pair null (RNA of one gene vs the protein of another, {int(s['n_null_pairs'])} pairs): median {fmt(s['null_mismatched_median'])}, 95th percentile {fmt(s['null_mismatched_q95'])}; "
                  f"{fmt(100 * s['frac_genes_above_null_q95'], 1)} % of genes exceed the null 95th percentile.",
                  f"- Same marker tissue (highest tissue mean) at RNA and protein: {fmt(100 * s['frac_same_marker_tissue'], 1)} % of genes (chance 1/7 = {fmt(100 / 7, 1)} %).",
                  f"- By protein cross-tissue range: median Spearman {fmt(s.get('spearman_median_Q1 smallest protein range'))} (Q1) → {fmt(s.get('spearman_median_Q4 largest protein range'))} (Q4) (`rna_protein_correlation_by_protein_range.csv`).", ""]
        if len(rp_panel_summary):
            q = rp_panel_summary.iloc[0]
            lines += ["### RNA panel genes at the protein level — `rna_protein_panel_genes.csv`, `rna_protein_panel_summary.csv`", "",
                      f"{int(q['n_rna_panel_genes'])} RNA panel genes (stable core, stability-51, k=20 all-animal panel); {int(q['n_with_protein_quantified'])} have a quantified protein, "
                      f"{int(q['n_with_protein_means_in_all_7'])} in all 7 tissues; {int(q['n_marker_in_prot7'])} have their RNA marker tissue among the 7 proteomics tissues.",
                      f"- Pre-registration (c): of the **{int(q['n_testable_c'])}** testable genes (marker tissue among the 7, protein quantified), the protein's highest tissue mean is the RNA marker tissue in "
                      f"**{int(q['n_same_marker'])}** ({fmt(100 * q['frac_same_marker'], 0)} %) → {'PASS' if q['prereg_c_pass'] else 'FAIL'} (threshold 70 %).",
                      f"- Median Spearman over panel genes with all-7 means {fmt(q['spearman_median_panel_genes'])} vs {fmt(q['spearman_median_all_genes'])} over all genes; {int(q['n_panel_genes_above_null_q95'])} panel genes above the null 95th percentile.", "",
                      report.df_to_md(rp_panel[["gene_symbol", "lists", "marker_tissue_19", "protein_quantified_any_tissue", "prot_marker_tissue", "rna_marker_tissue", "spearman", "testable_c", "protein_marker_equals_rna_marker19"]].head(70), floatfmt=".2f"), ""]
    lines += ["## What this phase does not show", "",
              "- Nothing about transfer: one plex holds one tissue, so tissue and plex are confounded exactly as in the audit; the tissue axis is real on this scale but its size cannot be separated from a per-plex processing offset without an external dataset (Phases 3 and 6).",
              "- Within-study accuracy is a ceiling task and is reported as context.",
              f"- The normalisation (channel total over the kept proteins) makes every value relative to that channel's quantified proteome; a tissue that quantifies fewer proteins has larger shares. The median-normalised sensitivity (`variance_partition_median_norm.csv`) checks that the PC1 result does not depend on this choice.",
              "", f"_Run time {(time.time() - t_start) / 60:.1f} min._"]
    (out / "README.md").write_text("\n".join(lines))
    json.dump({"args": vars(args), "minutes": (time.time() - t_start) / 60}, open(out / "run_info.json", "w"), indent=1)
    print(f"wrote {out} in {(time.time() - t_start) / 60:.1f} min")


if __name__ == "__main__":
    main()
