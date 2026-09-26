#!/usr/bin/env python
"""Phase 04b — what does a cross-tissue PROT classifier actually use? A diagnostic, not a result.

The EDA showed that tissue explains 0.001 of PC1 in the stacked `PROT_*_NORM_DATA` (ratios to a
per-tissue reference pool, median-centred per sample). If phase 04 still classifies tissue well
above chance, this script says what the model is using. On outer fold 0 (same folds and seed as
phase 04) it fits the logreg_l2 pipeline on the training animals and reports:
  1. the distribution of across-tissue R² over all shared features (do tissue means differ at all?);
  2. for the largest-|coefficient| features per tissue: their tissue R², within-tissue R² of sex and
     of TMT channel (data/raw/meta/PROT.csv), and per-tissue missingness;
  3. fold-0 test accuracy of (a) the model as fitted, (b) a model on missingness indicators only,
     (c) a model after per-tissue mean removal (uses the labels to centre — a diagnostic of what is
     left once mean differences are erased, never an evaluation).
Outputs (results/04_baselines/PROT/): diagnostic_tissue_r2.csv, diagnostic_top_coefs.csv,
diagnostic_accuracy.csv; one REPORT.md section.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

from motrpac import cli, config as C, io, models, report
from motrpac.splits import grouped_kfold


def r2_cat(x: np.ndarray, cat: np.ndarray) -> float:
    ok = ~np.isnan(x)
    x, cat = x[ok], cat[ok]
    if len(x) < 3 or x.var() == 0 or len(np.unique(cat)) < 2:
        return np.nan
    between = sum(((x[cat == c].mean() - x.mean()) ** 2) * (cat == c).sum() for c in np.unique(cat))
    return float(between / (x.var() * len(x)))


def main() -> None:
    ap = cli.common_parser("PROT cross-tissue diagnostic")
    ap.add_argument("--top", type=int, default=3, help="largest-|coef| features per tissue to inspect")
    args = ap.parse_args()
    cli.banner("04_prot_diagnostic", args)
    out = cli.outdir(f"04_baselines/{args.assay}", args.out)
    pheno = io.load_pheno()
    om = io.stack_tissues(args.assay, tissues=cli.parse_tissues(args.tissues), source="norm", join="inner", pheno=pheno,
                          complete=args.complete_features,
                          drop_incomplete_samples=args.drop_incomplete_samples)
    X = om.X.to_numpy(dtype=float)
    y = om.meta["tissue"].astype(str).to_numpy()
    g = om.groups()
    feat = np.asarray(om.X.columns)
    classes = sorted(np.unique(y))
    tr, te = next(iter(grouped_kfold(om.meta, "tissue", args.n_splits, args.seed)))

    # 1. tissue R² per protein, all samples (descriptive)
    tr2 = pd.Series([r2_cat(X[:, j], y) for j in range(X.shape[1])], index=feat, name="tissue_R2")
    q = tr2.quantile([0.5, 0.9, 0.99]).round(3).to_dict()
    tr2.sort_values(ascending=False).to_csv(out / "diagnostic_tissue_r2.csv", header=True)

    # 2. the fitted model and its top coefficients
    est = models.fit_tuned("logreg_l2", X[tr], y[tr], g[tr], prefilter=args.prefilter, quick=args.quick, seed=args.seed)
    clf = est.named_steps["clf"]
    kept = est.named_steps["prefilter"].idx_
    coef = pd.DataFrame(clf.coef_, index=list(clf.classes_), columns=feat[kept])
    meta_path = C.META_DIR / f"{args.assay}.csv"
    channel = None
    if meta_path.exists():
        meta = pd.read_csv(meta_path, dtype=str, low_memory=False).drop_duplicates("viallabel").set_index("viallabel")
        if "tmt11_channel" in meta.columns:
            channel = meta["tmt11_channel"].reindex(om.X.index).to_numpy()
    sex = om.meta["sex"].astype(str).to_numpy()
    rows = []
    for t in classes:
        top = coef.loc[t].abs().sort_values(ascending=False).head(args.top)
        for f in top.index:
            j = int(np.flatnonzero(feat == f)[0])
            mk = y == t
            nan_by_t = pd.Series(np.isnan(X[:, j])).groupby(y).mean()
            rows.append({"tissue": t, "feature_ID": f, "coef": float(coef.loc[t, f]),
                         "tissue_R2_all_samples": float(tr2[f]),
                         "mean_in_tissue": float(np.nanmean(X[mk, j])), "mean_other_tissues": float(np.nanmean(X[~mk, j])),
                         "R2_sex_within_tissue": r2_cat(X[mk, j], sex[mk]),
                         "R2_channel_within_tissue": r2_cat(X[mk, j], channel[mk]) if channel is not None else np.nan,
                         "frac_nan_in_tissue": float(nan_by_t[t]), "max_frac_nan_other_tissue": float(nan_by_t.drop(t).max())})
    top_df = pd.DataFrame(rows)
    top_df.to_csv(out / "diagnostic_top_coefs.csv", index=False)

    # 3. fold-0 accuracies: as fitted / missingness only / per-tissue means removed (label-using)
    acc = {"model_as_fitted": float(np.mean(est.predict(X[te]) == y[te]))}
    M = np.isnan(X).astype(float)
    if M[tr].std(axis=0).max() > 0:
        miss = LogisticRegression(C=0.1, max_iter=3000).fit(M[tr], y[tr])
        acc["missingness_indicators_only"] = float(np.mean(miss.predict(M[te]) == y[te]))
    else:
        acc["missingness_indicators_only"] = np.nan
    Xi = SimpleImputer(strategy="median").fit(X[tr]).transform(X)
    Xc = Xi.copy()
    for t in classes:  # centre every tissue on its TRAINING-fold feature means (uses labels: diagnostic only)
        mu = Xi[tr][y[tr] == t].mean(axis=0)
        Xc[y == t] -= mu
    sc = StandardScaler().fit(Xc[tr])
    cen = LogisticRegression(C=0.1, max_iter=3000).fit(sc.transform(Xc[tr]), y[tr])
    acc["per_tissue_means_removed"] = float(np.mean(cen.predict(sc.transform(Xc[te])) == y[te]))
    acc["chance_balanced"] = 1.0 / len(classes)
    acc_df = pd.DataFrame([acc])
    acc_df.to_csv(out / "diagnostic_accuracy.csv", index=False)

    n_hi = int((tr2 > 0.5).sum())
    body = [
        f"{om.notes[-1]}; fold 0 of {args.n_splits}: {len(np.unique(g[tr]))} training / {len(np.unique(g[te]))} test animals; "
        f"model `logreg_l2` fitted on the training animals only.",
        f"**Across-tissue R² of tissue per protein** (all samples, {len(tr2)} shared features): median {q[0.5]}, "
        f"90th percentile {q[0.9]}, 99th percentile {q[0.99]}; {n_hi} features with R² > 0.5. "
        "Even after per-tissue reference ratios and per-sample centring, a minority of features keep tissue-specific means "
        "(the reference pool is not the tissue mean, and the per-sample median shifts each feature by a tissue-dependent amount).",
        f"Largest-|coefficient| features per tissue ({args.top} each), with what they track:\n\n"
        + report.df_to_md(top_df, floatfmt=".2f"),
        "Fold-0 test accuracy (plain accuracy, 7 classes):\n\n" + report.df_to_md(acc_df, floatfmt=".3f"),
        "Read: `model_as_fitted` is what phase 04 measures; `missingness_indicators_only` is how far the NaN pattern alone "
        "identifies tissue (an imputation artefact, not biology); `per_tissue_means_removed` is what is left when every "
        "feature's tissue mean is erased (label-using diagnostic). If the first is high and the third collapses toward "
        f"chance ({1/len(classes):.2f}), the classifier lives on residual mean offsets of the per-tissue normalization — "
        "a signature of how each tissue was processed, **not a proteomic tissue fingerprint**. If sex or channel R² "
        "within tissue is high for the top features, the model is also leaning on sex- or channel-linked variation.",
    ]
    subset = f", {om.meta['tissue'].nunique()} tissues" if args.tissues else ""
    report.add_section(f"04 · {args.assay} cross-tissue diagnostic (not a result{subset})", "\n\n".join(body),
                       params={"assay": args.assay, "tissues": args.tissues, "fold": 0, "top": args.top})
    print(top_df.round(3).to_string(index=False))
    print(acc_df.round(3).to_string(index=False))
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
