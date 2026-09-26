#!/usr/bin/env python
"""Phase 03 — EDA: how much of the variance is tissue, sex, time point? Any batch structure?

Outputs (results/03_eda/):
  pca_<assay>_by_tissue.png, pca_<assay>_by_sex.png, pca_<assay>_by_group.png
  variance_partition_<assay>.csv   R² of tissue / sex / group / pid on each of the first 10 PCs
  batch_partition_<assay>.csv      R² of batch/QC covariates (data/raw/meta/) on the same PCs, next to tissue
                                   and the nested tissue|covariate factor
  within_tissue_pca_<assay>.csv    per tissue (own full feature table): R² of sex, group and each batch
                                   covariate on PCs 1–5, a permutation null, and design-residualized batch R²
  within_pca_<assay>_by_group.png  PC1 vs PC2 inside selected tissues, colored by time point
  tissue_by_sex_<assay>.csv        samples per tissue × sex (single-sex tissues = unseen classes under a sex shift)
  feature_platforms_<assay>.csv    METAB/IMMUNO: platform composition of the features kept by the join
  feature_stats_<assay>.csv        per-feature mean/var/missingness

EDA is descriptive only — nothing here is used to select features for later phases.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler

from motrpac import cli, config as C, io, plots, report

import matplotlib.pyplot as plt  # noqa: E402  (Agg backend is set by motrpac.plots)

# Sample-level batch covariates in data/raw/meta/<ASSAY>.csv (categorical) and continuous QC columns.
# METAB / IMMUNO have no sample-level meta object; platforms were merged per biospecimen, so the
# covariate there is how many platforms measured the sample (from meta/<ASSAY>_VIALS.csv).
BATCH_COLS = {
    "PROT": ["tmt_plex", "tmt11_channel"], "PHOSPHO": ["tmt_plex", "tmt11_channel"],
    "ACETYL": ["tmt_plex", "tmt11_channel"], "UBIQ": ["tmt_plex", "tmt11_channel"],
    "TRNSCRPT": ["GET_site", "RNA_extr_plate_ID", "Lib_batch_ID", "Seq_flowcell_ID"],
}
QC_COLS = {"TRNSCRPT": ["RIN", "pct_rRNA", "pct_globin", "pct_mrna", "pct_uniquely_mapped"]}


def r2_categorical(scores: np.ndarray, cat: pd.Series) -> np.ndarray:
    """Fraction of variance of each column of `scores` explained by group means of `cat`."""
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


def r2_of(scores: np.ndarray, cov) -> np.ndarray:
    """R² of each PC on a covariate: group means for a categorical one, squared correlation for a
    numeric one. Samples with a missing covariate are left out; all-NaN → NaN."""
    cov = pd.Series(np.asarray(cov)).reset_index(drop=True)
    numeric = pd.api.types.is_numeric_dtype(cov)
    ok = cov.notna().to_numpy()
    if ok.sum() < 3 or (not numeric and cov[ok].nunique() < 2):
        return np.full(scores.shape[1], np.nan)
    if numeric:
        c = cov[ok].to_numpy(dtype=float)
        out = []
        for j in range(scores.shape[1]):
            s = scores[ok, j]
            out.append(float(np.corrcoef(s, c)[0, 1] ** 2) if s.std() > 0 and c.std() > 0 else 0.0)
        return np.array(out)
    return r2_categorical(scores[ok], cov[ok])


def perm_null95(scores: np.ndarray, cat, rng: np.random.Generator, n_perm: int = 200) -> float:
    """95th percentile of the max-over-PCs R² when the labels are shuffled within the tissue:
    the chance level that an observed max R² has to beat."""
    cat = pd.Series(cat).astype(str).to_numpy()
    vals = [r2_categorical(scores, rng.permutation(cat)).max() for _ in range(n_perm)]
    return float(np.quantile(vals, 0.95))


def shared_with(cov, factor) -> float:
    """Association between a covariate and a design factor, on a 0–1 scale: R² of a numeric
    covariate on the factor's groups, or Cramér's V² for a categorical covariate. Near 1 means the
    covariate cannot be told apart from that factor in this tissue."""
    cov = pd.Series(np.asarray(cov)).reset_index(drop=True)
    fac = pd.Series(np.asarray(factor)).astype(str).reset_index(drop=True)
    ok = cov.notna().to_numpy()
    if ok.sum() < 3 or cov[ok].nunique() < 2 or fac[ok].nunique() < 2:
        return np.nan
    if pd.api.types.is_numeric_dtype(cov):
        return float(r2_categorical(cov[ok].to_numpy(dtype=float)[:, None], fac[ok])[0])
    ct = pd.crosstab(cov[ok].astype(str), fac[ok]).to_numpy(dtype=float)
    n = ct.sum()
    expected = np.outer(ct.sum(1), ct.sum(0)) / n
    chi2 = ((ct - expected) ** 2 / expected).sum()
    return float(chi2 / (n * (min(ct.shape) - 1)))


def residualize(scores: np.ndarray, cat) -> np.ndarray:
    """Remove the group means of `cat` (e.g. sex × group design cells) from each PC."""
    cat = pd.Series(cat).astype(str).to_numpy()
    out = scores.astype(float).copy()
    for c in np.unique(cat):
        m = cat == c
        out[m] -= scores[m].mean(axis=0)
    return out


def pca_of(om: io.OmicsMatrix, n_components: int = 10, top_var: int | None = 5000, seed: int = C.SEED):
    X = om.X.to_numpy(dtype=float)
    X = SimpleImputer(strategy="median").fit_transform(X)
    if top_var and X.shape[1] > top_var:
        v = X.var(axis=0)
        X = X[:, np.argsort(v)[::-1][:top_var]]
    X = StandardScaler().fit_transform(X)
    n = min(n_components, X.shape[0] - 1, X.shape[1])
    pca = PCA(n_components=n, random_state=seed)
    S = pca.fit_transform(X)
    return S, pca.explained_variance_ratio_


def batch_covariates(assay: str, viallabels) -> pd.DataFrame:
    """Sample-level batch/QC covariates for these viallabels (NaN where the meta table has no row)."""
    vl = pd.Index([str(v) for v in viallabels], name="viallabel")
    out = pd.DataFrame(index=vl)
    p = C.META_DIR / f"{assay}.csv"
    if p.exists():
        m = pd.read_csv(p, dtype=str, low_memory=False)
        if "viallabel" in m.columns:
            m = m.drop_duplicates("viallabel").set_index("viallabel")
            for c in BATCH_COLS.get(assay, []):
                if c in m.columns:
                    out[c] = m[c].reindex(vl).to_numpy()
            for c in QC_COLS.get(assay, []):
                if c in m.columns:
                    out[c] = pd.to_numeric(m[c].reindex(vl), errors="coerce").to_numpy()
    vials = C.META_DIR / f"{assay}_VIALS.csv"
    if vials.exists():
        v = pd.read_csv(vials, dtype=str)
        n = v.groupby("representative_viallabel")["dataset"].nunique().reindex(vl)
        out["n_platforms"] = [f"{int(x)}" if pd.notna(x) else np.nan for x in n]
    return out


def feature_platforms(assay: str, feature_ids) -> pd.Series | None:
    """METAB/IMMUNO: which platform(s) the features kept by the join come from."""
    p = C.META_DIR / f"{assay}_FEATURES.csv"
    if not p.exists():
        return None
    f = pd.read_csv(p, dtype=str)
    plat = f.groupby("feature_ID")["dataset"].agg(lambda s: "|".join(sorted(set(s))))
    ids = pd.Index([str(x).split("__dup")[0] for x in feature_ids])
    return plat.reindex(ids).fillna("unknown").value_counts()


def within_group_figure(panels, assay: str, source: str, path):
    """PC1 vs PC2 inside each selected tissue, colored by training group, marker by sex."""
    fig, axes = plt.subplots(1, len(panels), figsize=(4.4 * len(panels), 4.2), squeeze=False)
    cmap = plt.get_cmap("viridis", len(C.GROUP_ORDER))
    for ax, (t, S2, meta, ev) in zip(axes[0], panels):
        for gi, gname in enumerate(C.GROUP_ORDER):
            for sex, mk in (("female", "o"), ("male", "^")):
                m = ((meta["group"] == gname) & (meta["sex"] == sex)).to_numpy()
                if m.any():
                    ax.scatter(S2[m, 0], S2[m, 1], color=cmap(gi), marker=mk, s=30, edgecolor="k",
                               linewidth=0.3, label=f"{gname} {sex[0].upper()}")
        ax.set_title(f"{t} (n={len(meta)})")
        ax.set_xlabel(f"PC1 ({ev[0]:.0%} of within-tissue variance)")
        ax.set_ylabel(f"PC2 ({ev[1]:.0%})")
    axes[0][-1].legend(fontsize=6, ncol=2, frameon=False)
    fig.suptitle(f"{assay} ({source}): within-tissue PCA colored by training group (circle = female, triangle = male)",
                 fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def main() -> None:
    ap = cli.common_parser("Exploratory analysis")
    ap.add_argument("--assays", default="TRNSCRPT,PROT,METAB", help="assays to explore (stacked across tissues)")
    ap.add_argument("--join", default="outer", choices=["outer", "inner"],
                    help="feature join across tissues: outer keeps features present in >= --min-present of samples")
    ap.add_argument("--min-present", type=float, default=0.8)
    ap.add_argument("--within-plot-tissues", default="SKM-GN,HEART,WAT-SC",
                    help="tissues for the within-tissue PC1/PC2-by-group figure")
    ap.add_argument("--n-perm", type=int, default=200, help="label permutations for the within-tissue R² null")
    ap.add_argument("--within-features", default="full", choices=["full", "stacked"],
                    help="within-tissue PCA on the tissue's own full feature table (default) or on the cross-tissue join")
    args = ap.parse_args()
    n_perm = 50 if args.quick else args.n_perm
    rng = np.random.default_rng(args.seed)
    cli.banner("03_eda", args)
    out = cli.outdir("03_eda", args.out)
    pheno = io.load_pheno()
    tissues = cli.parse_tissues(args.tissues)
    plot_tissues = set(cli.parse_tissues(args.within_plot_tissues) or [])
    body = ["EDA is descriptive only — nothing here is used to select features for later phases."]

    for assay in [a.strip() for a in args.assays.split(",") if a.strip()]:
        source = cli.resolve_source(assay, args.source)
        try:
            om = io.stack_tissues(assay, tissues=tissues, source=source, join=args.join,
                                  min_present=args.min_present, pheno=pheno)
        except FileNotFoundError as e:
            body.append(f"**{assay}**: skipped ({e})")
            continue
        body.append(f"### {assay} ({source})\n\n{om.notes[-1]}")

        # tissue × sex: single-sex tissues are unseen classes under a held-out-sex shift (phase 08)
        tbs = pd.crosstab(om.meta["tissue"], om.meta["sex"])
        tbs.to_csv(out / f"tissue_by_sex_{assay}.csv")
        single = [t for t in tbs.index if int((tbs.loc[t] > 0).sum()) == 1]
        body.append("Samples per tissue × sex:\n\n" + report.df_to_md(tbs.reset_index())
                    + (f"\n\nSingle-sex tissues: **{', '.join(single)}**. Under a held-out-sex split the source never "
                       "contains them, so they are unseen classes on the target: not coverable by any set and "
                       "excluded from 'seen' coverage in phase 08." if single else "\n\nNo single-sex tissues."))

        # METAB / IMMUNO: which platforms the kept features come from
        fp = feature_platforms(assay, om.X.columns)
        if fp is not None:
            fp.rename_axis("platform").reset_index(name="n_features").to_csv(out / f"feature_platforms_{assay}.csv", index=False)
            body.append("Platform(s) of the features kept by the join (features measured on several platforms are "
                        "listed with `|`):\n\n" + report.df_to_md(fp.rename_axis("platform").reset_index(name="n_features").head(12)))

        # per-feature stats
        fs = pd.DataFrame({"mean": om.X.mean(axis=0), "var": om.X.var(axis=0), "frac_nan": om.X.isna().mean(axis=0)})
        fs.to_csv(out / f"feature_stats_{assay}.csv")

        # stacked PCA
        S, ev = pca_of(om, seed=args.seed)
        figs = []
        for color, shape in (("tissue", "sex"), ("sex", None), ("group", None)):
            if om.meta[color].nunique() > 1:
                figs.append(plots.pca_scatter(S, om.meta, color, shape, f"{assay} ({source}) PCA, stacked tissues",
                                              out / f"pca_{assay}_by_{color}.png", (ev[0], ev[1])))
        vp = pd.DataFrame({"PC": [f"PC{i+1}" for i in range(S.shape[1])], "explained": ev[: S.shape[1]]})
        for col in ("tissue", "sex", "group", "pid"):
            if om.meta[col].nunique() > 1:
                vp[f"R2_{col}"] = r2_categorical(S, om.meta[col])
        vp.to_csv(out / f"variance_partition_{assay}.csv", index=False)
        body.append("Variance explained by design factors on the leading PCs (R² of group means; `pid` is the "
                    "animal effect and is an upper bound because animals nest sex × group):\n\n"
                    + report.df_to_md(vp.head(6), floatfmt=".2f"))
        body.extend(report.figure_md(f, f"{assay}: PCA colored by {f.stem.split('_by_')[-1]}") for f in figs)

        # batch / QC covariates on the same PCs
        bc = batch_covariates(assay, om.X.index)
        if bc.shape[1]:
            bp = pd.DataFrame({"PC": vp["PC"], "explained": vp["explained"], "R2_tissue": r2_categorical(S, om.meta["tissue"])})
            tissue_str = om.meta["tissue"].astype(str).to_numpy()
            for c in bc.columns:
                bp[f"R2_{c}"] = r2_of(S, bc[c])
                if not pd.api.types.is_numeric_dtype(bc[c]):
                    nested = [f"{t}|{b}" if pd.notna(b) else np.nan for t, b in zip(tissue_str, bc[c])]
                    bp[f"R2_tissue|{c}"] = r2_of(S, nested)
            bp.to_csv(out / f"batch_partition_{assay}.csv", index=False)
            n_levels = {c: int(bc[c].nunique()) for c in bc.columns}
            body.append("Batch / QC covariates on the same PCs (`R2_tissue|x` is the nested tissue-within-covariate "
                        "factor, so its excess over `R2_tissue` is what the batch adds; levels: "
                        + ", ".join(f"{c}={n}" for c, n in n_levels.items()) + "; samples without a meta row: "
                        + f"{int(bc.isna().all(axis=1).sum())}):\n\n" + report.df_to_md(bp.head(3), floatfmt=".2f"))
        else:
            bp = None

        # within-tissue: is sex / time point / batch visible at all?
        rows, panels, skipped = [], [], []
        for t in sorted(om.meta["tissue"].unique()):
            if args.within_features == "full":  # the tissue's own table: all its features, not only the join
                sub = io.load_counts(t, pheno, assay=assay) if source == "counts" else io.load_norm(assay, t, pheno)
                sub = sub.subset(sub.meta.index.isin(om.meta.index[(om.meta["tissue"] == t).to_numpy()]))
            else:
                sub = om.subset((om.meta["tissue"] == t).to_numpy())
            keep = (sub.X.notna().mean(axis=0) >= args.min_present).to_numpy()
            n_dropped = int((~keep).sum())
            sub = sub.select_features(list(sub.X.columns[keep]))  # remaining NaNs are median-imputed in pca_of
            if sub.n_samples < 8 or sub.n_features < 5:
                print(f"  within-tissue PCA skipped for {t}: {sub.n_samples} samples, {sub.n_features} features "
                      f"present in >= {args.min_present:.0%} of samples ({n_dropped} dropped)")
                skipped.append(f"{t} ({sub.n_samples} samples, {sub.n_features} usable features)")
                continue
            St, evt = pca_of(sub, n_components=5, seed=args.seed)
            has_sex, has_grp = sub.meta["sex"].nunique() > 1, sub.meta["group"].nunique() > 1
            r_sex = r2_categorical(St, sub.meta["sex"]) if has_sex else np.zeros(St.shape[1])
            r_grp = r2_categorical(St, sub.meta["group"]) if has_grp else np.zeros(St.shape[1])
            n95_sex = perm_null95(St, sub.meta["sex"], rng, n_perm) if has_sex else np.nan
            n95_grp = perm_null95(St, sub.meta["group"], rng, n_perm) if has_grp else np.nan
            row = {"tissue": t, "n": sub.n_samples, "n_animals": sub.n_animals, "n_features": sub.n_features,
                   "n_features_dropped_missing": n_dropped, "PC1_explained": evt[0], "max_R2_sex_PC1-5": r_sex.max(), "null95_sex": n95_sex,
                   "max_R2_group_PC1-5": r_grp.max(), "null95_group": n95_grp,
                   "group_visible": bool(has_grp and r_grp.max() > n95_grp),
                   "PC_with_max_group_R2": int(np.argmax(r_grp)) + 1}
            if bc.shape[1]:
                bsub = bc.loc[sub.X.index]
                cell = (sub.meta["sex"].astype(str) + "/" + sub.meta["group"].astype(str)).to_numpy()
                St_res = residualize(St, cell)  # what is left after the sex × group design means
                for c in bsub.columns:
                    r, r_res = r2_of(St, bsub[c]), r2_of(St_res, bsub[c])
                    row[f"max_R2_{c}_PC1-5"] = float(np.nanmax(r)) if not np.all(np.isnan(r)) else np.nan
                    row[f"max_R2_{c}|design_PC1-5"] = float(np.nanmax(r_res)) if not np.all(np.isnan(r_res)) else np.nan
                    row[f"n_{c}"] = int(bsub[c].nunique())
                    row[f"R2_{c}~sex"] = shared_with(bsub[c], sub.meta["sex"]) if has_sex else np.nan
                    row[f"R2_{c}~group"] = shared_with(bsub[c], sub.meta["group"]) if has_grp else np.nan
            rows.append(row)
            if t in plot_tissues:
                panels.append((t, St[:, :2], sub.meta, evt))
        wt = pd.DataFrame(rows)
        wt.to_csv(out / f"within_tissue_pca_{assay}.csv", index=False)
        if len(wt):
            body.append("Within each tissue, the largest R² of sex, time point and each batch covariate on any of "
                        f"PCs 1–5. `null95_*` is the 95th percentile of the same maximum with labels shuffled "
                        f"({n_perm} permutations); `group_visible` = observed max above that. `R2_x~sex` / "
                        "`R2_x~group` is the association of a covariate with the design (near 1 = confounded, so "
                        "its PC R² is design, not batch); `max_R2_x|design` is the covariate's R² after removing the "
                        "sex × group cell means, i.e. the part that can only be technical. Within-tissue PCA uses "
                        f"the tissue's {'own full feature table' if args.within_features == 'full' else 'cross-tissue join features'}:\n\n"
                        + report.df_to_md(wt, floatfmt=".2f")
                        + (f"\n\nWithin-tissue PCA skipped (too few usable samples/features): {'; '.join(skipped)}." if skipped else ""))
        if panels:
            fig = within_group_figure(panels, assay, source, out / f"within_pca_{assay}_by_group.png")
            body.append(report.figure_md(fig, f"{assay}: within-tissue PC1 vs PC2 colored by training group"))

    default = args.join == "outer" and args.tissues is None and args.min_present == 0.8
    title = "03 · EDA" if default else f"03 · EDA ({args.assays}; join={args.join}; tissues={args.tissues or 'all'})"
    report.add_section(title, "\n\n".join(body), params={"assays": args.assays, "source": args.source, "join": args.join,
                                                          "min_present": args.min_present, "tissues": args.tissues})
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
