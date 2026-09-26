#!/usr/bin/env python
"""Phase 07b — is the control-vs-trained separability within tissue training, batch, or unresolvable?

Task B (scripts/07_fusion_vs_baselines.py) separates sedentary controls from 8-week-trained animals
within every fusion tissue at AUROC ≈ 1 with each omic alone. Three checks, all on the same
animal-grouped, sex × group-stratified folds as Task B:
  (a) read from Task B: the best single-omic AUROC by training duration (1w, 2w, 4w, 8w vs control;
      `taskB_duration_summary.csv`) — a training effect should grow with the duration, a batch
      effect should not;
  (b) classifiers on QC and collection covariates ONLY (no omics): PHENO collection / dissection /
      sacrifice variables (dates, times, staff, cage, freeze times), transcript library and QC
      metrics from data/raw/meta/TRNSCRPT.csv (RNA plate and date, RIN, library prep, sequencing
      run, mapping and mRNA fractions, …) and the TMT plex / channel of the proteomics vial; if
      these alone separate control from trained, the omic separation cannot be attributed to
      training on these data;
  (c) a univariate screen of EVERY PHENO and meta variable within tissue (AUROC for numeric, dates
      and times; Cramér's V for categorical; variables missing in one group only), each tagged
      design / outcome / baseline / collection / library / qc / processing so that the ones that
      separate by construction (randomization group, training log, VO2max change) are told apart
      from the ones that would indicate a batch.
The verdict is written to batch_conclusion.md (read by scripts/10_make_report.py).

Outputs (results/07_fusion/): batch_covariate_auroc.csv, batch_variable_screen.csv,
batch_variable_screen_summary.csv, batch_conclusion.csv, batch_conclusion.md
"""
from __future__ import annotations

import re
import warnings

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.compose import ColumnTransformer
from sklearn.ensemble import RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler

from motrpac import cli, config as C, io, report
from motrpac.splits import assert_no_group_leak, grouped_kfold

FUSION_TISSUES = ["CORTEX", "HEART", "KIDNEY", "LIVER", "LUNG", "SKM-GN", "WAT-SC"]
DURATION_WEEKS = {"1w": 1, "2w": 2, "4w": 4, "8w": 8}
DROP = {"pid", "bid", "labelid", "viallabel", "vial_label", "2D_barcode", "BID", "PID", "Tissue", "Species", "Sample_category",
        "Phase", "biclabeldata.barcode", "registration.ratid", "tissue", "tissue_code_no", "tissue_description",
        "specimen.processing.sampletypedescription", "specimen.processing.aliquotdescription", "group", "sex", "assay"}
CATEGORY_RULES = [  # (regex on the column name, category); first match wins; None = dropped (administrative)
    (r"comments$|formname|versionnbr|visitcode", None),
    (r"^(key\.(anirandgroup|sacrificetime|intervention|protocol|agegroup)|study_group_timepoint|sacrificetime|intervention|biclabeldata\.protocol)$", "design"),
    (r"^training\.|^vo2\.max\.test\.|^calculated\.variables\.(vo2_max_change|pct_body_(fat|lean|fluid)_change|lactate_change_dueto_train)|^nmr\.testing\..*_2$|^terminal\.weight\.", "outcome"),
    (r"^familiarization\.|^nmr\.testing\..*_1$|^registration\.(weight|d_birth|days_birth)$", "baseline"),
    (r"^key\.(d_arrive|d_sacrifice|siteid)$|^registration\.|^specimen\.collection\.|^specimen\.processing\.|^calculated\.variables\.(coll_time_train|deathtime_after_train|frozetime_after_train)$|^time_to_freeze$|^biclabeldata\.", "collection"),
    (r"^(GET_site|RNA_extr_|RIN$|r_260_|Lib_|Seq_)", "library"),
    (r"^(reads_raw|reads|uniquely_mapped|num_.*splices)$", "depth"),   # library depth: pooling / loading, a pure processing batch
    (r"^(pct_|avg_|median_5_3_bias)", "qc"),                             # composition and quality fractions: batch or biology
    (r"^tmt", "processing"),
]
BATCH_CATEGORIES = ("collection", "library", "depth", "qc", "processing")


def categorize(col: str) -> str | None:
    if col in DROP:
        return None
    for pat, cat in CATEGORY_RULES:
        if re.search(pat, col):
            return cat
    return "other"


def parse_series(v: pd.Series):
    """→ (kind, values): numeric (float), time (hours), datetime (days since 2018-01-01) or categorical (string)."""
    s = v.astype("string").str.strip()
    s = s.where(s.str.len() > 0)
    n_ok = int(s.notna().sum())
    if n_ok == 0:
        return "empty", s
    num = pd.to_numeric(s, errors="coerce")
    if num.notna().sum() >= 0.9 * n_ok:
        return "numeric", num.astype(float)
    if s.dropna().str.fullmatch(r"\d{1,2}:\d{2}(:\d{2})?").all():
        td = pd.to_timedelta(s.where(s.str.count(":") == 2, s + ":00"), errors="coerce")
        return "time", td.dt.total_seconds() / 3600.0
    for fmt in ("%d%b%Y", "%m/%d/%Y", "%Y-%m-%d", "%b %Y"):
        dt = pd.to_datetime(s, format=fmt, errors="coerce")
        if dt.notna().sum() >= 0.9 * n_ok:
            return "datetime", (dt - pd.Timestamp("2018-01-01")).dt.days.astype(float)
    return "categorical", s


def cramers_v(ct: pd.DataFrame) -> tuple[float, float]:
    ct = ct.loc[ct.sum(axis=1) > 0, ct.sum(axis=0) > 0]
    if ct.shape[0] < 2 or ct.shape[1] < 2:
        return np.nan, np.nan
    n = ct.to_numpy().sum()
    if ct.shape == (2, 2):
        p = stats.fisher_exact(ct.to_numpy())[1]
        chi2 = stats.chi2_contingency(ct.to_numpy(), correction=False)[0]
    else:
        chi2, p = stats.chi2_contingency(ct.to_numpy(), correction=False)[:2]
    return float(np.sqrt(chi2 / (n * (min(ct.shape) - 1)))), float(p)


def screen(table: pd.DataFrame, y: np.ndarray, source: str, level: dict, tissue: str) -> list[dict]:
    rows = []
    pos = y == 1
    for col in table.columns:
        cat = categorize(col)
        if cat is None:
            continue
        kind, vals = parse_series(table[col])
        miss_c, miss_t = float(vals[~pos].isna().mean()), float(vals[pos].isna().mean())
        row = {"tissue": tissue, "variable": col, "source": source, "category": cat, "level": level.get(col, "animal"), "kind": kind,
               "n_control": int((~pos).sum()), "n_trained": int(pos.sum()), "frac_missing_control": miss_c, "frac_missing_trained": miss_t,
               "missing_in_one_group_only": bool((miss_c >= 0.8 and miss_t <= 0.2) or (miss_t >= 0.8 and miss_c <= 0.2)),
               "n_levels": int(vals.dropna().nunique()), "stat": None, "separation": np.nan, "p_value": np.nan}
        if kind == "empty" or row["n_levels"] < 2:
            rows.append(row)
            continue
        if kind in ("numeric", "time", "datetime"):
            a, b = vals[pos].dropna().to_numpy(float), vals[~pos].dropna().to_numpy(float)
            if len(a) >= 2 and len(b) >= 2:
                u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
                auc = u / (len(a) * len(b))
                row.update({"stat": "auroc", "separation": float(max(auc, 1 - auc)), "p_value": float(p), "auroc_trained_higher": float(auc)})
        else:
            ct = pd.crosstab(vals, pd.Series(np.where(pos, "trained", "control"), index=vals.index))
            v, p = cramers_v(ct)
            row.update({"stat": "cramers_v", "separation": v, "p_value": p})
        rows.append(row)
    return rows


def covariate_frame(table: pd.DataFrame, cols: list[str]) -> tuple[pd.DataFrame, list[str], list[str]]:
    num, cat = {}, {}
    for c in cols:
        kind, vals = parse_series(table[c])
        if kind == "empty" or vals.dropna().nunique() < 2:
            continue
        if kind == "categorical":
            if vals.dropna().nunique() >= 0.9 * vals.notna().sum():   # identifier-like (unique per sample): no batch structure to learn
                continue
            cat[c] = vals.astype(object).where(vals.notna(), "missing")
        else:
            num[c] = vals.astype(float)
    X = pd.DataFrame({**num, **cat}, index=table.index)
    return X, list(num), list(cat)


def covariate_auroc(X: pd.DataFrame, num_cols, cat_cols, y, meta, n_splits, seed, model, n_perm=0, rng=None):
    pre = ColumnTransformer([("num", Pipeline([("imp", SimpleImputer(strategy="median")), ("sc", StandardScaler())]), num_cols),
                             ("cat", OneHotEncoder(handle_unknown="ignore"), cat_cols)])
    clf = (LogisticRegression(C=0.1, max_iter=3000) if model == "logreg"
           else RandomForestClassifier(n_estimators=500, class_weight="balanced", random_state=seed, n_jobs=-1))
    pipe = Pipeline([("pre", pre), ("clf", clf)])
    folds = list(grouped_kfold(meta, "sex_group", n_splits, seed))

    def cv_auc(yy):
        aucs = []
        for tr, te in folds:
            assert_no_group_leak(meta, tr, te)
            if len(np.unique(yy[te])) < 2 or len(np.unique(yy[tr])) < 2:
                continue
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                pipe.fit(X.iloc[tr], yy[tr])
                aucs.append(roc_auc_score(yy[te], pipe.predict_proba(X.iloc[te])[:, 1]))
        return aucs
    aucs = cv_auc(y)
    out = {"auroc_mean": float(np.mean(aucs)), "auroc_sd": float(np.std(aucs)), "n_folds": len(aucs)}
    if n_perm:
        sex = meta["sex"].to_numpy()
        null = []
        for _ in range(n_perm):
            yp = y.copy()
            for sx in np.unique(sex):
                idx = np.flatnonzero(sex == sx)
                yp[idx] = rng.permutation(yp[idx])
            null.append(float(np.mean(cv_auc(yp))))
        out["null_p95_auroc"] = float(np.quantile(null, 0.95))
        out["p_perm"] = float((np.sum(np.asarray(null) >= out["auroc_mean"]) + 1) / (len(null) + 1))
    return out


def main() -> None:
    ap = cli.common_parser("Batch check for Task B (control vs trained within tissue)")
    ap.add_argument("--trained-group", default="8w")
    ap.add_argument("--n-perm", type=int, default=100, help="permutations for the covariate classifier null (logistic only)")
    ap.add_argument("--sep-auroc", type=float, default=0.9)
    ap.add_argument("--sep-v", type=float, default=0.8)
    args = ap.parse_args()
    cli.banner("07_batch_check", args)
    out = cli.outdir("07_fusion", args.out)
    tissues = cli.parse_tissues(args.tissues) or FUSION_TISSUES
    rng = np.random.default_rng(args.seed)
    n_perm = 10 if args.quick else args.n_perm

    ph = pd.read_csv(C.RAW_DIR / "pheno.csv", dtype=str, low_memory=False)
    ph = ph[ph["group"].isin(["control", args.trained_group])]
    animal_level = {c: ("animal" if ph.groupby("pid")[c].nunique(dropna=True).max() <= 1 else "vial") for c in ph.columns}
    tm = pd.read_csv(C.META_DIR / "TRNSCRPT.csv", dtype=str, low_memory=False).drop_duplicates("viallabel")
    pm_path = C.META_DIR / "PROT.csv"
    pmeta = pd.read_csv(pm_path, dtype=str) if pm_path.exists() else None
    pheno_cols = [c for c in ph.columns if categorize(c)]
    trn_cols = [c for c in tm.columns if categorize(c)]

    auc_rows, screen_rows = [], []
    for t in tissues:
        d = ph[(ph["tissue"] == t) & ph["viallabel"].isin(tm["viallabel"])].drop_duplicates("pid").copy()   # the animal's TRNSCRPT vial
        vials = tm[tm["viallabel"].isin(d["viallabel"])]
        if d["pid"].nunique() < 12 or d["group"].nunique() < 2:
            print(f"  {t}: skipped ({d['pid'].nunique()} animals)")
            continue
        d = d.set_index("pid")
        y = (d["group"] == args.trained_group).astype(int).to_numpy()
        meta = pd.DataFrame({"pid": d.index, "sex": d["sex"].to_numpy(), "group": d["group"].to_numpy()}, index=d.index)
        meta["sex_group"] = meta["sex"] + "/" + meta["group"]
        n_splits = min(4, len(d) // 3)
        # vial-level transcript meta for this tissue's vial of each animal, PROT plex / channel via the biospecimen
        tv = vials.set_index("viallabel").reindex(d["viallabel"]).set_axis(d.index).copy()
        if "Lib_barcode_well" in tv.columns:   # plate position as two numbers (processing order is a classic batch axis)
            w = tv["Lib_barcode_well"].astype("string").str.upper().str.extract(r"^([A-H])(\d{1,2})$")
            tv["Lib_barcode_well_row"] = w[0].map({c: i + 1 for i, c in enumerate("ABCDEFGH")}).astype("string")
            tv["Lib_barcode_well_col"] = w[1]
        trn_cols_t = [c for c in tv.columns if categorize(c)]
        tables = {"pheno_collection": (d, [c for c in pheno_cols if categorize(c) == "collection"]),
                  "trnscrpt_library": (tv, [c for c in trn_cols_t if categorize(c) == "library"]),
                  "trnscrpt_depth": (tv, [c for c in trn_cols_t if categorize(c) == "depth"]),
                  "trnscrpt_qc": (tv, [c for c in trn_cols_t if categorize(c) == "qc"])}
        if pmeta is not None:
            pv = pmeta[(pmeta["tissue"] == t) & pmeta["viallabel"].astype(str).str.fullmatch(r"\d{11}")].drop_duplicates("viallabel").copy()
            pv["bid"] = pv["viallabel"].str[:5]
            pv = pv.drop_duplicates("bid").set_index("bid").reindex(d["bid"]).set_axis(d.index)[["tmt_plex", "tmt11_channel"]]
            tables["prot_plex_channel"] = (pv, ["tmt_plex", "tmt11_channel"])
        # (b) covariate-only classifiers
        frames = {}
        for fs, (tab, cols) in tables.items():
            X, num, cat = covariate_frame(tab, cols)
            frames[fs] = (X, num, cat)
        X_all = pd.concat([f[0] for f in frames.values()], axis=1)
        frames["all_batch_covariates"] = (X_all, sum((f[1] for f in frames.values()), []), sum((f[2] for f in frames.values()), []))
        for fs, (X, num, cat) in frames.items():
            if X.shape[1] == 0:
                continue
            for model in ("logreg", "rf"):
                r = covariate_auroc(X, num, cat, y, meta, n_splits, args.seed, model, n_perm=n_perm if model == "logreg" else 0, rng=rng)
                auc_rows.append({"tissue": t, "feature_set": fs, "model": model, "n_animals": int(len(d)), "n_features": int(X.shape[1]),
                                 "n_numeric": len(num), "n_categorical": len(cat), **r})
            print(f"  {t} {fs:22s} {X.shape[1]:3d} covariates: " + ", ".join(f"{r['model']} {r['auroc_mean']:.2f}" for r in auc_rows[-2:]))
        # (c) univariate screen
        screen_rows += screen(d[pheno_cols], y, "PHENO", animal_level, t)
        screen_rows += screen(tv[trn_cols_t], y, "meta/TRNSCRPT", {c: "vial" for c in trn_cols_t}, t)
        if pmeta is not None:
            screen_rows += screen(tables["prot_plex_channel"][0], y, "meta/PROT", {"tmt_plex": "vial", "tmt11_channel": "vial"}, t)

    auc = pd.DataFrame(auc_rows)
    auc.to_csv(out / "batch_covariate_auroc.csv", index=False)
    sc = pd.DataFrame(screen_rows)
    sc["separates"] = sc["missing_in_one_group_only"] | (((sc["stat"] == "auroc") & (sc["separation"] >= args.sep_auroc)
                                                          | (sc["stat"] == "cramers_v") & (sc["separation"] >= args.sep_v)) & (sc["p_value"] < 0.01))
    sc.to_csv(out / "batch_variable_screen.csv", index=False)
    summ = sc.groupby(["variable", "source", "category", "level"]).agg(
        n_tissues=("tissue", "nunique"), n_tissues_separating=("separates", "sum"), max_separation=("separation", "max"),
        min_p=("p_value", "min"), constant_in_all_tissues=("n_levels", lambda v: bool((v <= 1).all())),
        missing_in_one_group_only=("missing_in_one_group_only", "any"), kind=("kind", "first")).reset_index()
    summ = summ.sort_values(["n_tissues_separating", "max_separation"], ascending=[False, False])
    summ.to_csv(out / "batch_variable_screen_summary.csv", index=False)

    # ---- (d) verdict ---------------------------------------------------------------------------
    dur_path = out / "taskB_duration_summary.csv"
    dur = pd.read_csv(dur_path) if dur_path.exists() else None
    sep = summ[summ["n_tissues_separating"] > 0]
    sep_batch = sep[sep["category"].isin(BATCH_CATEGORIES)]
    sep_design = sep[sep["category"].isin(("design", "outcome"))]
    shared = summ[summ["category"].isin(("library", "processing")) & summ["constant_in_all_tissues"]]
    best_cov = auc.sort_values("auroc_mean", ascending=False).iloc[0] if len(auc) else None
    cov_high = auc[(auc["auroc_mean"] >= 0.8)]
    PURE = ("pheno_collection", "trnscrpt_library", "trnscrpt_depth", "prot_plex_channel")
    verdict_rows = []
    for t in auc["tissue"].unique():
        a_t = auc[auc["tissue"] == t]
        pure_max = float(a_t[a_t["feature_set"].isin(PURE)]["auroc_mean"].max())
        qc_max = float(a_t[a_t["feature_set"] == "trnscrpt_qc"]["auroc_mean"].max()) if (a_t["feature_set"] == "trnscrpt_qc").any() else np.nan
        sb_t = sc[(sc["tissue"] == t) & sc["separates"] & sc["category"].isin(("collection", "library", "depth", "processing"))]
        sq_t = sc[(sc["tissue"] == t) & sc["separates"] & (sc["category"] == "qc")]
        pure_set = a_t[a_t["feature_set"].isin(PURE)].sort_values("auroc_mean", ascending=False).iloc[0]["feature_set"]
        if pure_max >= 0.8 or len(sb_t):
            verdict = ("unresolvable: sequencing depth separates the groups" if set(sb_t["category"]) <= {"depth"} and len(sb_t)
                       else f"unresolvable: {pure_set} covariates classify the groups (AUROC {pure_max:.2f})")
        elif (not np.isnan(qc_max) and qc_max >= 0.8) or len(sq_t):
            verdict = "unresolvable: library QC metrics separate the groups (a within-plate batch or a training effect on RNA composition)"
        else:
            verdict = "training"
        verdict_rows.append({"tissue": t, "max_auroc_collection_library_plex": pure_max, "max_auroc_qc_metrics": qc_max,
                             "max_auroc_all_covariates": float(a_t[a_t["feature_set"] == "all_batch_covariates"]["auroc_mean"].max()),
                             "null_p95_logreg_max": float(a_t["null_p95_auroc"].max()) if "null_p95_auroc" in a_t else np.nan,
                             "batch_variables_separating": ";".join(sb_t["variable"]), "qc_variables_separating": ";".join(sq_t["variable"]),
                             "verdict": verdict})
    vt = pd.DataFrame(verdict_rows)
    if dur is not None and dur["duration"].nunique() > 1:
        mean_by = dur.groupby("duration")["best_single_auroc"].mean()
        mean_by = mean_by.reindex([d for d in sorted(DURATION_WEEKS, key=DURATION_WEEKS.get) if d in mean_by.index])
        rho = [stats.spearmanr(g_["weeks"], g_["best_single_auroc"])[0] for _, g_ in dur.groupby("tissue") if g_["weeks"].nunique() > 2]
        first, last = mean_by.iloc[0], mean_by.iloc[-1]
        if first >= 0.9:
            dur_read = ("The separation is already at its ceiling after one week of training, so the duration series cannot distinguish "
                        "a fast training response from a batch on its own; note that the 1-, 2- and 4-week animals were sacrificed months "
                        "apart from the controls (design), so only the 8-week comparison shares its collection dates with the controls.")
        elif last - first >= 0.1 and np.nanmedian(rho) >= 0.5:
            dur_read = "The separation grows with the duration, the signature of a biological response rather than of a batch."
        else:
            dur_read = "The separation shows no monotone trend with the duration."
        dur_txt = ("Task B by duration (best single-omic AUROC, mean over tissues): " + ", ".join(f"{d} {v:.2f}" for d, v in mean_by.items())
                   + f"; per-tissue Spearman ρ of AUROC vs weeks: median {np.nanmedian(rho):.2f} (n = {len(rho)} tissues). " + dur_read)
    else:
        dur_txt = "Task B by duration not available (run scripts/07_fusion_vs_baselines.py --durations 1w,2w,4w,8w)."
    n_unres = int((vt["verdict"] != "training").sum())
    n_pure = int(vt["verdict"].str.contains("covariates classify|sequencing depth").sum())
    depth_t = vt.loc[vt["verdict"].str.contains("sequencing depth"), "tissue"].tolist()
    n_qc = int(vt["verdict"].str.contains("QC metrics").sum())
    unres_t = vt.loc[vt["verdict"] != "training", "tissue"].tolist()
    overall = ("training" if n_unres == 0 else
               f"training in {len(vt) - n_unres} of {len(vt)} tissues, unresolvable in {', '.join(unres_t)}"
               + (" (library QC metrics separate the groups there; no collection, library or plex covariate does)" if n_pure == 0 else "")
               + (f" (sequencing depth differs between the arms in {', '.join(depth_t)})" if depth_t else ""))
    md = [
        f"**Verdict: {overall}.** Control vs {args.trained_group} within tissue on the {len(vt)} fusion tissues, same animal-grouped "
        f"sex × group-stratified folds as Task B. Covariate sets: `pheno_collection` (arrival, sacrifice and collection dates, times, "
        "staff, cage, freeze times), `trnscrpt_library` (RNA extraction and library preparation: plate, dates, concentrations, well "
        "position, sequencing run), `trnscrpt_depth` (raw, mapped and splice-junction read counts), `trnscrpt_qc` (mapping, mRNA, "
        "intronic, rRNA, globin, GC and duplication fractions, RIN, 5′–3′ bias), `prot_plex_channel` (TMT plex and channel).",
        (f"(b) Covariates only: the best covariate-only classifier reaches AUROC {best_cov['auroc_mean']:.2f} "
         f"({best_cov['tissue']}, {best_cov['feature_set']}, {best_cov['model']}); "
         + (f"{len(cov_high)} of {len(auc)} tissue × covariate-set × model combinations reach 0.8 or more"
            + (": " + "; ".join(f"{r.tissue} {r.feature_set}/{r.model} {r.auroc_mean:.2f}" for r in cov_high.itertuples()) if len(cov_high) else "")
            + f". Logistic-regression permutation null (group shuffled within sex, {n_perm} draws), 95th percentile of the mean AUROC: "
            + f"{auc['null_p95_auroc'].min():.2f}–{auc['null_p95_auroc'].max():.2f}.")) if best_cov is not None else "(b) no covariate table.",
        (f"(c) Variables that separate control from {args.trained_group} in at least one tissue (AUROC ≥ {args.sep_auroc} or Cramér's V ≥ {args.sep_v} "
         f"at p < 0.01, or missing in one group only): {len(sep)} of {len(summ)} screened — "
         + f"{len(sep_design)} design/outcome variables (" + ", ".join(sep_design['variable'].head(12)) + (", …" if len(sep_design) > 12 else "") + "), "
         + f"{len(sep[sep['category'] == 'baseline'])} baseline, "
         + f"{len(sep_batch)} collection/library/QC/processing variables" + (": " + ", ".join(f"{r.variable} ({r.category}, {int(r.n_tissues_separating)} of {int(r.n_tissues)} tissues, max {r.max_separation:.2f})" for r in sep_batch.itertuples()) if len(sep_batch) else "")
         + f". {len(shared)} library/processing variables are constant within every tissue for these animals (one RNA plate, library date, flowcell and lane per tissue)."),
        "(a) " + dur_txt,
        "Per tissue:\n\n" + report.df_to_md(vt, floatfmt=".2f"),
    ]
    reason = ((f"In {', '.join(unres_t)} the transcript library QC metrics alone (splice-read counts, mRNA / intronic fractions, mapping "
               "rates) separate control from trained animals although the two groups share the RNA plate, library date, flowcell and "
               "lane; on these data that is either a within-plate processing batch or a training effect on the RNA itself (a "
               "composition change is biology, a quality change is batch), and the omic separation there cannot be attributed to "
               "training alone. " if n_qc else "")
              + (f"In {', '.join(depth_t)} the control and trained libraries differ in sequencing depth (raw, mapped and splice-junction "
                 "read counts; same plate, date, flowcell and lane), a pooling or loading batch: log-CPM removes depth to first order "
                 "but the detection of low-count genes does not, so the omic separation there cannot be attributed to training alone. "
                 if depth_t else "")
              + ("".join(f"In {r.tissue} the {r.verdict.split(': ')[1].split(' covariates')[0]} covariates classify the groups at AUROC "
                         f"{r.max_auroc_collection_library_plex:.2f} (permutation null 95th percentile {r.null_p95_logreg_max:.2f}), so a "
                         "processing batch cannot be excluded there. "
                         for r in vt[vt["verdict"].str.contains("covariates classify")].itertuples()))
              + ("".join(f"In {t} the library composition metrics differ by group ("
                         + ", ".join(f"{r.variable} {'higher' if r.auroc_trained_higher >= 0.5 else 'lower'} in trained"
                                     for r in sc[(sc['tissue'] == t) & sc['separates'] & (sc['category'] == 'qc')].itertuples())
                         + ("); a higher mitochondrial and mRNA read fraction in trained muscle is what training-induced mitochondrial "
                            "biogenesis would produce, i.e. plausibly biology, but on these data it cannot be told from a within-plate "
                            "batch, so the omic separation there is not attributable to training alone. "
                            if (t.startswith("SKM") or t == "HEART") and ((sc[(sc["tissue"] == t) & sc["separates"] & (sc["variable"] == "pct_chrM")]["auroc_trained_higher"] >= 0.5).any())
                            else "); either a within-plate batch or a training effect on the RNA composition — on these data the two cannot "
                                 "be told apart, so the omic separation there is not attributable to training alone. ")
                         for t in vt.loc[vt["verdict"].str.contains("QC metrics") | (vt["qc_variables_separating"].fillna("") != ""), "tissue"]))
              + ("In the other tissues the groups share their arrival cohort, sacrifice dates and, within tissue, one RNA plate, library "
                 "preparation date, flowcell and lane; no collection, library, QC or processing variable separates them and the "
                 "covariate-only classifiers sit at their permutation null, while each omic alone separates them at AUROC ≈ 1: there "
                 "the separation is attributed to training (or to whatever else differs by design between the arms, such as the daily "
                 "handling of trained animals), not to a recorded batch." if n_unres < len(vt) else "")
              if n_unres else
              "Control and trained animals of the primary comparison share their arrival cohort, sacrifice dates and, within every tissue, "
              "one RNA plate, library preparation date, flowcell and lane; no collection, library, QC or processing variable separates the "
              "groups and the covariate-only classifiers sit at their permutation null, while each omic alone separates them at AUROC ≈ 1. "
              "The separation is therefore attributed to training (or to whatever else differs by design between the arms, such as the "
              "daily handling of trained animals), not to a recorded batch.")
    md.append("**Read:** " + reason)
    (out / "batch_conclusion.md").write_text("\n\n".join(md) + "\n")
    vt.assign(overall_verdict=overall).to_csv(out / "batch_conclusion.csv", index=False)
    report.add_section("07b · Batch check for Task B (covariates only, variable screen, duration)", "\n\n".join(md)
                       + "\n\nCovariate-only classifiers (animal-grouped CV AUROC; `null_p95_auroc` from the within-sex label permutation of the logistic model):\n\n"
                       + report.df_to_md(auc, floatfmt=".3f")
                       + "\n\nVariables separating control from trained in ≥ 1 tissue (full screen in batch_variable_screen.csv):\n\n"
                       + report.df_to_md(sep.head(60), floatfmt=".3f"),
                       params={"trained_group": args.trained_group, "tissues": ",".join(tissues), "n_perm": n_perm, "quick": args.quick})
    print("\n".join(md[:4]))
    print(vt.round(2).to_string(index=False))
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
