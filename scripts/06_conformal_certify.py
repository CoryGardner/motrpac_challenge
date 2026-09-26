#!/usr/bin/env python
"""Phase 06 — conformal prediction sets and the certified panel size.

Per outer animal-grouped fold: training animals → (fit, calibration) animals; test animals untouched.
  a) full-feature classifier → LAC and APS sets at α ∈ ALPHAS, calibrated two ways on the same
     calibration animals — pooled vials, and one vial per animal (exchangeability hygiene) — and
     evaluated on the test animals: coverage, set size, per-tissue coverage.
  b) certificate: for each k on the grid, the k-panel classifier (round-robin selection + the
     classifier, preprocessing fit on the fit animals) → calibration error → one-sided
     Clopper–Pearson upper bound at level δ → fixed-sequence testing from the largest k down →
     certified k. Two losses: the loss actually certified (one vial per animal with
     --one-per-animal, else pooled vials) and the per-animal "any tissue wrong" loss. The
     fit/calibration split is repeated --n-repeats times with different seeds so the certified k
     comes with a distribution; the test error at the certified k is the validity check.
Finally a deployment certificate on all animals (fit / calibration by animal), also repeated.

Outputs (results/06_conformal/<assay>/): coverage.csv, per_class_coverage.csv,
coverage_per_tissue_alpha0.1.png, certificate_by_fold.csv, certificate_validity.csv,
certificate_distribution.csv, certificate_final.csv, certificate_final_repeats.csv,
certificate_final.png, sizing_table.csv
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

from motrpac import cli, config as C, conformal as cp, io, models, plots, report
from motrpac.splits import assert_no_group_leak, fit_calibration_split, grouped_kfold

SIZING_PAIRS = [(0.05, 0.05), (0.10, 0.05), (0.10, 0.10), (0.20, 0.10)]


def proba_on(est, X, classes):
    p = est.predict_proba(X)
    order = [list(est.classes_).index(c) for c in classes]
    return p[:, order]


def split_train(meta, tr, cal_frac, cal_animals, seed):
    """(fit, calibration) by animal: an absolute number of calibration animals if given, else a fraction."""
    if cal_animals is None:
        return fit_calibration_split(meta, tr, cal_frac, seed)
    tr = np.asarray(tr)
    rng = np.random.default_rng(seed)
    pids = meta.iloc[tr]["pid"].astype(str).to_numpy()
    uniq = np.unique(pids)
    rng.shuffle(uniq)
    n_cal = min(max(1, int(cal_animals)), len(uniq) - 1)
    cal = set(uniq[:n_cal])
    is_cal = np.array([p in cal for p in pids])
    fit_idx, cal_idx = tr[~is_cal], tr[is_cal]
    assert_no_group_leak(meta, fit_idx, cal_idx)
    return fit_idx, cal_idx


def panel_family(X_fit, y_fit, grid, prefilter, model):
    """Preprocessing (impute → variance prefilter → scale → round-robin scores) fit once on the fit
    animals; then one k-panel classifier per k. Same model as the quick pipeline of phase 05
    (logreg_l2: C = 0.1). Returns predict(k, X)."""
    if model != "logreg_l2":
        raise ValueError("the certificate is defined for --model logreg_l2 (round-robin panel + L2 logistic regression)")
    imp = SimpleImputer(strategy="median").fit(X_fit)
    pre = models.VarianceTopK(prefilter).fit(imp.transform(X_fit))
    sc = StandardScaler().fit(pre.transform(imp.transform(X_fit)))
    Xs = sc.transform(pre.transform(imp.transform(X_fit)))
    sel = models.RoundRobinSelector(k=max(grid)).fit(Xs, y_fit)
    order = sel.order_
    fam = {}
    for k in grid:
        idx = np.sort(order[:k])
        fam[k] = (idx, LogisticRegression(solver="lbfgs", C=0.1, max_iter=3000).fit(Xs[:, idx], y_fit))

    def predict(k, X):
        Z = sc.transform(pre.transform(imp.transform(X)))
        idx, clf = fam[k]
        return clf.predict(Z[:, idx])

    def proba(k, X):
        """Class probabilities of the k-panel classifier and its class order (for --save-scores)."""
        Z = sc.transform(pre.transform(imp.transform(X)))
        idx, clf = fam[k]
        return clf.predict_proba(Z[:, idx]), list(clf.classes_)
    predict.proba = proba
    return predict


def reorder_columns(p, cls, classes):
    """Probability columns in `classes` order; a class the estimator never saw gets a zero column."""
    out = np.zeros((p.shape[0], len(classes)))
    for j, c in enumerate(cls):
        if c in classes:
            out[:, classes.index(c)] = p[:, j]
    return out


def score_tables(om, fold, fit_idx, cal_idx, cal_opa, te, est, p_te, fam, grid, classes, y_idx):
    """--save-scores rows for one outer fold: test-vial probabilities and LAC calibration scores (pooled
    vials and one vial per animal) for the full model and the k = 20 / 50 panels of the same fit animals."""
    meta = om.meta[["pid", "tissue", "sex", "group"]].copy()
    meta.index.name = "viallabel"
    prob_rows, cal_rows = [], []
    for mname in ["full"] + [f"k{k}" for k in (20, 50) if k in grid]:
        def proba_rows(idx):
            if mname == "full":
                return proba_on(est, om.X.to_numpy(dtype=float)[idx], classes)
            return reorder_columns(*fam.proba(int(mname[1:]), om.X.to_numpy(dtype=float)[idx]), classes)
        pt = p_te if mname == "full" else proba_rows(te)
        df = meta.iloc[te].reset_index().assign(fold=fold, model=mname)
        for j, c in enumerate(classes):
            df[f"p_{c}"] = pt[:, j]
        prob_rows.append(df)
        for mode, cidx in (("pooled", cal_idx), ("one_per_animal", cal_opa)):
            pc = proba_rows(cidx)
            dc = meta.iloc[cidx].reset_index().assign(fold=fold, model=mname, calibration=mode,
                                                      score_lac=cp.lac_scores(pc, y_idx[cidx]))
            cal_rows.append(dc)
    return pd.concat(prob_rows, ignore_index=True), pd.concat(cal_rows, ignore_index=True)


def certify_losses(predict, X, y, g, cal_idx, cal_opa, te, grid, alpha, delta, one_per_animal):
    """Certified k under the chosen loss and under the per-animal any-wrong loss; test error per k."""
    errors, errors_animal, test_err = {}, {}, {}
    m_opa = np.isin(cal_idx, cal_opa)
    for k in grid:
        wrong = predict(k, X[cal_idx]) != y[cal_idx]
        errors[k] = (int(wrong[m_opa].sum()), int(m_opa.sum())) if one_per_animal else (int(wrong.sum()), len(cal_idx))
        any_wrong = pd.Series(wrong).groupby(g[cal_idx]).any()
        errors_animal[k] = (int(any_wrong.sum()), int(len(any_wrong)))
        test_err[k] = float((predict(k, X[te]) != y[te]).mean()) if len(te) else np.nan
    k1, t1 = cp.certify_panel_size(errors, alpha, delta)
    k2, t2 = cp.certify_panel_size(errors_animal, alpha, delta)
    t1["loss"] = "one_vial_per_animal" if one_per_animal else "pooled_vials"
    t2["loss"] = "animal_any_tissue_wrong"
    return k1, k2, pd.concat([t1, t2], ignore_index=True), test_err, errors


def main() -> None:
    ap = cli.common_parser("Conformal sets and certified panel size")
    ap.add_argument("--label", default="tissue")
    ap.add_argument("--model", default="logreg_l2")
    ap.add_argument("--grid", default=",".join(map(str, C.PANEL_GRID)))
    ap.add_argument("--alpha", type=float, default=0.05, help="target error rate for the certificate")
    ap.add_argument("--delta", type=float, default=0.05, help="confidence level for the certificate (1-δ)")
    ap.add_argument("--cal-frac", type=float, default=0.3)
    ap.add_argument("--cal-animals", type=int, default=None, help="absolute number of calibration animals (overrides --cal-frac)")
    ap.add_argument("--n-repeats", type=int, default=1, help="repeat the fit/calibration split with different seeds")
    ap.add_argument("--one-per-animal", action="store_true", help="certify and calibrate with one vial per animal")
    ap.add_argument("--conditional", action="store_true",
                    help="also report Mondrian (class-conditional) and marginal-floor Mondrian sets for the pooled calibration")
    ap.add_argument("--alpha-delta-grid", default="0.05:0.05,0.10:0.05,0.10:0.10,0.20:0.05,0.20:0.10,0.05:0.10",
                    help="(alpha:delta) pairs re-evaluated on the same calibration errors")
    ap.add_argument("--k-of-interest", type=int, default=20)
    ap.add_argument("--save-scores", action="store_true",
                    help="also write per-vial test probabilities and calibration scores (full model and k = 20 / 50 "
                         "panels, repeat 0 of every fold): scores_test_probs.csv, scores_calibration.csv, classes.json")
    args = ap.parse_args()
    cli.banner("06_conformal_certify", args)
    source = cli.resolve_source(args.assay, args.source)
    out = cli.outdir(f"06_conformal/{args.assay}", args.out)
    pheno = io.load_pheno()
    om = io.stack_tissues(args.assay, tissues=cli.parse_tissues(args.tissues), source=source, join="inner", pheno=pheno,
                          complete=args.complete_features,
                          drop_incomplete_samples=args.drop_incomplete_samples)
    grid = [int(k) for k in args.grid.split(",")]
    n_splits = 3 if args.quick else args.n_splits
    prefilter = min(args.prefilter, 1000) if args.quick else args.prefilter
    n_repeats = max(1, args.n_repeats)

    X = om.X.to_numpy(dtype=float)
    y = om.meta[args.label].astype(str).to_numpy()
    classes = sorted(np.unique(y))
    y_idx = np.array([classes.index(v) for v in y])
    g = om.groups()
    grid = [k for k in grid if k <= X.shape[1]]

    cov_rows, pc_rows, cert_rows, valid_rows, err_store = [], [], [], [], []
    score_probs, score_cal = [], []
    ad_grid =[(float(a), float(d)) for a, d in (pair.split(":") for pair in args.alpha_delta_grid.split(",") if pair.strip())]
    for fold, (tr, te) in enumerate(grouped_kfold(om.meta, args.label, n_splits, args.seed)):
        for rep in range(n_repeats):
            seed_r = args.seed + 1000 * fold + rep
            fit_idx, cal_idx = split_train(om.meta, tr, args.cal_frac, args.cal_animals, seed_r)
            cal_opa = cp.one_per_group(cal_idx, g, seed_r)
            if rep == 0:
                # (a) full-feature model → conformal sets, calibrated pooled and one-per-animal on the same animals
                est = models.fit_tuned(args.model, X[fit_idx], y[fit_idx], g[fit_idx], k=None, prefilter=prefilter,
                                       quick=args.quick, seed=args.seed)
                p_te = proba_on(est, X[te], classes)
                for mode, cidx in (("pooled", cal_idx), ("one_per_animal", cal_opa)):
                    p_cal = proba_on(est, X[cidx], classes)
                    for method in ("lac", "aps"):
                        for alpha in C.ALPHAS:
                            variants = [("marginal", False)] + ([("mondrian", True), ("floored", "floored")] if (args.conditional and mode == "pooled") else [])
                            for cname, cond in variants:
                                r = cp.calibrate_and_evaluate(p_cal, y_idx[cidx], p_te, y_idx[te], classes, alpha, method,
                                                              extra_test=om.meta.iloc[te][["sex", "group"]], conditional=cond)
                                cov_rows.append({"fold": fold, "calibration": mode, "conformal": cname, "method": method, "alpha": alpha,
                                                 "qhat": r.qhat, "n_cal": r.n_cal, "n_cal_animals": int(om.meta.iloc[cidx]["pid"].nunique()),
                                                 "n_fit_animals": int(om.meta.iloc[fit_idx]["pid"].nunique()),
                                                 "n_test_animals": int(om.meta.iloc[te]["pid"].nunique()), **r.overall})
                                pc_rows.append(r.per_class.assign(fold=fold, calibration=mode, conformal=cname, method=method, alpha=alpha))
                lac = [c for c in cov_rows if c["fold"] == fold and c["method"] == "lac" and c["alpha"] == 0.1 and c["conformal"] == "marginal"]
                print(f"  fold {fold}: LAC α=0.1 coverage pooled={lac[0]['coverage']:.3f} one-per-animal={lac[1]['coverage']:.3f}")
                if args.save_scores:
                    fam0 = panel_family(X[fit_idx], y[fit_idx], grid, prefilter, args.model)
                    sp, sc_ = score_tables(om, fold, fit_idx, cal_idx, cal_opa, te, est, p_te, fam0, grid, classes, y_idx)
                    score_probs.append(sp)
                    score_cal.append(sc_)
            # (b) certificate on calibration animals; validity on test animals
            predict = panel_family(X[fit_idx], y[fit_idx], grid, prefilter, args.model)
            k1, k2, table, test_err, errs = certify_losses(predict, X, y, g, cal_idx, cal_opa, te, grid, args.alpha, args.delta,
                                                           args.one_per_animal)
            err_store.append(errs)
            table["fold"], table["repeat"] = fold, rep
            table["test_err"] = table["k"].map(test_err)
            cert_rows.append(table)
            valid_rows.append({"fold": fold, "repeat": rep, "n_fit_animals": int(om.meta.iloc[fit_idx]["pid"].nunique()),
                               "n_cal_animals": int(om.meta.iloc[cal_idx]["pid"].nunique()),
                               "n_cal_loss": int(len(cal_opa)) if args.one_per_animal else int(len(cal_idx)),
                               "certified_k": k1, "certified_k_animal_any_wrong": k2,
                               "test_err_at_certified_k": test_err.get(k1, np.nan) if k1 else np.nan,
                               "test_err_le_alpha": (test_err.get(k1, 1.0) <= args.alpha) if k1 else None})
        ks = [v["certified_k"] for v in valid_rows if v["fold"] == fold]
        print(f"  fold {fold}: certified k over {n_repeats} repeats: {ks}")

    cov = pd.DataFrame(cov_rows)
    cov.to_csv(out / "coverage.csv", index=False)
    if args.save_scores:
        pd.concat(score_probs, ignore_index=True).to_csv(out / "scores_test_probs.csv", index=False)
        pd.concat(score_cal, ignore_index=True).to_csv(out / "scores_calibration.csv", index=False)
        (out / "classes.json").write_text(json.dumps(classes))
        print(f"  wrote scores_test_probs.csv, scores_calibration.csv, classes.json to {out}")
    pc = pd.concat(pc_rows, ignore_index=True)
    pc.to_csv(out / "per_class_coverage.csv", index=False)
    cert = pd.concat(cert_rows, ignore_index=True)
    cert.to_csv(out / "certificate_by_fold.csv", index=False)
    valid = pd.DataFrame(valid_rows)
    valid.to_csv(out / "certificate_validity.csv", index=False)

    cov_all = cov
    cov = cov_all[cov_all["conformal"] == "marginal"]
    pc_all = pc
    pc = pc_all[pc_all["conformal"] == "marginal"]
    mond_md = ""
    if args.conditional:
        mc = cov_all[(cov_all["calibration"] == "pooled")].groupby(["conformal", "method", "alpha"]).agg(
            coverage=("coverage", "mean"), avg_set_size=("avg_set_size", "mean"), frac_empty=("frac_empty", "mean")).reset_index()
        mp = pc_all[(pc_all["calibration"] == "pooled") & (pc_all["method"] == "lac") & (pc_all["alpha"] == 0.1)]
        mpt = mp.groupby(["y_true", "conformal"]).agg(coverage=("coverage", "mean"), set_size=("avg_set_size", "mean")).unstack("conformal")
        mpt.columns = [f"{a}_{b}" for a, b in mpt.columns]
        mpt = mpt.reset_index()
        mpt.to_csv(out / "per_tissue_marginal_vs_mondrian_alpha0.1.csv", index=False)
        mc.to_csv(out / "coverage_marginal_vs_mondrian.csv", index=False)
        mond_md = ("Marginal vs Mondrian (class-conditional) vs marginal-floor Mondrian sets (each class quantile floored at the "
                   "marginal one), pooled calibration:\n\n" + report.df_to_md(mc, floatfmt=".3f")
                   + "\n\nPer-tissue coverage and set size at α = 0.1 (LAC), marginal vs Mondrian vs floored:\n\n" + report.df_to_md(mpt, floatfmt=".2f")
                   + "\n\nRead: Mondrian quantiles restore coverage inside every tissue by lowering the thresholds of the tissues the "
                   "model confuses (muscle, brain regions, vena cava): those tissues get a singleton set where the marginal threshold "
                   "left them empty, and the tissues that were easy lose coverage towards 1 − α. "
                   "The floored variant never tightens a tissue below the marginal threshold, so it contains both the marginal "
                   "and the Mondrian sets: it keeps the marginal guarantee that plain Mondrian can lose on the easy tissues, "
                   "at the price of slightly larger sets there.")
    cov_summary = cov.groupby(["calibration", "method", "alpha"]).agg(
        coverage=("coverage", "mean"), coverage_sd=("coverage", "std"), avg_set_size=("avg_set_size", "mean"),
        frac_singleton=("frac_singleton", "mean"), frac_empty=("frac_empty", "mean"), n_cal=("n_cal", "mean")).reset_index()
    pc10 = pc[(pc["method"] == "lac") & (pc["alpha"] == 0.1)]
    pc_summary = pc10.groupby(["y_true", "calibration"]).agg(coverage=("coverage", "mean"), avg_set_size=("avg_set_size", "mean"),
                                                             n=("n", "sum")).reset_index()
    pc_wide = pc_summary.pivot(index="y_true", columns="calibration", values="coverage").reset_index()
    under = pc_wide[pc_wide[["pooled", "one_per_animal"]].min(axis=1) < 0.9 - 0.05].sort_values("pooled")
    fig_pc = plots.coverage_bars(pc_summary[pc_summary["calibration"] == ("one_per_animal" if args.one_per_animal else "pooled")],
                                 "y_true", f"{args.assay}: per-tissue coverage (LAC, α=0.1)",
                                 out / "coverage_per_tissue_alpha0.1.png", alpha=0.1)

    # distribution of the certified k over folds × repeats
    dist = []
    for k in grid:
        dist.append({"k": k, "frac_certified": float(valid["certified_k"].apply(lambda v: v is not None and not pd.isna(v) and v <= k).mean()),
                     "frac_certified_animal_any_wrong": float(valid["certified_k_animal_any_wrong"].apply(lambda v: v is not None and not pd.isna(v) and v <= k).mean())})
    dist = pd.DataFrame(dist)
    dist.to_csv(out / "certificate_distribution.csv", index=False)
    k_half = dist.loc[dist["frac_certified"] >= 0.5, "k"].min() if (dist["frac_certified"] >= 0.5).any() else None
    ks = valid["certified_k"].dropna()
    n_none = int(valid["certified_k"].isna().sum())

    # deployment certificate on all animals, repeated
    all_idx = np.arange(len(y))
    final_rows, final_table = [], None
    for rep in range(n_repeats):
        fit_idx, cal_idx = split_train(om.meta, all_idx, args.cal_frac, args.cal_animals, args.seed + 999 + rep)
        cal_opa = cp.one_per_group(cal_idx, g, args.seed + 999 + rep)
        predict = panel_family(X[fit_idx], y[fit_idx], grid, prefilter, args.model)
        k1, k2, table, _, _ = certify_losses(predict, X, y, g, cal_idx, cal_opa, np.array([], dtype=int), grid,
                                             args.alpha, args.delta, args.one_per_animal)
        table["repeat"] = rep
        final_rows.append({"repeat": rep, "n_fit_animals": int(om.meta.iloc[fit_idx]["pid"].nunique()),
                           "n_cal_animals": int(om.meta.iloc[cal_idx]["pid"].nunique()), "certified_k": k1,
                           "certified_k_animal_any_wrong": k2})
        if rep == 0:
            final_table = table
            final_k = k1
    final = pd.DataFrame(final_rows)
    final.to_csv(out / "certificate_final_repeats.csv", index=False)
    final_table.to_csv(out / "certificate_final.csv", index=False)
    loss_name = "one_vial_per_animal" if args.one_per_animal else "pooled_vials"
    fig_cert = plots.certificate_plot(final_table[final_table["loss"] == loss_name].drop(columns=["loss", "repeat"]), args.alpha,
                                      f"{args.assay}: certified panel size (α={args.alpha}, δ={args.delta}, "
                                      f"n_cal={int(final_table.loc[final_table['loss'] == loss_name, 'n_cal'].iloc[0])})",
                                      out / "certificate_final.png")

    # (alpha, delta) grid re-evaluated on the stored calibration errors (same fits, same splits)
    ad_rows = []
    for a, d in ad_grid:
        ks_ad = [cp.certify_panel_size(e, a, d)[0] for e in err_store]
        ks_ok = [k for k in ks_ad if k is not None]
        ad_rows.append({"alpha": a, "delta": d, "n_splits": len(ks_ad),
                        f"frac_certified_k{args.k_of_interest}": float(np.mean([(k is not None and k <= args.k_of_interest) for k in ks_ad])),
                        "frac_any_certified": float(np.mean([k is not None for k in ks_ad])),
                        "median_certified_k": float(np.median(ks_ok)) if ks_ok else np.nan})
    ad = pd.DataFrame(ad_rows)
    ad.to_csv(out / "certificate_alpha_delta_grid.csv", index=False)

    # zero-error sizing: n >= ln δ / ln(1 − α)
    sizing = pd.DataFrame([{"alpha": a, "delta": d, "n_zero_error_needed": int(np.ceil(np.log(d) / np.log(1 - a)))}
                           for a, d in SIZING_PAIRS])
    sizing["n_cal_this_run"] = int(valid["n_cal_loss"].iloc[0])
    sizing.to_csv(out / "sizing_table.csv", index=False)
    n_0505 = int(sizing.loc[(sizing.alpha == 0.05) & (sizing.delta == 0.05), "n_zero_error_needed"].iloc[0])

    point_k = None
    pc_path = C.RESULTS_DIR / "05_panels" / args.assay / "panel_curve.csv"
    if pc_path.exists():
        curve = pd.read_csv(pc_path).groupby("k")["balanced_accuracy"].mean()
        point_k = int(curve[curve >= 1 - args.alpha].index.min()) if (curve >= 1 - args.alpha).any() else None

    body = [
        f"{om.notes[-1]}; model `{args.model}`; {n_splits} outer animal-grouped folds; calibration = "
        + (f"{args.cal_animals} animals" if args.cal_animals is not None else f"{args.cal_frac:.0%} of training animals")
        + f" ({int(valid['n_fit_animals'].iloc[0])} fit / {int(valid['n_cal_animals'].iloc[0])} calibration animals per fold), "
        f"{n_repeats} repeat(s) of the fit/calibration split; certified loss = **{loss_name}**; α = {args.alpha}, δ = {args.delta}.",
        "Coverage and set size on held-out test animals (mean over folds), calibrated pooled vs one vial per animal on "
        "the same calibration animals:\n\n" + report.df_to_md(cov_summary, floatfmt=".3f"),
        "Sanity: LAC coverage should sit near 1 − α. Far above → sets are trivially large; far below → leak or shift. "
        "Pooled calibration treats vials of one animal as exchangeable with vials of other animals; it is the larger-n "
        "but only approximately valid option, and its coverage on test *animals* is the check.",
        "Per-tissue (conditional) coverage at α = 0.1, LAC:\n\n" + report.df_to_md(pc_wide, floatfmt=".2f")
        + (f"\n\nUnder-covered tissues (coverage < 0.85 under either calibration): **{', '.join(under['y_true']) or 'none'}**."),
        report.figure_md(fig_pc, "Per-tissue coverage at α = 0.1 — under-covered tissues are the confusable ones"),
        f"Certified panel size per fold × repeat (α={args.alpha}, δ={args.delta}; Clopper–Pearson upper bound, fixed-sequence "
        "testing from large k down; `certified_k_animal_any_wrong` uses the per-animal loss 'any tissue of this animal wrong'):\n\n"
        + report.df_to_md(valid),
        "Fraction of fold × repeat splits certifying each k (the distribution of the certificate):\n\n" + report.df_to_md(dist, floatfmt=".2f"),
        f"Certified k: min {int(ks.min()) if len(ks) else 'n/a'}, median {ks.median() if len(ks) else 'n/a'}, max {int(ks.max()) if len(ks) else 'n/a'}"
        f"{f', no k certified in {n_none} split(s)' if n_none else ''}; smallest k certified in ≥ half of the splits: **{k_half}**.",
        f"Validity: fraction of splits where the test error at the certified k was ≤ α: "
        f"**{valid['test_err_le_alpha'].dropna().astype(bool).mean() if valid['test_err_le_alpha'].notna().any() else float('nan'):.2f}** "
        f"(should be ≥ 1 − δ = {1 - args.delta:.2f} in the long run).",
        f"Deployment certificate on all animals ({int(final['n_fit_animals'].iloc[0])} fit / {int(final['n_cal_animals'].iloc[0])} "
        f"calibration), {n_repeats} repeat(s): certified k = {sorted(final['certified_k'].dropna().astype(int).tolist())}"
        f"{' (none in ' + str(int(final['certified_k'].isna().sum())) + ')' if final['certified_k'].isna().any() else ''}; "
        f"per-animal any-wrong loss: {sorted(final['certified_k_animal_any_wrong'].dropna().astype(int).tolist())}. First repeat:\n\n"
        + report.df_to_md(final_table[final_table["loss"] == loss_name].drop(columns=["repeat"]), floatfmt=".3f"),
        report.figure_md(fig_cert, "Calibration error and upper confidence bound vs panel size (first repeat)"),
        mond_md,
        f"(α, δ) grid on the same {len(err_store)} calibration splits and fits (certified loss = {loss_name}): fraction certifying "
        f"k = {args.k_of_interest}, fraction certifying anything, median certified k:\n\n" + report.df_to_md(ad, floatfmt=".2f"),
        "Zero-error sizing, n ≥ ln δ / ln(1 − α): the number of calibration units needed to certify a panel that makes "
        "**no** calibration error:\n\n" + report.df_to_md(sizing),
        f"**Animals needed for (α, δ) = (0.05, 0.05): {n_0505} calibration animals with zero errors** — more animals than "
        f"the whole study has per tissue; with one vial per animal this study can certify at (0.10, 0.10) at best.",
        f"Why the certified k exceeds the point-estimate k ({point_k if point_k is not None else 'n/a'} from the phase-05 curve at "
        f"mean balanced accuracy ≥ {1 - args.alpha:.2f}): the curve reports an *average* error over folds, while the certificate "
        f"demands that the *upper confidence bound* of the error on the calibration animals be ≤ α, and with "
        f"{int(valid['n_cal_loss'].iloc[0])} calibration units the bound sits {'about ' + format(cp.clopper_pearson_upper(0, int(valid['n_cal_loss'].iloc[0]), args.delta), '.2f') + ' even with zero errors'} above the point estimate. "
        "Fixed-sequence testing from the largest k down stops at the first k that fails, so one extra calibration error at "
        "a small k ends the certificate there even if smaller panels happen to do well on other animals. The gap is the "
        "price of a guarantee rather than an estimate, and it shrinks only with more calibration animals.",
        "Caveat: the pooled-vial certificate treats the ~18 vials of one animal as independent draws; a wrong tissue call is "
        "rarely independent across an animal's vials, so the pooled bound is optimistic. The one-vial-per-animal loss and the "
        "per-animal any-wrong loss bracket the honest answer.",
    ]
    subset = f", {om.meta['tissue'].nunique()} tissues" if args.tissues else ""
    report.add_section(f"06 · Conformal sets and certified panel ({args.assay}/{source}{subset})", "\n\n".join(b for b in body if b),
                       params={"assay": args.assay, "model": args.model, "alpha": args.alpha, "delta": args.delta,
                               "one_per_animal": args.one_per_animal, "cal_animals": args.cal_animals, "cal_frac": args.cal_frac,
                               "n_repeats": n_repeats, "grid": args.grid, "quick": args.quick})
    print(cov_summary.to_string())
    print(dist.to_string())
    print(f"certified k distribution: {ks.tolist()}; k_half={k_half}; deployment: {final['certified_k'].tolist()}; wrote {out}")


if __name__ == "__main__":
    main()
