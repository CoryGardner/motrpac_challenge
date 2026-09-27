#!/usr/bin/env python
"""Phase 08 — does the fingerprint (and its coverage guarantee) survive a shift?

Splits: train on one sex → test on the other; hold out a time point (8w, 1w); train on
sedentary controls → test on trained animals. For each: full-feature and k-panel classifiers,
conformal sets calibrated on *source* animals (pooled vials, or one vial per animal with
--one-per-animal), evaluated on held-out source animals (in-distribution reference) and on the
target. With --recalibrate-target N, N target animals are held out as a calibration set and the
remaining target animals are the test set; their coverage under target calibration is compared
with the source-only calibration on the same remaining animals (repeated draws, mean reported).

Real-data wrinkle: OVARY and TESTES exist in only one sex, so under the sex shift some target
classes were never seen in training. Calibration samples whose class is absent from the fitted
estimator are dropped and counted (`n_cal_dropped_unseen`); test samples with such a label are
reported as unseen and count as not covered.

Outputs (results/08_shift/<assay>/): shift_table.csv, shift_per_class.csv, shift_recalibration.csv,
coverage_shift.png
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from tfp import cli, config as C, conformal as cp, io, models, plots, report
from tfp.splits import assert_no_group_leak, fit_calibration_split, leave_one_group_out, leave_one_sex_out, \
    train_controls_test_trained


def coverage_with(est, cls, X_te, y_te, qhat, method):
    sets = cp.predict_sets(est.predict_proba(X_te), qhat, method)
    cov, n_dropped = cp.coverage_seen(sets, y_te, cls)
    return cov, float(sets.sum(axis=1).mean()), n_dropped


def eval_split(name, tr, te, X, y, g, meta, k, alpha, model, prefilter, quick, seed, method="lac",
               one_per_animal=False, recal_ns=(), recal_repeats=10):
    fit_idx, cal_idx = fit_calibration_split(meta, tr, 0.3, seed)
    cal_used = cp.one_per_group(cal_idx, g, seed) if one_per_animal else cal_idx
    rows, per_class, recal_rows, variant_rows = [], [], [], []
    for arm, kk in (("full", None), (f"panel_k{k}", k)):
        est = models.fit_tuned(model, X[fit_idx], y[fit_idx], g[fit_idx], k=kk, prefilter=prefilter, quick=quick, seed=seed)
        cls = list(est.classes_)
        unseen = sorted(set(y[te]) - set(cls))
        scores, n_cal_dropped = cp.calibration_scores(est.predict_proba(X[cal_used]), y[cal_used], cls, method)
        qhat = cp.conformal_quantile(scores, alpha)
        p_te = est.predict_proba(X[te])
        sets_te = cp.predict_sets(p_te, qhat, method)
        yi_te, seen_mask = cp.class_indices(y[te], cls)
        covered = np.zeros(len(te), dtype=bool)
        covered[seen_mask] = sets_te[np.flatnonzero(seen_mask), yi_te[seen_mask]]
        size = sets_te.sum(axis=1)
        # class-conditional variants from the same calibration scores: Mondrian and marginal-floor Mondrian
        yi_cal, seen_cal = cp.class_indices(y[cal_used], cls)
        q_variants = {"mondrian": cp.conformal_quantile_per_class(scores, yi_cal[seen_cal], alpha, len(cls), fallback=qhat),
                      "floored": cp.conformal_quantile_per_class(scores, yi_cal[seen_cal], alpha, len(cls), fallback=qhat, floor=qhat)}
        variant_sets = {"marginal": sets_te}
        variant_cov = {}
        for vname, qv in q_variants.items():
            sv = cp.predict_sets_conditional(p_te, qv, method)
            cv_ = np.zeros(len(te), dtype=bool)
            cv_[seen_mask] = sv[np.flatnonzero(seen_mask), yi_te[seen_mask]]
            variant_sets[vname], variant_cov[vname] = sv, cv_
        variant_cov["marginal"] = covered
        y_pred = est.predict(X[te])
        acc_all = float(np.mean(y_pred == y[te]))
        bal_seen = models.balanced_accuracy_score(y[te][seen_mask], y_pred[seen_mask]) if seen_mask.any() else np.nan
        # source in-distribution reference: calibrate on half of the calibration animals, test on the other half
        s_fit, s_cal = fit_calibration_split(meta, cal_idx, 0.5, seed + 1)
        s_fit_used = cp.one_per_group(s_fit, g, seed + 1) if one_per_animal else s_fit
        src_scores, n_src_fit_dropped = cp.calibration_scores(est.predict_proba(X[s_fit_used]), y[s_fit_used], cls, method)
        q_src = cp.conformal_quantile(src_scores, alpha)
        cov_src, size_src, n_src_cal_dropped = coverage_with(est, cls, X[s_cal], y[s_cal], q_src, method)
        n_src_dropped = n_src_fit_dropped + n_src_cal_dropped
        # empty-set rate (LAC) and APS sets under the same shift, same calibration animals
        frac_empty_te = float((sets_te.sum(axis=1) == 0).mean())
        sets_src_lac = cp.predict_sets(est.predict_proba(X[s_cal]), q_src, method)
        frac_empty_src = float((sets_src_lac.sum(axis=1) == 0).mean())
        aps_scores_cal, _ = cp.calibration_scores(est.predict_proba(X[cal_used]), y[cal_used], cls, "aps")
        q_aps = cp.conformal_quantile(aps_scores_cal, alpha)
        cov_aps_te, size_aps_te, _ = coverage_with(est, cls, X[te], y[te], q_aps, "aps")
        aps_src_scores, _ = cp.calibration_scores(est.predict_proba(X[s_fit_used]), y[s_fit_used], cls, "aps")
        _, size_aps_src, _ = coverage_with(est, cls, X[s_cal], y[s_cal], cp.conformal_quantile(aps_src_scores, alpha), "aps")
        if n_cal_dropped or n_src_dropped:
            print(f"  {name} {arm}: dropped {n_cal_dropped} calibration and {n_src_dropped} source-check samples "
                  f"whose class is absent from the fitted estimator ({', '.join(sorted(set(y[cal_idx]) - set(cls))) or '-'})")
        row = {"split": name, "arm": arm, "n_train_animals": int(meta.iloc[tr]["pid"].nunique()),
               "n_test_animals": int(meta.iloc[te]["pid"].nunique()), "n_test": len(te),
               "n_cal": int(len(scores)), "n_cal_animals": int(meta.iloc[cal_idx]["pid"].nunique()),
               "n_cal_dropped_unseen": int(n_cal_dropped), "n_src_dropped_unseen": int(n_src_dropped),
               "unseen_classes": ",".join(unseen), "accuracy_all": acc_all, "bal_acc_seen": bal_seen,
               "coverage_source_id": cov_src, "avg_set_size_source": size_src,
               "coverage_target_all": float(covered.mean()),
               "coverage_target_seen": float(covered[seen_mask].mean()) if seen_mask.any() else np.nan,
               "avg_set_size_target": float(size.mean()), "qhat": qhat, "alpha": alpha,
               "lac_frac_empty_source": frac_empty_src, "lac_frac_empty_target": frac_empty_te,
               "aps_coverage_target_seen": cov_aps_te, "aps_set_size_source": size_aps_src, "aps_set_size_target": size_aps_te}
        for vname in ("mondrian", "floored"):
            row[f"coverage_target_seen_{vname}"] = float(variant_cov[vname][seen_mask].mean()) if seen_mask.any() else np.nan
            row[f"avg_set_size_target_{vname}"] = float(variant_sets[vname].sum(axis=1).mean())
        for vname, sv in variant_sets.items():
            variant_rows.append(pd.DataFrame({"split": name, "arm": arm, "conformal": vname, "y_true": y[te], "covered": variant_cov[vname],
                                              "size": sv.sum(axis=1), "seen": seen_mask})
                                .groupby(["split", "arm", "conformal", "y_true"]).agg(coverage=("covered", "mean"), avg_set_size=("size", "mean"),
                                                                                       n=("covered", "size"), seen=("seen", "all")).reset_index())
        # recalibration on N target animals, test on the remaining target animals (same animals for both quantiles)
        rng = np.random.default_rng(seed + 7)
        t_pids = np.unique(g[te])
        for n_recal in recal_ns:
            if n_recal >= len(t_pids):
                continue
            cov_t, cov_s, sz_t, sz_s, dropped = [], [], [], [], []
            for _ in range(recal_repeats):
                chosen = set(rng.choice(t_pids, size=n_recal, replace=False))
                m_cal = np.array([p in chosen for p in g[te]])
                tcal, ttest = te[m_cal], te[~m_cal]
                assert_no_group_leak(meta, tcal, ttest)
                t_scores, n_drop = cp.calibration_scores(est.predict_proba(X[tcal]), y[tcal], cls, method)
                q_t = cp.conformal_quantile(t_scores, alpha)
                c_t, s_t, _ = coverage_with(est, cls, X[ttest], y[ttest], q_t, method)
                c_s, s_s, _ = coverage_with(est, cls, X[ttest], y[ttest], qhat, method)
                cov_t.append(c_t); cov_s.append(c_s); sz_t.append(s_t); sz_s.append(s_s); dropped.append(n_drop)
            row[f"cov_target_recal_N{n_recal}"] = float(np.nanmean(cov_t))
            row[f"cov_target_srccal_N{n_recal}"] = float(np.nanmean(cov_s))
            recal_rows.append({"split": name, "arm": arm, "n_recal_animals": n_recal, "repeats": recal_repeats,
                               "coverage_target_recalibrated": float(np.nanmean(cov_t)), "coverage_target_recal_sd": float(np.nanstd(cov_t)),
                               "coverage_target_source_cal_same_test": float(np.nanmean(cov_s)),
                               "set_size_recalibrated": float(np.mean(sz_t)), "set_size_source_cal": float(np.mean(sz_s)),
                               "n_recal_dropped_unseen_mean": float(np.mean(dropped))})
        rows.append(row)
        per_class.append(pd.DataFrame({"split": name, "arm": arm, "y_true": y[te], "covered": covered, "size": size})
                         .groupby(["split", "arm", "y_true"]).agg(coverage=("covered", "mean"), avg_set_size=("size", "mean"),
                                                                   n=("covered", "size")).reset_index())
    return rows, pd.concat(per_class), recal_rows, pd.concat(variant_rows)


def main() -> None:
    ap = cli.common_parser("Shift tests")
    ap.add_argument("--label", default="tissue")
    ap.add_argument("--model", default="logreg_l2")
    ap.add_argument("--k", type=int, default=20, help="panel size to test (use the certified k from phase 06)")
    ap.add_argument("--alpha", type=float, default=0.1)
    ap.add_argument("--holdout-groups", default="8w,1w")
    ap.add_argument("--one-per-animal", action="store_true", help="calibrate on one vial per source animal")
    ap.add_argument("--recalibrate-target", default="", help="comma-separated numbers of target animals to hold out for recalibration, e.g. 3,5")
    ap.add_argument("--recal-repeats", type=int, default=10)
    ap.add_argument("--watch-tissues", default="SMLINT,BAT", help="tissues whose per-class coverage is quoted explicitly")
    args = ap.parse_args()
    cli.banner("08_shift_tests", args)
    source = cli.resolve_source(args.assay, args.source)
    out = cli.outdir(f"08_shift/{args.assay}", args.out)
    pheno = io.load_pheno()
    om = io.stack_tissues(args.assay, tissues=cli.parse_tissues(args.tissues), source=source, join="inner", pheno=pheno,
                          complete=args.complete_features, drop_incomplete_samples=args.drop_incomplete_samples)
    prefilter = min(args.prefilter, 1000) if args.quick else args.prefilter
    recal_ns = [int(v) for v in args.recalibrate_target.split(",") if v.strip()]
    X = om.X.to_numpy(dtype=float)
    y = om.meta[args.label].astype(str).to_numpy()
    g = om.groups()

    splits = list(leave_one_sex_out(om.meta))
    splits += list(leave_one_group_out(om.meta, "group", [h.strip() for h in args.holdout_groups.split(",")]))
    splits += list(train_controls_test_trained(om.meta))
    rows, pcs, recals, variants = [], [], [], []
    for name, (tr, te) in splits:
        r, pc, rc, vr = eval_split(name, tr, te, X, y, g, om.meta, args.k, args.alpha, args.model, prefilter, args.quick,
                                   args.seed, one_per_animal=args.one_per_animal, recal_ns=recal_ns, recal_repeats=args.recal_repeats)
        rows.extend(r); pcs.append(pc); recals.extend(rc); variants.append(vr)
        for rr in r:
            print(f"  {name:28s} {rr['arm']:10s} acc={rr['accuracy_all']:.3f} cov_src={rr['coverage_source_id']:.3f} "
                  f"cov_target={rr['coverage_target_all']:.3f} (seen: {rr['coverage_target_seen']:.3f}) "
                  f"size={rr['avg_set_size_target']:.2f} unseen={rr['unseen_classes'] or '-'}")
    table = pd.DataFrame(rows)
    table["coverage_drop"] = table["coverage_source_id"] - table["coverage_target_seen"]
    table["set_size_change"] = table["avg_set_size_target"] - table["avg_set_size_source"]
    table.to_csv(out / "shift_table.csv", index=False)
    pc = pd.concat(pcs, ignore_index=True)
    pc.to_csv(out / "shift_per_class.csv", index=False)
    vpc = pd.concat(variants, ignore_index=True)
    vpc.to_csv(out / "shift_per_class_variants.csv", index=False)
    sex_splits = [sp for sp in vpc["split"].unique() if sp.startswith("train_male") or sp.startswith("train_female")]
    vsex = vpc[vpc["split"].isin(sex_splits) & (vpc["arm"] == f"panel_k{args.k}") & vpc["seen"]]
    vsex_w = vsex.pivot_table(index="y_true", columns=["split", "conformal"], values="coverage").reset_index()
    vsex_w.columns = ["_".join(c).rstrip("_") if isinstance(c, tuple) else c for c in vsex_w.columns]
    recal = pd.DataFrame(recals)
    if len(recal):
        recal.to_csv(out / "shift_recalibration.csv", index=False)
    plot_df = table.assign(x=table["split"] + "\n" + table["arm"]).rename(columns={"coverage_target_seen": "coverage"})
    fig = plots.coverage_bars(plot_df, "x", f"{args.assay}: target coverage under shift (α={args.alpha}, seen classes)",
                              out / "coverage_shift.png", alpha=args.alpha)
    cols = ["split", "arm", "n_train_animals", "n_test_animals", "n_cal", "n_cal_dropped_unseen", "unseen_classes",
            "accuracy_all", "bal_acc_seen", "coverage_source_id", "coverage_target_seen", "coverage_drop",
            "coverage_target_seen_mondrian", "coverage_target_seen_floored",
            "avg_set_size_source", "avg_set_size_target", "set_size_change"] + [c for c in table.columns if c.startswith("cov_target_recal_N")]
    watch = [t.strip() for t in args.watch_tissues.split(",") if t.strip()]
    watch_md = ""
    ctt = pc[(pc["split"] == "train_control_test_trained") & (pc["y_true"].isin(watch))]
    if len(ctt):
        watch_md = ("Per-tissue coverage under controls→trained for the tissues whose panel markers are training-regulated:\n\n"
                    + report.df_to_md(ctt, floatfmt=".2f"))
    cal_txt = "one vial per source animal" if args.one_per_animal else "pooled vials of the source calibration animals"
    body = [
        f"{om.notes[-1]}; model `{args.model}`; panel k={args.k}; α={args.alpha}; calibration on 30% of source animals "
        f"({cal_txt}); 'coverage_source_id' is in-distribution coverage on held-out *source* animals for comparison; "
        f"`cov_target_recal_N` is coverage on the remaining target animals after calibrating on N target animals "
        f"(pooled vials, mean of {args.recal_repeats} draws), `cov_target_srccal_N` the source-calibrated coverage on the same animals.",
        report.df_to_md(table[cols], floatfmt=".3f"),
        (report.df_to_md(recal, floatfmt=".3f") if len(recal) else ""),
        f"Per-tissue coverage on the held-out sex (panel k={args.k}, seen classes only): marginal vs Mondrian vs marginal-floor "
        "Mondrian, all calibrated on the source sex (`coverage_target_seen_mondrian` / `_floored` in shift_table.csv give the pooled "
        "numbers for every split):\n\n" + report.df_to_md(vsex_w, floatfmt=".2f"),
        "Empty-set rate (LAC) and APS set size under the same shifts:\n\n"
        + report.df_to_md(table[["split", "arm", "lac_frac_empty_source", "lac_frac_empty_target", "aps_coverage_target_seen",
                                 "aps_set_size_source", "aps_set_size_target"]], floatfmt=".3f"),
        report.figure_md(fig, "Coverage on the shifted target vs the 1−α line"),
        watch_md,
        "Calibration samples whose class is absent from the fitted estimator are dropped and counted "
        "(`n_cal_dropped_unseen`), never scored against another class's column.",
        "Read: a coverage drop with a stable accuracy means the *guarantee* broke before the classifier did — "
        "that is the headline. Unseen classes (sex-specific tissues) can never be covered without a "
        "sex-conditional or open-set design; they are excluded from 'seen' columns and listed explicitly. "
        "Under the sex shift the source holds 25 animals, so an animal-level certificate at (α, δ) = (0.10, 0.10) "
        "(22 zero-error calibration animals) is impossible there; the sex-shift rows are pooled-vial calibrations.",
    ]
    subset = f", {om.meta['tissue'].nunique()} tissues" if args.tissues else ""
    report.add_section(f"08 · Shift tests ({args.assay}/{source}{subset})", "\n\n".join(b for b in body if b),
                       params={"assay": args.assay, "k": args.k, "alpha": args.alpha, "one_per_animal": args.one_per_animal,
                               "recalibrate_target": args.recalibrate_target, "tissues": args.tissues, "quick": args.quick})
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
