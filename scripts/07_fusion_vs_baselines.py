#!/usr/bin/env python
"""Phase 07 — does multi-omic fusion beat the best tuned single-omic baseline?

Task A (easy): tissue identity on the tissues that have every requested assay. PROT is excluded
by default: its plex-structured missingness is a tissue label in disguise (phase 04 diagnostic),
so any cross-tissue arm containing it would be measuring the plex layout.
Task B (hard): control vs 8w within each fusion tissue (n ≈ 20 animals). TRNSCRPT from counts
(per-tissue filter), PROT restricted to proteins with no missing value in that tissue, METAB to
complete features. Folds are animal-grouped and stratified on sex × group; sex is an explicit
covariate column in every arm. A permutation null (group shuffled within sex, --n-perm times)
records the maximum AUROC over all arms each time, so each tissue's best arm is judged against
the null of "best of 13 arms at n ≈ 20", not against 0.5. The null uses the untuned (quick)
versions of the same arms with 100-tree forests to keep it affordable; that is stated in the report.
--durations 1w,2w,4w,8w repeats the tuned Task B arms for every training duration vs control
(the permutation null only for --primary-duration): a training effect should grow with the
duration, a batch effect should not (scripts/07_batch_check.py reads the result).

Outputs (results/07_fusion/): taskA_per_fold.csv, taskA_summary.csv, taskB_per_fold.csv,
taskB_summary.csv, taskB_best_vs_null.csv, taskB_null.csv, taskB_batch_balance.csv,
taskB_by_duration.csv, taskB_duration_summary.csv, taskB_auroc_by_duration.png
"""
from __future__ import annotations

import time
import warnings

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

from motrpac import cli, config as C, io, models, report
from motrpac.splits import assert_no_group_leak, grouped_kfold

DURATION_WEEKS = {"1w": 1, "2w": 2, "4w": 4, "8w": 8}


def duration_figure(ds: pd.DataFrame, path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6, 4))
    for t, d in ds.groupby("tissue"):
        d = d.assign(w=d["duration"].map(DURATION_WEEKS)).sort_values("w")
        ax.plot(d["w"], d["best_single_auroc"], marker="o", label=t)
    ax.set_xscale("log", base=2)
    ax.set_xticks(sorted(DURATION_WEEKS.values()))
    ax.set_xticklabels([f"{w} wk" for w in sorted(DURATION_WEEKS.values())])
    ax.set_xlabel("training duration vs sedentary control")
    ax.set_ylabel("best single-omic AUROC (animal-grouped CV)")
    ax.set_ylim(0.4, 1.02)
    ax.axhline(0.5, color="grey", ls=":", lw=1)
    ax.set_title("Task B: control vs trained, by duration")
    ax.legend(frameon=False, fontsize=7, ncol=2)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def load_block(assay: str, tissue: str, pheno, source_trn: str, complete: bool):
    om = io.load_counts(tissue, pheno) if (assay == "TRNSCRPT" and source_trn == "counts") else io.load_norm(assay, tissue, pheno)
    n0, s0 = om.n_features, om.n_samples
    if complete and assay != "TRNSCRPT":
        # samples missing a whole platform / plex block first (else no feature is complete), then complete features
        keep_s = (om.X.isna().mean(axis=1) <= 0.2).to_numpy()
        om = om.subset(keep_s)
        keep = om.X.notna().all(axis=0)
        om = om.select_features(list(om.X.columns[keep.to_numpy()]))
        om.notes.append(f"{assay}/{tissue}: dropped {s0 - om.n_samples} samples missing > 20% of features; "
                        f"{om.n_features}/{n0} features with no missing value in the remaining samples")
    return om


def build_blocks(tissues, assays, pheno, source_trn, complete=True, complete_stage="tissue"):
    """Per-assay matrices indexed by bid, aligned across assays within tissue and stacked across tissues.
    complete_stage='tissue': drop incomplete samples / keep complete features inside each tissue (Task B,
    one tissue). complete_stage='stacked': join features across tissues first, then drop samples missing
    > 20% of the joined features and keep complete features (Task A) — a sample missing one platform
    would otherwise remove every feature of that platform from the cross-tissue intersection."""
    blocks = {a: [] for a in assays}
    metas, notes = [], []
    for t in tissues:
        mats = {a: load_block(a, t, pheno, source_trn, complete and complete_stage == "tissue") for a in assays}
        n_before = {a: m.n_samples for a, m in mats.items()}
        al = io.align_by_animal(mats, key="bid")
        notes.append(f"{t}: {n_before} → {al[assays[0]].n_samples} shared biospecimens; features "
                     + ", ".join(f"{a}={al[a].n_features}" for a in assays)
                     + "; " + "; ".join(n for m in mats.values() for n in m.notes if "dropped" in n))
        for a in assays:
            blocks[a].append(al[a].X)
        metas.append(al[assays[0]].meta)
    out = {}
    for a in assays:
        common = set(blocks[a][0].columns)
        for b in blocks[a][1:]:
            common &= set(b.columns)
        common = [c for c in blocks[a][0].columns if c in common]
        out[a] = pd.concat([b[common] for b in blocks[a]], axis=0)
    meta = pd.concat(metas, axis=0)
    if complete and complete_stage == "stacked":
        keep_s = pd.Series(True, index=meta.index)
        for a in assays:
            if a != "TRNSCRPT":
                keep_s &= (out[a].isna().mean(axis=1) <= 0.2)
        n_drop = int((~keep_s).sum())
        for a in assays:
            out[a] = out[a].loc[keep_s.to_numpy()]
            if a != "TRNSCRPT":
                n0 = out[a].shape[1]
                out[a] = out[a].loc[:, out[a].notna().all(axis=0).to_numpy()]
                notes.append(f"{a}: {out[a].shape[1]}/{n0} joined features with no missing value after dropping "
                             f"{n_drop} samples missing > 20% of the joined features")
        meta = meta.loc[keep_s.to_numpy()]
    return out, meta, notes


def with_sex(X: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    return X.assign(sex_male=(meta["sex"].astype(str) == "male").astype(float).to_numpy())


def run_task(blocks, meta, label, strat_label, n_splits, prefilter, quick, seed, kinds, add_sex=False,
             rf_trees=None, verbose=True):
    assays = list(blocks)
    y = meta[label].astype(str).to_numpy()
    g = meta["pid"].astype(str).to_numpy()
    classes = sorted(np.unique(y))
    Xb = {a: (with_sex(blocks[a], meta) if add_sex else blocks[a]) for a in assays}
    Xs = {a: Xb[a].to_numpy(dtype=float) for a in assays}
    Xf, sizes = models.early_fusion_matrix(Xb)
    Xf = Xf.to_numpy(dtype=float)
    rows = []

    def fit(kind, X, ytr, gtr, fusion=False):
        if fusion:
            pipe = models.make_fusion_pipeline(kind, sizes, seed=seed)
        else:
            pipe = models.make_pipeline(kind, k=None, prefilter=prefilter, seed=seed)
        if rf_trees and kind == "rf":
            pipe.set_params(clf__n_estimators=rf_trees)
        if quick or fusion:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                return pipe.fit(X, ytr)
        return models.fit_tuned(kind, X, ytr, gtr, prefilter=prefilter, quick=quick, seed=seed)

    for fold, (tr, te) in enumerate(grouped_kfold(meta, strat_label, n_splits, seed)):
        assert_no_group_leak(meta, tr, te)
        n_te_animals = int(meta.iloc[te]["pid"].nunique())
        per_block_est = {}
        for a in assays:
            for kind in kinds:
                est = fit(kind, Xs[a][tr], y[tr], g[tr])
                proba = est.predict_proba(Xs[a][te]) if models.has_proba(est) else None
                m = models.classification_metrics(y[te], est.predict(Xs[a][te]), proba, list(est.classes_) if proba is not None else classes)
                rows.append({"fold": fold, "arm": f"single:{a}", "model": kind, "n_test_animals": n_te_animals, **m})
                if kind == "logreg_l2":
                    per_block_est[a] = est
        for kind in kinds:
            pipe = fit(kind, Xf[tr], y[tr], g[tr], fusion=True)
            proba = pipe.predict_proba(Xf[te]) if models.has_proba(pipe) else None
            m = models.classification_metrics(y[te], pipe.predict(Xf[te]), proba, list(pipe.classes_) if proba is not None else classes)
            rows.append({"fold": fold, "arm": "early_fusion", "model": kind, "n_test_animals": n_te_animals, **m})
        if len(per_block_est) == len(assays):
            P, cls = models.late_fusion_proba([per_block_est[a] for a in assays], [Xs[a][te] for a in assays])
            y_pred = np.asarray(cls)[P.argmax(axis=1)]
            m = models.classification_metrics(y[te], y_pred, P, list(cls))
            rows.append({"fold": fold, "arm": "late_fusion", "model": "logreg_l2", "n_test_animals": n_te_animals, **m})
        if verbose:
            last = [r for r in rows if r["fold"] == fold]
            best = max(last, key=lambda r: r["balanced_accuracy"])
            print(f"  fold {fold}: best arm {best['arm']}/{best['model']} bal.acc={best['balanced_accuracy']:.3f}")
    return pd.DataFrame(rows)


def summarize(per_fold, metric, by=("arm", "model")):
    s = per_fold.groupby(list(by))[metric].agg(["mean", "std"]).reset_index()
    s.columns = list(by) + [f"{metric}_mean", f"{metric}_sd"]
    return s.sort_values(f"{metric}_mean", ascending=False)


def permuted_max_auroc(blocks, meta, n_splits, prefilter, seed, kinds, perm_seed):
    """One permutation: shuffle group within sex (animal level), rerun all arms quick, return the max mean AUROC."""
    rng = np.random.default_rng(perm_seed)
    m = meta.copy()
    grp = m["group"].astype(str).to_numpy().copy()
    for sx in np.unique(m["sex"]):
        idx = np.flatnonzero((m["sex"] == sx).to_numpy())
        grp[idx] = rng.permutation(grp[idx])
    m["group"] = grp
    m["sex_group"] = m["sex"].astype(str) + "/" + m["group"].astype(str)
    pf = run_task(blocks, m, "group", "sex_group", n_splits, prefilter, True, seed, kinds, add_sex=True, rf_trees=100, verbose=False)
    metric = "auroc" if "auroc" in pf.columns else "balanced_accuracy"
    return float(pf.groupby(["arm", "model"])[metric].mean().max())


def batch_balance(tissue: str, meta_sub: pd.DataFrame) -> dict:
    """Group balance across TMT plex and channel for the animals of a Task B subset (from meta/PROT.csv)."""
    p = C.META_DIR / "PROT.csv"
    out = {"tissue": tissue, "n_control": int((meta_sub["group"] == "control").sum()), "n_8w": int((meta_sub["group"] == "8w").sum())}
    if not p.exists():
        return out
    pm = pd.read_csv(p, dtype=str)
    pm = pm[pm["tissue"] == tissue].drop_duplicates("viallabel")
    pm["bid"] = pm["viallabel"].str[:5]
    j = pm.merge(meta_sub.reset_index()[["bid", "group"]], on="bid", how="inner")
    for col in ("tmt_plex", "tmt11_channel"):
        ct = pd.crosstab(j[col], j["group"])
        n = ct.to_numpy().sum()
        if n == 0 or ct.shape[1] < 2:
            out[f"V_{col}"] = np.nan
            continue
        exp = np.outer(ct.sum(1), ct.sum(0)) / n
        chi2 = ((ct.to_numpy() - exp) ** 2 / exp).sum()
        out[f"V_{col}"] = float(np.sqrt(chi2 / (n * (min(ct.shape) - 1))))
        out[f"n_{col}_levels"] = int(ct.shape[0])
        out[f"max_frac_one_group_in_a_{col}"] = float((ct.max(axis=1) / ct.sum(axis=1)).max())
    return out


def main() -> None:
    ap = cli.common_parser("Fusion vs single-omic baselines")
    ap.add_argument("--assays", default="TRNSCRPT,METAB", help="Task A blocks (PROT excluded on purpose)")
    ap.add_argument("--assays-b", default="TRNSCRPT,PROT,METAB", help="Task B blocks")
    ap.add_argument("--models", default="centroid,logreg_l2,rf")
    ap.add_argument("--task", default="both", choices=["A", "B", "both"])
    ap.add_argument("--n-perm", type=int, default=200)
    ap.add_argument("--perm-jobs", type=int, default=-1)
    ap.add_argument("--time-one-tissue", action="store_true", help="time Task B and one permutation on the first tissue, then stop")
    ap.add_argument("--durations", default="8w", help="trained groups compared with control in Task B, e.g. 1w,2w,4w,8w")
    ap.add_argument("--primary-duration", default="8w", help="the duration that gets the permutation null and the taskB_* files")
    args = ap.parse_args()
    cli.banner("07_fusion_vs_baselines", args)
    out = cli.outdir("07_fusion", args.out)
    pheno = io.load_pheno()
    kinds = [k.strip() for k in args.models.split(",")]
    source_trn = cli.resolve_source("TRNSCRPT", args.source)
    n_splits = 3 if args.quick else args.n_splits
    prefilter = min(args.prefilter, 1000) if args.quick else args.prefilter
    n_perm = 20 if args.quick else args.n_perm
    files = io.list_sample_files("norm")
    body = []

    def tissues_for(assays):
        return cli.parse_tissues(args.tissues) or sorted(set.intersection(*[set(files[files["assay"] == a]["tissue"]) for a in assays]))

    # ---- Task A ---------------------------------------------------------------------------
    if args.task in ("A", "both") and not args.time_one_tissue:
        assays = [a.strip() for a in args.assays.split(",")]
        tissues = tissues_for(assays)
        blocks, meta, notes = build_blocks(tissues, assays, pheno, source_trn, complete=True, complete_stage="stacked")
        pfA = run_task(blocks, meta, "tissue", "tissue", n_splits, prefilter, args.quick, args.seed, kinds)
        pfA.to_csv(out / "taskA_per_fold.csv", index=False)
        sA = summarize(pfA, "balanced_accuracy")
        sA.to_csv(out / "taskA_summary.csv", index=False)
        body.append(f"**Task A — tissue identity on {', '.join(assays)} ({meta['tissue'].nunique()} tissues, {meta['pid'].nunique()} "
                    f"animals):** PROT is excluded because its plex-structured missingness is a tissue label (phase 04 diagnostic: "
                    "missingness alone classifies tissue with accuracy 1.00), so a cross-tissue arm containing it would score the "
                    "plex layout, not biology.\n\nAlignment: " + "; ".join(notes) + "\n\n" + report.df_to_md(sA, floatfmt=".3f"))
        print(sA.to_string())

    # ---- Task B ---------------------------------------------------------------------------
    if args.task in ("B", "both"):
        assays_b = [a.strip() for a in args.assays_b.split(",")]
        tissues = tissues_for(assays_b)
        durations = [d.strip() for d in args.durations.split(",") if d.strip()]
        primary = args.primary_duration
        if primary not in durations:
            durations.append(primary)
        rowsB, nulls, balance, per_t, rows_dur = [], [], [], {}, []
        t_iter = tissues[:1] if args.time_one_tissue else tissues
        for t in t_iter:
            blocks_t, meta_t, notes_t = build_blocks([t], assays_b, pheno, source_trn, complete=True)
            for dur in (durations if not args.time_one_tissue else [primary]):
                m = meta_t["group"].isin(["control", dur]).to_numpy()
                sub_blocks = {a: b.loc[m] for a, b in blocks_t.items()}
                sub_meta = meta_t.loc[m].copy()
                sub_meta["sex_group"] = sub_meta["sex"].astype(str) + "/" + sub_meta["group"].astype(str)
                if sub_meta["pid"].nunique() < 12 or sub_meta["group"].nunique() < 2:
                    print(f"  task B {t} {dur}: skipped ({sub_meta['pid'].nunique()} animals)")
                    continue
                n_splits_b = min(4, sub_meta["pid"].nunique() // 3)
                t0 = time.time()
                pf = run_task(sub_blocks, sub_meta, "group", "sex_group", n_splits_b, min(prefilter, 2000), args.quick, args.seed,
                              kinds, add_sex=True, verbose=False)
                t_task = time.time() - t0
                pf["tissue"], pf["duration"] = t, dur
                pf["n_animals"] = int(sub_meta["pid"].nunique())
                rows_dur.append(pf)
                print(f"  task B {t} control vs {dur}: {sub_meta['pid'].nunique()} animals, {n_splits_b} folds, features "
                      + ", ".join(f"{a}={b.shape[1]}" for a, b in sub_blocks.items()) + f"; tuned arms {t_task:.0f}s")
                if dur != primary:
                    continue
                rowsB.append(pf)
                balance.append(batch_balance(t, sub_meta))
            t0 = time.time()
            if n_perm > 0 and rowsB and rowsB[-1]["tissue"].iloc[0] == t:
                null = Parallel(n_jobs=1 if args.time_one_tissue else args.perm_jobs)(
                    delayed(permuted_max_auroc)(sub_blocks, sub_meta, n_splits_b, min(prefilter, 2000), args.seed, kinds, args.seed + 100 * i)
                    for i in range(1 if args.time_one_tissue else n_perm))
                t_perm = time.time() - t0
                print(f"    permutation null: {len(null)} draw(s) in {t_perm:.0f}s" + (f" → projected {t_perm * n_perm / 60:.1f} min per tissue sequential" if args.time_one_tissue else ""))
                nulls.append(pd.DataFrame({"tissue": t, "perm": range(len(null)), "max_auroc": null}))
            if args.time_one_tissue:
                print(f"projected Task B for {len(tissues)} tissues: tuned {t_task * len(tissues) / 60:.1f} min + null "
                      f"{t_perm * n_perm * len(tissues) / 60 / max(1, 12 if args.perm_jobs == -1 else args.perm_jobs):.1f} min with parallel draws")
                return
        dur_md = ""
        if rows_dur:
            pfD = pd.concat(rows_dur, ignore_index=True)
            pfD.to_csv(out / "taskB_by_duration.csv", index=False)
            metric_d = "auroc" if "auroc" in pfD.columns else "balanced_accuracy"
            sD = pfD.groupby(["tissue", "duration", "arm", "model"]).agg(auroc_mean=(metric_d, "mean"), auroc_sd=(metric_d, "std"),
                                                                        n_animals=("n_animals", "first")).reset_index()
            drows = []
            for (t, dur), st in sD.groupby(["tissue", "duration"]):
                single = st[st["arm"].str.startswith("single:")].sort_values("auroc_mean", ascending=False)
                best = st.sort_values("auroc_mean", ascending=False).iloc[0]
                row = {"tissue": t, "duration": dur, "weeks": DURATION_WEEKS.get(dur, np.nan), "n_animals": int(st["n_animals"].iloc[0]),
                       "best_arm": f"{best['arm']}/{best['model']}", "best_auroc": best["auroc_mean"],
                       "best_single_arm": f"{single.iloc[0]['arm']}/{single.iloc[0]['model']}", "best_single_auroc": single.iloc[0]["auroc_mean"],
                       "best_single_sd": single.iloc[0]["auroc_sd"]}
                for a in assays_b:
                    sa = single[single["arm"] == f"single:{a}"]
                    row[f"best_{a}_auroc"] = float(sa["auroc_mean"].max()) if len(sa) else np.nan
                drows.append(row)
            dsum = pd.DataFrame(drows).sort_values(["tissue", "weeks"])
            dsum.to_csv(out / "taskB_duration_summary.csv", index=False)
            if dsum["duration"].nunique() > 1:
                duration_figure(dsum, out / "taskB_auroc_by_duration.png")
                wide = dsum.pivot(index="tissue", columns="duration", values="best_single_auroc").reset_index()
                wide = wide[["tissue"] + [d for d in sorted(DURATION_WEEKS, key=DURATION_WEEKS.get) if d in wide.columns]]
                dur_md = ("**Task B by training duration** (control vs 1w / 2w / 4w / 8w within tissue; best single-omic arm's AUROC, "
                          "same tuned protocol and fold construction as the primary comparison; a training effect should grow with "
                          "the duration, a collection or processing batch should not):\n\n" + report.df_to_md(wide, floatfmt=".3f")
                          + "\n\n" + report.figure_md(out / "taskB_auroc_by_duration.png", "Best single-omic AUROC vs training duration, per tissue")
                          + "\n\nPer-assay best AUROC by duration:\n\n"
                          + report.df_to_md(dsum[["tissue", "duration", "n_animals"] + [f"best_{a}_auroc" for a in assays_b]], floatfmt=".3f"))
        if rowsB:
            pfB = pd.concat(rowsB, ignore_index=True)
            pfB.to_csv(out / "taskB_per_fold.csv", index=False)
            metric = "auroc" if "auroc" in pfB.columns else "balanced_accuracy"
            sB = pfB.groupby(["tissue", "arm", "model"])[metric].agg(["mean", "std"]).reset_index()
            sB.columns = ["tissue", "arm", "model", f"{metric}_mean", f"{metric}_sd"]
            sB.to_csv(out / "taskB_summary.csv", index=False)
            nl = pd.concat(nulls, ignore_index=True) if nulls else pd.DataFrame(columns=["tissue", "perm", "max_auroc"])
            nl.to_csv(out / "taskB_null.csv", index=False)
            bal = pd.DataFrame(balance)
            bal.to_csv(out / "taskB_batch_balance.csv", index=False)
            rows = []
            for t, st in sB.groupby("tissue"):
                best = st.sort_values(f"{metric}_mean", ascending=False).iloc[0]
                single = st[st["arm"].str.startswith("single:")].sort_values(f"{metric}_mean", ascending=False).iloc[0]
                fus = st[~st["arm"].str.startswith("single:")].sort_values(f"{metric}_mean", ascending=False).iloc[0]
                q = nl.loc[nl["tissue"] == t, "max_auroc"]
                rows.append({"tissue": t, "best_arm": f"{best['arm']}/{best['model']}", "best_auroc": best[f"{metric}_mean"],
                             "best_sd": best[f"{metric}_sd"], "best_single": f"{single['arm']}/{single['model']}",
                             "single_auroc": single[f"{metric}_mean"], "single_sd": single[f"{metric}_sd"],
                             "best_fusion": f"{fus['arm']}/{fus['model']}", "fusion_auroc": fus[f"{metric}_mean"],
                             "fusion_minus_single": fus[f"{metric}_mean"] - single[f"{metric}_mean"],
                             "fusion_beats_single_by_gt_sd": bool(fus[f"{metric}_mean"] - single[f"{metric}_mean"] > single[f"{metric}_sd"]),
                             "null_median_max_auroc": float(q.median()) if len(q) else np.nan,
                             "null_p95_max_auroc": float(q.quantile(0.95)) if len(q) else np.nan,
                             "p_perm": float((np.sum(q >= best[f"{metric}_mean"]) + 1) / (len(q) + 1)) if len(q) else np.nan,
                             "best_beats_null_p95": bool(len(q) and best[f"{metric}_mean"] > q.quantile(0.95))})
            bv = pd.DataFrame(rows)
            bv.to_csv(out / "taskB_best_vs_null.csv", index=False)
            n_fus = int(bv["fusion_beats_single_by_gt_sd"].sum())
            n_null = int(bv["best_beats_null_p95"].sum())
            body.append(f"**Task B — control vs {primary} within tissue ({metric}; {', '.join(assays_b)}; folds animal-grouped and stratified on "
                        "sex × group; `sex_male` is a covariate column in every arm; PROT and METAB restricted to features with no "
                        "missing value in the tissue):**\n\n" + report.df_to_md(bv, floatfmt=".3f"))
            body.append("Full Task B table:\n\n" + report.df_to_md(sB, floatfmt=".3f"))
            body.append("Group balance across TMT plex and channel for the Task B animals (Cramér's V of group vs plex/channel; "
                        "V near 1 means group is confounded with the batch):\n\n" + report.df_to_md(bal, floatfmt=".2f"))
            body.append(f"Permutation null: group shuffled within sex {n_perm} times per tissue, all arms rerun (untuned quick "
                        "versions, 100-tree forests), the maximum mean AUROC over arms recorded each time. `p_perm` is the fraction of "
                        "draws whose best arm matched or beat the observed best arm. The null is what 'best of 13 arms at n ≈ 20' "
                        "produces from noise alone.")
            body.append(f"**Conclusion.** Fusion beats the best single-omic baseline by more than that baseline's fold sd in "
                        f"{n_fus} of {len(bv)} tissues; the best arm beats the permutation null's 95th percentile in {n_null} of {len(bv)} "
                        "tissues. " + ("Nothing here supports fusion as a reliable improvement, and most 'signal' is what the best of "
                                      "13 arms produces from shuffled labels." if n_fus == 0 and n_null <= 1 else
                                      "Read the tissues that pass both tests as the only candidates for a real fusion effect, and "
                                      "check their plex balance before believing them."))
            print(bv.round(3).to_string(index=False))
        if dur_md:
            body.append(dur_md)
            print(dsum.round(3).to_string(index=False))
    report.add_section("07 · Fusion vs single-omic baselines", "\n\n".join(body),
                       params={"assays_A": args.assays, "assays_B": args.assays_b, "models": args.models, "n_perm": n_perm,
                               "durations": args.durations, "quick": args.quick})
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
