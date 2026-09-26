"""Cross-dataset transfer of the tissue fingerprint (rat BodyMap, human GTEx): shared code path.

Panels are fit on the source (MoTrPAC) raw log matrix: median imputation → variance prefilter on the
raw values → per-gene z-score within the source → round-robin panel order (optionally with a
class × gene exclusion mask for transferability) → one classifier per panel size on one of three
representations: per-gene z-scores (within each dataset), within-sample ranks of the panel genes,
or top-scoring pairs (which of two panel genes is higher within the sample). Targets are scored
through the same imputation and prefilter, z-scored within the target. Target tissues map to source
classes, some as super-classes; unmapped ones are out-of-distribution. Conformal sets calibrated on
held-out source animals are compared with sets recalibrated on a few target individuals.
"""
from __future__ import annotations

import itertools
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

from . import config as C, conformal as cp, models

REPRESENTATIONS = ("zscore", "rank", "pairs")


def zscore(df: pd.DataFrame) -> np.ndarray:
    return StandardScaler().fit_transform(SimpleImputer(strategy="median").fit_transform(np.asarray(df, dtype=float)))


def pair_features(X: np.ndarray) -> np.ndarray:
    """Top-scoring-pair indicators: for every pair i < j of the k panel genes, 1 if gene i is higher
    than gene j within the sample (k(k−1)/2 binary features)."""
    i, j = np.triu_indices(X.shape[1], k=1)
    return (X[:, i] > X[:, j]).astype(float)


def represent(kind: str, Lp_rows: np.ndarray, Zp_rows: np.ndarray, idx: np.ndarray) -> np.ndarray:
    if kind == "zscore":
        return Zp_rows[:, idx]
    if kind == "rank":
        return rankdata(Lp_rows[:, idx], axis=1) / len(idx)
    if kind == "pairs":
        return pair_features(Lp_rows[:, idx])
    raise ValueError(f"unknown representation {kind}")


class PanelModels:
    """Fit on the source raw log matrix L (samples × genes). See the module docstring.
    `exclude` is a (classes × genes) boolean mask aligned to `genes` and np.unique(y)."""

    def __init__(self, L, y, genes, grid, prefilter, representation="zscore", exclude=None, C_=0.1, full=True):
        L = np.asarray(L, dtype=float)
        self.genes = np.asarray(genes)
        self.representation = representation
        self.imp = SimpleImputer(strategy="median").fit(L)
        Li = self.imp.transform(L)
        self.pre = models.VarianceTopK(prefilter).fit(Li)
        Lp = self.pre.transform(Li)
        self.genes_pre = self.genes[self.pre.idx_]
        self.scaler = StandardScaler().fit(Lp)
        Zp = self.scaler.transform(Lp)
        ex = None if exclude is None else np.asarray(exclude, dtype=bool)[:, self.pre.idx_]
        sel = models.RoundRobinSelector(k=max(grid), exclude=ex).fit(Zp, y)
        self.order, self.scores, self.classes = sel.order_, sel.scores_, list(sel.classes_)
        self.fits, self.rep_scalers = {}, {}
        for k in grid:
            idx = np.sort(self.order[:k])
            F = represent(representation, Lp, Zp, idx)
            sc = StandardScaler().fit(F) if representation == "rank" else None
            self.rep_scalers[f"k{k}"] = sc
            self.fits[f"k{k}"] = (idx, LogisticRegression(solver="lbfgs", C=C_, max_iter=3000).fit(sc.transform(F) if sc else F, y))
        if full:
            self.fits["full"] = (np.arange(Zp.shape[1]), LogisticRegression(solver="lbfgs", C=C_, max_iter=3000).fit(Zp, y))

    def panel_genes(self, k):
        return list(self.genes_pre[np.sort(self.order[:k])])

    def source_matrices(self, L):
        """Prefiltered raw values and z-scores (source scaler) for source rows, e.g. calibration animals."""
        Lp = self.pre.transform(self.imp.transform(np.asarray(L, dtype=float)))
        return Lp, self.scaler.transform(Lp)

    def target_matrices(self, L_target):
        """Prefiltered raw values and z-scores computed within the target (all its samples)."""
        Lp = self.pre.transform(self.imp.transform(np.asarray(L_target, dtype=float)))
        return Lp, StandardScaler().fit_transform(Lp)

    def proba(self, name, Lp_rows, Zp_rows):
        idx, clf = self.fits[name]
        if name == "full":
            F = Zp_rows
        else:
            F = represent(self.representation, Lp_rows, Zp_rows, idx)
            sc = self.rep_scalers.get(name)
            F = sc.transform(F) if sc else F
        return clf.predict_proba(F), list(clf.classes_)


def superclass_correct(pred, organ, organ_map):
    ts = organ_map.get(organ)
    return bool(ts) and pred in ts


def score_block(pm, name, Lp_rows, Zp_rows, mt_rows, organ_map):
    p, cls = pm.proba(name, Lp_rows, Zp_rows)
    pred = np.asarray(cls)[p.argmax(axis=1)]
    ok = np.array([superclass_correct(pr, o, organ_map) for pr, o in zip(pred, mt_rows["organ"])])
    return pred, ok, p, cls


CONFORMAL_VARIANTS = ("marginal", "mondrian", "floored")


def calibrate_models(pm_c, Lp_cal, Zp_cal, y_cal, classes, alpha, model_names):
    """Per model: marginal LAC quantile (`q`), Mondrian per-class quantiles (`q_class`), the
    marginal-floor Mondrian quantiles (`q_floor` = max(q_class, q)) from held-out source animals,
    and a `proba(Lp, Zp)` closure."""
    out = {}
    y_idx = np.array([classes.index(v) for v in y_cal])
    for name in model_names:
        p_cal, cls = pm_c.proba(name, Lp_cal, Zp_cal)
        assert cls == classes
        scores = cp.lac_scores(p_cal, y_idx)
        q = cp.conformal_quantile(scores, alpha)
        q_c = cp.conformal_quantile_per_class(scores, y_idx, alpha, len(classes), fallback=q)
        q_f = cp.conformal_quantile_per_class(scores, y_idx, alpha, len(classes), fallback=q, floor=q)
        out[name] = {"q": q, "q_class": q_c, "q_floor": q_f, "n_cal": int(len(scores)),
                     "scores": scores, "y_idx": y_idx,
                     "proba": (lambda Lr, Zr, n=name: pm_c.proba(n, Lr, Zr)[0])}
    return out


def sets_for(cal: dict, p: np.ndarray, variant: str) -> np.ndarray:
    """Prediction sets (LAC) for one calibrated model under one conformal variant."""
    if variant == "marginal":
        return cp.predict_sets(p, cal["q"], "lac")
    if variant == "mondrian":
        return cp.predict_sets_conditional(p, cal["q_class"], "lac")
    if variant == "floored":
        return cp.predict_sets_conditional(p, cal["q_floor"], "lac")
    raise ValueError(f"unknown conformal variant {variant}")


def per_organ_coverage(calib, classes, Lp_t, Zp_t, mt, organ_map, model_names, variants=CONFORMAL_VARIANTS):
    """Per target organ: coverage (super-class), empty-set rate and set size for every model ×
    conformal variant on the rows given; unmapped organs get NaN coverage."""
    rows = []
    for name in model_names:
        cal = calib[name]
        p = cal["proba"](Lp_t, Zp_t)
        for variant in variants:
            sets = sets_for(cal, p, variant)
            size = sets.sum(axis=1)
            for organ, idx in mt.reset_index(drop=True).groupby("organ").indices.items():
                idx = np.asarray(idx)
                ts = [t for t in (organ_map.get(organ) or ()) if t in classes]
                cov = float(np.mean([any(sets[i, classes.index(t)] for t in ts) for i in idx])) if ts else np.nan
                rows.append({"model": name, "conformal": variant, "organ": organ, "n": int(len(idx)), "coverage": cov,
                             "frac_empty": float((size[idx] == 0).mean()), "avg_set_size": float(size[idx].mean())})
    return pd.DataFrame(rows)


def conformal_transfer(calib, classes, Lp_t, Zp_t, mt, organ_map, alpha, stage_col, primary_stage, recal_ns, recal_repeats, rng,
                       model_names, ood_organs, collect: list | None = None):
    """Source-calibrated sets (marginal, Mondrian, floored Mondrian) on the target per stage, what the sets hold for the
    out-of-distribution organs at the primary stage, and recalibration on N target individuals
    (pooled tissues; the calibration label is the organ's best-scoring mapped class), tested on the
    remaining individuals of the primary stage. `mt` needs columns organ, group_id and `stage_col`.
    `collect` (a list, optional) receives one dict per recalibration draw: model, n_recal, draw, the chosen
    individuals, the recalibrated threshold q_t and the number of calibration scores behind it."""
    conf_rows, ood_rows, recal_rows = [], [], []
    primary = (mt[stage_col] == primary_stage).to_numpy()
    for name in model_names:
        cal = calib[name]
        for stage in sorted(mt[stage_col].unique()):
            m_s = (mt[stage_col] == stage).to_numpy()
            p_b = cal["proba"](Lp_t[m_s], Zp_t[m_s])
            for cname in CONFORMAL_VARIANTS:
                sets = sets_for(cal, p_b, cname)
                sub = mt[m_s]
                in_dist = np.array([bool(organ_map.get(o)) for o in sub["organ"]])
                covered = np.array([bool(organ_map.get(o)) and any(sets[i, classes.index(t)] for t in organ_map[o] if t in classes)
                                    for i, o in enumerate(sub["organ"])])
                size = sets.sum(axis=1)
                conf_rows.append({stage_col: stage, "model": name, "conformal": cname, "n_mapped": int(in_dist.sum()),
                                  "coverage_mapped": float(covered[in_dist].mean()) if in_dist.any() else np.nan,
                                  "frac_empty_mapped": float((size[in_dist] == 0).mean()) if in_dist.any() else np.nan,
                                  "avg_set_size_mapped": float(size[in_dist].mean()) if in_dist.any() else np.nan,
                                  "ood_frac_empty": float((size[~in_dist] == 0).mean()) if (~in_dist).any() else np.nan,
                                  "ood_avg_set_size": float(size[~in_dist].mean()) if (~in_dist).any() else np.nan})
                if stage == primary_stage:
                    for organ in ood_organs:
                        mo = (sub["organ"] == organ).to_numpy()
                        if not mo.any():
                            continue
                        counts = pd.Series(sets[mo].sum(axis=0), index=classes).sort_values(ascending=False)
                        ood_rows.append({"model": name, "conformal": cname, "organ": organ, "n": int(mo.sum()),
                                         "frac_empty": float((size[mo] == 0).mean()), "avg_set_size": float(size[mo].mean()),
                                         "max_set_size": int(size[mo].max()),
                                         "most_frequent_members": "; ".join(f"{t} ({int(c)})" for t, c in counts.head(4).items() if c > 0)})
        ids = sorted(mt.loc[primary, "group_id"].unique())
        p_ad = cal["proba"](Lp_t[primary], Zp_t[primary])
        sub_ad = mt[primary].reset_index(drop=True)
        mapped_ad = np.array([bool(organ_map.get(o)) for o in sub_ad["organ"]])
        for n_recal in recal_ns:
            if n_recal >= len(ids):
                continue
            combos = list(itertools.combinations(ids, n_recal)) if len(ids) <= 12 else None
            if combos is None or len(combos) > recal_repeats:
                combos = [tuple(rng.choice(ids, size=n_recal, replace=False)) for _ in range(recal_repeats)]
            cov_r, cov_s, emp_r, emp_s, sz_r = [], [], [], [], []
            for draw, chosen in enumerate(combos):
                m_c = sub_ad["group_id"].isin(chosen).to_numpy() & mapped_ad
                m_t = ~sub_ad["group_id"].isin(chosen).to_numpy() & mapped_ad
                if not m_c.any() or not m_t.any():
                    continue
                y_c = [max([t for t in organ_map[sub_ad.loc[i, "organ"]] if t in classes], key=lambda t: p_ad[i, classes.index(t)])
                       for i in np.flatnonzero(m_c)]
                sc, _ = cp.calibration_scores(p_ad[m_c], np.array(y_c), classes, "lac")
                q_t = cp.conformal_quantile(sc, alpha)
                if collect is not None:
                    collect.append({"model": name, "n_recal": n_recal, "draw": draw, "chosen": ";".join(str(c) for c in chosen),
                                    "q_t": q_t, "n_cal_scores": int(len(sc))})
                for q_use, cov_list, emp_list in ((q_t, cov_r, emp_r), (cal["q"], cov_s, emp_s)):
                    sets = cp.predict_sets(p_ad[m_t], q_use, "lac")
                    covered = [any(sets[j, classes.index(t)] for t in organ_map[sub_ad.loc[i, "organ"]] if t in classes)
                               for j, i in enumerate(np.flatnonzero(m_t))]
                    cov_list.append(float(np.mean(covered)))
                    emp_list.append(float((sets.sum(axis=1) == 0).mean()))
                    if q_use is q_t:
                        sz_r.append(float(sets.sum(axis=1).mean()))
            recal_rows.append({"model": name, "n_recal": n_recal, "n_test_individuals": len(ids) - n_recal, "draws": len(cov_r),
                               "coverage_recalibrated": float(np.mean(cov_r)), "coverage_source_cal_same_test": float(np.mean(cov_s)),
                               "frac_empty_recalibrated": float(np.mean(emp_r)), "frac_empty_source_cal": float(np.mean(emp_s)),
                               "set_size_recalibrated": float(np.mean(sz_r))})
    return pd.DataFrame(conf_rows), pd.DataFrame(ood_rows), pd.DataFrame(recal_rows)


def save_transfer_scores(out_dir, calib, classes, mt, keep_cols, Lp_c, Zp_c, cal_meta, collect, organ_map, model_names,
                         pm_all=None, Lp_all=None, Zp_all=None):
    """--save-scores for the transfer phases: per target sample, the probabilities of the calibrated models
    (`calib`, fit without the calibration animals) under every model name; per model, the LAC calibration
    scores of the held-out source animals (`cal_meta` rows, same order as `calib[name]["scores"]`); the
    per-draw recalibration thresholds collected by `conformal_transfer`; the class list and the organ map.
    If the all-animal models `pm_all` are given (with their target matrices), each sample also carries
    `pred_all_animals`, the call of the model the accuracy tables use. Everything a viewer needs to rebuild
    the prediction sets outside the pipeline."""
    import json
    out_dir = Path(out_dir)
    base = mt[keep_cols].copy()
    base.insert(0, "sample", mt.index.astype(str))
    rows = []
    for name in model_names:
        p = calib[name]["proba"](Lp_c, Zp_c)
        df = base.copy()
        df["model"] = name
        for j, c in enumerate(classes):
            df[f"p_{c}"] = p[:, j]
        if pm_all is not None:
            pa, cls = pm_all.proba(name, Lp_all, Zp_all)
            df["pred_all_animals"] = np.asarray(cls)[pa.argmax(axis=1)]
        rows.append(df)
    pd.concat(rows, ignore_index=True).to_csv(out_dir / "scores_target_probs.csv", index=False)
    cal_rows = []
    for name in model_names:
        df = cal_meta[["pid", "tissue"]].copy()
        df.insert(0, "viallabel", cal_meta.index.astype(str))
        df["model"] = name
        df["score_lac"] = calib[name]["scores"]
        cal_rows.append(df)
    pd.concat(cal_rows, ignore_index=True).to_csv(out_dir / "scores_calibration.csv", index=False)
    pd.DataFrame(collect or [], columns=["model", "n_recal", "draw", "chosen", "q_t", "n_cal_scores"]).to_csv(
        out_dir / "recal_thresholds.csv", index=False)
    (out_dir / "classes.json").write_text(json.dumps(list(classes)))
    (out_dir / "organ_map.json").write_text(json.dumps({o: (sorted(ts) if ts else None) for o, ts in organ_map.items()}))


def gene_check(pm, genes, sym, Zm_df, y, Zt_df, mt_primary, tissue_to_organ, flags):
    """Per panel gene: marker tissue on the source, effect there vs the target organ (mean z in the
    marker organ minus the highest other-organ mean); 'fails' = effect <= 0 on the target."""
    rows = []
    for gg in genes:
        j = list(pm.genes_pre).index(gg)
        marker = pm.classes[int(np.argmax(pm.scores[:, j]))]
        organ = tissue_to_organ.get(marker)
        mu_m = Zm_df[gg].groupby(y).mean()
        eff_m = float(mu_m[marker] - mu_m.drop(marker).max())
        row = {"feature_ID": gg, "gene_symbol": sym.get(gg, gg), "marker_tissue": marker, "target_organ": organ or "not in target",
               "source_effect_z": eff_m}
        if organ and organ in set(mt_primary["organ"]):
            mu_b = Zt_df[gg].groupby(mt_primary["organ"]).mean()
            eff_b = float(mu_b[organ] - mu_b.drop(organ).max())
            row.update({"target_effect_z": eff_b, "target_top_organ": str(mu_b.idxmax()), "fails_in_target": bool(eff_b <= 0),
                        "weakened": bool(0 < eff_b < 0.5 * eff_m)})
        else:
            row.update({"target_effect_z": np.nan, "target_top_organ": None, "fails_in_target": None, "weakened": None})
        row.update(flags.get(gg, {"risk_T7_regulated": None, "risk_qc_correlated": None}))
        rows.append(row)
    return pd.DataFrame(rows)


def one_to_one_orthologs(raw_dir=None) -> pd.DataFrame:
    """Strict 1:1 rat–human pairs from data/raw/rat_to_human_gene.csv (versionless human Ensembl ids)."""
    r = pd.read_csv((raw_dir or C.RAW_DIR) / "rat_to_human_gene.csv", dtype=str)
    r = r.dropna(subset=["RAT_ENSEMBL_ID"])
    r["HUMAN_ORTHOLOG_ENSEMBL_ID"] = r["HUMAN_ORTHOLOG_ENSEMBL_ID"].str.split(".").str[0]
    r = r[r["HUMAN_ORTHOLOG_ENSEMBL_ID"].notna() | r["HUMAN_ORTHOLOG_SYMBOL"].notna()]
    r["human_key"] = r["HUMAN_ORTHOLOG_ENSEMBL_ID"].fillna("SYM:" + r["HUMAN_ORTHOLOG_SYMBOL"].astype(str))
    n_rat = r.groupby("RAT_ENSEMBL_ID")["human_key"].nunique()
    n_hum = r.groupby("human_key")["RAT_ENSEMBL_ID"].nunique()
    keep = r["RAT_ENSEMBL_ID"].map(n_rat).eq(1) & r["human_key"].map(n_hum).eq(1)
    return r[keep].drop_duplicates("RAT_ENSEMBL_ID")[["RAT_ENSEMBL_ID", "RAT_SYMBOL", "HUMAN_ORTHOLOG_ENSEMBL_ID", "HUMAN_ORTHOLOG_SYMBOL"]]


def qc_correlation(L, y, classes, qc) -> np.ndarray:
    """(classes × genes) Pearson r between each gene and a sample-level QC covariate within each class."""
    L = np.asarray(L, dtype=float)
    qc = np.asarray(qc, dtype=float)
    out = np.full((len(classes), L.shape[1]), np.nan)
    for i, t in enumerate(classes):
        m = (y == t) & ~np.isnan(qc)
        if m.sum() < 8:
            continue
        X = L[m] - L[m].mean(axis=0)
        q = qc[m] - qc[m].mean()
        denom = np.sqrt((X ** 2).sum(axis=0)) * np.sqrt((q ** 2).sum())
        with np.errstate(invalid="ignore", divide="ignore"):
            out[i] = (X.T @ q) / denom
    return out


def exclusion_mask(L, y, genes, classes, qc=None, regulated=None, r_thresh=0.5) -> np.ndarray:
    """(classes × genes) True where a gene may not be picked for that class: |r| > r_thresh with the
    QC covariate inside the class, or training-regulated in that tissue (regulated: tissue → set of ids)."""
    genes = np.asarray(genes)
    mask = np.zeros((len(classes), len(genes)), dtype=bool)
    if qc is not None:
        r = qc_correlation(L, y, classes, qc)
        mask |= np.abs(np.nan_to_num(r)) > r_thresh
    if regulated:
        for i, t in enumerate(classes):
            mask[i] |= np.isin(genes, list(regulated.get(t, ())))
    return mask


def match_gtex_orthologs(orth: pd.DataFrame, gtex_genes, gtex_symbols: pd.Series | None = None):
    """Match 1:1 rat–human pairs to the GTEx columns by human Ensembl id, then by symbol for the
    rest (when a GTEx id → symbol map is given). Returns (pairs, n_by_ensembl, n_by_symbol) with
    pairs strictly 1:1 on both sides."""
    hum = set(gtex_genes)
    by_ens = orth[orth["HUMAN_ORTHOLOG_ENSEMBL_ID"].isin(hum)]
    n_by_ensembl, used_symbol = len(by_ens), 0
    if gtex_symbols is not None:
        sym_to_ens = {sym: e for e, sym in gtex_symbols.items() if isinstance(sym, str)}
        rest = orth[~orth["HUMAN_ORTHOLOG_ENSEMBL_ID"].isin(hum) & orth["HUMAN_ORTHOLOG_SYMBOL"].isin(sym_to_ens)].copy()
        rest["HUMAN_ORTHOLOG_ENSEMBL_ID"] = rest["HUMAN_ORTHOLOG_SYMBOL"].map(sym_to_ens)
        used_symbol = len(rest)
        by_ens = pd.concat([by_ens, rest])
    by_ens = by_ens.drop_duplicates("RAT_ENSEMBL_ID").drop_duplicates("HUMAN_ORTHOLOG_ENSEMBL_ID")
    return by_ens, n_by_ensembl, used_symbol


def gtex_symbols(external_dir=None) -> pd.Series | None:
    p = (external_dir or C.EXTERNAL_DIR) / "gtex_gene_symbols.csv"
    return pd.read_csv(p, index_col=0).iloc[:, 0] if p.exists() else None
