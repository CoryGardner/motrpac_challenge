# %% [markdown]
# ## 11. Multi-omic fusion, the permutation null, and the batch check
#
# **Question.** Does fusing omics (transcriptome + metabolome + proteome) beat the best *single* omic,
# (A) for tissue identity and (B) for control vs 8-week-trained animals within one tissue? And when
# Task B separates the groups almost perfectly, is that training — or a processing batch?
#
# **Why it matters.** "Fusion helps" is the default expectation of a multi-omic hackathon project. With a
# few dozen animals per tissue (counts printed below), 13 model arms and a best-of-13 pick, a
# near-perfect AUROC is cheap; the
# honest comparator is the *null of the best arm* under shuffled labels, not 0.5. And a perfect
# separation proves nothing about training if QC or collection covariates alone also separate the
# groups.
#
# **What to look for.** (1) Task A: is every arm at the ceiling, so fusion has no room to help?
# (2) Task B: does fusion beat the best single omic by more than one fold sd anywhere, and does the
# best arm beat the permutation null's 95th percentile? (3) Batch check: which tissues does the verdict
# rule call "unresolvable", and why? (4) The two different SKM-GN covariate-only AUROCs that circulate in the findings —
# which row of the table each one is, and which statement each supports.
#
# **Mode.** Fusion is the slowest phase of the pipeline (tuned arms plus a 200-draw permutation null per
# tissue), and the covariate-only classifiers carry 100-draw nulls. In the default mode these are
# **loaded** from `results/07_fusion/`; the cheap derived steps (best-vs-null table, verdict rule, and
# the two SKM-GN covariate classifiers without their null) are recomputed live and compared with the
# published files. `RECOMPUTE = True` reruns the copied script logic into `_outputs/s11/`.

# %%
section("11 fusion")
# ---- helpers copied from scripts/07_fusion_vs_baselines.py (used by the RECOMPUTE branch) ----
import warnings as _s11_warn
from joblib import Parallel, delayed
from scipy import stats
from sklearn.compose import ColumnTransformer
from sklearn.ensemble import RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler

S11_OUT = OUT / "s11"
S11_OUT.mkdir(parents=True, exist_ok=True)
S11_DURATION_WEEKS = {"1w": 1, "2w": 2, "4w": 4, "8w": 8}


# from scripts/07_fusion_vs_baselines.py::load_block (verbatim)
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


# from scripts/07_fusion_vs_baselines.py::build_blocks (verbatim)
def build_blocks(tissues, assays, pheno, source_trn, complete=True, complete_stage="tissue"):
    """Per-assay matrices indexed by bid, aligned across assays within tissue and stacked across tissues.
    complete_stage='tissue': drop incomplete samples / keep complete features inside each tissue (Task B,
    one tissue). complete_stage='stacked': join features across tissues first, then drop samples missing
    > 20% of the joined features and keep complete features (Task A)."""
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


# from scripts/07_fusion_vs_baselines.py::with_sex (verbatim)
def with_sex(X: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    return X.assign(sex_male=(meta["sex"].astype(str) == "male").astype(float).to_numpy())


# from scripts/07_fusion_vs_baselines.py::run_task (verbatim)
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
            with _s11_warn.catch_warnings():
                _s11_warn.simplefilter("ignore")
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


# from scripts/07_fusion_vs_baselines.py::summarize (verbatim)
def summarize(per_fold, metric, by=("arm", "model")):
    s = per_fold.groupby(list(by))[metric].agg(["mean", "std"]).reset_index()
    s.columns = list(by) + [f"{metric}_mean", f"{metric}_sd"]
    return s.sort_values(f"{metric}_mean", ascending=False)


# from scripts/07_fusion_vs_baselines.py::permuted_max_auroc (verbatim)
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


# from scripts/07_fusion_vs_baselines.py::batch_balance (verbatim)
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


# from src/motrpac/cli.py::resolve_source (verbatim; cli.py is not pasted into the notebook)
def resolve_source(assay: str, source: str) -> str:
    if source != "auto":
        return source
    if assay == "TRNSCRPT" and any((C.COUNTS_DIR).glob(f"{assay}__*.csv")):
        return "counts"
    return "norm"


# The `make fusion` settings: --models centroid,logreg_l2,rf --durations 1w,2w,4w,8w; script defaults otherwise
S11_ARGS = SimpleNamespace(assays="TRNSCRPT,METAB", assays_b="TRNSCRPT,PROT,METAB", kinds=["centroid", "logreg_l2", "rf"],
                           n_splits=5, prefilter=5000, n_perm=200, perm_jobs=-1, durations=["1w", "2w", "4w", "8w"],
                           primary="8w", seed=C.SEED)
print("fusion helpers defined; settings:", vars(S11_ARGS))

# %%
# ---- Task A: tissue identity on the tissues that have both TRNSCRPT and METAB (PROT excluded on purpose) ----
if RECOMPUTE:
    # from scripts/07_fusion_vs_baselines.py::main, Task A block (trimmed: report text dropped, output → OUT/s11)
    _files = io.list_sample_files("norm")
    _assays = S11_ARGS.assays.split(",")
    _tissues = sorted(set.intersection(*[set(_files[_files["assay"] == a]["tissue"]) for a in _assays]))
    _blocks, _meta, _notes = build_blocks(_tissues, _assays, cached("pheno", io.load_pheno), resolve_source("TRNSCRPT", "auto"),
                                          complete=True, complete_stage="stacked")
    s11_pfA = run_task(_blocks, _meta, "tissue", "tissue", S11_ARGS.n_splits, S11_ARGS.prefilter, False, S11_ARGS.seed, S11_ARGS.kinds)
    s11_pfA.to_csv(S11_OUT / "taskA_per_fold.csv", index=False)
    s11_taskA = summarize(s11_pfA, "balanced_accuracy")
    s11_taskA.to_csv(S11_OUT / "taskA_summary.csv", index=False)
else:
    s11_pfA = pd.read_csv(res("07_fusion", "taskA_per_fold.csv"))
    s11_taskA = pd.read_csv(res("07_fusion", "taskA_summary.csv"))

# n animals per test fold (evaluation rule: report it next to the mean)
print("Task A test animals per fold:", s11_pfA.groupby("fold")["n_test_animals"].first().tolist())
print(s11_taskA.to_string(index=False, float_format="%.4f"))
_single = s11_taskA[s11_taskA["arm"].str.startswith("single:")]["balanced_accuracy_mean"].max()
_fusion = s11_taskA[~s11_taskA["arm"].str.startswith("single:")]["balanced_accuracy_mean"].max()
print(f"\nbest single-omic arm {_single:.4f}  vs  best fusion arm {_fusion:.4f};  lowest arm {s11_taskA['balanced_accuracy_mean'].min():.4f}")
record("s11.taskA_best_single", _single, "11")
record("s11.taskA_best_fusion", _fusion, "11")
record("s11.taskA_min_arm", s11_taskA["balanced_accuracy_mean"].min(), "11")

# %% [markdown]
# **What Task A shows.** Every arm — single omic or fused, any model — sits at or within a hair of
# perfect balanced accuracy, and a single omic (METAB) is already at the top. Tissue identity is a
# ceiling task here: fusion has nothing left to add.
#
# **What it does not show.** That fusion is useless in general. A ceiling task cannot rank methods;
# it only says that tissue identity is not the place to demonstrate fusion.

# %%
# ---- Task B: control vs 8w within each fusion tissue; per-arm AUROC and the best-of-13-arms permutation null ----
if RECOMPUTE:
    # from scripts/07_fusion_vs_baselines.py::main, Task B block (trimmed: --time-one-tissue branch and report text dropped)
    _pheno = cached("pheno", io.load_pheno)
    _files = io.list_sample_files("norm")
    _assays_b = S11_ARGS.assays_b.split(",")
    _tissues = sorted(set.intersection(*[set(_files[_files["assay"] == a]["tissue"]) for a in _assays_b]))
    _rowsB, _nulls, _balance, _rows_dur = [], [], [], []
    for t in _tissues:
        blocks_t, meta_t, notes_t = build_blocks([t], _assays_b, _pheno, resolve_source("TRNSCRPT", "auto"), complete=True)
        for dur in S11_ARGS.durations:
            m = meta_t["group"].isin(["control", dur]).to_numpy()
            sub_blocks = {a: b.loc[m] for a, b in blocks_t.items()}
            sub_meta = meta_t.loc[m].copy()
            sub_meta["sex_group"] = sub_meta["sex"].astype(str) + "/" + sub_meta["group"].astype(str)
            if sub_meta["pid"].nunique() < 12 or sub_meta["group"].nunique() < 2:
                continue
            n_splits_b = min(4, sub_meta["pid"].nunique() // 3)
            pf = run_task(sub_blocks, sub_meta, "group", "sex_group", n_splits_b, min(S11_ARGS.prefilter, 2000), False,
                          S11_ARGS.seed, S11_ARGS.kinds, add_sex=True, verbose=False)
            pf["tissue"], pf["duration"], pf["n_animals"] = t, dur, int(sub_meta["pid"].nunique())
            _rows_dur.append(pf)
            if dur != S11_ARGS.primary:
                continue
            _rowsB.append(pf)
            _balance.append(batch_balance(t, sub_meta))
            null = Parallel(n_jobs=S11_ARGS.perm_jobs)(
                delayed(permuted_max_auroc)(sub_blocks, sub_meta, n_splits_b, min(S11_ARGS.prefilter, 2000), S11_ARGS.seed,
                                            S11_ARGS.kinds, S11_ARGS.seed + 100 * i) for i in range(S11_ARGS.n_perm))
            _nulls.append(pd.DataFrame({"tissue": t, "perm": range(len(null)), "max_auroc": null}))
            print(f"  {t}: tuned arms + {len(null)} permutations done")
    _pfD = pd.concat(_rows_dur, ignore_index=True)
    _pfD.to_csv(S11_OUT / "taskB_by_duration.csv", index=False)
    _sD = _pfD.groupby(["tissue", "duration", "arm", "model"]).agg(auroc_mean=("auroc", "mean"), auroc_sd=("auroc", "std"),
                                                                    n_animals=("n_animals", "first")).reset_index()
    _drows = []
    for (t, dur), st in _sD.groupby(["tissue", "duration"]):
        single = st[st["arm"].str.startswith("single:")].sort_values("auroc_mean", ascending=False)
        _drows.append({"tissue": t, "duration": dur, "weeks": S11_DURATION_WEEKS.get(dur, np.nan), "n_animals": int(st["n_animals"].iloc[0]),
                       "best_single_arm": f"{single.iloc[0]['arm']}/{single.iloc[0]['model']}",
                       "best_single_auroc": single.iloc[0]["auroc_mean"], "best_single_sd": single.iloc[0]["auroc_sd"]})
    s11_dur = pd.DataFrame(_drows).sort_values(["tissue", "weeks"])
    s11_dur.to_csv(S11_OUT / "taskB_duration_summary.csv", index=False)
    _pfB = pd.concat(_rowsB, ignore_index=True)
    _pfB.to_csv(S11_OUT / "taskB_per_fold.csv", index=False)
    s11_sB = _pfB.groupby(["tissue", "arm", "model"])["auroc"].agg(["mean", "std"]).reset_index()
    s11_sB.columns = ["tissue", "arm", "model", "auroc_mean", "auroc_sd"]
    s11_sB.to_csv(S11_OUT / "taskB_summary.csv", index=False)
    s11_null = pd.concat(_nulls, ignore_index=True)
    s11_null.to_csv(S11_OUT / "taskB_null.csv", index=False)
    s11_bal = pd.DataFrame(_balance)
    s11_bal.to_csv(S11_OUT / "taskB_batch_balance.csv", index=False)
else:
    s11_sB = pd.read_csv(res("07_fusion", "taskB_summary.csv"))
    s11_null = pd.read_csv(res("07_fusion", "taskB_null.csv"))
    s11_bal = pd.read_csv(res("07_fusion", "taskB_batch_balance.csv"))
    s11_dur = pd.read_csv(res("07_fusion", "taskB_duration_summary.csv"))
print("Task B arms per tissue:", s11_sB.groupby("tissue").size().to_dict(),
      "| permutation draws per tissue:", s11_null.groupby("tissue").size().to_dict())

# %%
# ---- best arm vs best single vs best fusion, and vs the null — recomputed live from the summary + null tables ----
# from scripts/07_fusion_vs_baselines.py::main, the best-vs-null loop (verbatim, metric = "auroc")
metric = "auroc"
_rows = []
for t, st in s11_sB.groupby("tissue"):
    best = st.sort_values(f"{metric}_mean", ascending=False).iloc[0]
    single = st[st["arm"].str.startswith("single:")].sort_values(f"{metric}_mean", ascending=False).iloc[0]
    fus = st[~st["arm"].str.startswith("single:")].sort_values(f"{metric}_mean", ascending=False).iloc[0]
    q = s11_null.loc[s11_null["tissue"] == t, "max_auroc"]
    _rows.append({"tissue": t, "best_arm": f"{best['arm']}/{best['model']}", "best_auroc": best[f"{metric}_mean"],
                  "best_sd": best[f"{metric}_sd"], "best_single": f"{single['arm']}/{single['model']}",
                  "single_auroc": single[f"{metric}_mean"], "single_sd": single[f"{metric}_sd"],
                  "best_fusion": f"{fus['arm']}/{fus['model']}", "fusion_auroc": fus[f"{metric}_mean"],
                  "fusion_minus_single": fus[f"{metric}_mean"] - single[f"{metric}_mean"],
                  "fusion_beats_single_by_gt_sd": bool(fus[f"{metric}_mean"] - single[f"{metric}_mean"] > single[f"{metric}_sd"]),
                  "null_median_max_auroc": float(q.median()) if len(q) else np.nan,
                  "null_p95_max_auroc": float(q.quantile(0.95)) if len(q) else np.nan,
                  "p_perm": float((np.sum(q >= best[f"{metric}_mean"]) + 1) / (len(q) + 1)) if len(q) else np.nan,
                  "best_beats_null_p95": bool(len(q) and best[f"{metric}_mean"] > q.quantile(0.95))})
s11_bv = pd.DataFrame(_rows)
_cols = ["tissue", "best_arm", "best_auroc", "single_auroc", "single_sd", "fusion_auroc", "fusion_minus_single",
         "null_median_max_auroc", "null_p95_max_auroc", "p_perm", "fusion_beats_single_by_gt_sd", "best_beats_null_p95"]
print(s11_bv[_cols].to_string(index=False, float_format="%.3f"))
# recomputed vs published (default mode: the same inputs, so this checks the copied logic, not the models)
_pub = pd.read_csv(res("07_fusion", "taskB_best_vs_null.csv"))
_num = ["best_auroc", "fusion_auroc", "null_p95_max_auroc", "p_perm"]
print("\nmax |recomputed − published| over", _num, ":",
      float((s11_bv.set_index("tissue")[_num] - _pub.set_index("tissue")[_num]).abs().max().max()))
_nf, _nn = int(s11_bv["fusion_beats_single_by_gt_sd"].sum()), int(s11_bv["best_beats_null_p95"].sum())
print(f"fusion beats best single by > 1 sd: {_nf} of {len(s11_bv)} tissues (published {int(_pub['fusion_beats_single_by_gt_sd'].sum())})")
print(f"best arm beats the null's 95th percentile: {_nn} of {len(s11_bv)} tissues (published {int(_pub['best_beats_null_p95'].sum())})")
print(f"null 95th percentile of the best-of-arms AUROC, range over tissues: "
      f"{s11_bv['null_p95_max_auroc'].min():.3f}–{s11_bv['null_p95_max_auroc'].max():.3f}")
record("s11.taskB_n_fusion_beats_single", _nf, "11")
record("s11.taskB_n_beats_null", _nn, "11")
record("s11.taskB_null_p95_max", s11_bv["null_p95_max_auroc"].max(), "11")

# %%
# ---- Task B by training duration (best single-omic AUROC): a training effect should grow with duration, a batch should not ----
_wide = s11_dur.pivot(index="tissue", columns="duration", values="best_single_auroc")
_wide = _wide[[d for d in sorted(S11_DURATION_WEEKS, key=S11_DURATION_WEEKS.get) if d in _wide.columns]]
print("best single-omic AUROC, control vs each duration (animal-grouped CV):")
print(_wide.to_string(float_format="%.3f"))
print("\nmean over tissues:", _wide.mean().round(3).to_dict())
# TMT plex / channel balance of the Task B animals (Cramér's V of group vs plex or channel; near 1 = confounded)
print("\nTMT plex/channel balance (Cramér's V, group vs batch):")
print(s11_bal[[c for c in s11_bal.columns if c in ("tissue", "n_control", "n_8w", "V_tmt_plex", "V_tmt11_channel")]]
      .to_string(index=False, float_format="%.2f"))

# %% [markdown]
# **What Task B shows.** In every fusion tissue the best arm separates control from 8-week animals
# almost perfectly and beats the 95th percentile of the best-of-arms permutation null — but in no tissue
# does fusion beat the best single omic by more than that omic's fold sd. So: there is a real group
# difference, and one omic already captures it. The duration table shows the separation is already high
# after one week of training, so the duration series cannot, by itself, tell a fast training response
# from a batch (and by design only the 8-week animals share their collection dates with the controls).
#
# **What it does not show.** That the separation is *training*. The permutation null tests "is there any
# difference between the two groups", not "what caused it". That is what the batch check below is for.
# Note also that the null used the untuned (quick) arms with 100-tree forests, a cheaper version of
# the observed arms; the pipeline states this in its report.

# %%
# ---- batch check: helpers copied from scripts/07_batch_check.py ----
FUSION_TISSUES = ["CORTEX", "HEART", "KIDNEY", "LIVER", "LUNG", "SKM-GN", "WAT-SC"]
# from scripts/07_batch_check.py (verbatim constants)
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


# from scripts/07_batch_check.py::categorize (verbatim)
def categorize(col: str) -> str | None:
    if col in DROP:
        return None
    for pat, cat in CATEGORY_RULES:
        if re.search(pat, col):
            return cat
    return "other"


# from scripts/07_batch_check.py::parse_series (verbatim)
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


# from scripts/07_batch_check.py::cramers_v (verbatim)
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


# from scripts/07_batch_check.py::screen (verbatim)
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


# from scripts/07_batch_check.py::covariate_frame (verbatim)
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


# from scripts/07_batch_check.py::covariate_auroc (verbatim)
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
            with _s11_warn.catch_warnings():
                _s11_warn.simplefilter("ignore")
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


# from scripts/07_batch_check.py::main, the per-tissue table construction (restructured: the loop body up to the
# covariate frames is a function here, so the live SKM-GN check below and the RECOMPUTE branch share it)
def batch_tables(t, ph, tm, pmeta, trained_group="8w"):
    d = ph[(ph["tissue"] == t) & ph["viallabel"].isin(tm["viallabel"])].drop_duplicates("pid").copy()   # the animal's TRNSCRPT vial
    vials = tm[tm["viallabel"].isin(d["viallabel"])]
    if d["pid"].nunique() < 12 or d["group"].nunique() < 2:
        return None
    d = d.set_index("pid")
    y = (d["group"] == trained_group).astype(int).to_numpy()
    meta = pd.DataFrame({"pid": d.index, "sex": d["sex"].to_numpy(), "group": d["group"].to_numpy()}, index=d.index)
    meta["sex_group"] = meta["sex"] + "/" + meta["group"]
    n_splits = min(4, len(d) // 3)
    tv = vials.set_index("viallabel").reindex(d["viallabel"]).set_axis(d.index).copy()
    if "Lib_barcode_well" in tv.columns:   # plate position as two numbers (processing order is a classic batch axis)
        w = tv["Lib_barcode_well"].astype("string").str.upper().str.extract(r"^([A-H])(\d{1,2})$")
        tv["Lib_barcode_well_row"] = w[0].map({c: i + 1 for i, c in enumerate("ABCDEFGH")}).astype("string")
        tv["Lib_barcode_well_col"] = w[1]
    pheno_cols = [c for c in ph.columns if categorize(c)]
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
    frames = {fs: covariate_frame(tab, cols) for fs, (tab, cols) in tables.items()}
    X_all = pd.concat([f[0] for f in frames.values()], axis=1)
    frames["all_batch_covariates"] = (X_all, sum((f[1] for f in frames.values()), []), sum((f[2] for f in frames.values()), []))
    return SimpleNamespace(d=d, y=y, meta=meta, n_splits=n_splits, tv=tv, trn_cols_t=trn_cols_t, pheno_cols=pheno_cols,
                           tables=tables, frames=frames)


# the batch-check inputs (read-only): pheno.csv (control + 8w rows), transcript library meta, proteomics plex meta
s11_ph = pd.read_csv(C.RAW_DIR / "pheno.csv", dtype=str, low_memory=False)
s11_ph = s11_ph[s11_ph["group"].isin(["control", "8w"])]
s11_tm = pd.read_csv(C.META_DIR / "TRNSCRPT.csv", dtype=str, low_memory=False).drop_duplicates("viallabel")
s11_pm = pd.read_csv(C.META_DIR / "PROT.csv", dtype=str) if (C.META_DIR / "PROT.csv").exists() else None
print("batch-check inputs:", len(s11_ph), "control/8w pheno rows;", len(s11_tm), "TRNSCRPT library-meta vials")

# %%
# ---- (b) covariate-only classifiers and (c) the univariate screen: loaded (100-draw nulls), or recomputed ----
if RECOMPUTE:
    # from scripts/07_batch_check.py::main (trimmed: the markdown report and the duration read-out dropped; output → OUT/s11)
    _rng = np.random.default_rng(C.SEED)
    _animal_level = {c: ("animal" if s11_ph.groupby("pid")[c].nunique(dropna=True).max() <= 1 else "vial") for c in s11_ph.columns}
    _auc_rows, _screen_rows = [], []
    for t in FUSION_TISSUES:
        bt = batch_tables(t, s11_ph, s11_tm, s11_pm)
        if bt is None:
            continue
        for fs, (X, num, cat) in bt.frames.items():
            if X.shape[1] == 0:
                continue
            for model in ("logreg", "rf"):
                r = covariate_auroc(X, num, cat, bt.y, bt.meta, bt.n_splits, C.SEED, model,
                                    n_perm=100 if model == "logreg" else 0, rng=_rng)
                _auc_rows.append({"tissue": t, "feature_set": fs, "model": model, "n_animals": int(len(bt.d)), "n_features": int(X.shape[1]),
                                  "n_numeric": len(num), "n_categorical": len(cat), **r})
        _screen_rows += screen(bt.d[bt.pheno_cols], bt.y, "PHENO", _animal_level, t)
        _screen_rows += screen(bt.tv[bt.trn_cols_t], bt.y, "meta/TRNSCRPT", {c: "vial" for c in bt.trn_cols_t}, t)
        if s11_pm is not None:
            _screen_rows += screen(bt.tables["prot_plex_channel"][0], bt.y, "meta/PROT", {"tmt_plex": "vial", "tmt11_channel": "vial"}, t)
    s11_auc = pd.DataFrame(_auc_rows)
    s11_auc.to_csv(S11_OUT / "batch_covariate_auroc.csv", index=False)
    s11_sc = pd.DataFrame(_screen_rows)
    s11_sc.to_csv(S11_OUT / "batch_variable_screen.csv", index=False)
else:
    s11_auc = pd.read_csv(res("07_fusion", "batch_covariate_auroc.csv"))
    s11_sc = pd.read_csv(res("07_fusion", "batch_variable_screen.csv"), low_memory=False)

# "separates" flag, recomputed from the screen's statistics (script defaults --sep-auroc 0.9 --sep-v 0.8, p < 0.01)
_sep = s11_sc["missing_in_one_group_only"].astype(bool) | (((s11_sc["stat"] == "auroc") & (s11_sc["separation"] >= 0.9)
                                                             | (s11_sc["stat"] == "cramers_v") & (s11_sc["separation"] >= 0.8)) & (s11_sc["p_value"] < 0.01))
if "separates" in s11_sc.columns:
    print("recomputed 'separates' flag equals the published one:", bool((_sep == s11_sc["separates"].astype(bool)).all()))
s11_sc["separates"] = _sep
print(f"{len(s11_auc)} covariate-classifier rows ({s11_auc['tissue'].nunique()} tissues × feature sets × 2 models); "
      f"{len(s11_sc)} screened tissue × variable rows")

# %%
# ---- (d) the verdict rule, recomputed live from the two tables ----
# from scripts/07_batch_check.py::main, the per-tissue verdict loop (verbatim)
PURE = ("pheno_collection", "trnscrpt_library", "trnscrpt_depth", "prot_plex_channel")
auc, sc = s11_auc, s11_sc
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
s11_vt = pd.DataFrame(verdict_rows)
print(s11_vt[["tissue", "max_auroc_collection_library_plex", "max_auroc_qc_metrics", "max_auroc_all_covariates", "verdict"]]
      .to_string(index=False, float_format="%.3f"))
_pubv = pd.read_csv(res("07_fusion", "batch_conclusion.csv"))
print("\nverdicts equal to published batch_conclusion.csv:",
      bool((s11_vt.set_index("tissue")["verdict"] == _pubv.set_index("tissue")["verdict"]).all()))
print("published overall verdict:", _pubv["overall_verdict"].iloc[0])
_n_unres = int((s11_vt["verdict"] != "training").sum())
print(f"unresolvable tissues: {_n_unres} of {len(s11_vt)}")
record("s11.batch_n_unresolvable", _n_unres, "11")
del auc, sc, metric  # generic names from the verbatim loops; keep s11_auc / s11_sc / s11_vt

# %%
# ---- SKM-GN: "0.917" vs "0.958" — both are real rows of batch_covariate_auroc.csv ----
_sk = s11_auc[(s11_auc["tissue"] == "SKM-GN") & s11_auc["feature_set"].isin(["trnscrpt_depth", "trnscrpt_qc"])]
print(_sk[["tissue", "feature_set", "model", "n_features", "auroc_mean", "auroc_sd", "null_p95_auroc", "p_perm"]]
      .to_string(index=False, float_format="%.4f"))
_depth_rf = _sk.query("feature_set == 'trnscrpt_depth' and model == 'rf'").iloc[0]
_qc_lr = _sk.query("feature_set == 'trnscrpt_qc' and model == 'logreg'").iloc[0]
_qc_rf = _sk.query("feature_set == 'trnscrpt_qc' and model == 'rf'").iloc[0]
_vsk = s11_vt.set_index("tissue").loc["SKM-GN"]
print(f"\n(1) sequencing-depth covariates, RF: AUROC {_depth_rf['auroc_mean']:.3f} = the verdict rule's max over the 'pure' "
      f"(collection/library/depth/plex) sets for SKM-GN: {_vsk['max_auroc_collection_library_plex']:.3f}; its own null: "
      f"{'none (RF rows have no permutation null)' if pd.isna(_depth_rf['p_perm']) else _depth_rf['p_perm']}")
print(f"(2) library-QC fractions: logreg {_qc_lr['auroc_mean']:.3f}, RF {_qc_rf['auroc_mean']:.3f}; logreg permutation "
      f"p = {_qc_lr['p_perm']:.3f} (null 95th pct {_qc_lr['null_p95_auroc']:.3f})")
_dl = _sk.query("feature_set == 'trnscrpt_depth' and model == 'logreg'").iloc[0]
print(f"    for comparison, depth covariates with logreg: {_dl['auroc_mean']:.3f}, p = {_dl['p_perm']:.3f}")
_other = s11_auc[(s11_auc["tissue"] == "SKM-GN") & (s11_auc["p_perm"] <= 0.05)]
print("    SKM-GN logreg rows with p_perm ≤ 0.05:", [(r.feature_set, round(r.auroc_mean, 3), round(r.p_perm, 4)) for r in _other.itertuples()])
record("s11.skmgn_depth_rf", _depth_rf["auroc_mean"], "11")
record("s11.skmgn_qc_logreg", _qc_lr["auroc_mean"], "11")
record("s11.skmgn_qc_logreg_p", _qc_lr["p_perm"], "11")

# Live recompute of those three classifiers (no permutation null: the null's RNG stream runs through every
# tissue and set in order, so a single row's null cannot be replayed on its own).
_bt = batch_tables("SKM-GN", s11_ph, s11_tm, s11_pm)
for fs, model, pub in (("trnscrpt_depth", "rf", _depth_rf), ("trnscrpt_qc", "logreg", _qc_lr), ("trnscrpt_qc", "rf", _qc_rf)):
    X, num, cat = _bt.frames[fs]
    r = covariate_auroc(X, num, cat, _bt.y, _bt.meta, _bt.n_splits, C.SEED, model)
    print(f"  live {fs:15s} {model:6s}: AUROC {r['auroc_mean']:.4f}  (published {pub['auroc_mean']:.4f}); "
          f"{X.shape[1]} covariates, {len(_bt.d)} animals, {r['n_folds']} folds")
    record(f"s11.skmgn_{fs.split('_')[1]}_{model}_live", r["auroc_mean"], "11")
for _n in ("X", "num", "cat", "r", "fs", "model", "pub", "t", "m", "q", "best", "single", "fus", "verdict", "verdict_rows"):
    globals().pop(_n, None)   # generic loop names from the verbatim code

# %% [markdown]
# **What the batch check shows.** The verdict rule, rerun here on the loaded tables, reproduces the
# published verdicts: several tissues are "unresolvable" — their control-vs-8w separation cannot be
# attributed to training on these data — and the rest are attributed to training (or to anything else
# that differs by design between the arms, such as daily handling).
#
# **The two SKM-GN numbers.** Both are real rows of `batch_covariate_auroc.csv`, and they support different statements:
# - **The lower number** (row 1 of the printout) is the *sequencing-depth* covariates (`trnscrpt_depth`: raw, mapped and splice-junction
#   read counts) with a **random forest**. It is the maximum over the "pure" batch sets (collection,
#   library, depth, TMT plex) — the number the verdict rule uses — so it supports *"in SKM-GN a pure
#   processing covariate set classifies control vs trained, so the verdict is unresolvable"*. It has no
#   permutation null of its own (the script runs nulls for the logistic model only), and the logistic
#   model on the same depth covariates is not significant — so read it as a flag, not as a test.
# - **The higher number** (row 2) is the *library-QC fractions* (`trnscrpt_qc`: mapping, mRNA, intronic, mitochondrial
#   fractions, …) with **both** logistic regression and random forest. The logistic row has its own
#   within-sex permutation null, and it clears it (the p-value printed above). It supports *"library
#   composition differs between control and trained muscle, beyond chance"*. The all-covariates
#   logistic row clears its null too, but it contains the QC set. Whether a composition difference is
#   batch or biology (a higher mitochondrial read fraction is what training-induced mitochondrial
#   biogenesis would produce) cannot be told from these data.
#
# **What it does not show.** That the other tissues are free of batch: "training" here means "no
# *recorded* covariate separates the groups", at the small per-tissue n printed above, where the logistic nulls'
# 95th percentiles are themselves high (the table above). Unrecorded batches are not tested.
