#!/usr/bin/env python
"""Phase 40 — validation numbers for the "Check samples" product page (site/index.html). No model is fit here.

The product scores uploads with the exported 20-gene transfer model (site/data/panel_model.json, phase 34) exactly as
the browser does: log2 CPM → per-gene z-score → softmax(z · Wᵀ + b) → LAC sets from the 15-animal calibration scores.
This script measures, from existing inputs only:

  scaling.csv        rat BodyMap scored three ways: the whole set z-scored together (the pipeline, "within_all"),
                     each organ z-scored alone ("within_organ_alone": what a single-tissue upload would get), and
                     each sample z-scored with the MoTrPAC means and scales ("reference"); per organ and pooled, all
                     ages and 21-week adults; accuracy, coverage, set size and empty sets at α = 0.10.
  flag_rates.csv     the claimed-label flag: false-Mismatch rate on correctly labelled samples and detection of
                     simulated label swaps (fixed seed, 1,000 swaps between samples of different tissues), per α in
                     {0.05, 0.10, 0.20}; in-study from the held-out k20 scores of the phase-06 regeneration
                     (results_frozen/31_site_regen/06_conformal/TRNSCRPT/scores_*.csv, each fold's own calibration),
                     and on BodyMap from the product model's own scores.
  venacv_cases.csv   the held-out vena cava vials the 20-gene model calls brown fat, with their sets and flag status
                     under the claimed label "vena cava", and whether the consortium flagged them.
  pca_*.csv          a 2-D PCA of the 20-gene z-space (MoTrPAC z-scoring of the model) fitted on the MoTrPAC
                     reference vials of site/data/expr_motrpac.json: centre, loadings and reference coordinates.
  fixtures           tests/fixtures/bodymap_counts_subset.csv.gz (real BodyMap counts, every gene, a few samples)
                     and tests/fixtures/product_parity.json (the pipeline's log2 CPM of the panel genes for them), for
                     the browser's count → log2 CPM parity test.

Default output: results_product/40_product/ (committed, like results_multiomic/: results/ is git-ignored and
results_frozen/ is a fixed snapshot). Usage: python scripts/40_product_validation.py [--out DIR]
"""
from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd

from tfp import config as C, io
from tfp.conformal import conformal_quantile

ROOT = C.ROOT
SITE = ROOT / "site" / "data"
FZ = ROOT / "results_frozen"
ALPHAS = (0.05, 0.10, 0.20)
SEED = 20260927
N_SWAPS = 1000


def softmax(z):
    z = z - z.max(axis=1, keepdims=True)
    e = np.exp(z)
    return e / e.sum(axis=1, keepdims=True)


def probs(Z, model):
    return softmax(Z @ np.asarray(model["coef"]).T + np.asarray(model["intercept"]))


def z_within(X):
    mu, sd = X.mean(axis=0), X.std(axis=0)          # population sd, as StandardScaler / score.js
    return np.where(sd > 0, (X - mu) / np.where(sd > 0, sd, 1), 0.0)


def z_reference(X, model):
    mu, sc = np.asarray(model["source_z"]["mean"]), np.asarray(model["source_z"]["scale"])
    return (X - mu) / sc


def claim_status(in_set: np.ndarray, empty: np.ndarray) -> np.ndarray:
    """Claim status per sample: 'consistent' (claimed ∈ set), 'mismatch' (set non-empty, claimed ∉ set), 'cant_confirm' (empty)."""
    return np.where(empty, "cant_confirm", np.where(in_set, "consistent", "mismatch"))


def flag_block(setting, P, true_sets, q_per_sample, alpha, rng_seed, extra):
    """P: n × C probabilities; true_sets: list of index lists (the correct claim, 1 class or a super-class);
    q_per_sample: threshold per sample. Returns one row with false-flag and swap-detection rates."""
    n = len(P)
    S = P >= (1 - q_per_sample)[:, None]
    empty = ~S.any(axis=1)
    in_true = np.array([S[i, t].any() for i, t in enumerate(true_sets)])
    st = claim_status(in_true, empty)
    rng = np.random.default_rng(rng_seed)
    key = np.array([",".join(map(str, t)) for t in true_sets])
    det_vial = det_pair = both = miss = cc = 0
    k = 0
    while k < N_SWAPS:
        i, j = rng.integers(0, n, size=2)
        if key[i] == key[j] or set(true_sets[i]) & set(true_sets[j]):
            continue
        k += 1
        # vial i now claims j's tissue and vice versa
        si = "cant_confirm" if empty[i] else ("consistent" if S[i, true_sets[j]].any() else "mismatch")
        sj = "cant_confirm" if empty[j] else ("consistent" if S[j, true_sets[i]].any() else "mismatch")
        det_vial += (si == "mismatch") + (sj == "mismatch")
        det_pair += (si == "mismatch") or (sj == "mismatch")
        both += (si == "mismatch") and (sj == "mismatch")
        miss += (si == "consistent") + (sj == "consistent")
        cc += (si == "cant_confirm") + (sj == "cant_confirm")
    return {"setting": setting, "alpha": alpha, "n_samples": n, **extra,
            "consistent_rate": float((st == "consistent").mean()), "false_mismatch_rate": float((st == "mismatch").mean()),
            "cant_confirm_rate": float((st == "cant_confirm").mean()), "n_false_mismatch": int((st == "mismatch").sum()),
            "n_swaps": N_SWAPS, "seed": rng_seed, "swap_vial_detection_rate": det_vial / (2 * N_SWAPS),
            "swap_pair_detection_rate": det_pair / N_SWAPS, "swap_pair_both_flagged_rate": both / N_SWAPS,
            "swap_vial_missed_rate": miss / (2 * N_SWAPS), "swap_vial_cant_confirm_rate": cc / (2 * N_SWAPS)}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", default=str(ROOT / "results_product" / "40_product"))
    ap.add_argument("--fixtures", default=str(ROOT / "tests" / "fixtures"))
    args = ap.parse_args()
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    fx = Path(args.fixtures); fx.mkdir(parents=True, exist_ok=True)
    model = json.loads((SITE / "panel_model.json").read_text())
    classes = model["classes"]; ci = {c: i for i, c in enumerate(classes)}
    gids = [g["id"] for g in model["genes"]]
    cal = np.asarray(model["calibration"]["scores"], dtype=float)

    # ---- BodyMap: counts → log2 CPM (io.log_cpm, total library over every gene) → the 20 genes ---------------------
    counts = pd.read_csv(C.EXTERNAL_DIR / "bodymap_counts.csv", index_col=0)
    counts.index = counts.index.astype(str)
    meta = pd.read_csv(C.EXTERNAL_DIR / "bodymap_meta.csv", dtype=str).set_index("sample").loc[counts.columns]
    meta["stage_weeks"] = meta["stage_weeks"].astype(int)
    lcpm = io.log_cpm(counts.T, log=True)
    X = lcpm[gids].to_numpy(dtype=float)
    organ_map = json.loads((FZ / "31_site_regen" / "12_bodymap" / "organ_map.json").read_text())
    organs = meta["organ"].to_numpy()
    # sanity: the whole-set within scaling reproduces the pipeline's BodyMap probabilities
    ref = pd.read_csv(FZ / "31_site_regen" / "12_bodymap" / "scores_target_probs.csv")
    ref = ref[ref["model"] == "k20"].set_index("sample").loc[counts.columns]
    P_all = probs(z_within(X), model)
    dmax = float(np.abs(P_all - ref[[f"p_{c}" for c in classes]].to_numpy()).max())
    assert dmax < 1e-9, f"within_all does not reproduce the pipeline (max |Δp| {dmax})"

    q10 = conformal_quantile(cal, 0.10)
    rows = []
    subsets = {"all_ages": np.ones(len(meta), bool), "adult_21wk": (meta["stage_weeks"] == 21).to_numpy()}
    for sub, msk in subsets.items():
        P_modes = {"within_all": P_all, "reference": probs(z_reference(X, model), model)}
        P_alone = np.zeros_like(P_all)
        for o in np.unique(organs[msk]):
            idx = np.where(msk & (organs == o))[0]
            P_alone[idx] = probs(z_within(X[idx]), model)
        P_modes["within_organ_alone"] = P_alone
        for mode, P in P_modes.items():
            S = P >= 1 - q10
            call = P.argmax(axis=1)
            mapped_rows = []
            for o in sorted(np.unique(organs[msk])):
                idx = np.where(msk & (organs == o))[0]
                tgt = organ_map.get(o)
                r = {"subset": sub, "mode": mode, "organ": o, "n": int(len(idx)), "mapped": tgt is not None,
                     "rat_classes": ";".join(tgt) if tgt else "", "alpha": 0.10, "q": q10,
                     "avg_set_size": float(S[idx].sum(axis=1).mean()), "frac_empty": float((~S[idx].any(axis=1)).mean())}
                if tgt:
                    ti = [ci[c] for c in tgt]
                    r["accuracy"] = float(np.isin(call[idx], ti).mean())
                    r["coverage"] = float(S[np.ix_(idx, ti)].any(axis=1).mean())
                    mapped_rows.append(idx)
                rows.append(r)
            idx = np.concatenate(mapped_rows)
            tsets = [[ci[c] for c in organ_map[o]] for o in organs[idx]]
            rows.append({"subset": sub, "mode": mode, "organ": "ALL_MAPPED", "n": int(len(idx)), "mapped": True, "rat_classes": "", "alpha": 0.10, "q": q10,
                         "accuracy": float(np.mean([call[i] in t for i, t in zip(idx, tsets)])),
                         "coverage": float(np.mean([S[i, t].any() for i, t in zip(idx, tsets)])),
                         "avg_set_size": float(S[idx].sum(axis=1).mean()), "frac_empty": float((~S[idx].any(axis=1)).mean())})
    scaling = pd.DataFrame(rows)
    scaling.to_csv(out / "scaling.csv", index=False)

    # ---- flag rates ---------------------------------------------------------------------------------------------------
    fr = []
    sp = pd.read_csv(FZ / "31_site_regen" / "06_conformal" / "TRNSCRPT" / "scores_test_probs.csv", dtype={"viallabel": str, "pid": str})
    sp = sp[sp["model"] == "k20"].reset_index(drop=True)
    sc = pd.read_csv(FZ / "31_site_regen" / "06_conformal" / "TRNSCRPT" / "scores_calibration.csv", dtype={"viallabel": str})
    sc = sc[(sc["model"] == "k20") & (sc["calibration"] == "pooled")]
    Pm = sp[[f"p_{c}" for c in classes]].to_numpy(dtype=float)
    tm = [[ci[t]] for t in sp["tissue"]]
    for a in ALPHAS:
        qf = {f: conformal_quantile(g["score_lac"].to_numpy(dtype=float), a) for f, g in sc.groupby("fold")}
        qs = sp["fold"].map(qf).to_numpy(dtype=float)
        fr.append(flag_block("motrpac_heldout", Pm, tm, qs, a, SEED,
                             {"n_individuals": int(sp["pid"].nunique()), "calibration": "each fold's pooled k20 calibration scores (phase-06 design)",
                              "source": "results_frozen/31_site_regen/06_conformal/TRNSCRPT/scores_test_probs.csv, scores_calibration.csv"}))
    for sub, msk in subsets.items():
        m2 = msk & np.array([organ_map.get(o) is not None for o in organs])
        Pb = P_all[m2]
        tb = [[ci[c] for c in organ_map[o]] for o in organs[m2]]
        for a in ALPHAS:
            q = conformal_quantile(cal, a)
            fr.append(flag_block(f"bodymap_{sub}", Pb, tb, np.full(len(Pb), q), a, SEED,
                                 {"n_individuals": int(meta.loc[m2, "replicate"].str.rsplit("_", n=1).str[-1].add(meta.loc[m2, "sex"] + meta.loc[m2, "stage_weeks"].astype(str)).nunique()),
                                  "calibration": "the product model's 15-animal calibration scores (site/data/panel_model.json)",
                                  "source": "data/external/bodymap_counts.csv scored with site/data/panel_model.json (equals results_frozen/31_site_regen/12_bodymap/scores_target_probs.csv)"}))
    pd.DataFrame(fr).to_csv(out / "flag_rates.csv", index=False)

    # ---- the real case: vena cava vials called brown fat -------------------------------------------------------------
    flagged = set(pd.read_csv(FZ / "15_time_course" / "design" / "flagged_vials.csv", dtype=str)["viallabel"])
    q10f = {f: conformal_quantile(g["score_lac"].to_numpy(dtype=float), 0.10) for f, g in sc.groupby("fold")}
    vb = sp[(sp["tissue"] == "VENACV") & (Pm.argmax(axis=1) == ci["BAT"])]
    vrows = []
    for _, r in vb.iterrows():
        p = r[[f"p_{c}" for c in classes]].to_numpy(dtype=float)
        q = q10f[r["fold"]]
        s = [c for c, v in zip(classes, p) if v >= 1 - q]
        st = "cant_confirm" if not s else ("consistent" if "VENACV" in s else "mismatch")
        vrows.append({"viallabel": r["viallabel"], "pid": r["pid"], "sex": r["sex"], "group": r["group"], "fold": int(r["fold"]),
                      "p_BAT": float(p[ci["BAT"]]), "p_VENACV": float(p[ci["VENACV"]]), "alpha": 0.10, "q": q, "set": ";".join(s),
                      "status_claimed_venacv": st, "consortium_flagged": r["viallabel"] in flagged})
    pd.DataFrame(vrows).to_csv(out / "venacv_cases.csv", index=False)

    # ---- reference PCA of the 20-gene z-space ------------------------------------------------------------------------
    em = json.loads((SITE / "expr_motrpac.json").read_text())
    gi = {g: i for i, g in enumerate(em["genes"])}
    Xm = np.array([em["values"][gi[g]] for g in gids], dtype=float).T          # vials × 20 genes (log2 CPM, 4 d.p.)
    Zm = z_reference(Xm, model)
    centre = Zm.mean(axis=0)
    U, s, Vt = np.linalg.svd(Zm - centre, full_matrices=False)
    L = Vt[:2]
    for k in range(2):                                                         # sign convention: largest |loading| positive
        if L[k, np.argmax(np.abs(L[k]))] < 0:
            L[k] = -L[k]
    coords = (Zm - centre) @ L.T
    ev = (s ** 2) / (s ** 2).sum()
    pd.DataFrame({"feature_ID": gids, "symbol": [g["symbol"] for g in model["genes"]], "centre": centre, "PC1": L[0], "PC2": L[1]}).to_csv(out / "pca_loadings.csv", index=False)
    pd.DataFrame({"viallabel": [s_["id"] for s_ in em["samples"]], "tissue": [s_["tissue"] for s_ in em["samples"]], "PC1": coords[:, 0], "PC2": coords[:, 1]}).to_csv(out / "pca_reference_coords.csv", index=False)
    pd.DataFrame({"PC": ["PC1", "PC2"], "explained": ev[:2], "n_reference_vials": len(Zm), "z_space": "MoTrPAC z-scoring of the model (panel_model.json source_z)"}).to_csv(out / "pca_summary.csv", index=False)

    # ---- parity fixtures: real counts, every gene, a few samples ------------------------------------------------------
    pick = [meta.index[(meta["organ"] == o) & (meta["stage_weeks"] == 21)][0] for o in ("Liver", "Brain", "Testes", "Thymus")]
    sub = counts[pick]
    with gzip.open(fx / "bodymap_counts_subset.csv.gz", "wt", compresslevel=9) as f:
        sub.to_csv(f)
    lsub = io.log_cpm(sub.T, log=True)
    (fx / "product_parity.json").write_text(json.dumps({
        "note": "tests/fixtures/bodymap_counts_subset.csv.gz (genes × samples, real rat BodyMap counts, every gene) → io.log_cpm (total library, log2(CPM + 1)); expected values for the 20 panel genes. Written by scripts/40_product_validation.py.",
        "samples": pick, "genes": gids, "library_size": [float(sub[c].sum()) for c in pick],
        "log2cpm": [[float(lsub.loc[c, g]) for g in gids] for c in pick]}, indent=1))

    (out / "NOTES.md").write_text(
        "# Phase 40 — product validation\n\nWritten by `scripts/40_product_validation.py`; no model is fit. The within_all BodyMap "
        f"probabilities reproduce `results_frozen/31_site_regen/12_bodymap/scores_target_probs.csv` (max |Δp| {dmax:.1e}).\n\n"
        "- `scaling.csv`: BodyMap under three scalings (whole set z-scored together / each organ alone / MoTrPAC reference z-scoring).\n"
        "- `flag_rates.csv`: false-Mismatch rate on correct labels and detection of 1,000 simulated label swaps (seed "
        f"{SEED}), α ∈ {{0.05, 0.10, 0.20}}; in-study (held-out k20 scores, phase-06 design) and BodyMap (the product model).\n"
        "- `venacv_cases.csv`: held-out vena cava vials called brown fat, their sets and flag status under the claim 'vena cava'.\n"
        "- `pca_*.csv`: the reference map (MoTrPAC z-space of the model, 2 components).\n")
    print(f"wrote {out} (scaling {len(scaling)} rows, flag {len(fr)} rows, venacv {len(vrows)} rows, PCA {len(Zm)} vials; within_all max |Δp| {dmax:.1e})")


if __name__ == "__main__":
    main()
