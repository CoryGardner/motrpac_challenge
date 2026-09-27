#!/usr/bin/env python
"""Phase 34 — export the 20-gene transfer model so a browser can score new samples ("Score your own samples").

The scoring model is the one that carries the guarantee in the rat BodyMap transfer (scripts/12_bodymap_validate.py):
a PanelModels fit on the 35 MoTrPAC animals outside the calibration split (k = 20, z-score representation, multinomial
logistic regression with C = 0.1), whose prediction sets are calibrated on the other 15 animals. The "all-animal call"
is the same panel fit on all 50 animals. This script re-fits both with the same code path and seed, checks that the
scoring model's probabilities for every BodyMap sample equal the --save-scores regeneration
(results/31_site_regen/12_bodymap/scores_target_probs.csv) and that a plain per-gene z-score + softmax on the 20 genes
alone reproduces them (which is what the browser does), and writes results/34_panel_model/:
  genes.csv                  feature_ID, symbol, position (the model's feature order), in_all_animal_panel
  coefficients.csv           one row per class: intercept and one column per feature_ID
  source_z.csv               per gene: the mean and scale of the MoTrPAC z-scoring (fallback for uploads too small to z-score)
  call_model_genes.csv, call_model_coefficients.csv, call_model_source_z.csv   the all-animal panel (the "all-animal call")
  validation.csv             max |Δp| and call agreement against the regeneration, per model
  NOTES.md
Usage: python scripts/34_panel_model.py [--out DIR]   (needs data/raw, data/external and the regeneration run)
"""
from __future__ import annotations

import importlib.util
import json

import numpy as np
import pandas as pd

from tfp import cli, config as C, io, report
from tfp.transfer import PanelModels

K = 20


def _bm12():
    spec = importlib.util.spec_from_file_location("bm12", C.ROOT / "scripts" / "12_bodymap_validate.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def dump_model(m: PanelModels, name: str, out, sym: dict, prefix: str, all_panel: set) -> pd.DataFrame:
    idx, clf = m.fits[name]
    genes = list(m.genes_pre[idx])
    classes = list(clf.classes_)
    g = pd.DataFrame({"feature_ID": genes, "symbol": [sym.get(x, x) for x in genes], "position": range(len(genes)),
                      "in_all_animal_panel": [x in all_panel for x in genes]})
    g.to_csv(out / f"{prefix}genes.csv", index=False)
    coef = pd.DataFrame(clf.coef_, columns=genes)
    coef.insert(0, "intercept", clf.intercept_)
    coef.insert(0, "class", classes)
    coef.to_csv(out / f"{prefix}coefficients.csv", index=False)
    pd.DataFrame({"feature_ID": genes, "symbol": g["symbol"], "mean": m.scaler.mean_[idx], "scale": m.scaler.scale_[idx]}).to_csv(out / f"{prefix}source_z.csv", index=False)
    return g


def softmax(z):
    z = z - z.max(axis=1, keepdims=True)
    e = np.exp(z)
    return e / e.sum(axis=1, keepdims=True)


def main() -> None:
    ap = cli.common_parser("Export the 20-gene transfer model (scoring tool)")
    ap.add_argument("--grid", default="20,50")
    ap.add_argument("--cal-frac", type=float, default=0.3)
    args = ap.parse_args()
    cli.banner("34_panel_model", args)
    out = cli.outdir("34_panel_model", args.out)
    grid = [int(k) for k in args.grid.split(",")]
    bm12 = _bm12()

    pheno = io.load_pheno()
    om = io.stack_tissues("TRNSCRPT", source="counts", pheno=pheno, verbose=False)
    Xb_counts, mb = bm12.load_bodymap()
    shared = [g for g in om.X.columns if g in set(Xb_counts.columns)]
    f2g = io.load_feature_to_gene()
    sym = io.map_to_gene_symbols(om.X.columns, f2g)
    lb = io.log_cpm(Xb_counts, log=True)
    y = om.meta["tissue"].astype(str).to_numpy()
    g = om.groups()
    L_src = om.X[shared].to_numpy(dtype=float)
    print(f"  shared genes {len(shared)}; BodyMap {len(mb)} samples")

    # the same split as the transfer phase (seeded), then the two fits
    rng = np.random.default_rng(args.seed)
    uniq = np.unique(g)
    rng.shuffle(uniq)
    n_cal = max(1, int(round(args.cal_frac * len(uniq))))
    cal_animals = set(uniq[:n_cal])
    is_cal = np.array([a in cal_animals for a in g])
    fit_idx, cal_idx = np.flatnonzero(~is_cal), np.flatnonzero(is_cal)
    pm = PanelModels(L_src, y, shared, grid, args.prefilter, full=False)
    pm_c = PanelModels(L_src[fit_idx], y[fit_idx], shared, grid, args.prefilter, full=False)
    name = f"k{K}"
    all_panel = set(pm.panel_genes(K))
    genes_c = dump_model(pm_c, name, out, sym, "", all_panel)
    dump_model(pm, name, out, sym, "call_model_", all_panel)
    print(f"  scoring model: {int(len(np.unique(g[fit_idx])))} fit animals, {int(len(np.unique(g[cal_idx])))} calibration animals; "
          f"{int(genes_c['in_all_animal_panel'].sum())} of {K} genes shared with the all-animal panel")

    # validation against the regeneration and against the browser's arithmetic
    Lt_c, Zt_c = pm_c.target_matrices(lb[shared].to_numpy(dtype=float))
    p_ref, cls = pm_c.proba(name, Lt_c, Zt_c)
    Lb_pre, Zb_pre = pm.target_matrices(lb[shared].to_numpy(dtype=float))
    p_all, _ = pm.proba(name, Lb_pre, Zb_pre)
    regen = C.RESULTS_DIR / "31_site_regen" / "12_bodymap" / "scores_target_probs.csv"
    rows = []
    if regen.exists():
        r = pd.read_csv(regen)
        r = r[r["model"] == name].set_index("sample").loc[mb.index]
        pcols = [f"p_{c}" for c in cls]
        d_regen = float(np.abs(r[pcols].to_numpy(dtype=float) - p_ref).max())
        calls_agree = float(np.mean(r[pcols].to_numpy().argmax(axis=1) == p_ref.argmax(axis=1)))
        call_all_agree = float(np.mean(r["pred_all_animals"].to_numpy() == np.array(cls)[p_all.argmax(axis=1)])) if "pred_all_animals" in r else np.nan
        rows.append({"model": "scoring (35 fit / 15 calibration animals)", "compared_to": "results/31_site_regen/12_bodymap/scores_target_probs.csv",
                     "n_samples": int(len(mb)), "max_abs_p_diff": d_regen, "call_agreement": calls_agree})
        rows.append({"model": "all-animal call", "compared_to": "results/31_site_regen/12_bodymap/scores_target_probs.csv (pred_all_animals)",
                     "n_samples": int(len(mb)), "max_abs_p_diff": np.nan, "call_agreement": call_all_agree})
        assert d_regen < 1e-9, d_regen
    # the browser's arithmetic: z-score the 20 genes within the upload (population sd), then softmax(z W^T + b)
    idx, clf = pm_c.fits[name]
    genes = list(pm_c.genes_pre[idx])
    X = lb[genes].to_numpy(dtype=float)
    Z = (X - X.mean(axis=0)) / X.std(axis=0)
    p_js = softmax(Z @ clf.coef_.T + clf.intercept_)
    d_js = float(np.abs(p_js - p_ref).max())
    rows.append({"model": "scoring, per-gene z-score + softmax on the 20 genes only (the browser's arithmetic)", "compared_to": "PanelModels.proba",
                 "n_samples": int(len(mb)), "max_abs_p_diff": d_js, "call_agreement": float(np.mean(p_js.argmax(axis=1) == p_ref.argmax(axis=1)))})
    assert d_js < 1e-9, d_js
    pd.DataFrame(rows).to_csv(out / "validation.csv", index=False)
    print(pd.DataFrame(rows).to_string(index=False))
    notes = (f"# Phase 34 — the 20-gene transfer model, exported for the browser\n\n"
             f"Scoring model: PanelModels (z-score representation, multinomial logistic regression C = 0.1) fit on the {int(len(np.unique(g[fit_idx])))} "
             f"MoTrPAC animals outside the transfer phase's calibration split (seed {args.seed}, cal_frac {args.cal_frac}); its prediction sets are "
             f"calibrated on the other {int(len(np.unique(g[cal_idx])))} animals (results/31_site_regen/12_bodymap/scores_calibration.csv). "
             f"All-animal call: the same panel fit on all 50 animals. A new sample is scored by z-scoring each of the {K} genes within the uploaded "
             f"set (population sd, as StandardScaler does; the source mean/scale in source_z.csv is the fallback for small uploads), then "
             f"softmax(z · coefficients^T + intercept). validation.csv shows the re-fit reproduces the regeneration and that the 20-gene "
             f"arithmetic reproduces PanelModels exactly.\n")
    (out / "NOTES.md").write_text(notes)
    report.add_section("34 · The 20-gene transfer model for the scoring tool", notes.split("\n", 1)[1], params={"k": K, "seed": args.seed})
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
