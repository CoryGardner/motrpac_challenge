#!/usr/bin/env python
"""Part 6 robustness: repeat the four held-out-8w designs over N random fit/calibration animal splits
(fit_calibration_split seeds SEED+1..SEED+N; model seed fixed), so the (a) vs (b) comparison includes
calibration-set randomness, which the single-split bootstrap in tc56 does not.
Output: results/15_time_course/5_6_fingerprint/part6_seeds.csv, part6_seeds_summary.csv
Run from code/pipeline:  PYTHONPATH=src python investigations/time_course/tc6_fingerprint_seeds.py [N]
"""
from __future__ import annotations

import os

os.environ.setdefault("LOKY_MAX_CPU_COUNT", "8")
os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from tc56_fingerprint_shift import K, OUT, SEED, fit_and_calibrate, metrics, score_vials  # noqa: E402
from tfp import cli, io  # noqa: E402
from tfp.splits import assert_no_group_leak, fit_calibration_split  # noqa: E402

DESIGNS = {"a_control+1w": ["control", "1w"], "b_control+1w+2w+4w": ["control", "1w", "2w", "4w"],
           "c_control+4w": ["control", "4w"], "d_control_only": ["control"]}


def main() -> None:
    n_seeds = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    t0 = time.time()
    om = io.stack_tissues("TRNSCRPT", source=cli.resolve_source("TRNSCRPT", "auto"), join="inner", pheno=io.load_pheno(),
                          verbose=False)
    X, y, g, meta = om.X.to_numpy(dtype=float), om.meta["tissue"].astype(str).to_numpy(), om.groups(), om.meta
    grp = meta["group"].astype(str).to_numpy()
    te8 = np.flatnonzero(grp == "8w")
    rows = []
    for s in range(1, n_seeds + 1):
        for dname, groups in DESIGNS.items():
            trg = np.flatnonzero(np.isin(grp, groups))
            assert_no_group_leak(meta, trg, te8)
            f_idx, c_idx = fit_calibration_split(meta, trg, 0.3, SEED + s)
            for arm, k in (("full", None), (f"panel_k{K}", K)):
                est, cls, qhat, info = fit_and_calibrate(X, y, g, meta, f_idx, c_idx, k)
                dv = score_vials(est, cls, qhat, X, y, meta, te8)
                rows.append({"split_seed": SEED + s, "design": dname, "arm": arm, "n_fit_animals": info["n_fit_animals"],
                             "n_cal_animals": info["n_cal_animals"], "n_cal_vials": info["n_cal_vials"], "qhat": qhat,
                             "n_test_animals": dv["pid"].nunique(), "n_test_vials": len(dv), **metrics(dv)})
        print(f"  seed {s}/{n_seeds} ({time.time() - t0:.0f}s)", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "part6_seeds.csv", index=False)
    summ = []
    for (dname, arm), d in df.groupby(["design", "arm"]):
        r = {"design": dname, "arm": arm, "n_splits": len(d), "n_test_animals": int(d["n_test_animals"].iloc[0]),
             "n_cal_animals": int(d["n_cal_animals"].iloc[0])}
        for key in ("accuracy", "coverage", "empty_rate", "mean_set_size", "qhat"):
            r[f"{key}_mean"], r[f"{key}_sd"] = float(d[key].mean()), float(d[key].std(ddof=1))
            r[f"{key}_min"], r[f"{key}_max"] = float(d[key].min()), float(d[key].max())
        summ.append(r)
    # paired over seeds: design minus (b)
    for arm in df["arm"].unique():
        piv = df[df["arm"] == arm].pivot(index="split_seed", columns="design", values="coverage")
        for dname in ("a_control+1w", "c_control+4w", "d_control_only"):
            dd = piv[dname] - piv["b_control+1w+2w+4w"]
            summ.append({"design": f"coverage[{dname}] - coverage[b]", "arm": arm, "n_splits": len(dd),
                         "coverage_mean": float(dd.mean()), "coverage_sd": float(dd.std(ddof=1)),
                         "coverage_min": float(dd.min()), "coverage_max": float(dd.max()),
                         "n_splits_positive": int((dd > 0).sum())})
    pd.DataFrame(summ).to_csv(OUT / "part6_seeds_summary.csv", index=False)
    print(f"done ({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()
