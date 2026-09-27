#!/usr/bin/env python
"""Phase 15a — study-design dates and physiology of the training groups (PHENO only, no omics).

Design. Per group × sex: number of animals, arrival date, first and last sacrifice date, number of distinct
sacrifice days and the span in days, days from arrival to first sacrifice, median and range of the time of death.
Per duration × sex (and pooled over sex): the gap in days between that group's median sacrifice date and the
controls' median, the same for the arrival date, and whether the group shares the controls' arrival cohort.
Two control-vs-8w tables: animals per sacrifice day and the time of death (Mann–Whitney per sex).

Physiology. Per animal: VO2max pre/post and change, NMR fat % pre/post and change, lean and weight change and
the study days of the tests (columns `vo2.max.test.*`, `nmr.testing.*`, `training.day1_days` of PHENO).
Per group × sex: n, median and quartiles. Group tests: 8w vs control (cohort-matched) and 4w vs control
(different arrival cohort and sacrifice season), within sex, two-sided Mann–Whitney on the per-animal values.

Identities checked (the script raises if any fails): 147 animals; the group parsed by `io.load_pheno()` equals
`key.anirandgroup`; `calculated.variables.vo2_max_change` = `vo2_max_2 − vo2_max_1`, `pct_body_fat_change` =
`nmr_fat_2 − nmr_fat_1` and `pct_body_lean_change` = `nmr_lean_2 − nmr_lean_1` (to 1e-6); `days_vo2_1 <
training.day1_days < days_vo2_2` for every trained animal with a post test; 1w and 2w animals have no
post-training VO2max or NMR measurement.

Outputs: results/15_time_course/design/{design_by_group_sex, design_contrasts, control_8w_by_sacrifice_day,
control_8w_time_of_death}.csv, results/15_time_course/physiology/{physiology_per_animal, physiology_by_group,
physiology_group_tests}.csv, results/15_time_course/design/NOTES.md. Seconds; `--quick` changes nothing.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

from tfp import cli, config as C, io, report

GROUP_MAP = {"Eight-week program Control Group": "control", "One-week program": "1w", "Two-week program": "2w",
             "Four-week program": "4w", "Eight-week program Training Group": "8w"}
DURATIONS = ["1w", "2w", "4w", "8w"]
PHYS_VARS = ["vo2max_pre", "vo2max_post", "vo2max_change", "fat_pre", "fat_post", "fat_change", "lean_change",
             "weight_change", "days_vo2_pre", "days_vo2_post", "days_train_start"]
TEST_VARS = ["vo2max_change", "fat_change", "lean_change", "weight_change", "vo2max_pre", "fat_pre"]
TEST_CONTRASTS = [("8w", True), ("4w", False)]      # (trained group, cohort-matched to the controls?)


def animals() -> pd.DataFrame:
    """One row per animal with the parsed dates, times and physiology columns; identities asserted."""
    ph = io.load_pheno().reset_index()
    a = ph.drop_duplicates("pid").copy()
    assert len(a) == 147, f"expected 147 animals, found {len(a)}"
    assert (a["key.anirandgroup"].map(GROUP_MAP) == a["group"]).all(), "group parse disagrees with key.anirandgroup"
    a["d_arrive"] = pd.to_datetime(a["key.d_arrive"], format="%d%b%Y")
    a["d_sac"] = pd.to_datetime(a["key.d_sacrifice"], format="%d%b%Y")
    a["death_hour"] = pd.to_timedelta(a["specimen.collection.t_death"], errors="coerce").dt.total_seconds() / 3600
    num = lambda c: pd.to_numeric(a[c], errors="coerce")
    a["days_vo2_pre"], a["days_vo2_post"] = num("vo2.max.test.days_vo2_1"), num("vo2.max.test.days_vo2_2")
    a["days_train_start"] = num("training.day1_days")
    a["days_nmr_pre"], a["days_nmr_post"] = num("nmr.testing.days_nmr_1"), num("nmr.testing.days_nmr_2")
    a["vo2max_pre"], a["vo2max_post"] = num("vo2.max.test.vo2_max_1"), num("vo2.max.test.vo2_max_2")
    a["vo2max_change"] = a["vo2max_post"] - a["vo2max_pre"]
    a["fat_pre"], a["fat_post"] = num("nmr.testing.nmr_fat_1"), num("nmr.testing.nmr_fat_2")
    a["fat_change"] = a["fat_post"] - a["fat_pre"]
    a["lean_change"] = num("nmr.testing.nmr_lean_2") - num("nmr.testing.nmr_lean_1")
    a["weight_change"] = num("nmr.testing.nmr_weight_2") - num("nmr.testing.nmr_weight_1")
    checks = []
    for mine, theirs in [("vo2max_change", "calculated.variables.vo2_max_change"),
                         ("fat_change", "calculated.variables.pct_body_fat_change"),
                         ("lean_change", "calculated.variables.pct_body_lean_change")]:
        diff = float((a[mine] - num(theirs)).abs().max())
        assert not (diff > 1e-6), f"{mine} disagrees with {theirs} (max |diff| {diff})"
        checks.append(f"`{mine}` equals `{theirs}` (max |difference| {diff:.1e})")
    post = a.dropna(subset=["days_vo2_post"])
    assert (post["days_vo2_pre"] < post["days_vo2_post"]).all()
    tr = a.dropna(subset=["days_train_start", "days_vo2_post"])
    assert ((tr["days_vo2_pre"] < tr["days_train_start"]) & (tr["days_train_start"] < tr["days_vo2_post"])).all()
    checks.append(f"`days_vo2_1 < training.day1_days < days_vo2_2` for all {len(tr)} animals with a post test "
                  f"(controls, 4w, 8w; `training.day1_days` is the start of the controls' handling period too)")
    for g in ("1w", "2w"):
        n_post = int(a.loc[a["group"] == g, ["vo2max_post", "fat_post"]].notna().any(axis=1).sum())
        assert n_post == 0, f"{g} animals have post-training measurements ({n_post})"
    checks.append("no 1w or 2w animal has a post-training VO2max or NMR value (0 of 30 in each group)")
    a.attrs["checks"] = checks
    return a


def design_tables(a: pd.DataFrame, out) -> tuple[pd.DataFrame, pd.DataFrame]:
    tab = a.groupby(["group", "sex"]).agg(
        n_animals=("pid", "size"), arrived=("d_arrive", lambda s: str(s.min().date())),
        sac_first=("d_sac", lambda s: str(s.min().date())), sac_last=("d_sac", lambda s: str(s.max().date())),
        sac_days=("d_sac", "nunique"), sac_span_days=("d_sac", lambda s: (s.max() - s.min()).days),
        death_time_median_h=("death_hour", "median"),
        death_time_range_h=("death_hour", lambda s: f"{s.min():.1f}-{s.max():.1f}")).reindex(C.GROUP_ORDER, level=0)
    tab["days_arrival_to_first_sacrifice"] = (pd.to_datetime(tab["sac_first"]) - pd.to_datetime(tab["arrived"])).dt.days
    tab = tab.reset_index()[["group", "sex", "n_animals", "arrived", "sac_first", "sac_last", "sac_days", "sac_span_days",
                             "days_arrival_to_first_sacrifice", "death_time_median_h", "death_time_range_h"]]
    tab.to_csv(out / "design_by_group_sex.csv", index=False)

    # contrasts against the controls: median dates, per sex and pooled over sex
    def med_day(s: pd.Series) -> float:
        return float(np.median(s.map(pd.Timestamp.toordinal).to_numpy()))

    rows = []
    for dur in DURATIONS:
        for sex in ("female", "male", "pooled"):
            d = a[(a["group"] == dur) & ((a["sex"] == sex) | (sex == "pooled"))]
            c = a[(a["group"] == "control") & ((a["sex"] == sex) | (sex == "pooled"))]
            arr_gap = med_day(d["d_arrive"]) - med_day(c["d_arrive"])
            rows.append({"duration": dur, "sex": sex, "n_animals": len(d), "n_control": len(c),
                         "sacrifice_gap_days": med_day(d["d_sac"]) - med_day(c["d_sac"]), "arrival_gap_days": arr_gap,
                         "same_cohort": bool(set(d["d_arrive"].dt.date) == set(c["d_arrive"].dt.date)),
                         "sac_median": str(pd.Timestamp.fromordinal(int(round(med_day(d["d_sac"])))).date()),
                         "sac_median_control": str(pd.Timestamp.fromordinal(int(round(med_day(c["d_sac"])))).date()),
                         "arrival_dates": ";".join(sorted(str(x) for x in set(d["d_arrive"].dt.date))),
                         "arrival_dates_control": ";".join(sorted(str(x) for x in set(c["d_arrive"].dt.date)))})
    contrasts = pd.DataFrame(rows)
    contrasts.to_csv(out / "design_contrasts.csv", index=False)

    cw = a[a["group"].isin(["control", "8w"])]
    pd.crosstab([cw["sex"], cw["d_sac"].dt.date], cw["group"]).to_csv(out / "control_8w_by_sacrifice_day.csv")
    tod = []
    for s in ("male", "female"):
        d = cw[cw["sex"] == s]
        c, e = d[d["group"] == "control"]["death_hour"].dropna(), d[d["group"] == "8w"]["death_hour"].dropna()
        tod.append({"sex": s, "n_control": len(c), "n_8w": len(e), "death_time_median_control_h": c.median(),
                    "death_time_median_8w_h": e.median(), "mannwhitney_p": stats.mannwhitneyu(c, e).pvalue})
    pd.DataFrame(tod).to_csv(out / "control_8w_time_of_death.csv", index=False)
    return tab, contrasts


def physiology_tables(a: pd.DataFrame, out) -> tuple[pd.DataFrame, pd.DataFrame]:
    per = a[["pid", "bid", "sex", "group", "days_vo2_pre", "days_vo2_post", "days_train_start", "days_nmr_pre",
             "days_nmr_post", "vo2max_pre", "vo2max_post", "vo2max_change", "fat_pre", "fat_post", "fat_change",
             "lean_change", "weight_change"]]
    per.to_csv(out / "physiology_per_animal.csv", index=False)
    rows = []
    for var in PHYS_VARS:
        for grp in C.GROUP_ORDER:
            for sex in ["all", "female", "male"]:
                s = per[(per["group"] == grp) & ((per["sex"] == sex) | (sex == "all"))][var]
                x = s.dropna()
                rows.append({"variable": var, "group": grp, "sex": sex, "n_animals": int(len(x)), "n_in_group": int(len(s)),
                             "median": x.median() if len(x) else np.nan, "q25": x.quantile(0.25) if len(x) else np.nan,
                             "q75": x.quantile(0.75) if len(x) else np.nan})
    by_group = pd.DataFrame(rows)
    by_group.to_csv(out / "physiology_by_group.csv", index=False)
    tests = []
    for var in TEST_VARS:
        for grp, matched in TEST_CONTRASTS:
            for sex in ["female", "male"]:
                t = per[(per["group"] == grp) & (per["sex"] == sex)][var].dropna()
                c = per[(per["group"] == "control") & (per["sex"] == sex)][var].dropna()
                u = stats.mannwhitneyu(t, c, alternative="two-sided")
                tests.append({"variable": var, "contrast": f"{grp} vs control", "sex": sex, "n_trained": len(t),
                              "n_control": len(c), "median_trained": t.median(), "median_control": c.median(),
                              "median_diff": t.median() - c.median(), "mw_p": u.pvalue, "cohort_matched": matched})
    tests = pd.DataFrame(tests)
    tests.to_csv(out / "physiology_group_tests.csv", index=False)
    return by_group, tests


def notes(a: pd.DataFrame, tab: pd.DataFrame, contrasts: pd.DataFrame, tests: pd.DataFrame, out_design, out_phys) -> str:
    key = tests[tests["variable"].isin(["vo2max_change", "fat_change"])]
    body = f"""# Phase 15a: study-design dates and physiology

Script: `scripts/15_time_course_design.py` (PHENO only: `data/raw/pheno.csv` via `io.load_pheno()`, c1.0 dotted
column names). Rerun with `make time-course`. [measured] = read from the CSVs here; [interpreted] = inference.

## Identities checked [measured]

- {len(a)} animals (one `pid` each); the group parsed by `io.load_pheno()` equals `key.anirandgroup` for every animal.
- Index `_1` is the pre-training test and `_2` the post-training test:
""" + "".join(f"  - {c}.\n" for c in a.attrs["checks"]) + f"""
## `design/design_by_group_sex.csv`

One row per group × sex. `n_animals`; `arrived` (arrival date, `key.d_arrive`); `sac_first`, `sac_last` (first and
last sacrifice date, `key.d_sacrifice`); `sac_days` (number of distinct sacrifice days); `sac_span_days`
(`sac_last − sac_first`, days); `days_arrival_to_first_sacrifice`; `death_time_median_h`, `death_time_range_h`
(time of death, `specimen.collection.t_death`, hours after midnight).

{report.df_to_md(tab, floatfmt=".2f")}

## `design/design_contrasts.csv`

One row per duration × sex, plus a `pooled` row per duration (both sexes together). `sacrifice_gap_days` = median
sacrifice date of the group minus the controls' median (same sex, or pooled), in days; `arrival_gap_days` likewise
for the arrival date; `same_cohort` = the group's arrival dates are exactly the controls' arrival dates;
`sac_median`, `sac_median_control` (the medians, rounded to a day); `arrival_dates`, `arrival_dates_control`.

{report.df_to_md(contrasts[["duration", "sex", "n_animals", "n_control", "sacrifice_gap_days", "arrival_gap_days", "same_cohort"]], floatfmt=".1f")}

Only 8w shares the controls' arrival cohort and sacrifice dates; within every group the two sexes arrived about a
month apart [measured]. Any 1w/2w/4w-vs-control contrast is also a contrast of cohort, season and date [interpreted].

## `design/control_8w_by_sacrifice_day.csv`, `design/control_8w_time_of_death.csv`

Animals per sacrifice day (sex × date × group) and the time of death per sex (`n_control`, `n_8w`, medians in hours,
two-sided Mann–Whitney `mannwhitney_p`).

## `physiology/physiology_per_animal.csv`

One row per animal: `pid`, `bid`, `sex`, `group`; `days_vo2_pre/post`, `days_train_start`, `days_nmr_pre/post`
(study day of each test and of the first training day); `vo2max_pre/post/change` (`vo2.max.test.vo2_max_1/2`);
`fat_pre/post/change` (`nmr.testing.nmr_fat_1/2`, % body fat, change in percentage points); `lean_change`,
`weight_change` (NMR lean %, NMR weight in g). 1w and 2w animals have no post values.

## `physiology/physiology_by_group.csv`

Long table: `variable`, `group`, `sex` (`all`, `female`, `male`), `n_animals` (with a value), `n_in_group`, `median`,
`q25`, `q75`.

## `physiology/physiology_group_tests.csv`

`variable`, `contrast` (`8w vs control`, `4w vs control`), `sex`, `n_trained`, `n_control`, `median_trained`,
`median_control`, `median_diff` (trained minus control), `mw_p` (two-sided Mann–Whitney on per-animal values),
`cohort_matched` (True only for 8w). Variables: the four changes and the two baselines.

{report.df_to_md(key, floatfmt=".4g")}

Training worked physiologically in the cohort-matched contrast [measured]. The 4w row is not a clean dose step:
different arrival cohort, tests and sacrifice in Nov–Dec 2018 against Sep–Oct for the controls, and a shorter
pre-to-post interval [interpreted].
"""
    (out_design / "NOTES.md").write_text(body)
    return body


def main() -> None:
    ap = cli.common_parser("Study-design dates and physiology of the training groups")
    args = ap.parse_args()
    cli.banner("15_time_course_design", args)
    out_design = cli.outdir("15_time_course/design", args.out)
    out_phys = cli.outdir("15_time_course/physiology", None if args.out is None else str(out_design.parent / "physiology"))
    a = animals()
    tab, contrasts = design_tables(a, out_design)
    by_group, tests = physiology_tables(a, out_phys)
    print(tab.to_string(index=False))
    print(contrasts[["duration", "sex", "n_animals", "sacrifice_gap_days", "arrival_gap_days", "same_cohort"]].to_string(index=False))
    print(tests[tests["variable"].isin(["vo2max_change", "fat_change"])].to_string(index=False))
    body = notes(a, tab, contrasts, tests, out_design, out_phys)
    report.add_section("15a · Study-design dates and physiology", body.split("\n", 1)[1],
                       params={"n_animals": len(a)})
    # the consortium's outlier flags (metadata): the vena-cava vials marked brown-fat contaminated, for the confusion note
    ol = io.load_outliers()
    fl = ol[(ol["tissue"].astype(str) == "VENACV") & ol["reason"].astype(str).str.contains("BAT contamination")].copy()
    fl = fl[["viallabel", "pid", "group", "reason"]].drop_duplicates("viallabel")
    fl.insert(2, "sex", fl["group"].astype(str).str.split("_").str[0])
    fl.to_csv(out_design / "flagged_vials.csv", index=False)
    print(f"flagged vena-cava vials (brown-fat contamination): {len(fl)} vials, {fl['pid'].nunique()} animals")
    print(f"wrote {out_design} and {out_phys}")


if __name__ == "__main__":
    main()
