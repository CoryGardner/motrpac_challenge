"""Part 1 (design) and the proteomics plex/channel balance check (part 2 addendum).
Run from code/pipeline:  python investigations/time_course/tc1_design_and_plex.py
Writes results/15_time_course/1_design/*.csv and 2_3_gradients/prot_{plex_by_duration,channel_vs_group}.csv."""
from pathlib import Path
import pandas as pd
from scipy import stats

R = Path("results/15_time_course")
(R / "1_design").mkdir(parents=True, exist_ok=True)
G = {"Eight-week program Control Group": "control", "One-week program": "1w", "Two-week program": "2w",
     "Four-week program": "4w", "Eight-week program Training Group": "8w"}
p = pd.read_csv("data/raw/pheno.csv", dtype=str)
a = p.drop_duplicates("pid").copy()
a["group"] = a["key.anirandgroup"].map(G)
a["d_arrive"] = pd.to_datetime(a["key.d_arrive"], format="%d%b%Y")
a["d_sac"] = pd.to_datetime(a["key.d_sacrifice"], format="%d%b%Y")
a["death_hour"] = pd.to_timedelta(a["specimen.collection.t_death"], errors="coerce").dt.total_seconds() / 3600
order = ["control", "1w", "2w", "4w", "8w"]
tab = a.groupby(["group", "sex"]).agg(
    n_animals=("pid", "size"), arrived=("d_arrive", lambda s: str(s.min().date())),
    sac_first=("d_sac", lambda s: str(s.min().date())), sac_last=("d_sac", lambda s: str(s.max().date())),
    sac_days=("d_sac", "nunique"), death_time_median_h=("death_hour", "median"),
    death_time_range_h=("death_hour", lambda s: f"{s.min():.1f}-{s.max():.1f}")).reindex(order, level=0)
tab["days_arrival_to_first_sacrifice"] = (pd.to_datetime(tab.sac_first) - pd.to_datetime(tab.arrived)).dt.days
tab.to_csv(R / "1_design" / "design_by_group_sex.csv")
cw = a[a.group.isin(["control", "8w"])]
pd.crosstab([cw.sex, cw.d_sac.dt.date], cw.group).to_csv(R / "1_design" / "control_8w_by_sacrifice_day.csv")
rows = []
for s in ("male", "female"):
    d = cw[cw.sex == s]
    c, e = d[d.group == "control"].death_hour.dropna(), d[d.group == "8w"].death_hour.dropna()
    rows.append({"sex": s, "n_control": len(c), "n_8w": len(e), "death_time_median_control_h": c.median(),
                 "death_time_median_8w_h": e.median(), "mannwhitney_p": stats.mannwhitneyu(c, e).pvalue})
pd.DataFrame(rows).to_csv(R / "1_design" / "control_8w_time_of_death.csv", index=False)
print(tab.to_string())

# proteomics: is TMT plex or channel associated with training group?
m = pd.read_csv("data/raw/meta/PROT.csv", dtype=str)
m = m[m.viallabel.str.startswith("9")].merge(p[["viallabel", "key.anirandgroup"]].drop_duplicates("viallabel"), on="viallabel")
m["group"] = m["key.anirandgroup"].map(G)
plex, chan = [], []
for t, d in m.groupby("tissue"):
    ct = pd.crosstab(d.tmt_plex, d.group)
    for dur in ("1w", "2w", "4w", "8w"):
        sub = ct[["control", dur]]
        sub = sub[sub.sum(axis=1) > 0]
        plex.append({"tissue": t, "duration": dur, "n_plexes": len(sub),
                     "plexes_with_both": int(((sub["control"] > 0) & (sub[dur] > 0)).sum()),
                     "chi2_p_plex_vs_group": stats.chi2_contingency(sub)[1] if len(sub) > 1 else float("nan")})
    chan.append({"tissue": t, "chi2_p_channel_vs_group": stats.chi2_contingency(pd.crosstab(d.tmt11_channel, d.group))[1],
                 "channels": d.tmt11_channel.nunique()})
pd.DataFrame(plex).to_csv(R / "2_3_gradients" / "prot_plex_by_duration.csv", index=False)
pd.DataFrame(chan).to_csv(R / "2_3_gradients" / "prot_channel_vs_group.csv", index=False)
print(pd.DataFrame(chan).round(3).to_string(index=False))
