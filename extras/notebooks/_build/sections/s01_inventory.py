# %% [markdown]
# ## 1. Setup and data inventory
#
# **Question.** What is actually on disk in the pipeline's export (`code/pipeline/data/raw/`): how many
# animals, vials, tissues and assays; how the animals split by sex × training group; which tissues have
# transcriptomics, proteomics *and* metabolomics (the "fusion" tissues); and whether every sample column
# can be tied to an animal.
#
# **Why it matters.** Every later split is by animal (`pid`), and every cross-assay join is by `bid`. If
# `bid` and `pid` were not one-to-one, or if sample columns silently failed to match the phenotype table,
# grouped cross-validation would leak and joins would be wrong.
#
# **What to look for.** The design table (animals per sex × group); the tissue × assay table (which
# tissues are single-assay); that `bid ↔ pid` is one-to-one; and — a detail the pipeline's own inventory
# hides — how many sample columns were matched to the phenotype table only through the `bid` fallback in
# `io.sample_meta_for` rather than by their exact `viallabel`. The pipeline's `n_unmatched` column counts
# what is still unmatched *after* that fallback, so it is zero by construction.

# %%
section("1 Setup and data inventory")
# ---- the phenotype table: one row per vial; standardized pid / bid / sex / group columns ----
pheno = cached("pheno", io.load_pheno)
s01_manifest = io.load_manifest()
print(f"data package: {s01_manifest.get('package', '?')}  version {s01_manifest.get('version', '?')}")

s01_animals = pheno.drop_duplicates("pid")
print(f"vials in PHENO: {pheno.shape[0]}   animals (unique pid): {s01_animals.shape[0]}   unique bid: {pheno['bid'].nunique()}")
record("s01.n_animals", s01_animals.shape[0], "01")
record("s01.n_vials", pheno.shape[0], "01")

# ---- bid <-> pid must be one-to-one: every bid maps to one pid and every pid to one bid ----
s01_pid_per_bid = pheno.groupby("bid")["pid"].nunique()
s01_bid_per_pid = pheno.groupby("pid")["bid"].nunique()
print(f"max pids per bid: {s01_pid_per_bid.max()}   max bids per pid: {s01_bid_per_pid.max()}   "
      f"-> one-to-one: {bool(s01_pid_per_bid.max() == 1 and s01_bid_per_pid.max() == 1)}")
# and bid is the first 5 characters of the viallabel (the rule the loader's fallback relies on)
s01_prefix_ok = (pheno.index.str[:5] == pheno["bid"]).mean()
print(f"fraction of vials whose viallabel starts with their bid: {s01_prefix_ok:.4f}")
record("s01.max_pid_per_bid", s01_pid_per_bid.max(), "01")
record("s01.max_bid_per_pid", s01_bid_per_pid.max(), "01")

# ---- animals by sex x training group (from scripts/02_inventory.py::main, verbatim) ----
design = s01_animals.groupby(["sex", "group"]).size().unstack(fill_value=0)
design = design.reindex(columns=[g for g in C.GROUP_ORDER if g in design.columns])
display(design)
s01_design_pub = pd.read_csv(res("02_inventory", "animals_by_sex_group.csv"), index_col=0)
print("identical to results/02_inventory/animals_by_sex_group.csv:",
      bool((design.to_numpy() == s01_design_pub[design.columns].to_numpy()).all()))

# %%
# ---- one row per sample-level table (norm and raw counts) ----
# from scripts/02_inventory.py::main (trimmed: dropped argparse/report, and added the explicit
# exact-match / bid-fallback / unresolved counts that the script's n_unmatched hides)
S01_OUT = OUT / "s01"
S01_OUT.mkdir(parents=True, exist_ok=True)
rows = []
s01_pids = {}  # [notebook] distinct animals per (kind, assay), for the per-assay summary below
for kind in ("norm", "counts"):
    files = io.list_sample_files(kind)
    for _, r in files.iterrows():
        feat, vals = io.read_sample_table(r["path"])
        meta = io.sample_meta_for(vals.columns, pheno, r["tissue"], r["assay"])
        # [notebook] how each sample column was matched: exact viallabel in PHENO, or only via its bid prefix
        vl = pd.Index([str(v) for v in vals.columns])
        exact = vl.isin(pheno.index)
        s01_pids.setdefault((kind, r["assay"]), set()).update(meta["pid"].dropna())
        rows.append({"kind": kind, "assay": r["assay"], "tissue": r["tissue"], "n_features": vals.shape[0],
                     "n_samples": vals.shape[1], "n_animals": int(meta["pid"].nunique()),
                     "n_unmatched": int(meta["pid"].isna().sum()),
                     "n_exact_viallabel": int(exact.sum()),                                   # [notebook]
                     "n_bid_fallback": int((~exact & meta["pid"].notna().to_numpy()).sum()),  # [notebook]
                     "n_female": int((meta["sex"] == "female").sum()), "n_male": int((meta["sex"] == "male").sum()),
                     "groups": ",".join(f"{g}:{n}" for g, n in meta["group"].value_counts().sort_index().items()),
                     "frac_nan": float(vals.isna().to_numpy().mean())})
inv = pd.DataFrame(rows).sort_values(["kind", "assay", "tissue"]).reset_index(drop=True)
inv.to_csv(S01_OUT / "inventory.csv", index=False)
print(f"sample-level tables: {len(inv)}")
record("s01.n_tables", len(inv), "01")

# ---- recomputed vs published (results/02_inventory/inventory.csv) on the shared columns ----
s01_inv_pub = pd.read_csv(res("02_inventory", "inventory.csv"))
s01_cmp = inv.merge(s01_inv_pub, on=["kind", "assay", "tissue"], suffixes=("", "_pub"), how="outer", indicator=True)
for c in ("n_features", "n_samples", "n_animals", "n_unmatched"):
    print(f"  {c:12s} identical in all {len(s01_cmp)} tables: {bool((s01_cmp[c] == s01_cmp[c + '_pub']).all())}")

# ---- the matching detail: n_unmatched is 0 because unmatched vials fall back to their bid ----
print(f"\nn_unmatched (pipeline column, after fallback): {int(inv['n_unmatched'].sum())}")
print(f"sample columns matched by exact viallabel    : {int(inv['n_exact_viallabel'].sum())}")
print(f"sample columns matched only via bid fallback : {int(inv['n_bid_fallback'].sum())}")
record("s01.n_unmatched_after_fallback", inv["n_unmatched"].sum(), "01")
record("s01.n_bid_fallback", inv["n_bid_fallback"].sum(), "01")
s01_fb = inv[inv["n_bid_fallback"] > 0][["kind", "assay", "tissue", "n_samples", "n_exact_viallabel", "n_bid_fallback"]]
display(s01_fb if len(s01_fb) else Markdown("_no table needed the bid fallback_"))

# %%
# ---- per-assay summary of the normalized tables (the Data table of results/SUMMARY.md) ----
# from scripts/10_make_report.py (verbatim aggregation). Note: `animals` there is the MAX over tissues
# of the per-table animal count, not the number of distinct animals; the distinct count is added next to it.
norm = inv[inv["kind"] == "norm"]
s01_by_assay = norm.groupby("assay").agg(tissues=("tissue", "nunique"), samples=("n_samples", "sum"),
                                         animals=("n_animals", "max"), features_median=("n_features", "median"))
s01_by_assay["distinct_animals"] = [len(s01_pids[("norm", a)]) for a in s01_by_assay.index]  # [notebook]
display(s01_by_assay)

# ---- tissue x assay: samples per tissue in each normalized assay (from 02_inventory, verbatim) ----
present = norm.pivot_table(index="tissue", columns="assay", values="n_samples", aggfunc="sum", fill_value=0)
display(present)
show_png(plots.missingness_heatmap(present, "samples per tissue × assay (normalized data)",
                                   S01_OUT / "tissue_by_assay.png"), width=520)

# %%
# ---- the fusion tissues: TRNSCRPT, PROT and METAB all present; animal overlap joined on bid ----
# from scripts/02_inventory.py::main (verbatim)
core = [a for a in ("TRNSCRPT", "PROT", "METAB") if a in present.columns]
fusion_tissues = list(present.index[(present[core] > 0).all(axis=1)]) if core else []
print("fusion tissues:", fusion_tissues)
record("s01.n_fusion_tissues", len(fusion_tissues), "01")
overlap_rows = []
for t in fusion_tissues:
    mats = {a: io.load_norm(a, t, pheno) for a in core}
    bids = {a: set(m.meta["bid"]) for a, m in mats.items()}
    common = set.intersection(*bids.values())
    overlap_rows.append({"tissue": t, **{f"n_{a}": len(b) for a, b in bids.items()}, "n_common_bid": len(common)})
overlap = pd.DataFrame(overlap_rows)
s01_ov_pub = pd.read_csv(res("02_inventory", "fusion_overlap.csv"))
display(overlap.merge(s01_ov_pub[["tissue", "n_common_bid"]], on="tissue", suffixes=("", "_published")))
record("s01.fusion_common_bid_total", overlap["n_common_bid"].sum(), "01")

# %% [markdown]
# **What this shows.** The design table above is the full set of animals in the export, split by sex and
# training group; groups are roughly balanced by sex. `bid` and `pid` are one-to-one, and `bid` is the
# first five characters of every viallabel, so joining assays on `bid` and splitting on `pid` refer to
# the same animals. Transcriptomics and metabolomics cover most tissues; proteomics (and its PTM
# assays) covers a small subset, and the fusion tissues are the tissues where all three core assays
# exist, with nearly all animals shared across the three assays (last table). The matching printout
# separates what the pipeline's `n_unmatched` column merges: sample columns found by their exact
# viallabel, and those found only through the `bid` fallback. The pipeline's inventory reports zero
# unmatched samples; the fallback count above says how much of that zero is the fallback at work.
#
# **What it does not show.** A fallback match is correct at the *animal* level (bid is one-to-one with
# pid) but it cannot confirm that the vial itself is the one described in PHENO. The design is also not
# balanced in *time*: arrival cohort and sacrifice date differ between groups and between sexes
# (`docs/DATA_GUIDE.md`, design confounds), which no count table can reveal. The `animals` column of the
# pipeline's summary is the largest per-tissue count, not the number of distinct animals per assay;
# both are printed above.
