# %% [markdown]
# ## 2. PCA and variance partition: what the leading axes of each assay are
#
# **Question.** When all tissues of one assay are stacked into one matrix, what do the first principal
# components track — tissue, sex, training group, or the individual animal?
#
# **Why it matters.** A tissue-identity classifier (sections 4 onward) is only meaningful if tissue is a
# real axis of the data. If the leading axes are tissue, the task is easy and the interesting question
# becomes *how few features* suffice; if they are not, a high accuracy has to come from somewhere else.
#
# **What to look for.** Two different quantities, printed as separate columns:
# - `explained` — the fraction of the *total* (scaled) variance that the PC carries;
# - `R2_<factor>` — the fraction of *that PC's own* variance explained by the group means of a factor.
#
# A PC can have R² near 1 for tissue while carrying only a modest share of the total variance; the two
# should not be multiplied or confused. `R2_pid` is the animal effect; animals nest sex × group, so it
# is an upper bound on anything animal-level. Compare the proteomics rows with the other two assays.
#
# Settings are the pipeline defaults of `scripts/03_eda.py` (`make eda`): outer feature join keeping
# features present in ≥ 80 % of samples, TRNSCRPT from raw counts (log2 CPM, stacked gene filter),
# median imputation, top 5,000 features by variance, standardized, 10 PCs.

# %%
section("2 PCA and variance partition")
from sklearn.decomposition import PCA
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler


# from scripts/03_eda.py::r2_categorical (verbatim)
def r2_categorical(scores: np.ndarray, cat: pd.Series) -> np.ndarray:
    """Fraction of variance of each column of `scores` explained by group means of `cat`."""
    cat = pd.Series(cat).astype(str).to_numpy()
    out = []
    for j in range(scores.shape[1]):
        s = scores[:, j]
        tot = np.var(s) * len(s)
        if tot == 0:
            out.append(0.0)
            continue
        between = sum(len(s[cat == c]) * (s[cat == c].mean() - s.mean()) ** 2 for c in np.unique(cat))
        out.append(float(between / tot))
    return np.array(out)


# from scripts/03_eda.py::pca_of (verbatim)
def pca_of(om: io.OmicsMatrix, n_components: int = 10, top_var: int | None = 5000, seed: int = C.SEED):
    X = om.X.to_numpy(dtype=float)
    X = SimpleImputer(strategy="median").fit_transform(X)
    if top_var and X.shape[1] > top_var:
        v = X.var(axis=0)
        X = X[:, np.argsort(v)[::-1][:top_var]]
    X = StandardScaler().fit_transform(X)
    n = min(n_components, X.shape[0] - 1, X.shape[1])
    pca = PCA(n_components=n, random_state=seed)
    S = pca.fit_transform(X)
    return S, pca.explained_variance_ratio_


# from src/tfp/cli.py::resolve_source (verbatim; cli.py is not part of the pasted library)
def resolve_source(assay: str, source: str) -> str:
    if source != "auto":
        return source
    if assay == "TRNSCRPT" and any((C.COUNTS_DIR).glob(f"{assay}__*.csv")):
        return "counts"
    return "norm"

# %%
# ---- stacked PCA per assay, outer join (03_eda default) and inner join (03_eda --join inner) ----
# from scripts/03_eda.py::main (trimmed: only the stacked PCA + variance-partition block; plots done below)
S02_OUT = OUT / "s02"
S02_OUT.mkdir(parents=True, exist_ok=True)
s02_scores, s02_rows, s02_oms = {}, [], {}
for assay in ("TRNSCRPT", "PROT", "METAB"):
    source = resolve_source(assay, "auto")
    for join in ("outer", "inner"):
        if source == "counts" and join == "inner":  # same matrix as the outer pass (see below): reuse its PCs
            s02_rows.append(s02_rows[-1].assign(join="inner"))
            continue
        if source == "counts":
            # the counts path ignores `join` (stacked gene filter, absent = 0 counts), so the matrix is the
            # same one section 4 uses — share it through the cache
            om = cached("TRNSCRPT_counts", lambda: io.stack_tissues(
                "TRNSCRPT", tissues=None, source="counts", join="inner", pheno=pheno,
                complete=False, drop_incomplete_samples=None))
        else:
            om = io.stack_tissues(assay, tissues=None, source=source, join=join, min_present=0.8,
                                  pheno=pheno, verbose=False)
        print(f"{assay:8s} join={join:5s}: {om.notes[-1] if source == 'norm' else om.notes[0][:120]}")
        S, ev = pca_of(om, seed=C.SEED)
        vp = pd.DataFrame({"assay": assay, "join": join, "PC": [f"PC{i+1}" for i in range(S.shape[1])],
                           "explained": ev[: S.shape[1]]})
        for col in ("tissue", "sex", "group", "pid"):
            if om.meta[col].nunique() > 1:
                vp[f"R2_{col}"] = r2_categorical(S, om.meta[col])
        s02_rows.append(vp)
        if join == "outer":
            s02_scores[assay] = (S, ev, om.meta)
            s02_oms[assay] = om            # notebook addition: kept for the UMAP in 2b
s02_vp = pd.concat(s02_rows, ignore_index=True)
s02_vp.to_csv(S02_OUT / "variance_partition_all.csv", index=False)

# %%
# ---- the headline table: PCs 1–3 per assay (outer join), recomputed next to the published values ----
s02_pub = pd.concat([pd.read_csv(res("03_eda", f"variance_partition_{a}.csv")).assign(assay=a)
                     for a in ("TRNSCRPT", "PROT", "METAB")])
s02_pub_inner = pd.concat([pd.read_csv(res("03_eda_inner", f"variance_partition_{a}.csv")).assign(assay=a)
                           for a in ("TRNSCRPT", "PROT", "METAB")])
s02_top = s02_vp[(s02_vp["join"] == "outer") & s02_vp["PC"].isin(["PC1", "PC2", "PC3"])].drop(columns="join")
s02_top = s02_top.merge(s02_vp[s02_vp["join"] == "inner"][["assay", "PC", "R2_tissue"]]
                        .rename(columns={"R2_tissue": "R2_tissue_inner_join"}), on=["assay", "PC"])
cols = ["explained", "R2_tissue", "R2_sex", "R2_group", "R2_pid"]
s02_cmp = s02_top.merge(s02_pub[["assay", "PC"] + cols], on=["assay", "PC"], suffixes=("", "_pub"))
print("PCs 1–3, stacked tissues. `explained` = share of total variance; `R2_*` = share of that PC's variance.")
display(s02_top.set_index(["assay", "PC"]).round(4))
s02_dev = max(float((s02_cmp[c] - s02_cmp[c + "_pub"]).abs().max()) for c in cols)
s02_inner_pub = s02_top.merge(s02_pub_inner[["assay", "PC", "R2_tissue"]], on=["assay", "PC"], suffixes=("", "_pub"))
s02_dev_inner = float((s02_inner_pub["R2_tissue_inner_join"] - s02_inner_pub["R2_tissue_pub"]).abs().max())
print(f"max |recomputed - published| over these cells: outer join {s02_dev:.2e}; inner-join R2_tissue {s02_dev_inner:.2e}")

for a in ("TRNSCRPT", "PROT", "METAB"):
    r = s02_top[(s02_top["assay"] == a) & (s02_top["PC"] == "PC1")].iloc[0]
    p = s02_pub[(s02_pub["assay"] == a) & (s02_pub["PC"] == "PC1")].iloc[0]
    print(f"{a:8s} PC1: explained {r['explained']:.4f} (published {p['explained']:.4f})   "
          f"R2_tissue {r['R2_tissue']:.4f} (published {p['R2_tissue']:.4f})   R2_sex {r['R2_sex']:.4f} (published {p['R2_sex']:.4f})")
    record(f"s02.{a}_PC1_explained", r["explained"], "02")
    record(f"s02.{a}_PC1_R2_tissue", r["R2_tissue"], "02")
    record(f"s02.{a}_PC1_R2_sex", r["R2_sex"], "02")

# %%
# ---- PC1 vs PC2 per assay, coloured by tissue (outer join; same scores as the table) ----
fig, axes = plt.subplots(1, 3, figsize=(16, 5))
cmap = plt.get_cmap("tab20")
for ax, (a, (S, ev, meta)) in zip(axes, s02_scores.items()):
    tissues = sorted(meta["tissue"].unique())
    for i, t in enumerate(tissues):
        m = (meta["tissue"] == t).to_numpy()
        ax.scatter(S[m, 0], S[m, 1], s=8, color=cmap(i % 20), label=t, alpha=0.8)
    ax.set_title(a)
    ax.set_xlabel(f"PC1 ({ev[0]:.1%} of variance)")
    ax.set_ylabel(f"PC2 ({ev[1]:.1%})")
axes[0].legend(fontsize=6, ncol=2, frameon=False, markerscale=1.5)
axes[1].legend(fontsize=7, frameon=False, markerscale=1.5)
fig.suptitle("Stacked-tissue PCA, coloured by tissue")
plt.tight_layout()
plt.show()

# PROT again, coloured by sex: the axis PC1 does track
S, ev, meta = s02_scores["PROT"]
fig, ax = plt.subplots(figsize=(5.5, 4.5))
for sex, c in (("female", "tab:red"), ("male", "tab:blue")):
    m = (meta["sex"] == sex).to_numpy()
    ax.scatter(S[m, 0], S[m, 1], s=8, color=c, label=sex, alpha=0.7)
ax.set_xlabel(f"PC1 ({ev[0]:.1%})")
ax.set_ylabel(f"PC2 ({ev[1]:.1%})")
ax.set_title("PROT, coloured by sex")
ax.legend(frameon=False)
plt.tight_layout()
plt.show()

# %% [markdown]
# **What this shows.** For transcripts and metabolites, the first three PCs are almost entirely tissue:
# R² of tissue on each PC is close to 1, and the scatter plots separate into tissue clusters. Proteomics
# is the exception: its PC1 has essentially no tissue component, and sex and the individual animal
# explain more of it than tissue does; tissue appears only weakly on a later PC. The inner-join column
# shows the same picture, so the choice of feature join does not drive it. This is what the
# normalization of `PROT_*_NORM_DATA` predicts: values are log2 ratios to a *per-tissue* reference pool,
# median-centred per sample, so every tissue sits near zero and between-tissue abundance differences are
# removed by design.
#
# **What it does not show.** R² of a PC is not a share of the total variance: a tissue R² near 1 on a PC
# that carries a minority of the total variance does not mean tissue explains the data "almost entirely" — read
# `explained` and `R2_tissue` separately. This is descriptive, on all samples, and nothing here selects
# features for later sections. It also does not say whether training (the `group` factor) is visible
# *within* a tissue: stacked across tissues, group R² is negligible for every assay, which is expected
# when tissue dominates, and says nothing about within-tissue effects.


# %% [markdown]
# ### 2b. UMAP: the same matrices, viewed non-linearly
#
# **Question.** PCA shows only linear axes, and only a few of them at a time. Does a non-linear
# embedding of the same data show the same grouping: tissue for transcripts and metabolites, and
# something other than tissue for proteomics?
#
# **Why it matters.** Several tissues that PCA stacks on top of each other may still be separable in the
# full space. That is the question the classifiers in sections 4–5 answer formally. UMAP is a quick visual
# check of which tissues sit next to each other: the muscle pair, the brain regions, vena cava and fat.
#
# **How it is computed.** For each assay: the same matrix as the PCA above, with the same preprocessing
# (median imputation, 5,000 highest-variance features, z-scoring), then its first 50 principal
# components, then UMAP (`umap-learn`, 15 neighbours, minimum distance 0.1, fixed seed). The
# first run in a kernel compiles UMAP's numba code, which takes a few extra seconds.
#
# **What to look for.** Whether each tissue forms its own island, and which tissues share one. The table
# makes that quantitative: the share of each sample's 15 nearest neighbours that come from the same
# tissue, in the UMAP map and in the 50-PC space it was built from.

# %%
# notebook addition (not in the pipeline): UMAP of the section-2 matrices; umap-learn is optional → import-guarded
try:
    import umap
except ImportError:
    umap = None
    print("umap-learn is not installed in this kernel: 2b is skipped (conda install -c conda-forge umap-learn)")

from sklearn.neighbors import NearestNeighbors


def s02_neighbour_purity(Z: np.ndarray, labels: pd.Series, k: int = 15) -> float:
    """Mean share of each sample's k nearest neighbours (itself excluded) that carry the same label."""
    idx = NearestNeighbors(n_neighbors=k + 1).fit(Z).kneighbors(Z, return_distance=False)[:, 1:]
    lab = labels.to_numpy()
    return float((lab[idx] == lab[:, None]).mean())


s02_umap, s02_purity = {}, []
if umap is not None:
    for a, om in s02_oms.items():
        P50, _ = pca_of(om, n_components=50, seed=C.SEED)          # same preprocessing as the PCA above
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")                        # "n_jobs overridden by random_state"
            U = umap.UMAP(n_neighbors=15, min_dist=0.1, random_state=C.SEED).fit_transform(P50)
        s02_umap[a] = (U, om.meta)
        for col in ("tissue", "sex"):
            s02_purity.append({"assay": a, "label": col, "n_samples": len(U),
                               "same-label share of 15 neighbours, 50-PC space": s02_neighbour_purity(P50, om.meta[col]),
                               "same-label share of 15 neighbours, UMAP": s02_neighbour_purity(U, om.meta[col]),
                               "chance (share of the commonest label)": float(om.meta[col].value_counts(normalize=True).iloc[0])})
    s02_purity = pd.DataFrame(s02_purity)
    s02_purity.to_csv(S02_OUT / "umap_neighbour_purity.csv", index=False)
    display(s02_purity.round(3))
    for r in s02_purity.itertuples(index=False):
        record(f"s02.umap_purity_{r.assay}_{r.label}", r[4], "02", note="UMAP neighbour purity")
        record(f"s02.pc50_purity_{r.assay}_{r.label}", r[3], "02", note="50-PC neighbour purity")

# %%
# UMAP per assay, coloured by tissue (same colours as the PCA figure), and PROT again coloured by sex
if s02_umap:
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    cmap = plt.get_cmap("tab20")
    for ax, (a, (U, meta)) in zip(axes, s02_umap.items()):
        for i, t in enumerate(sorted(meta["tissue"].unique())):
            m = (meta["tissue"] == t).to_numpy()
            ax.scatter(U[m, 0], U[m, 1], s=8, color=cmap(i % 20), label=t, alpha=0.8)
        ax.set_title(a)
        ax.set_xlabel("UMAP 1")
        ax.set_ylabel("UMAP 2")
        ax.set_xticks([]); ax.set_yticks([])
    axes[0].legend(fontsize=6, ncol=2, frameon=False, markerscale=1.5)
    axes[1].legend(fontsize=7, frameon=False, markerscale=1.5)
    fig.suptitle("UMAP of the first 50 PCs, coloured by tissue (15 neighbours, min_dist 0.1)")
    plt.tight_layout()
    plt.show()

    U, meta = s02_umap["PROT"]
    fig, ax = plt.subplots(figsize=(5.5, 4.5))
    for sex, c in (("female", "tab:red"), ("male", "tab:blue")):
        m = (meta["sex"] == sex).to_numpy()
        ax.scatter(U[m, 0], U[m, 1], s=8, color=c, label=sex, alpha=0.7)
    ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2"); ax.set_xticks([]); ax.set_yticks([])
    ax.set_title("PROT UMAP, coloured by sex")
    ax.legend(frameon=False)
    plt.tight_layout()
    plt.show()

# %% [markdown]
# **What this shows.** The neighbour table gives the reading in numbers. For each assay and label it
# compares the share of a sample's nearest neighbours that share its tissue (or sex) in the UMAP map
# and in the 50-PC space the map was built from, against the share expected from the commonest label
# alone. Where the tissue share is near 1 in both spaces, tissues form their own islands. The figure
# then shows which tissues, if any, share an island, and whether a tissue splits into two nearby islands
# (compare the sex rows of the table). For proteomics, a tissue grouping that PC1 did not show is **not**
# evidence of tissue biology. This outer-join matrix median-imputes proteins that are missing from
# whole tissues, and section 3 shows that the missingness pattern alone identifies the tissue. The
# per-tissue reference normalization has already removed abundance differences between tissues.
#
# **What it does not show.** UMAP keeps local neighbourhoods, not global geometry. Distances *between*
# islands, their sizes and their positions mean little, and they change with the number of neighbours
# and the seed. Only who-is-next-to-whom is interpretable. Like the PCA, this is descriptive, on all
# samples, and selects nothing for later sections. A tissue island is not a classifier's accuracy:
# sections 4–5 measure that on held-out animals.
