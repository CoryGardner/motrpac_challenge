"""I/O round-trip on a tiny synthetic export (independent of data/raw)."""
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from tfp import config as C


@pytest.fixture(scope="module")
def tiny_raw(tmp_path_factory):
    out = tmp_path_factory.mktemp("raw")
    subprocess.run([sys.executable, str(C.ROOT / "scripts" / "make_synthetic_data.py"), "--out", str(out),
                    "--per-cell", "2", "--n-genes", "200"], check=True, capture_output=True)
    return out


def test_pheno_and_norm(tiny_raw):
    from tfp import io
    ph = io.load_pheno(tiny_raw)
    assert {"pid", "bid", "sex", "group"} <= set(ph.columns)
    assert set(ph["sex"]) == {"female", "male"}
    om = io.load_norm("TRNSCRPT", "SKM-GN", ph, raw_dir=tiny_raw)
    assert om.n_samples == 20 and om.n_features == 200
    assert (om.meta["bid"] == om.X.index.str[:5]).all()
    assert om.meta["pid"].notna().all()


def test_counts_are_log_cpm(tiny_raw):
    from tfp import io
    ph = io.load_pheno(tiny_raw)
    om = io.load_counts("HEART", ph, raw_dir=tiny_raw)
    assert om.X.min().min() >= 0
    assert om.n_features <= 200


def test_stack_and_align(tiny_raw):
    from tfp import io
    ph = io.load_pheno(tiny_raw)
    st = io.stack_tissues("TRNSCRPT", pheno=ph, raw_dir=tiny_raw, verbose=False)
    assert st.meta["tissue"].nunique() == 8 and st.n_samples == 8 * 20
    a = io.load_norm("TRNSCRPT", "HEART", ph, raw_dir=tiny_raw)
    b = io.load_norm("PROT", "HEART", ph, raw_dir=tiny_raw)
    al = io.align_by_animal({"a": a, "b": b}, key="bid")
    assert al["a"].X.index.equals(al["b"].X.index)
    assert al["a"].n_samples == 20


def test_gene_mapping(tiny_raw):
    from tfp import io
    s = io.map_to_gene_symbols(["ENSRNOG00000000001", "NP_000001", "nope"], raw_dir=tiny_raw)
    assert s.iloc[0] == "GENE0001" and isinstance(s.iloc[1], str) and s.isna().iloc[2]


def test_conformal_quantile_and_certificate():
    from tfp import conformal as cp
    rng = np.random.default_rng(0)
    scores = rng.uniform(size=200)
    q = cp.conformal_quantile(scores, 0.1)
    assert 0.85 <= q <= 0.95
    k, table = cp.certify_panel_size({5: (30, 100), 10: (2, 100), 20: (0, 100), 50: (0, 100)}, alpha=0.05, delta=0.05)
    assert k == 20 and table.loc[table["k"] == 10, "certified"].item() is np.False_ or not table.loc[table["k"] == 10, "certified"].item()


def test_stacked_counts_filter_keeps_tissue_restricted_genes(tiny_raw):
    """A gene expressed in one tissue only must survive the stacked filter (absent = 0 counts),
    while the old per-tissue inner join drops it."""
    import pandas as pd
    from tfp import io
    ph = io.load_pheno(tiny_raw)
    gene = "ENSRNOG00000000005"
    for p in sorted((tiny_raw / "counts").glob("TRNSCRPT__*.csv")):
        if p.stem.endswith("HEART"):
            continue
        df = pd.read_csv(p, dtype={"feature": str, "feature_ID": str, "tissue": str, "assay": str})
        df.loc[df["feature_ID"] == gene, df.columns[4:]] = 0
        df.to_csv(p, index=False)
    stacked = io.stack_tissues("TRNSCRPT", source="counts", pheno=ph, raw_dir=tiny_raw, verbose=False)
    per_tissue = io.stack_tissues("TRNSCRPT", source="counts", pheno=ph, raw_dir=tiny_raw, verbose=False,
                                  counts_filter="per_tissue", join="inner")
    assert gene in stacked.X.columns and gene not in per_tissue.X.columns
    assert not stacked.X.isna().any().any()
    assert stacked.n_features >= per_tissue.n_features
    heart = (stacked.meta["tissue"] == "HEART").to_numpy()
    assert stacked.X.loc[~heart, gene].max() == 0.0 and stacked.X.loc[heart, gene].max() > 0
    assert "recovered" in stacked.notes[0]
    with pytest.raises(ValueError):
        io.stack_tissues("TRNSCRPT", source="counts", pheno=ph, raw_dir=tiny_raw, verbose=False, counts_filter="nope")
