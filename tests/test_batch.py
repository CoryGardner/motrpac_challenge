"""Batch nesting: which tissue pairs are estimable (share a level of every processing variable)?"""
import numpy as np
import pandas as pd
import pytest

from tfp import batch


def _meta():
    """Four tissues. `plate` is nested in tissue except that C and D share P3; `lib` is fully nested;
    `flowcell` is one level for everyone; `noise` is balanced within every tissue (independent of it)."""
    rows = []
    plates = {"A": "P1", "B": "P2", "C": "P3", "D": "P3"}
    for t in "ABCD":
        for i in range(12):
            rows.append({"viallabel": f"9{t}{i:03d}", "tissue": t, "plate": plates[t], "lib": f"L_{t}",
                         "flowcell": "F1", "noise": ["x", "y", "z"][i % 3]})
    return pd.DataFrame(rows)


def test_estimable_pairs_need_a_shared_level_of_every_variable():
    meta = _meta()
    pairs, total = batch.estimable_pairs(meta, "tissue", ["plate", "lib", "flowcell"])
    assert total == 6
    assert pairs == []                     # lib separates every pair
    pairs2, _ = batch.estimable_pairs(meta, "tissue", ["plate", "flowcell"])
    assert pairs2 == [("C", "D")]          # the only pair sharing a plate


def test_cramers_v_is_one_when_nested_and_near_zero_when_independent():
    meta = _meta()
    assert batch.cramers_v(meta["tissue"], meta["lib"]) == pytest.approx(1.0)
    assert batch.cramers_v(meta["tissue"], meta["noise"]) == pytest.approx(0.0, abs=1e-12)
    assert batch.cramers_v(meta["tissue"], meta["flowcell"]) == pytest.approx(0.0)   # one level: no association


def test_nesting_table_counts_levels_and_sharing():
    meta = _meta()
    tab = batch.nesting_table(meta, "tissue", ["plate", "lib", "flowcell"]).set_index("variable")
    assert list(tab.columns) == ["n_levels", "n_tissues", "median_levels_per_tissue", "max_levels_per_tissue",
                                 "max_tissues_per_level", "n_levels_shared", "tissues_in_one_level", "cramers_v",
                                 "n_pairs_sharing_level", "n_pairs_total"]
    p = tab.loc["plate"]
    assert (p["n_levels"], p["n_tissues"], p["max_tissues_per_level"], p["n_levels_shared"]) == (3, 4, 2, 1)
    assert (p["tissues_in_one_level"], p["n_pairs_sharing_level"], p["n_pairs_total"]) == (4, 1, 6)
    assert p["median_levels_per_tissue"] == 1 and p["max_levels_per_tissue"] == 1
    f = tab.loc["flowcell"]
    assert f["n_levels"] == 1 and f["n_pairs_sharing_level"] == 6
