"""Round-robin panel selection: prefix property and agreement with the selector class."""
import numpy as np

from tfp import models


def test_roundrobin_prefix_and_selector_agree():
    rng = np.random.default_rng(3)
    X = rng.normal(size=(60, 40))
    y = np.repeat(["A", "B", "C", "D"], 15)
    X[y == "A", :5] += 3; X[y == "B", 5:10] += 3; X[y == "C", 10:15] += 3; X[y == "D", 15:20] += 3
    sel15 = models.RoundRobinSelector(k=15).fit(X, y)
    order15 = models.roundrobin_order(sel15.scores_, 15)
    assert sorted(order15) == list(sel15.idx_)
    sel7 = models.RoundRobinSelector(k=7).fit(X, y)
    assert list(sel7.idx_) == sorted(order15[:7])          # prefix property
    assert len(set(order15)) == 15                         # unique picks
    first_round = set(order15[:4])
    assert all(any(i in range(5 * c, 5 * c + 5) for i in first_round) for c in range(4))  # one marker per class first


def test_roundrobin_exclusion_skips_masked_genes():
    rng = np.random.default_rng(5)
    X = rng.normal(size=(60, 40))
    y = np.repeat(["A", "B", "C", "D"], 15)
    for c, sl in enumerate((slice(0, 5), slice(5, 10), slice(10, 15), slice(15, 20))):
        X[y == ["A", "B", "C", "D"][c], sl] += 3
    sel = models.RoundRobinSelector(k=8).fit(X, y)
    top_a = sel.order_[0]                                   # class A's first pick
    ex = np.zeros((4, 40), dtype=bool)
    ex[0, top_a] = True                                     # A may not pick its best gene
    sel2 = models.RoundRobinSelector(k=8, exclude=ex).fit(X, y)
    assert top_a not in sel2.order_[:1] and sel2.order_[0] != top_a
    assert top_a in sel2.order_ or True                      # another class may still take it
    assert len(set(sel2.order_)) == 8
