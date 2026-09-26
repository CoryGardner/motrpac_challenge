#!/usr/bin/env python3
"""List every number that moved when the conformal quantile was fixed (2026-09-25).

Compares each CSV in the pre-fix snapshot (backup/pipeline_history/results_pre_quantile_fix_2026-09-25/)
with the same file in results/ after the rerun. Rows are matched on their non-numeric columns (plus row
order within duplicates); every numeric cell whose value changed is written out, before → after.

    PYTHONPATH=src python scripts/quantile_fix_diff.py
writes results/QUANTILE_FIX_CHANGES.csv (all changed cells) and prints a per-file summary.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

PIPE = Path(__file__).resolve().parents[1]
SNAP = PIPE.parents[1] / "backup" / "pipeline_history" / "results_pre_quantile_fix_2026-09-25"
RES = PIPE / "results"


def keyed(df: pd.DataFrame) -> pd.DataFrame:
    """Index a table by its non-numeric columns, plus an occurrence counter so duplicate keys stay distinct."""
    keys = [c for c in df.columns if not pd.api.types.is_numeric_dtype(df[c]) or pd.api.types.is_bool_dtype(df[c])]
    if not keys:
        return df.rename_axis("row").assign(_occ=0).set_index("_occ", append=True)
    df = df.copy()
    df["_occ"] = df.groupby(keys, dropna=False).cumcount()
    for k in keys:
        df[k] = df[k].astype(str)
    return df.set_index(keys + ["_occ"])


def main():
    rows, summary = [], []
    for old_path in sorted(SNAP.rglob("*.csv")):
        rel = old_path.relative_to(SNAP)
        new_path = RES / rel
        if not new_path.exists():
            summary.append((str(rel), "missing after rerun", 0, 0))
            continue
        old, new = keyed(pd.read_csv(old_path)), keyed(pd.read_csv(new_path))
        num = [c for c in old.columns if c in new.columns and pd.api.types.is_numeric_dtype(old[c])
               and pd.api.types.is_numeric_dtype(new[c])]
        common = old.index.intersection(new.index)
        n_changed = 0
        for c in num:
            a, b = old.loc[common, c].astype(float), new.loc[common, c].astype(float)
            same = np.isclose(a, b, rtol=0, atol=1e-12, equal_nan=True) | ((a == b) & np.isinf(a))
            for idx in a.index[~same]:
                key = idx if isinstance(idx, tuple) else (idx,)
                rows.append({"file": str(rel), "row": " | ".join(map(str, key[:-1])), "column": c,
                             "before": a[idx], "after": b[idx], "delta": b[idx] - a[idx]})
                n_changed += 1
        extra = f"rows only before {len(old.index.difference(new.index))}, only after {len(new.index.difference(old.index))}"
        summary.append((str(rel), extra, len(common) * len(num), n_changed))
    out = pd.DataFrame(rows)
    out.to_csv(RES / "QUANTILE_FIX_CHANGES.csv", index=False)
    s = pd.DataFrame(summary, columns=["file", "row_alignment", "numeric_cells_compared", "cells_changed"])
    s.to_csv(RES / "QUANTILE_FIX_FILES.csv", index=False)
    with pd.option_context("display.max_rows", 500, "display.width", 200):
        print(s.to_string(index=False))
    print(f"\n{len(out)} changed cells -> {RES / 'QUANTILE_FIX_CHANGES.csv'}")


if __name__ == "__main__":
    main()
