#!/usr/bin/env python
"""Check the Python environment, R availability, and whether data are present."""
from __future__ import annotations

import importlib
import shutil
import sys

from tfp import config as C


def main() -> int:
    ok = True
    print(f"python {sys.version.split()[0]}")
    for mod in ("numpy", "pandas", "sklearn", "scipy", "matplotlib"):
        try:
            m = importlib.import_module(mod)
            print(f"  {mod:12s} {getattr(m, '__version__', '?')}")
        except ImportError:
            print(f"  {mod:12s} MISSING")
            ok = False
    for mod in ("umap", "pyarrow", "mapie", "statsmodels"):
        try:
            importlib.import_module(mod)
            print(f"  {mod:12s} (optional) present")
        except ImportError:
            print(f"  {mod:12s} (optional) absent")
    r = shutil.which("Rscript")
    print(f"Rscript: {r or 'not found (needed only for R/export_motrpac.R)'}")
    print(f"RAW_DIR: {C.RAW_DIR} {'(exists)' if C.RAW_DIR.exists() else '(missing)'}")
    n_norm = len(list((C.RAW_DIR / 'norm').glob('*.csv'))) if (C.RAW_DIR / 'norm').exists() else 0
    n_counts = len(list((C.RAW_DIR / 'counts').glob('*.csv'))) if (C.RAW_DIR / 'counts').exists() else 0
    print(f"  norm files: {n_norm}, count files: {n_counts}, pheno: {(C.RAW_DIR / 'pheno.csv').exists()}")
    if (C.RAW_DIR / "_SYNTHETIC").exists():
        print("  NOTE: data/raw holds SYNTHETIC data (delete it before the real export)")
    print("environment OK" if ok else "environment INCOMPLETE")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
