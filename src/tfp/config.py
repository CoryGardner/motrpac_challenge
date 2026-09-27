"""Paths, constants, and small helpers shared by every phase."""
from __future__ import annotations

import os
from pathlib import Path

ROOT = Path(os.environ.get("MOTRPAC_ROOT", Path(__file__).resolve().parents[2]))
DATA_DIR = ROOT / "data"
RAW_DIR = Path(os.environ.get("MOTRPAC_RAW", DATA_DIR / "raw"))
NORM_DIR = RAW_DIR / "norm"
COUNTS_DIR = RAW_DIR / "counts"
DA_DIR = RAW_DIR / "da"
META_DIR = RAW_DIR / "meta"
CODES_DIR = RAW_DIR / "codes"
EXTERNAL_DIR = DATA_DIR / "external"
RESULTS_DIR = Path(os.environ.get("MOTRPAC_RESULTS", ROOT / "results"))
FIG_DIR = RESULTS_DIR / "figures"
# The committed snapshot of every result file the site export reads (scripts/33_freeze_results.py).
FROZEN_DIR = ROOT / "results_frozen"


def results_root(explicit=None) -> Path:
    """Where the READ side (site export, tests, figures) finds the result tables.
    Order: an explicit path (--results) > $TFP_RESULTS > results/ when it holds a complete run
    (06_conformal/TRNSCRPT/coverage.csv) > results_frozen/. Writers (the numbered phases) always use RESULTS_DIR;
    nothing ever writes into results_frozen/."""
    if explicit:
        return Path(explicit)
    env = os.environ.get("TFP_RESULTS")
    if env:
        return Path(env)
    if (RESULTS_DIR / "06_conformal" / "TRNSCRPT" / "coverage.csv").exists():
        return RESULTS_DIR
    return FROZEN_DIR

SEED = 20260925

TISSUES = [
    "ADRNL", "BAT", "BLOOD", "COLON", "CORTEX", "HEART", "HIPPOC", "HYPOTH", "KIDNEY", "LIVER",
    "LUNG", "OVARY", "PLASMA", "SKM-GN", "SKM-VL", "SMLINT", "SPLEEN", "TESTES", "VENACV", "WAT-SC",
]
ASSAYS = ["TRNSCRPT", "PROT", "PHOSPHO", "ACETYL", "UBIQ", "METAB", "IMMUNO", "ATAC", "METHYL"]
GROUP_ORDER = ["control", "1w", "2w", "4w", "8w"]
SEXES = ["female", "male"]

# Sample-level tables start with these four columns, then one column per viallabel.
LEADING_COLS = ["feature", "feature_ID", "tissue", "assay"]

PANEL_GRID = [1, 2, 3, 5, 8, 10, 15, 20, 30, 50, 100]
ALPHAS = [0.05, 0.10, 0.20]


def tissue_token(tissue: str) -> str:
    """'SKM-GN' -> 'SKMGN' (R object / file-name form)."""
    return tissue.replace("-", "").upper()


_TOKEN_TO_TISSUE = {tissue_token(t): t for t in TISSUES}


def tissue_from_token(token: str) -> str:
    """'SKMGN' -> 'SKM-GN'; unknown tokens pass through unchanged."""
    return _TOKEN_TO_TISSUE.get(token.upper(), token)


def ensure_dirs() -> None:
    for d in (RESULTS_DIR, FIG_DIR):
        d.mkdir(parents=True, exist_ok=True)
