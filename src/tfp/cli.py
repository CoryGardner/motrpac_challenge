"""Shared argparse plumbing for the numbered phase scripts."""
from __future__ import annotations

import argparse
from pathlib import Path

from . import config as C


def common_parser(description: str) -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=description)
    ap.add_argument("--assay", default="TRNSCRPT", help="assay code (TRNSCRPT, PROT, METAB, ...)")
    ap.add_argument("--source", default="auto", choices=["auto", "norm", "counts"],
                    help="sample-level source; auto = counts for TRNSCRPT if present, else norm")
    ap.add_argument("--tissues", default=None, help="comma-separated subset of tissues (default: all available)")
    ap.add_argument("--n-splits", type=int, default=5)
    ap.add_argument("--prefilter", type=int, default=5000, help="variance top-k prefilter inside each fold")
    ap.add_argument("--quick", action="store_true", help="small/fast run (fewer folds, smaller grids)")
    ap.add_argument("--seed", type=int, default=C.SEED)
    ap.add_argument("--out", default=None, help="results subdirectory (default: results/<phase>)")
    ap.add_argument("--drop-incomplete-samples", type=float, default=None, metavar="FRAC",
                    help="before --complete-features: drop samples missing more than FRAC of the joined features "
                         "(platform-level gaps), logged")
    ap.add_argument("--complete-features", action="store_true",
                    help="after the cross-tissue join, keep only features with no missing value in any sample "
                         "(use when the missingness pattern itself encodes tissue)")
    return ap


def resolve_source(assay: str, source: str) -> str:
    if source != "auto":
        return source
    if assay == "TRNSCRPT" and any((C.COUNTS_DIR).glob(f"{assay}__*.csv")):
        return "counts"
    return "norm"


def parse_tissues(arg: str | None) -> list[str] | None:
    if not arg:
        return None
    return [C.tissue_from_token(C.tissue_token(t.strip())) for t in arg.split(",") if t.strip()]


def is_synthetic() -> bool:
    return (C.RAW_DIR / "_SYNTHETIC").exists()


def outdir(phase: str, override: str | None = None) -> Path:
    d = Path(override) if override else C.RESULTS_DIR / phase
    d.mkdir(parents=True, exist_ok=True)
    return d


def banner(phase: str, args) -> None:
    tag = " [SYNTHETIC DATA]" if is_synthetic() else ""
    print(f"=== {phase}{tag} ===")
    print("   " + ", ".join(f"{k}={v}" for k, v in vars(args).items()))
