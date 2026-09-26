"""Append-only markdown report at results/REPORT.md, plus small table helpers."""
from __future__ import annotations

import subprocess
import sys
from datetime import datetime
from pathlib import Path

import pandas as pd

from . import config as C


def _git_hash() -> str:
    try:
        return subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], cwd=C.ROOT,
                                       stderr=subprocess.DEVNULL, text=True).strip()
    except Exception:
        return "n/a"


def df_to_md(df: pd.DataFrame, floatfmt: str = ".3f", max_rows: int = 60) -> str:
    """Markdown table without the tabulate dependency."""
    d = df.copy()
    if len(d) > max_rows:
        d = d.head(max_rows)
        truncated = f"\n\n_(first {max_rows} of {len(df)} rows)_"
    else:
        truncated = ""
    for c in d.columns:
        if pd.api.types.is_float_dtype(d[c]):
            d[c] = d[c].map(lambda v: "" if pd.isna(v) else format(v, floatfmt))
    cols = [str(c) for c in d.columns]
    lines = ["| " + " | ".join(cols) + " |", "|" + "|".join(["---"] * len(cols)) + "|"]
    for _, r in d.iterrows():
        lines.append("| " + " | ".join(str(v) for v in r.tolist()) + " |")
    return "\n".join(lines) + truncated


import os


def add_section(title: str, body: str, script: str | None = None, params: dict | None = None,
                path: Path | None = None) -> Path:
    if os.environ.get("MOTRPAC_NO_REPORT"):  # verification runs: write results, skip the log
        return
    path = Path(path) if path else C.RESULTS_DIR / "REPORT.md"
    path.parent.mkdir(parents=True, exist_ok=True)
    script = script or Path(sys.argv[0]).name
    header = f"\n\n## {title}\n\n_{datetime.now():%Y-%m-%d %H:%M} · `{script}` · git {_git_hash()}"
    if params:
        header += " · " + ", ".join(f"{k}={v}" for k, v in params.items())
    header += "_\n\n"
    if not path.exists():
        path.write_text("# MoTrPAC hackathon — running report\n\nSections are appended by each phase script; "
                        "re-running a script appends a new section (keep the latest).\n")
    with path.open("a") as f:
        f.write(header + body.rstrip() + "\n")
    return path


def figure_md(path: Path, caption: str) -> str:
    rel = Path(path).relative_to(C.RESULTS_DIR) if str(path).startswith(str(C.RESULTS_DIR)) else Path(path)
    return f"![{caption}]({rel.as_posix()})\n\n_{caption}_"
