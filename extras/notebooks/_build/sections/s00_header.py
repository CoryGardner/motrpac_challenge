# %% [markdown]
# # MoTrPAC pre-hackathon findings — replication notebook
#
# This notebook re-derives the pre-hackathon findings of the evaluation pipeline in `code/pipeline/`
# as one readable narrative. It is a presentation and a check, not new analysis (except the last,
# clearly-marked section on a covariates-only baseline).
#
# **Two notebooks.** The replication is split at the section 8/9 boundary: `01_replication.ipynb` holds
# sections 1–8 (inventory to shift tests) and `02_transfer.ipynb` holds sections 9–14 (BodyMap, GTEx,
# fusion, discordance, the QC-only baseline, limitations). Each runs on its own from the same header and
# library, and each ends with the self-check of its own rows of `expected_values.csv`.
#
# **How to read it.** Each section has three parts: *(a)* the question, why it matters and what to look
# for; *(b)* the code, commented; *(c)* what the output shows — and what it does not show. The markdown
# never states a result number: every number is printed by a code cell, computed in this run or read
# from a `results/` CSV in this run. The last cell compares the printed numbers with
# `expected_values.csv` and prints a pass/fail table.
#
# **Self-contained.** The notebook does not import the pipeline package. The pipeline's library
# (`src/tfp/`) is pasted in below, verbatim, and the logic of each pipeline script is copied into
# the section that needs it, with a comment naming its source file and function.
#
# **Two modes.**
# - `RECOMPUTE = False` (default): the expensive legs (the 100-split certificate, stability selection,
#   fusion tasks, batch-covariate permutation nulls, BodyMap/GTEx transfer) are *loaded* from
#   `results/`; everything cheap is recomputed live. Target: under 10 minutes.
# - `RECOMPUTE = True`: everything is recomputed from the exported data (about an hour); outputs go to
#   `notebooks/_outputs/`, never over `results/`.
#
# Set the mode in the next cell, or with the environment variable `NB_RECOMPUTE=1`.

# %%
# ---- imports: numpy / pandas / scikit-learn / scipy / matplotlib only (the pipeline's core-path rule) ----
from __future__ import annotations

import itertools
import json
import os
import re
import time
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from types import SimpleNamespace
from typing import Iterable, Iterator

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from IPython.display import Image, Markdown, display

warnings.filterwarnings("ignore", category=FutureWarning)
pd.set_option("display.width", 160, "display.max_columns", 30, "display.precision", 4)

# ---- mode ----
RECOMPUTE = os.environ.get("NB_RECOMPUTE", "0") == "1"


# ---- paths: found from this notebook's location (code/pipeline/notebooks/), overridable ----
def _find_pipeline_root() -> Path:
    """The pipeline directory is the nearest ancestor of the working directory holding src/tfp and
    results/. Override with MOTRPAC_PIPELINE_ROOT (the pipeline's own MOTRPAC_ROOT is not used here,
    because in the rest of the repo MOTRPAC_ROOT means the repository root)."""
    if os.environ.get("MOTRPAC_PIPELINE_ROOT"):
        return Path(os.environ["MOTRPAC_PIPELINE_ROOT"]).resolve()
    for p in [Path.cwd(), *Path.cwd().parents]:
        if (p / "src" / "tfp").is_dir() and (p / "results").is_dir():
            return p
    raise FileNotFoundError("run this notebook from code/pipeline/notebooks/ or set MOTRPAC_PIPELINE_ROOT")


PIPE = _find_pipeline_root()                 # code/pipeline
REPO = PIPE.parents[1]                       # the private workspace above the repository: only the pre-fix comparison cells (s06, s07) and the portal QC table (s14) read from it
RES = PIPE / "results"                       # the pipeline's published results (read-only here)
NB = PIPE / "notebooks"
OUT = NB / "_outputs"                        # everything this notebook writes goes here
OUT.mkdir(parents=True, exist_ok=True)
print("pipeline root :", PIPE)
print("results       :", RES)
print("outputs       :", OUT)
print("RECOMPUTE     :", RECOMPUTE)


# ---- section timer: section("3 …") closes the previous section's clock and starts a new one ----
_TIMES: dict[str, float] = {}
_CURRENT = {"name": None, "t0": None}


def section(name: str | None) -> None:
    now = time.perf_counter()
    if _CURRENT["name"] is not None:
        _TIMES[_CURRENT["name"]] = _TIMES.get(_CURRENT["name"], 0.0) + now - _CURRENT["t0"]
    _CURRENT.update(name=name, t0=now)


def timing_table() -> pd.DataFrame:
    section(None)
    t = pd.Series(_TIMES, name="seconds").round(1).to_frame()
    t.loc["TOTAL"] = t["seconds"].sum()
    return t


section("0 setup + library")


# ---- self-check registry: record(key, value) stores a number this run printed ----
CHECKS: dict[str, dict] = {}


def record(key: str, value, section_id: str, note: str = "") -> None:
    """Keep a value for the final self-check against expected_values.csv. Returns nothing; print the
    value in the section itself so the reader sees it where it is computed."""
    CHECKS[key] = {"value": float(value), "section": section_id, "note": note}


def res(*parts) -> Path:
    """Path of a published result file; fails loudly if the pipeline has not produced it."""
    p = RES.joinpath(*parts)
    if not p.exists():
        raise FileNotFoundError(f"missing pipeline result {p} (run the pipeline phase that writes it)")
    return p




# ---- small cache so sections that need the same big matrix load it once per run ----
_CACHE: dict[str, object] = {}


def cached(key: str, fn):
    if key not in _CACHE:
        _CACHE[key] = fn()
    return _CACHE[key]


# %% [markdown]
# **Every figure and table is also written to a run folder.** The next cell installs small hooks, so the
# sections need no extra code. Each run of a notebook gets a folder `notebooks/runs/<YYYY-MM-DD_HHMMSS>/`,
# which holds:
# - `figures/`: every figure shown with `plt.show()` or `show_png()`, as PNG;
# - `tables/`: every table shown with `display()` or printed with `.to_string()`, as CSV;
# - `logs/`: everything printed, one text file per section;
# - `index.csv`: one row per file (notebook, section, kind, shape, columns).
#
# To collect both notebooks in one folder, set the same `NB_RUN_STAMP` for both runs.

# %%
# ---- capture: figures / tables / printed output of this run → runs/<date_time>/ (notebook-only helper) ----
import shutil
import sys
from datetime import datetime

RUN_STAMP = os.environ.get("NB_RUN_STAMP") or datetime.now().strftime("%Y-%m-%d_%H%M%S")
RUN_DIR = NB / "runs" / RUN_STAMP
for _d in ("figures", "tables", "logs"):
    (RUN_DIR / _d).mkdir(parents=True, exist_ok=True)
_CAP = {"n": {}}


def _cap_name(kind: str, ext: str) -> Path:
    nb = globals().get("NOTEBOOK_NAME", "notebook.ipynb").replace(".ipynb", "")
    sec = re.sub(r"[^0-9A-Za-z]+", "-", str(_CURRENT["name"] or "setup")).strip("-")
    k = _CAP["n"][(nb, sec, kind)] = _CAP["n"].get((nb, sec, kind), 0) + 1
    return RUN_DIR / f"{kind}s" / f"{nb}__{sec}__{kind}{k:02d}.{ext}"


def _cap_index(path: Path, kind: str, shape="", columns="") -> None:
    row = pd.DataFrame([{"file": str(path.relative_to(RUN_DIR)), "kind": kind,
                         "notebook": globals().get("NOTEBOOK_NAME", ""), "section": _CURRENT["name"],
                         "shape": shape, "columns": columns}])
    idx = RUN_DIR / "index.csv"            # appended, so two notebooks can share one run folder
    row.to_csv(idx, mode="a", header=not idx.exists(), index=False)


def _cap_table(obj) -> None:
    if hasattr(obj, "data") and type(obj).__name__ == "Styler":
        obj = obj.data
    if isinstance(obj, pd.Series):
        obj = obj.to_frame()
    if isinstance(obj, pd.DataFrame):
        p = _cap_name("table", "csv")
        obj.to_csv(p)
        _cap_index(p, "table", "x".join(map(str, obj.shape)), ";".join(map(str, obj.columns))[:300])


# display(): save tables, then show them as usual
_display_orig = display


def display(*objs, **kw):
    for o in objs:
        try:
            _cap_table(o)
        except Exception as e:                      # capturing must never break the analysis
            print(f"[capture] table not saved: {e}")
    return _display_orig(*objs, **kw)


# DataFrame/Series.to_string(): save when called from notebook code (not from pandas' own repr)
def _wrap_to_string(cls):
    orig = cls.to_string

    def to_string(self, *a, **kw):
        caller = sys._getframe(1).f_code.co_filename
        if "site-packages" not in caller and "/pandas/" not in caller:
            try:
                _cap_table(self)
            except Exception as e:
                print(f"[capture] table not saved: {e}")
        return orig(self, *a, **kw)
    to_string.__wrapped__ = orig
    cls.to_string = to_string


for _cls in (pd.DataFrame, pd.Series):
    if not hasattr(_cls.to_string, "__wrapped__"):   # idempotent if the header cell is re-run
        _wrap_to_string(_cls)

# plt.show(): save every open figure first
_plt_show_orig = getattr(plt.show, "__wrapped__", plt.show)


def _show(*a, **kw):
    for num in plt.get_fignums():
        p = _cap_name("figure", "png")
        plt.figure(num).savefig(p, dpi=150, bbox_inches="tight")
        _cap_index(p, "figure")
    return _plt_show_orig(*a, **kw)


_show.__wrapped__ = _plt_show_orig
plt.show = _show


# show_png(): the library's plot functions save a PNG themselves; copy it into the run folder
def show_png(path: Path, width: int = 700) -> None:
    p = _cap_name("figure", "png")
    shutil.copyfile(path, p)
    _cap_index(p, "figure")
    _display_orig(Image(filename=str(path), width=width))


# printed output: tee stdout into one log file per section
class _Tee:
    def __init__(self, stream):
        self._s, self._f, self._name = stream, None, None

    def write(self, t):
        name = _CURRENT["name"]
        if name != self._name:
            if self._f:
                self._f.close()
            nb = globals().get("NOTEBOOK_NAME", "notebook.ipynb").replace(".ipynb", "")
            sec = re.sub(r"[^0-9A-Za-z]+", "-", str(name or "setup")).strip("-")
            self._f, self._name = open(RUN_DIR / "logs" / f"{nb}__{sec}.txt", "a"), name
        self._f.write(t)
        self._f.flush()
        return self._s.write(t)

    def flush(self):
        return self._s.flush()

    def __getattr__(self, k):
        return getattr(self._s, k)


if type(sys.stdout).__name__ != "_Tee":     # idempotent if the cell is re-run
    sys.stdout = _Tee(sys.stdout)
print("figures, tables and printed output of this run are saved under", RUN_DIR)
