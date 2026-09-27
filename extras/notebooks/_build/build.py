#!/usr/bin/env python3
"""Assemble the replication notebooks from section files, with the pipeline library copied in verbatim.

Why a builder: the notebooks must contain every line of code they run (no `import tfp`), and the
library code they need must stay identical to src/tfp/. Rather than retyping ~1,600 library lines,
this script copies each module's source into notebook cells (dropping only the package-relative
imports and matplotlib's forced "Agg" backend), then appends the hand-written section files.

Section files live in _build/sections/ in "percent" format:
    # %% [markdown]
    # ## A heading
    # %%
    x = 1
Each "# %%" starts a new cell; "# %% [markdown]" cells have their leading "# " stripped.

Usage (from code/pipeline/notebooks):
    python _build/build.py                       # build 01_replication.ipynb and 02_transfer.ipynb
    python _build/build.py --only s01,s02 --out /tmp/test.ipynb --execute   # quick test of sections
"""
from __future__ import annotations

import argparse
import ast
import re
import subprocess
import sys
from pathlib import Path

import nbformat as nbf

HERE = Path(__file__).resolve().parent            # notebooks/_build
NB_DIR = HERE.parent                              # notebooks/
PIPE = NB_DIR.parent                              # code/pipeline
SRC = PIPE / "src" / "tfp"
SECTIONS = HERE / "sections"

# Library modules in dependency order, and the namespace name each is exposed under in the notebook
# (the scripts refer to them as C., io., splits., models., cp., discordance., transfer., plots.).
MODULES = [("config", "C"), ("io", "io"), ("splits", "splits"), ("models", "models"),
           ("conformal", "cp"), ("discordance", "discordance"), ("transfer", "transfer"), ("plots", "plots")]

NOTEBOOKS = {
    "01_replication.ipynb": ["s00_header", "LIB", "s01_inventory", "s02_pca", "s03_prot", "s04_baselines",
                             "s05_panels", "s06_conformal", "s07_certificate", "s08_shift", "s99_selfcheck"],
    "02_transfer.ipynb": ["s00_header", "LIB", "s09_bodymap", "s10_gtex", "s11_fusion", "s12_discordance",
                          "s14_qc_baseline", "s13_limitations", "s99_selfcheck"],
}


def md(text: str):
    return nbf.v4.new_markdown_cell(text.strip("\n"))


def code(text: str):
    return nbf.v4.new_code_cell(text.strip("\n"))


def transform_module(mod: str, src: str) -> str:
    """The three mechanical edits applied to a src/tfp module before it is pasted into the notebook:
    (1) package-relative imports removed (indented ones become `pass`), (2) matplotlib.use("Agg") removed,
    (3) in config, ROOT and RESULTS_DIR computed from the notebook's location. Also used verbatim by the
    notebook's drift check, so the check and the build can never disagree about what "the same" means."""
    lines = []
    for line in src.splitlines():
        m = re.match(r"^(\s*)from \.(\w+)? import ", line)
        if m or re.match(r"^from __future__ import", line):
            ind = m.group(1) if m else ""
            # indented (in-function) relative imports become `pass` so the block stays valid
            lines.append(f"{ind}pass  # [notebook] removed: {line.strip()}" if ind else f"# [notebook] removed: {line}")
            continue
        if line.strip() == 'matplotlib.use("Agg")':
            lines.append(f"# [notebook] removed: {line}  (keep the inline backend)")
            continue
        lines.append(line)
    body = "\n".join(lines)
    if mod == "config":
        body = re.sub(r'^ROOT = .*$', 'ROOT = PIPE  # [notebook] was: Path(os.environ.get("MOTRPAC_ROOT", Path(__file__).resolve().parents[2]))',
                      body, count=1, flags=re.M)
        body = re.sub(r'^RESULTS_DIR = .*$', 'RESULTS_DIR = PIPE / "results"  # [notebook] read-only here; notebook outputs go to OUT',
                      body, count=1, flags=re.M)
    return f"# ===== copied verbatim from src/tfp/{mod}.py (see the note above for the 3 mechanical edits) =====\n" + body


def library_cells() -> list:
    """One markdown + one code cell per library module, copied verbatim from src/tfp/*.py."""
    cells = [md("""## Appendix-in-place: the pipeline library, copied verbatim

The cells below are the source of `src/tfp/` (the pipeline's library), pasted in so this notebook runs
without importing it. Three mechanical edits only: package-relative imports (`from . import …`) are
removed because everything shares one namespace; `matplotlib.use("Agg")` is removed so figures can show;
and in `config`, the path constants are computed from this notebook's location instead of the file's.
After each module, a namespace object (`C`, `io`, `splits`, `models`, `cp`, `discordance`, `transfer`,
`plots`) is created so code copied from the scripts can keep calling e.g. `cp.predict_sets(...)`.

**Drift check.** The cell right after the library re-reads each `src/tfp/*.py` from disk (when the
notebook sits in its repository), applies the same three edits with the same function the notebook
builder uses (`transform_module`, pasted in), and compares the result with the library code that actually
executed in this kernel. It prints an identical / differs table and **raises** if any module differs. So
a run that gets past it used exactly the library in `src/`. If `src/tfp` is not reachable (the
notebook was copied elsewhere), it says the inline copy could not be verified, and continues.""")]
    for mod, alias in MODULES:
        path = SRC / f"{mod}.py"
        src = path.read_text()
        tree = ast.parse(src)
        names = []
        for node in tree.body:
            if isinstance(node, (ast.FunctionDef, ast.ClassDef)):
                names.append(node.name)
            elif isinstance(node, ast.Assign):
                names += [t.id for t in node.targets if isinstance(t, ast.Name)]
            elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name):
                names.append(node.target.id)
            elif isinstance(node, (ast.Import, ast.ImportFrom)) and getattr(node, "level", 0) == 0 \
                    and getattr(node, "module", None) != "__future__":
                # the module's own imports are attributes of the module too (scripts call e.g.
                # models.balanced_accuracy_score), so expose them on the namespace as well
                names += [(a.asname or a.name).split(".")[0] for a in node.names]
        body = transform_module(mod, src)
        cells.append(code(body))
        ns = ", ".join(f"{n}={n}" for n in dict.fromkeys(names) if not n.startswith("__"))
        cells.append(code(f"# expose module `{mod}` as `{alias}.<name>`, the way the pipeline scripts call it\n"
                          f"{alias} = SimpleNamespace({ns})"))
    import inspect
    cells.append(code(DRIFT_CHECK.replace("__TRANSFORM_SOURCE__", inspect.getsource(transform_module).rstrip())))
    return cells


DRIFT_CHECK = '''# ---- drift check: is the inline library (the cells above) still identical to src/tfp/*.py? ----
# notebook-only helper (not in the pipeline). transform_module is pasted from notebooks/_build/build.py, the
# function that generated the cells above, so "identical" means exactly what the build did.
__TRANSFORM_SOURCE__


def _executed_library_cells() -> dict:
    """The library cells as they actually ran in this kernel (IPython's input history), falling back to
    the saved notebook file. Keyed by module name, from each cell's first line."""
    hist = list(globals().get("In", []))
    if not any(h.startswith("# ===== copied verbatim from src/tfp/") for h in hist):
        import nbformat as _nbf
        hist = [c.source for c in _nbf.read(NB / NOTEBOOK_NAME, 4).cells if c.cell_type == "code"]
    out = {}
    for h in hist:
        m = re.match(r"# ===== copied verbatim from src/tfp/(\\w+)\\.py", h)
        if m:
            out[m.group(1)] = h.rstrip("\\n")      # the last execution of each module's cell wins
    return out


_src_dir = PIPE / "src" / "tfp"
_mods = ["config", "io", "splits", "models", "conformal", "discordance", "transfer", "plots"]
if not _src_dir.is_dir():
    print(f"src/tfp not found at {_src_dir}: the inline copy of the library could NOT be verified; continuing.")
else:
    _ran = _executed_library_cells()
    _rows = []
    for _m in _mods:
        _expected = transform_module(_m, (_src_dir / f"{_m}.py").read_text()).rstrip("\\n")
        _got = _ran.get(_m)
        _rows.append({"module": _m, "status": "missing" if _got is None else ("identical" if _got == _expected else "DIFFERS"),
                      "inline_lines": None if _got is None else _got.count("\\n") + 1,
                      "src_lines": _expected.count("\\n") + 1})
    _drift = pd.DataFrame(_rows)
    print(_drift.to_string(index=False))
    _bad = _drift[_drift["status"] != "identical"]
    if len(_bad):
        raise RuntimeError(f"inline library differs from src/tfp for: {', '.join(_bad['module'])}. Rebuild the "
                           "notebook (python notebooks/_build/build.py) or revert the edit; results below would not "
                           "come from the pipeline's library.")
    print("inline library identical to src/tfp (after the 3 documented edits)")
'''


def section_cells(name: str) -> list:
    text = (SECTIONS / f"{name}.py").read_text()
    chunks = re.split(r"^# %%(.*)$", text, flags=re.M)
    cells = []
    # chunks: [preamble, tag1, body1, tag2, body2, ...]
    for tag, body in zip(chunks[1::2], chunks[2::2]):
        if "[markdown]" in tag:
            body = "\n".join(re.sub(r"^# ?", "", ln) for ln in body.strip("\n").splitlines())
            cells.append(md(body))
        else:
            if body.strip():
                cells.append(code(body))
    return cells


def build(order: list[str], out: Path) -> nbf.NotebookNode:
    nb = nbf.v4.new_notebook()
    nb.metadata["kernelspec"] = {"name": "motrpac-py", "display_name": "Python (motrpac-py)", "language": "python"}
    for part in order:
        if part == "LIB":
            nb.cells += library_cells()
        elif (SECTIONS / f"{part}.py").exists():
            cells = section_cells(part)
            if part == "s00_header":   # which notebook this is (self-check rows, output file names): set first
                cells.insert(1, code(f'NOTEBOOK_NAME = "{out.name}"'))
            nb.cells += cells
        else:
            print(f"  (skipping missing section {part})", file=sys.stderr)
    nbf.write(nb, out)
    return nb


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", help="comma-separated section names to include after header+library (testing)")
    ap.add_argument("--out", type=Path, help="output path (with --only)")
    ap.add_argument("--execute", action="store_true", help="execute with nbconvert on the motrpac-py kernel")
    ap.add_argument("--timeout", type=int, default=3600)
    a = ap.parse_args()
    targets = {}
    if a.only:
        parts = ["s00_header", "LIB"] + [p for p in a.only.split(",") if p not in ("s00_header", "LIB")]
        targets[a.out or Path("/tmp/section_test.ipynb")] = parts
    else:
        targets = {NB_DIR / k: v for k, v in NOTEBOOKS.items()}
    for out, order in targets.items():
        nb = build(order, out)
        print(f"wrote {out}  ({len(nb.cells)} cells)")
        if a.execute:
            cmd = [sys.executable, "-m", "jupyter", "nbconvert", "--to", "notebook", "--execute", "--inplace",
                   f"--ExecutePreprocessor.timeout={a.timeout}", "--ExecutePreprocessor.kernel_name=motrpac-py", str(out)]
            print(" ".join(cmd))
            import os
            env = {**os.environ, "MOTRPAC_PIPELINE_ROOT": str(PIPE)}  # test notebooks may live outside notebooks/
            r = subprocess.run(cmd, cwd=NB_DIR, env=env)
            if r.returncode:
                sys.exit(r.returncode)


if __name__ == "__main__":
    main()
