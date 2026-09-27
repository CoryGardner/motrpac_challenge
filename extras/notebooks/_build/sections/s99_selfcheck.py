# %% [markdown]
# ## Self-check: did this run reproduce the expected values?
#
# Every headline number printed above was also stored with `record(...)`. This cell compares each one
# with `expected_values.csv`:
# - `published_value` is what the pipeline reported before the 2026-09-25 conformal-quantile fix;
# - `reference_value` is what this notebook should reproduce now, after the fix (the same as published
#   where the fix did not touch the number).
#
# A row passes when |value − reference| ≤ tolerance. Tolerances were set per row before running the
# notebook, and are never widened to make a row pass. A failure is a finding to explain, not a number to
# tune. Rows the expected table lists for the *other* notebook are skipped.

# %%
section("99 self-check")
EXPECTED = pd.read_csv(NB / "expected_values.csv", dtype={"key": str, "section": str})
_got = pd.DataFrame.from_dict(CHECKS, orient="index")
_chk = EXPECTED.set_index("key").join(_got[["value"]] if len(_got) else pd.DataFrame(columns=["value"]), how="left")
_chk = _chk[_chk["notebook"] == NOTEBOOK_NAME] if "notebook" in _chk.columns else _chk
_chk["abs_diff"] = (_chk["value"] - _chk["reference_value"]).abs()
_chk["status"] = np.where(_chk["value"].isna(), "NOT RUN",
                          np.where(_chk["abs_diff"] <= _chk["tolerance"] + 1e-12, "pass", "FAIL"))
_chk["moved_by_fix"] = (_chk["published_value"] - _chk["reference_value"]).abs() > 1e-9
with pd.option_context("display.max_rows", 500, "display.precision", 4):
    display(_chk[["section", "published_value", "reference_value", "value", "tolerance", "abs_diff",
                  "status", "moved_by_fix", "source"]])
_unlisted = sorted(set(CHECKS) - set(EXPECTED["key"]))
print(_chk["status"].value_counts().to_string())
if _unlisted:
    print("recorded but not in expected_values.csv:", _unlisted)
_chk.reset_index().to_csv(OUT / f"selfcheck_{NOTEBOOK_NAME.replace('.ipynb', '')}.csv", index=False)

# %%
# Run time, per section and in total (this run, this machine).
_t = timing_table()
display(_t)
_t.to_csv(OUT / f"timings_{NOTEBOOK_NAME.replace('.ipynb', '')}_{'recompute' if RECOMPUTE else 'default'}.csv")

# %% [markdown]
# ### Where the code came from
# The cell below reads this notebook's own file and lists every provenance comment. These are the lines
# that name the pipeline source a block was copied from: `src/tfp/*.py` for the library cells, and
# `scripts/*.py::function` for the script logic. It also lists every helper marked notebook-only (not in
# the pipeline).

# %%
import nbformat as _nbf
_nbpath = NB / NOTEBOOK_NAME
if _nbpath.exists():
    _rows = []
    for _i, _c in enumerate(_nbf.read(_nbpath, 4).cells):
        if _c.cell_type != "code":
            continue
        for _ln in _c.source.splitlines():
            _m = re.search(r"(copied verbatim from src/tfp/\S+|from (?:scripts|code/probes)/[\w./]+(?:::[\w, .]+)?|notebook-only helper)", _ln)
            if _m and _ln.lstrip().startswith("#"):
                _rows.append({"cell": _i, "provenance": _m.group(1).strip(" ,.")})
    _prov = pd.DataFrame(_rows).drop_duplicates("provenance")
    with pd.option_context("display.max_rows", 300, "display.max_colwidth", 120):
        display(_prov.reset_index(drop=True))
else:
    print("notebook file not found next to the kernel's working directory:", _nbpath)
