# %% [markdown]
# ## 14. Limitations
#
# **Question.** What should a reader *not* conclude from the sections above?
#
# The printed facts behind each point are computed in the next cell, so the list stays tied to this run.
#
# 1. **Small n.** Each tissue has a few dozen animals, so every accuracy and coverage above has
#    wide uncertainty. Per-fold spread is printed in each section; read it before the mean.
# 2. **Batch is nested in tissue.** Each tissue was extracted, library-prepped and sequenced as
#    one batch. Section 13 shows that QC covariates alone classify tissue, so within-study accuracy cannot, on
#    its own, separate biology from processing. The external transfers (sections 9–10, notebook 02) carry
#    that argument, not the within-study numbers.
# 3. **Pre-split gene filter (known deviation, not fixed).** When the transcript counts are stacked,
#    a gene is kept if it has enough counts in enough samples of *at least one tissue*. That filter
#    uses all samples and their tissue labels *before* the animal split, contrary to the pipeline's own
#    rule 2 ("everything data-dependent happens inside the training fold"). The leak is presumably
#    mild (it selects genes, not a model), but its size was never measured. It is left as is, so the
#    notebook reproduces the published numbers; see `docs/findings/FINDINGS_REPORT.md` §14, item 13(b).
#    The code that does it is printed below.
# 4. **Single-sex tissues.** Ovary and testes exist for one sex only, so the sex-shift tests
#    (section 8) either drop them or score them as unseen classes. Sex is also confounded with arrival
#    cohort within every group.
# 5. **Design confounds in the time course.** Only 8 w vs control is matched on arrival cohort and
#    sacrifice date. 1 w / 2 w / 4 w shift results mix training duration with sacrifice date.
# 6. **Super-class scoring in external transfer.** BodyMap scores muscle ({SKM-GN, SKM-VL}) and brain
#    ({CORTEX, HIPPOC, HYPOTH}) as super-classes, and GTEx scores only muscle that way. Accuracy there
#    answers "right organ" for those tissues, not "right tissue". Seven tissues have no BodyMap organ.
# 7. **Assumed BodyMap animal IDs.** BodyMap samples carry no public animal ID. Its calibration and
#    recalibration splits group on *assumed* animal IDs, so its conformal guarantees rest on that
#    assumption.
# 8. **Conformal quantile (fixed 2026-09-25).** The pipeline's quantile was one rank more conservative than
#    the textbook rule when there were 19 or more calibration points. With 8 or fewer, it used the largest
#    score instead of "all classes", which is less conservative. It is now fixed and the affected phases
#    rerun; `results/QUANTILE_FIX_CHANGES.md` lists what moved. Section 7 checks that the certificate,
#    which counts errors, did not move.
# 9. **One release.** Everything here is portal release c1.0 (rn6), the R package. The rn7 c2.0
#    reprocessing changes training-regulated sets substantially. Any claim about *which* genes carry the
#    signal is release-specific.

# %%
section("14 limitations")
# The facts behind points 1, 3 and 4, computed from this run's data rather than stated from memory.
# Point 1: animals per tissue in the transcript classifier's matrix (section 4), else in all study vials.
if "S04_OM" in globals() or "s14_om" in globals():
    _m, _what = (S04_OM if "S04_OM" in globals() else s14_om).meta, "TRNSCRPT classifier matrix (sections 4 / 13)"
else:
    _m = cached("pheno", io.load_pheno)
    _m, _what = _m[_m["viallabel"].astype(str).str.startswith("9")], "all study vials, all assays (PHENO)"
_by = _m.groupby("tissue")["pid"].nunique()
print(f"animals per tissue — {_what}: min / median / max =", int(_by.min()), "/", int(_by.median()), "/", int(_by.max()))
# Point 4: tissues sampled from one sex only.
_sx = _m.groupby("tissue")["sex"].nunique()
print("single-sex tissues:", sorted(_sx[_sx == 1].index))
# Point 3: show the filter as it sits in the (pasted) library, so the deviation is visible, not just asserted.
import inspect
print("\nio._stack_counts — the per-tissue expression mask, computed on all samples before any split:")
for _ln in inspect.getsource(io._stack_counts).splitlines():
    if re.search(r"masks|keep_any", _ln):
        print("   ", _ln.strip())
