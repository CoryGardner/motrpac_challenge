#!/usr/bin/env python
"""The README's numbered text, generated from site/data (whose values carry provenance entries in
site/data/provenance.json). Every README block that contains a result number sits between
`<!-- BEGIN generated:NAME -->` and `<!-- END generated:NAME -->` and is written here; nothing inside is typed by hand.

  python tools/readme_blocks.py --write   # rewrite the blocks in README.md
  python tools/readme_blocks.py --check   # exit 1 if README.md differs from the generator (tests/test_readme.py)

Blocks: how (how the submission meets the track brief), useit (the Check samples summary and flag table), answer (the
four-paragraph answer), outputs (track outputs → where they are), keyresults (scripts/30_export_site_data.py
--readme-table, the same function the site export uses).
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SITE = ROOT / "site" / "data"
README = ROOT / "README.md"
LIVE = "https://corygardner.github.io/motrpac_challenge/"
TEAM = "https://github.com/Stanford-Bioinformatics-Center/multiomics-hackathon-2026-track-3"


def load(name):
    return json.loads((SITE / name).read_text())


def f3(v, d=3):
    return f"{v:.{d}f}"


def pct1(v):
    """97.6 %, 100 %: one decimal where needed (as on the site)."""
    s = f"{100 * v:.1f}"
    return (s[:-2] if s.endswith(".0") else s) + " %"


def js_pct(v, d=1):
    """The site's pct(v, 1): JavaScript toFixed on 100 × v (so 0.6265 → 62.6 %, as the page shows it)."""
    x = 100 * v
    s = f"{x:.{d}f}"
    # toFixed rounds the binary value; Python's format does the same (round-half-even on exact ties is not reached here)
    return s + " %"


def ids():
    return {e["id"]: e.get("value") for e in load("provenance.json")["entries"]}


def blocks() -> dict[str, str]:
    H, PR, MO, SC, PC, BM = load("headline.json"), load("product.json"), load("multiomic.json"), load("stable_core.json"), load("panel_curve.json"), load("bodymap.json")
    P = ids()
    ex, d = H["extras"], H["design"]
    tile = {t["id"]: t for t in H["tiles"]}
    lad = lambda rung, model="k20": next(r for r in H["ladder"] if r["rung_id"] == rung and r["model"] == model and r["variant"] == "marginal" and r.get("calibration", "pooled") in ("pooled", None))
    bl = {b["model"]: b["balanced_accuracy_mean"] for b in PC["baselines"]}
    acc = H["accuracy"]
    out = {}

    # ---- how we meet the brief ----------------------------------------------------------------------------------
    out["how"] = "\n".join([
        f"- **No animal-level leakage.** Whole animals are held out in {d['n_outer_folds']} animal-grouped folds "
        f"({d['n_train_animals']} training and {d['n_test_animals']} test animals per fold); every data-dependent step "
        "(imputation, variance prefilter, gene selection, scaling, tuning) is fit inside the fold, and the "
        f"{d['n_cal_animals']} calibration animals of the prediction sets are disjoint from the fit and test animals "
        "(`docs/EVALUATION_RULES.md`, `src/tfp/splits.py`).",
        f"- **Simple baselines on the same folds.** All-gene nearest centroid {f3(bl['centroid'])}, L2 logistic regression "
        f"{f3(bl['logreg_l2'])} and random forest {f3(bl['rf'])} balanced accuracy, against {f3(acc['k20']['mean'])} for the "
        f"{d['k_panel']}-gene panel; a univariate F-test selector at the same {d['k_panel']} genes reaches "
        f"{f3(ex['acc_fclassif_k20'])}, which is why the class-aware selector matters.",
        "- **The four outputs.** Each is mapped to where it lives in [Track outputs](#track-outputs--where-they-are).",
    ])

    # ---- use it ---------------------------------------------------------------------------------------------------
    F = PR["flag_key"]
    n_ex = sum(1 for s in load("expr_bodymap.json")["samples"] if s["age_weeks"] == 21)
    row = lambda lab, s: f"| {lab} | {F[s + '.0.1.n_samples']} | {js_pct(F[s + '.0.1.false_mismatch_rate'])} | {js_pct(F[s + '.0.1.cant_confirm_rate'])} | {js_pct(F[s + '.0.1.swap_vial_detection_rate'])} | {js_pct(F[s + '.0.1.swap_pair_detection_rate'])} |"
    K = PR["scaling_key"]
    out["useit"] = "\n".join([
        f"**[Check samples]({LIVE})**: *is this sample the tissue you think it is?* Drop in a rat RNA-seq table (the "
        f"{d['k_panel']} panel genes in log2 CPM, or a full raw-count matrix) with an optional claimed-tissue column; for "
        f"every sample the page returns a tissue call, its {pct1(1 - d['alpha'])} prediction set (the tissues the model "
        "cannot rule out), the genes behind the call, and a flag when the label does not fit. It runs in the browser; "
        f"nothing is uploaded. “Try the example” loads {n_ex} rat BodyMap samples from another laboratory with two labels "
        "deliberately swapped: both are flagged.",
        "",
        f"How good the flag is, at α = {d['alpha']:.2f}, on existing held-out scores (`results_product/40_product/flag_rates.csv`):",
        "",
        "| data | samples | correct labels flagged Mismatch | correct labels Can't confirm | swapped labels flagged Mismatch | swaps with ≥ 1 of the pair flagged |",
        "|---|---|---|---|---|---|",
        row("MoTrPAC held-out animals", "motrpac_heldout"),
        row("rat BodyMap 21-week adults (another lab)", "bodymap_adult_21wk"),
        "",
        f"A small or single-tissue upload is never z-scored within itself: on the same BodyMap adults, scaling each organ "
        f"alone names {pct1(K['adult_21wk.within_organ_alone.accuracy'])} of mapped organs, against "
        f"{pct1(K['adult_21wk.within_all.accuracy'])} for the whole mixed set and {pct1(K['adult_21wk.reference.accuracy'])} "
        "with MoTrPAC reference scaling (`results_product/40_product/scaling.csv`).",
    ])

    # ---- the answer ---------------------------------------------------------------------------------------------------
    bm, gt, tr = lad("different_lab"), lad("different_species"), lad("train_control_test_trained")
    om = BM["organ_map"]
    RD = PR["recal_draws"]
    b3, g5, g3 = RD["bodymap.3"], RD["gtex.5"], RD["gtex.3"]
    pc1 = next(r for r in MO["scales"]["rows"] if r["PC"] == "PC1")
    jr = next(r for r in MO["ladder_species"] if r["target"].startswith("protein → Jiang 2020 (cleaned") and r["model"] == "k20")
    p42 = next(r for r in MO["ladder_species"] if r["layer"] == "protein" and r["target"].startswith("protein, 7-class") and r["model"] == "k20")
    r42 = next(r for r in MO["ladder_species"] if r["layer"] == "RNA" and r["target"].startswith("RNA, 7-class") and r["model"] == "k20")
    late = next(r for r in MO["fusion"] if r["model"] == "k20" and r["layer"] == "late_mean")
    mw = next(m for m in MO["metabolites"] if m["leg"] == "deep_mw")
    out["answer"] = "\n\n".join([
        f"**{d['k_panel']} genes are enough.** Selected inside each animal-grouped fold by a class-aware round-robin rule, a "
        f"{d['k_panel']}-gene panel identifies {ex['n_tissues']} rat tissues at {f3(acc['k20']['mean'])} ± {f3(acc['k20']['sd'])} "
        f"balanced accuracy ({int('k50'[1:])} genes "
        f"{f3(acc['k50']['mean'])}, all genes {f3(acc['full']['mean'])}). The panel fit on all animals names "
        f"{'every mapped adult organ' if bm['accuracy'] == 1 else f3(bm['accuracy']) + ' of mapped adult organs'} correctly in another "
        f"laboratory's rats (rat BodyMap: {sum(1 for v in om.values() if v)} of {len(om)} organs have a MoTrPAC counterpart; "
        f"{bm['n_samples']} samples from {bm['n_individuals']} animals).",
        f"**The guarantee holds in the study; coverage travels honestly.** The {pct1(1 - d['alpha'])} conformal guarantee "
        f"holds within the study and on trained animals (fit on the {tr['n_source_animals']} sedentary controls alone, the "
        f"panel names the tissue of all {tr['n_individuals']} trained animals at {f3(tr['accuracy'])} with coverage "
        f"{f3(tr['coverage'])}). Beyond the study it abstains rather than errs: calibrated on MoTrPAC it covers "
        f"{f3(bm['coverage'])} of the BodyMap adults and {f3(gt['coverage'])} of human GTEx samples, and the shortfall is "
        f"empty sets, not confident error. {b3['n_recal']} animals from the new laboratory restore observed coverage "
        f"({f3(b3['mean_coverage_all_draws'])}; {b3['min_coverage']:.2f}–{b3['max_coverage']:.2f} per draw); across species {g5['n_recal']} "
        f"donors restore observed coverage of {f3(g5['mean_coverage_all_draws'])} at {g5['mean_set_size_finite']:.1f} tissues per "
        f"set, while with {g3['n_recal']} donors {g3['n_infinite']} of {g3['draws']} draws have no finite threshold: the number returns "
        "before the information does. Recalibrated coverage is observed across draws, not a guarantee for new animals.",
        f"**It is biology, not processing.** Like every large multi-tissue design, this study processed each tissue as a "
        "unit, so within-study accuracy alone cannot say how much of a fingerprint is biology. Two things can: the external "
        f"replicate, and MoTrPAC's bridging reference pools, on which batch measured directly is "
        f"{pct1(ex['bridge_pools_min_sum_ratio_all_genes'])}–{pct1(ex['bridge_pools_max_sum_ratio_all_genes'])} of the variance "
        f"that separates tissues across the {ex['bridge_n_pools']} pools. Training itself barely moves the fingerprint "
        "(see the Exercise page).",
        f"**Proteins and metabolites carry the tissue axis too, more weakly.** On the portal's reporter-ion intensities tissue "
        f"explains R² {f3(pc1['r2_rii'])} of the first proteomics component ({pc1['r2_ratio']:.4f} on the distributed ratios, "
        f"which are built to remove it); a {d['k_panel']}-protein panel names the tissue of {f3(jr['accuracy'])} of "
        f"{jr['n_samples']} human samples from another laboratory (chance {1 / jr['n_classes']:.3f}), below RNA on the same "
        f"{p42['n_samples']} samples ({f3(p42['accuracy'])} vs {f3(r42['accuracy'])}), and fusing the two does not help "
        f"({f3(late['accuracy'])}). A {d['k_panel']}-metabolite panel names the organ of {f3(mw['acc_k20'])} of mouse samples "
        f"(chance {f3(mw['chance'])}). The core fingerprint is RNA ([Multiomic page]({LIVE}multiomic.html)).",
    ])

    # ---- track outputs --------------------------------------------------------------------------------------------------
    k, ncore = d["k_panel"], len(SC["core"])
    out["outputs"] = "\n".join([
        "| output the track names | what it is | where |",
        "|---|---|---|",
        f"| A classifier | the {k}-gene logistic regression with {pct1(1 - d['alpha'])} conformal prediction sets, applied to your samples | [Check samples]({LIVE}) (calls, sets, claim checks, CSV and report); `site/data/panel_model.json` (the model the browser runs); [`src/tfp/models.py`](src/tfp/models.py) |",
        f"| Minimal tissue-signature panel | the {k} genes and the {ncore}-gene stable core | [Panel page]({LIVE}fingerprint.html); `site/data/panel_card.csv`, `site/data/panel_card.json` |",
        f"| Feature-selection workflow | the class-aware round-robin selector, fitted inside animal-grouped folds | [Methods: the selector]({LIVE}methods.html#selector); [Reference atlas: panel builder]({LIVE}explore.html#panel-builder); [`RoundRobinSelector`](src/tfp/models.py#L58) in `src/tfp/models.py` |",
        f"| Interactive model-explanation tool | per sample: probabilities against the calibrated threshold, per-gene contributions (“why X, not Y”), gene values against the reference tissues, a reference map | [Check samples]({LIVE}) (the sample drawer); [Reference atlas]({LIVE}explore.html): [tissue card]({LIVE}explore.html#tissue-card), [gene explorer]({LIVE}explore.html#gene-explorer), [panel builder]({LIVE}explore.html#panel-builder) |",
    ])

    # ---- key results: the site export's own generator -----------------------------------------------------------------------
    spec = importlib.util.spec_from_file_location("export30", ROOT / "scripts" / "30_export_site_data.py")
    mod = importlib.util.module_from_spec(spec)
    sys.path.insert(0, str(ROOT / "src"))
    spec.loader.exec_module(mod)
    out["keyresults"] = mod.readme_table(load("provenance.json")["entries"])
    return out


PAT = re.compile(r"(<!-- BEGIN generated:(\w+) -->\n)(.*?)\n?(<!-- END generated:\2 -->)", re.S)


def render(text: str, B: dict[str, str]) -> str:
    missing = set(B) - {m.group(2) for m in PAT.finditer(text)}
    if missing:
        raise SystemExit(f"README.md lacks the generated blocks: {sorted(missing)}")
    return PAT.sub(lambda m: m.group(1) + B[m.group(2)] + "\n" + m.group(4), text)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--write", action="store_true")
    g.add_argument("--check", action="store_true")
    a = ap.parse_args()
    text = README.read_text()
    new = render(text, blocks())
    if a.write:
        README.write_text(new)
        print("README.md generated blocks written")
        return 0
    if new != text:
        print("README.md differs from tools/readme_blocks.py: run `python tools/readme_blocks.py --write`")
        return 1
    print("README.md generated blocks up to date")
    return 0


if __name__ == "__main__":
    sys.exit(main())
