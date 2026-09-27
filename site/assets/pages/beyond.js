import { mountChrome, loadJSON, el, fmt, callout, pendingBlock, tableFrom } from "../site.js";
import { figure, bar, line, tokens, palette } from "../charts.js";

async function main() {
  await mountChrome("beyond.html");
  const [BEY, BV, M] = await Promise.all([loadJSON("data/beyond.json"), loadJSON("data/batch_verdict.json"), loadJSON("data/manifest.json")]);

  // ---- status ------------------------------------------------------------------------------------------
  const blocks = document.getElementById("status-blocks");
  const row = (name, present, what, where) => el("div", { class: "callout " + (present ? "note" : "pending") }, [
    el("div", { class: "kind" }, present ? "present in results/" : "not in this copy of results/"),
    el("p", {}, [el("b", {}, name + ". "), what, " ", el("span", { class: "muted" }, present ? `(${where})` : `(expected under ${where})`)]),
  ]);
  blocks.append(
    row("Training transfer (phase 17)", BEY.training_transfer.present, "Leave-one-tissue-out prediction of training state: does a training signature learned in 18 tissues detect training in the 19th, after the tissue axis is removed? Shown only once its split function is confirmed to block on both animal and tissue and its permutation null is present.", BEY.training_transfer.present ? BEY.training_transfer.dirs.join(", ") : "results/17_*"),
    row("Identifiability audit (phase 21)", BEY.identifiability_audit.present, "The estimable-pairs framing across every omic layer and the bridge-sample variance measurement. The Identifiability page recomputes the nesting from metadata in the meantime.", BEY.identifiability_audit.present ? BEY.identifiability_audit.dirs.join(", ") : "results/21_*"),
    row("Exercise decomposition (phases 22–26)", BEY.decomposition.present, "Per-gene tissue-agnostic training effects: forest plots of the top genes (the heat-shock genes lead), the I² distribution across tissues, the cross-tissue response correlation map, the leave-one-tissue-out AUROC with a QC-covariate control arm, and the RNA-vs-protein shared-effect scatter.", BEY.decomposition.present ? BEY.decomposition.dirs.join(", ") : "results/22_* … results/26_*"),
    row("Pre-registration", BEY.preregistration.present, BEY.preregistration.note + ". The planned analyses are described from the build brief, not from a registered file.", "docs/"),
    row("Time-course investigation (phase 15)", BEY.time_course.present, BEY.time_course.note, BEY.time_course.present ? "results/15_time_course" : "results/15_*"),
  );

  // ---- invariance: phase 15 parts 5–6 --------------------------------------------------------------------
  const p5 = BV.part5_by_duration;
  if (p5 && p5.length) {
    const prim = p5.filter((r) => String(r.design).toLowerCase().includes("control") || r.design === "primary" || p5.every((x) => x.design === p5[0].design));
    const use = prim.length ? prim : p5;
    const designs = [...new Set(use.map((r) => r.design))];
    const d0 = designs.includes("control_only_fit7_cal3") ? "control_only_fit7_cal3" : designs[0];
    const rows = use.filter((r) => r.design === d0);
    const k20 = rows.filter((r) => r.arm === "k20" || r.arm === "panel_k20"), full = rows.filter((r) => r.arm === "full");
    document.getElementById("p-invariance").replaceChildren(
      `A model fit and calibrated on sedentary controls only (${d0}) was tested separately on animals trained for 1, 2, 4 and 8 weeks. Accuracy of the 20-gene panel stays between ${fmt(Math.min(...k20.map((r) => r.accuracy)))} and ${fmt(Math.max(...k20.map((r) => r.accuracy)))} across durations, coverage of the α = 0.10 sets between ${fmt(Math.min(...k20.map((r) => r.coverage)))} and ${fmt(Math.max(...k20.map((r) => r.coverage)))}: training does not move the tissue fingerprint. `,
      "That invariance is what makes the fingerprint usable as an instrument: the tissue axis can be fixed and removed before the training axis is examined.",
    );
    await figure(document.getElementById("fig-invariance"), {
      title: "Trained animals are identified as well as controls, at every duration",
      subtitle: "Accuracy and α = 0.10 coverage of a control-trained model (fit on 7, calibrated on 3 control animals) on the 10 animals of each training duration; whiskers: 95 % cluster-bootstrap over animals.",
      build: () => {
        const groups = ["1w", "2w", "4w", "8w"].filter((g) => k20.some((r) => r.test_group === g));
        const g = (rowsArm, key) => groups.map((gg) => rowsArm.find((r) => r.test_group === gg)?.[key] ?? null);
        const err = (rowsArm, key) => ({ type: "data", symmetric: false, array: groups.map((gg) => { const r = rowsArm.find((x) => x.test_group === gg); return r ? r[key + "_hi"] - r[key] : 0; }), arrayminus: groups.map((gg) => { const r = rowsArm.find((x) => x.test_group === gg); return r ? r[key] - r[key + "_lo"] : 0; }), visible: true, thickness: 1.5, width: 4, color: tokens().ink2 });
        return { traces: [
          { ...bar(groups, g(k20, "accuracy"), { name: "accuracy, 20 genes", slot: 1, hover: "%{x}: accuracy %{y:.3f}<extra>20 genes</extra>" }), error_y: err(k20, "accuracy") },
          { ...bar(groups, g(k20, "coverage"), { name: "coverage, 20 genes", slot: 2, hover: "%{x}: coverage %{y:.3f}<extra>20 genes</extra>" }), error_y: err(k20, "coverage") },
          { ...bar(groups, g(full, "accuracy"), { name: "accuracy, all genes", slot: 4, hover: "%{x}: accuracy %{y:.3f}<extra>all genes</extra>" }), error_y: err(full, "accuracy") },
          { ...bar(groups, g(full, "coverage"), { name: "coverage, all genes", slot: 5, hover: "%{x}: coverage %{y:.3f}<extra>all genes</extra>" }), error_y: err(full, "coverage") },
        ], layout: { barmode: "group", yaxis: { range: [0.7, 1.05], title: { text: "fraction" } }, xaxis: { title: { text: "training duration of the test animals" } }, margin: { t: 40 }, legend: { y: 1.14 },
                     shapes: [{ type: "line", xref: "paper", x0: 0, x1: 1, yref: "y", y0: 0.9, y1: 0.9, line: { color: tokens().ink2, width: 1, dash: "dash" } }] },
          table: { columns: ["design", "arm", "test_group", "n_test_animals", "n_test_vials", "accuracy", "accuracy_lo", "accuracy_hi", "coverage", "coverage_lo", "coverage_hi", "empty_rate", "mean_set_size"], rows } };
      },
      source: "results/15_time_course/5_6_fingerprint/part5_by_duration.csv",
      notShow: "the cohort confound: 1-, 2- and 4-week animals were collected months apart from the controls (only 8 weeks is date-matched), so those contrasts also compare arrival cohorts.",
    });
    document.getElementById("invariance-caveat").replaceChildren(callout("caveat", "Read the 1-, 2- and 4-week bars as cohort-confounded",
      "Only the 8-week animals share arrival cohort and sacrifice dates with the controls. Most of the few errors at 1 and 2 weeks are vena cava vials from female animals that the consortium flags as brown-fat contaminated (results/15_time_course/5_6_fingerprint/NOTES.md)."));
  } else {
    document.getElementById("p-invariance").textContent = "";
    document.getElementById("fig-invariance").replaceChildren(pendingBlock("Fingerprint invariance to training duration", "results/15_time_course/5_6_fingerprint/part5_by_duration.csv is not present in this copy."));
  }

  // ---- planned ------------------------------------------------------------------------------------------
  const planned = document.getElementById("planned-blocks");
  const plan = (title, files, items) => el("div", { class: "pending-block" }, [el("div", { class: "kind" }, "Pending · " + files), el("b", {}, title), el("ul", {}, items.map((i) => el("li", {}, i)))]);
  planned.append(
    plan("Training transfer across tissues", "results/17_*", [
      "Per-tissue AUROC of the leave-one-tissue-out training classifier against its permutation null, badged preliminary until the null is present.",
      "The split function is read and confirmed to block on both pid and tissue, with a passing unit test, before any number is shown.",
    ]),
    plan("Exercise decomposition", "results/22_* to results/26_*", [
      "Forest-plot picker for the top tissue-agnostic genes, leading with the heat-shock genes.",
      "The I² distribution of the per-gene training effect across tissues.",
      "The cross-tissue response correlation map.",
      "The leave-one-tissue-out AUROC with the QC-covariate control arm.",
      "The RNA-vs-protein shared-effect scatter.",
    ]),
  );
  planned.append(el("p", { class: "small" }, `Manifest: ${M.git_hash}, generated ${M.generated}; absent phases recorded in site/data/manifest.json: ${Object.keys(M.absent_phases || {}).join(", ") || "none"}.`));
}

main().catch((e) => { console.error(e); document.querySelector(".lede").textContent = "Failed to load: " + e.message; });
