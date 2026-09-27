// The transfer ladder — the site's signature figure, shared by the home and transfer pages.
import { el, fmt, segmented, control, tableFrom } from "./site.js";
import { figure, bar, refLine, tokens, template, CONFIG } from "./charts.js";

export const RUNGS = [
  { id: "in_distribution", label: "In-distribution<br>(held-out animals)", short: "In-distribution" },
  { id: "train_control_test_trained", label: "Trained animals<br>(fit on controls only)", short: "Training state (controls → trained)" },
  { id: "train_male_test_female", label: "Held-out sex<br>males → females", short: "Held-out sex (M → F)" },
  { id: "train_female_test_male", label: "Held-out sex<br>females → males", short: "Held-out sex (F → M)" },
  { id: "different_lab", label: "Other laboratory<br>(rat BodyMap adults)", short: "Other laboratory (BodyMap)" },
  { id: "different_species", label: "Other species<br>(human GTEx)", short: "Other species (GTEx)" },
];

export function pick(ladder, id, model, variant, calib = "pooled") {
  return ladder.find((r) => r.rung_id === id && r.model === model && r.variant === variant && (r.calibration || "pooled") === calib);
}

export function ladderBuild(H, state) {
  const t = tokens();
  const x = RUNGS.map((r) => r.label);
  const rows = RUNGS.map((r) => pick(H.ladder, r.id, state.model, state.variant, r.id === "in_distribution" ? state.calib : "pooled"));
  const acc = rows.map((r) => (r && !r.pending ? r.accuracy : null));
  const cov = rows.map((r) => (r && !r.pending ? r.coverage : null));
  // whiskers: sd over folds where the design has folds, else the 95 % cluster-bootstrap interval over individuals
  const err = (key, ciKey, sdKey) => {
    const plus = rows.map((r) => (r && !r.pending ? (r[ciKey] ? r[ciKey][1] - r[key] : r[sdKey] || 0) : 0));
    const minus = rows.map((r) => (r && !r.pending ? (r[ciKey] ? r[key] - r[ciKey][0] : r[sdKey] || 0) : 0));
    return plus.some((v) => v) || minus.some((v) => v) ? { type: "data", symmetric: false, array: plus, arrayminus: minus, visible: true, color: tokens().ink2, thickness: 1.5, width: 4 } : undefined;
  };
  const ci = (c) => (c ? ` [${fmt(c[0])}, ${fmt(c[1])}]` : "");
  const hoverAcc = rows.map((r) => (r && !r.pending ? `accuracy ${fmt(r.accuracy)}${r.accuracy_sd ? " ± " + fmt(r.accuracy_sd) + " (sd over folds)" : ""}${r.accuracy_ci ? " 95 % CI" + ci(r.accuracy_ci) : ""}<br>n = ${r.n_samples} samples, ${r.n_individuals} individuals${r.n_source_animals ? "<br>source: " + r.n_source_animals + " animals, " + r.n_calibration_animals + " of them calibration" : ""}${r.accuracy_cv40 ? "<br>same panel size under the 40-animal CV of the Fingerprint page: " + fmt(r.accuracy_cv40) : ""}${r.accuracy_refit !== null && r.accuracy_refit !== undefined ? "<br>the refit model that carries the sets: " + fmt(r.accuracy_refit) : ""}${r.accuracy_seen !== null && r.accuracy_seen !== undefined ? "<br>the unseen sex-specific tissue's vials count as wrong; balanced accuracy over seen classes " + fmt(r.accuracy_seen) : ""}` : "pending"));
  const hoverCov = rows.map((r) => (r && !r.pending ? `coverage ${fmt(r.coverage)}${r.coverage_sd ? " ± " + fmt(r.coverage_sd) + " (sd over folds)" : ""}${r.coverage_ci ? " 95 % CI" + ci(r.coverage_ci) : ""}${r.n_samples_coverage ? " (n = " + r.n_samples_coverage + " seen-class vials)" : ""}<br>empty sets ${r.empty === null || r.empty === undefined ? "—" : fmt(r.empty)}${r.wrong_non_empty !== null && r.wrong_non_empty !== undefined ? ", non-empty but wrong " + fmt(r.wrong_non_empty) : ""}<br>mean set size ${fmt(r.set_size, 2)}${r.n_classes_seen && r.set_size >= r.n_classes_seen ? " — every seen tissue in every set (per-class thresholds +∞)" : ""}${r.recal_n3 !== null && r.recal_n3 !== undefined ? "<br>recalibrated on 3 target individuals: " + fmt(r.recal_n3) : ""}` : "pending"));
  const traces = [
    { ...bar(x, acc, { name: "accuracy", slot: 1, hover: "%{customdata}<extra>accuracy</extra>" }), customdata: hoverAcc, error_y: err("accuracy", "accuracy_ci", "accuracy_sd") },
    { ...bar(x, cov, { name: `coverage of the 90 % set (${state.variant})`, slot: 2, hover: "%{customdata}<extra>coverage</extra>" }), customdata: hoverCov, error_y: err("coverage", "coverage_ci", "coverage_sd") },
  ];
  const ref = refLine(0.9, "1 − α = 0.90");
  const annotations = [...ref.annotations];
  // value labels above the whisker, not on it: grouped bars sit at ±0.2 category units from the tick
  const top = (key, ciKey, sdKey) => rows.map((r) => (r && !r.pending ? (r[ciKey] ? r[ciKey][1] : r[key] + (r[sdKey] || 0)) : null));
  [[acc, top("accuracy", "accuracy_ci", "accuracy_sd"), -0.2], [cov, top("coverage", "coverage_ci", "coverage_sd"), 0.2]].forEach(([vals, tops, dx]) => {
    vals.forEach((v, i) => { if (v !== null) annotations.push({ x: i + dx, y: Math.max(v, tops[i] ?? v) + 0.01, xref: "x", yref: "y", text: fmt(v), showarrow: false, yanchor: "bottom", font: { color: t.ink2, size: 11 } }); });
  });
  rows.forEach((r, i) => {
    if (!r || r.pending) annotations.push({ x: x[i], y: 0.5, xref: "x", yref: "y", text: "pending:<br>" + ((r && r.reason) || "not run"), showarrow: false, font: { color: t.muted, size: 11 } });
  });
  const table = { columns: ["shift", "model", "variant", "calibration", "accuracy", "accuracy_ci", "accuracy_sd", "accuracy_refit", "accuracy_seen", "coverage", "coverage_ci", "coverage_sd", "empty", "wrong_non_empty", "set_size", "recal_n3", "n_samples", "n_samples_coverage", "n_individuals", "source"],
                  rows: rows.map((r, i) => (r ? { shift: RUNGS[i].short, model: r.model, variant: r.variant, calibration: r.calibration, accuracy: r.accuracy, accuracy_ci: r.accuracy_ci ? r.accuracy_ci.map((v) => v.toFixed(3)).join(" – ") : "", accuracy_sd: r.accuracy_sd, accuracy_refit: r.accuracy_refit, accuracy_seen: r.accuracy_seen, coverage: r.coverage,
                                                  coverage_ci: r.coverage_ci ? r.coverage_ci.map((v) => v.toFixed(3)).join(" – ") : "", coverage_sd: r.coverage_sd, empty: r.empty, wrong_non_empty: r.wrong_non_empty, set_size: r.set_size, recal_n3: r.recal_n3, n_samples: r.n_samples, n_samples_coverage: r.n_samples_coverage ?? r.n_samples, n_individuals: r.n_individuals,
                                                  source: r.pending ? r.reason : (r.source || []).join("; ") } : { shift: RUNGS[i].short })) };
  return { traces, layout: { yaxis: { range: [0, 1.12], title: { text: "fraction" }, tickformat: ".1f" }, xaxis: { tickfont: { size: 11 } }, shapes: ref.shapes, annotations, barmode: "group",
                            legend: { y: 1.14 }, margin: { t: 40, b: 60 } }, table };
}

/** The title states what the data show: which models keep the 90 % coverage on the training-state rung (point estimate ≥ 0.90, or interval reaching it); the rest of the sentence is the recalibration result. */
export function ladderTitle(H) {
  const keep = ["k20", "k50", "full"].filter((m) => { const r = pick(H.ladder, "train_control_test_trained", m, "marginal"); return r && !r.pending && (r.coverage >= 0.9 || (r.coverage_ci && r.coverage_ci[1] >= 0.9)); });
  const lose = ["k20", "k50", "full"].filter((m) => !keep.includes(m));
  const name = { k20: "the 20-gene panel", k50: "the 50-gene panel", full: "the all-gene model" };
  const list = (ms) => ms.map((m) => name[m]).join(ms.length === 2 ? " and " : ", ").replace(/, ([^,]*)$/, " and $1");
  // "well above chance" is literally true on every rung (chance is 1/19; the lowest rung is the human one)
  const head = "Accuracy stays well above chance at every rung; the 90 % guarantee ";
  const tail = " and is restored beyond it by recalibrating on three target animals (within species)";
  if (keep.length === 0) return head + "is restored beyond the calibration data by recalibrating on three target animals (within species)";
  if (keep.length === 3) return head + "holds within the study for every model" + tail;
  return head + `holds within the study for ${list(keep)}` + tail;
}

/** The same title split at the semicolon: a one-line head and a tail sentence (capitalised) for the subtitle. */
export function ladderTitleParts(H) {
  const full = ladderTitle(H);
  const i = full.indexOf("; ");
  const tail = full.slice(i + 2);
  return { head: full.slice(0, i), tail: tail.charAt(0).toUpperCase() + tail.slice(1) + "." };
}

/** Which models keep the 90 % coverage on the training-state rung and which do not (for the collapsed notes). */
export function ladderKeepNote(H) {
  const rows = ["k20", "k50", "full"].map((m) => [m, pick(H.ladder, "train_control_test_trained", m, "marginal")]).filter(([, r]) => r && !r.pending);
  const name = { k20: "the 20-gene panel", k50: "the 50-gene panel", full: "the all-gene model" };
  const ci = (c) => (c ? ` [${fmt(c[0])}, ${fmt(c[1])}]` : "");
  return rows.length ? " On the training-state rung the marginal coverage is " + rows.map(([m, r]) => `${fmt(r.coverage)}${ci(r.coverage_ci)} for ${name[m]}`).join(", ") + "; the title counts a model as keeping the guarantee when the point estimate or the upper end of its interval reaches 0.90." : "";
}

/** Mount the ladder with its controls into `container`. opts.full adds the calibration control (in-distribution rung). */
export async function mountLadder(container, H, opts = {}) {
  const state = { model: opts.model || "k20", variant: opts.variant || "marginal", calib: opts.calib || "pooled" };
  const ctl = el("div", { class: "controls" });
  let fig;
  async function rerender() {
    if (!fig) return;
    const b = ladderBuild(H, state);
    fig.traces = b.traces; fig.table = b.table;
    await window.Plotly.react(fig.chart, b.traces, { ...template(), ...b.layout }, CONFIG);
    if (!fig.tableWrap.hidden) fig.tableWrap.replaceChildren(tableFrom(b.table));
  }
  ctl.append(
    control("Model", segmented([["k20", "20 genes"], ["k50", "50 genes"], ["full", "all genes"]], state.model, (v) => { state.model = v; rerender(); }, "model")),
    control("Prediction set", segmented([["marginal", "marginal"], ["mondrian", "Mondrian"], ["floored", "floored Mondrian"]], state.variant, (v) => { state.variant = v; rerender(); }, "conformal variant")),
  );
  if (opts.full) ctl.append(control("In-distribution calibration", segmented([["pooled", "pooled vials"], ["one_per_animal", "one vial per animal"]], state.calib, (v) => { state.calib = v; rerender(); }, "calibration")));
  const opa = pick(H.ladder, "in_distribution", "full", "marginal", "one_per_animal");
  const tr20 = pick(H.ladder, "train_control_test_trained", "k20", "marginal");
  const sx20 = pick(H.ladder, "train_male_test_female", "k20", "marginal");
  fig = await figure(container, {
    title: opts.title || (opts.shortTitle ? ladderTitleParts(H).head : ladderTitle(H)),
    subtitle: (opts.shortTitle ? ladderTitleParts(H).tail + " " : "") + "Per rung: accuracy and coverage of α = 0.10 sets calibrated on the source; whiskers are fold sd in-distribution and 95 % cluster-bootstrap intervals elsewhere.",
    build: () => ladderBuild(H, state),
    source: "results/05_panels/TRNSCRPT/panel_curve.csv, results/04_baselines/TRNSCRPT/summary.csv, results/06_conformal/TRNSCRPT/coverage.csv (full model) and results/31_site_regen/06_conformal/TRNSCRPT/scores_*.csv (k20/k50, recomputed with the same design), results/08_shift/TRNSCRPT/shift_table.csv, results/12_bodymap/{age_shift_accuracy,conformal_transfer}.csv, results/13_gtex/{accuracy_overall,conformal_transfer}.csv",
    notShow: "why coverage falls: it falls through empty sets, not through wrong confident sets." + ladderKeepNote(H) + " Accuracy is balanced in-distribution and sample-weighted on mapped tissues elsewhere; the whiskers are sd over 5 folds in-distribution and 95 % cluster-bootstrap intervals over animals or donors on the shifts. Each rung pairs the accuracy and the coverage of the same models: in-distribution, the models fit on 18 animals per fold whose sets are calibrated on 22 (the headline 0.976 on the Fingerprint page is the 40-animal CV of the same panel size); on the BodyMap and GTEx rungs the accuracy bar is the panel fit on all 50 animals and the sets come from its 35-animal refit (the refit's own accuracy is in the hover). The held-out-sex accuracy counts the unseen sex-specific tissue's vials as wrong, while coverage is over seen-class vials (both n in the hover); in-distribution calibration is pooled vials unless switched (one vial per animal gives " + fmt(opa?.coverage) + " for the full model). The training-state rung fits the model on " + (tr20 ? tr20.n_source_animals - tr20.n_calibration_animals : "—") + " of the " + (tr20?.n_source_animals ?? "—") + " sedentary control animals, calibrates on the other " + (tr20?.n_calibration_animals ?? "—") + " and tests every vial of the " + (tr20?.n_individuals ?? "—") + " trained animals. On the held-out-sex and training-state rungs the Mondrian and floored sets are full sets (every tissue the model knows), which is why they read 1.0: the source calibration holds " + (sx20?.n_calibration_animals ?? "—") + " (sex) or at most " + (tr20?.n_calibration_animals ?? "—") + " (training state) vials per class, and with n ≤ " + (sx20?.n_calibration_animals ?? "—") + " scores the α = 0.10 rank ⌈(n + 1) × 0.9⌉ exceeds n, so every per-class threshold is +∞.",
    toolbar: ctl,
    height: "tall",
  });
  return { state, fig, rerender };
}
