// The transfer ladder — the site's signature figure, shared by the home and transfer pages.
import { el, fmt, segmented, control, tableFrom } from "./site.js";
import { figure, bar, refLine, tokens, template, CONFIG } from "./charts.js";

export const RUNGS = [
  { id: "in_distribution", label: "In-distribution<br>(held-out animals)", short: "In-distribution" },
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
  const accSd = rows.map((r) => (r && !r.pending && r.accuracy_sd ? r.accuracy_sd : 0));
  const covSd = rows.map((r) => (r && !r.pending && r.coverage_sd ? r.coverage_sd : 0));
  const hoverAcc = rows.map((r) => (r && !r.pending ? `accuracy ${fmt(r.accuracy)}${r.accuracy_sd ? " ± " + fmt(r.accuracy_sd) : ""}<br>n = ${r.n_samples} samples, ${r.n_individuals} individuals` : "pending"));
  const hoverCov = rows.map((r) => (r && !r.pending ? `coverage ${fmt(r.coverage)}${r.coverage_sd ? " ± " + fmt(r.coverage_sd) : ""}<br>empty sets ${r.empty === null || r.empty === undefined ? "—" : fmt(r.empty)}<br>mean set size ${fmt(r.set_size, 2)}${r.recal_n3 !== null && r.recal_n3 !== undefined ? "<br>recalibrated on 3 target individuals: " + fmt(r.recal_n3) : ""}` : "pending"));
  const traces = [
    { ...bar(x, acc, { name: "accuracy", slot: 1, sd: accSd.some((v) => v) ? accSd : null, text: acc.map((v) => (v === null ? "" : fmt(v))), hover: "%{customdata}<extra>accuracy</extra>" }), customdata: hoverAcc },
    { ...bar(x, cov, { name: `coverage of the 90 % set (${state.variant})`, slot: 2, sd: covSd.some((v) => v) ? covSd : null, text: cov.map((v) => (v === null ? "" : fmt(v))), hover: "%{customdata}<extra>coverage</extra>" }), customdata: hoverCov },
  ];
  const ref = refLine(0.9, "1 − α = 0.90");
  const annotations = [...ref.annotations];
  rows.forEach((r, i) => {
    if (!r || r.pending) annotations.push({ x: x[i], y: 0.5, xref: "x", yref: "y", text: "pending:<br>" + ((r && r.reason) || "not run"), showarrow: false, font: { color: t.muted, size: 11 } });
  });
  const table = { columns: ["shift", "model", "variant", "calibration", "accuracy", "accuracy_sd", "coverage", "coverage_sd", "empty", "set_size", "recal_n3", "n_samples", "n_individuals", "source"],
                  rows: rows.map((r, i) => (r ? { shift: RUNGS[i].short, model: r.model, variant: r.variant, calibration: r.calibration, accuracy: r.accuracy, accuracy_sd: r.accuracy_sd, coverage: r.coverage,
                                                  coverage_sd: r.coverage_sd, empty: r.empty, set_size: r.set_size, recal_n3: r.recal_n3, n_samples: r.n_samples, n_individuals: r.n_individuals,
                                                  source: r.pending ? r.reason : (r.source || []).join("; ") } : { shift: RUNGS[i].short })) };
  return { traces, layout: { yaxis: { range: [0, 1.08], title: { text: "fraction" }, tickformat: ".1f" }, xaxis: { tickfont: { size: 11 } }, shapes: ref.shapes, annotations, barmode: "group",
                            legend: { y: 1.14 }, margin: { t: 40, b: 60 } }, table };
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
  fig = await figure(container, {
    title: opts.title || "Accuracy survives every shift; the 90 % guarantee survives none of them",
    subtitle: "Balanced accuracy (in-distribution) or accuracy on mapped tissues (shifts), and the coverage of α = 0.10 prediction sets calibrated on the source, per shift. Bars: mean; whiskers: sd over 5 folds where the design has folds.",
    build: () => ladderBuild(H, state),
    source: "results/05_panels/TRNSCRPT/panel_curve.csv, results/04_baselines/TRNSCRPT/summary.csv, results/06_conformal/TRNSCRPT/coverage.csv (full model) and results/31_site_regen/06_conformal/TRNSCRPT/scores_*.csv (k20/k50, recomputed with the same design), results/08_shift/TRNSCRPT/shift_table.csv, results/12_bodymap/{age_shift_accuracy,conformal_transfer}.csv, results/13_gtex/{accuracy_overall,conformal_transfer}.csv",
    notShow: "why coverage falls: it falls through empty sets, not through wrong confident sets. The held-out-sex rungs exist for the 20-gene panel and the full model only; the in-distribution rung uses pooled-vial calibration unless switched (one vial per animal gives " + fmt(opa?.coverage) + " for the full model). Mondrian and floored sets on the held-out sex are full 18-tissue sets: the source calibration has no vials of the unseen sex-specific tissue.",
    toolbar: ctl,
    height: "tall",
  });
  return { state, fig, rerender };
}
