import { mountChrome, loadJSON, el, fmt, segmented, control, tableFrom, callout } from "../site.js";
import { figure, bar, line, band, heatmap, strip, tokens, palette, organSystem, tissueLabel, template, CONFIG, jitter } from "../charts.js";

async function main() {
  await mountChrome("fingerprint.html");
  const [H, PC, SC, CM, GENES] = await Promise.all([loadJSON("data/headline.json"), loadJSON("data/panel_curve.json"), loadJSON("data/stable_core.json"), loadJSON("data/confusion_motrpac.json"), loadJSON("data/genes.json")]);
  const ex = H.extras;
  const k20 = PC.curve.find((r) => r.k === 20), k15 = PC.curve.find((r) => r.k === 15), k10 = PC.curve.find((r) => r.k === 10), k50 = PC.curve.find((r) => r.k === 50);
  const base = Object.fromEntries(PC.baselines.map((r) => [r.model, r]));
  document.getElementById("lede").replaceChildren(
    `A class-aware selector and a plain logistic regression identify 19 rat tissues at ${fmt(k20.roundrobin_mean)} ± ${fmt(k20.roundrobin_sd)} balanced accuracy with 20 genes, ${fmt(k50.roundrobin_mean)} with 50, and ${fmt(base.logreg_l2.balanced_accuracy_mean)} with all ${ex.motrpac_genes} genes. `,
    `Ten of the genes are stable across animal bootstraps; most are textbook markers.`,
  );

  // ---- panel curve ------------------------------------------------------------------------------------
  document.getElementById("p-curve").replaceChildren(
    `The accuracy-vs-size curve is flat from 15 genes on (${fmt(k15.roundrobin_mean)} at 15, ${fmt(k20.roundrobin_mean)} at 20, ${fmt(k50.roundrobin_mean)} at 50) and steep below: ${fmt(k10.roundrobin_mean)} at 10. `,
    "The step is the selector's construction: it picks one marker per tissue in turn, so at k = 10 nine tissues have no marker yet. ",
    `The univariate F-test ranks genes by a one-vs-rest statistic dominated by the easy tissues; its top 20 leave most tissues unrepresented and reach ${fmt(k20.fclassif_mean)}. `,
    `For reference, the tuned baselines on the same folds with all genes: nearest centroid ${fmt(base.centroid.balanced_accuracy_mean)}, L2 logistic regression ${fmt(base.logreg_l2.balanced_accuracy_mean)}, random forest ${fmt(base.rf.balanced_accuracy_mean)}.`,
  );
  await figure(document.getElementById("fig-curve"), {
    title: "Fifteen to twenty genes are enough; the F-test never gets there",
    subtitle: "Balanced accuracy vs panel size on a log axis: mean (line) ± sd (band) over 5 animal-grouped folds, with every fold's value as a dot; round-robin vs F-test selection, logreg_l2 on the selected genes.",
    build: () => {
      const c = PC.curve;
      const ks = c.map((r) => r.k);
      const rr = c.map((r) => r.roundrobin_mean), rs = c.map((r) => r.roundrobin_sd || 0);
      const fc = c.map((r) => r.fclassif_mean), fs = c.map((r) => r.fclassif_sd || 0);
      const p = palette();
      const dots = (key, slot, name) => ({ ...strip(c.flatMap((r) => r[key].map((_, i) => r.k * (1 + (i - 2) * 0.012))), c.flatMap((r) => r[key]), { color: p[slot - 1], size: 6, opacity: 0.6, name: name + " (per fold)", hover: "k = %{customdata}: %{y:.3f}<extra>" + name + ", one fold</extra>" }), customdata: c.flatMap((r) => r[key].map(() => r.k)), legendgroup: name, showlegend: false });
      return { traces: [
        band(ks, rr.map((v, i) => v - rs[i]), rr.map((v, i) => v + rs[i]), { slot: 1, legendgroup: "rr" }),
        band(ks, fc.map((v, i) => v - fs[i]), fc.map((v, i) => v + fs[i]), { slot: 2, legendgroup: "fc" }),
        dots("roundrobin_folds", 1, "round-robin"), dots("fclassif_folds", 2, "F-test"),
        line(ks, rr, { name: "round-robin selector", slot: 1, legendgroup: "rr", hover: "k = %{x}: %{y:.3f}<extra>round-robin, mean</extra>" }),
        line(ks, fc, { name: "F-test selector", slot: 2, legendgroup: "fc", hover: "k = %{x}: %{y:.3f}<extra>F-test, mean</extra>" }),
      ], layout: { xaxis: { type: "log", title: { text: "panel size k (genes)" }, tickvals: ks, ticktext: ks.map(String) }, yaxis: { range: [0, 1.05], title: { text: "balanced accuracy" } }, margin: { t: 40 } },
        table: { columns: ["k", "roundrobin_mean", "roundrobin_sd", "fclassif_mean", "fclassif_sd", "n_folds", "n_train_animals", "n_test_animals"], rows: c } };
    },
    source: "results/05_panels/TRNSCRPT/panel_curve.csv, results/05_panels/TRNSCRPT/panel_curve_fclassif.csv",
    notShow: `the certified panel size, which is larger (Methods): with 22 calibration animals a 10 %/90 % certificate certifies k = 20 in ${(100 * ex.cert_k20_frac).toFixed(0)} % of repeated splits and no size at all in ${(100 * ex.cert_none_frac).toFixed(0)} % (results/06_conformal/TRNSCRPT/certificate_distribution.csv, certificate_alpha_delta_grid.csv).`,
  });

  // ---- stable core ------------------------------------------------------------------------------------
  const gi = Object.fromEntries(GENES.genes.map((g) => [g.id, g]));
  const core = SC.core;
  document.getElementById("p-core").replaceChildren(
    `Selection was repeated on 50 bootstrap resamples of animals; ${core.length} genes were chosen in at least 80 % of them (Source: results/05_panels/TRNSCRPT/candidate_panel_annotated.csv). `,
    `They cover ${new Set(core.map((r) => r.marker_tissue)).size} tissues with single markers such as Umod (kidney), Cyp21a1 (adrenal) and Pmch (hypothalamus). ${core.filter((r) => r.risk_T7_regulated).length} carry a training-regulated flag and ${core.filter((r) => r.risk_qc_correlated).length} a QC-correlation flag: annotations, not exclusions.`,
  );
  const coreRows = core.map((r) => { const g = gi[r.feature_ID] || {}; return {
    gene: r.gene_symbol, tissue: tissueLabel(r.marker_tissue).split(" · ").pop(), freq: r.selection_frequency, effect: r.effect_size, "runner-up": r.next_highest_tissue, OvR: r.ovr_score,
    "r mRNA": r.r_pct_mrna_in_marker_tissue, regulated: r.risk_T7_regulated ? "yes" : "no", "QC flag": r.risk_qc_correlated ? "yes" : "no",
    BodyMap: g.fails_bodymap === null || g.fails_bodymap === undefined ? "no organ" : g.fails_bodymap ? "fails" : g.weakened_bodymap ? "weakened" : "holds",
    GTEx: g.fails_gtex === null || g.fails_gtex === undefined ? "not testable" : g.fails_gtex ? "fails" : g.weakened_gtex ? "weakened" : "holds" }; });
  const f2 = (v) => v.toFixed(2);
  document.getElementById("tbl-core").replaceChildren(tableFrom({ columns: Object.keys(coreRows[0]), rows: coreRows, format: { freq: f2, effect: f2, OvR: f2, "r mRNA": f2 } }),
    el("p", { class: "small" }, "Columns: freq = selection frequency over the 50 bootstraps; effect = log2 CPM above the runner-up tissue; OvR = one-vs-rest score; r mRNA = correlation with the library mRNA fraction in the marker tissue; regulated = training-regulated in the marker tissue; QC flag = QC-correlated. Sources: results/05_panels/TRNSCRPT/candidate_panel_annotated.csv (selection frequency, effect size, one-vs-rest score, QC correlation, flags); results/12_bodymap/panel_gene_check.csv and results/13_gtex/panel_gene_check.csv (holds / weakened / fails on the external targets; the GTEx panel is re-selected in ortholog space, so a core gene can be untested there)."));

  // ---- per-tissue accuracy and confusion --------------------------------------------------------------
  const pt20 = CM.confusion.k20.per_tissue_accuracy, ptFull = CM.confusion.full.per_tissue_accuracy;
  const tissues = Object.keys(pt20).sort((a, b) => pt20[a] - pt20[b]);
  document.getElementById("p-tissue").replaceChildren(
    `Pooled over the five folds, the 20-gene panel is right on every vial of ${Object.values(pt20).filter((v) => v === 1).length} tissues and drops below 0.9 on ${tissues.filter((t) => pt20[t] < 0.9).map((t) => tissueLabel(t)).join(", ") || "none"}. `,
    `Its ${CM.confusable_k20.reduce((a, r) => a + r.count, 0)} errors are mostly ${CM.confusable_k20[0].true} → ${CM.confusable_k20[0].predicted} (${CM.confusable_k20[0].count}) and the two skeletal muscles for each other. With all genes the same folds leave ${CM.confusable_full.reduce((a, r) => a + r.count, 0)} errors.`,
  );
  await figure(document.getElementById("fig-tissue"), {
    title: "Accuracy per tissue: brown fat, the two skeletal muscles and vena cava are the hard ones",
    subtitle: "Fraction of each tissue's 899 pooled out-of-fold vials called correctly, 20-gene panel vs all genes (logreg_l2), 5 animal-grouped folds.",
    build: () => ({
      traces: [bar(tissues.map((t) => pt20[t]), tissues.map((t) => tissueLabel(t)), { horizontal: true, name: "20-gene panel", slot: 1, hover: "%{y}: %{x:.3f}<extra>20 genes</extra>" }),
               bar(tissues.map((t) => ptFull[t] ?? null), tissues.map((t) => tissueLabel(t)), { horizontal: true, name: "all genes", slot: 4, hover: "%{y}: %{x:.3f}<extra>all genes</extra>" })],
      layout: { barmode: "group", xaxis: { range: [0, 1.02], title: { text: "per-tissue accuracy (pooled out-of-fold vials)" } }, yaxis: { automargin: true, tickfont: { size: 11 } }, margin: { t: 40, l: 10 }, bargap: 0.25, legend: { y: 1.1 } },
      table: { columns: ["tissue", "k20", "full"], rows: tissues.map((t) => ({ tissue: t, k20: pt20[t], full: ptFull[t] })) },
    }),
    source: "results/05_panels/TRNSCRPT/confusion_k20.csv, results/04_baselines/TRNSCRPT/confusion_logreg_l2.csv (diagonal / row sum)",
    notShow: "per-fold spread per tissue (10 test animals per fold, ~50 vials per tissue in total).", height: "tall",
  });
  const cState = { key: "k20" };
  const ctl = el("div", { class: "controls" }, [control("Panel", segmented([["k20", "20 genes"], ["k30", "30 genes"], ["full", "all genes"]], cState.key, (v) => { cState.key = v; cf.rerender(); }, "panel"))]);
  const confBuild = () => {
    const c = CM.confusion[cState.key];
    const rowsum = c.counts.map((r) => r.reduce((a, b) => a + b, 0));
    const z = c.counts.map((r, i) => r.map((v) => (rowsum[i] ? v / rowsum[i] : 0)));
    return { traces: [heatmap(z, c.cols, c.rows, { zmin: 0, zmax: 1, ztitle: "row fraction", hover: "%{y} → %{x}: %{z:.2f}<extra></extra>", text: z.map((r) => r.map((v) => (v > 0 && v < 1 ? v.toFixed(2) : ""))), texttemplate: "%{text}" })],
             layout: { xaxis: { title: { text: "predicted" }, tickangle: -60, tickfont: { size: 10 } }, yaxis: { title: { text: "true tissue" }, autorange: "reversed", tickfont: { size: 10 }, automargin: true }, margin: { t: 20, b: 90 } },
             table: { columns: ["true", ...c.cols], rows: c.rows.map((r, i) => ({ true: r, ...Object.fromEntries(c.cols.map((cc, j) => [cc, c.counts[i][j]])) })) } };
  };
  const cf = await figure(document.getElementById("fig-confusion"), {
    title: "The confusions are anatomical neighbours: vena cava → brown fat, gastrocnemius ↔ vastus lateralis",
    subtitle: "Row-normalised confusion matrix over the pooled out-of-fold vials; off-diagonal fractions labelled.",
    build: confBuild, toolbar: ctl, source: "results/05_panels/TRNSCRPT/confusion_k20.csv, confusion_k30.csv, results/04_baselines/TRNSCRPT/confusion_logreg_l2.csv", notShow: "counts (the data table holds them).", height: "tall",
  });
  cf.rerender = async () => { const b = confBuild(); cf.traces = b.traces; cf.table = b.table; await window.Plotly.react(cf.chart, b.traces, { ...template(), ...b.layout }, CONFIG); };

  // ---- hard tissues -------------------------------------------------------------------------------------
  const stab = SC.stability_k20;
  const best = {};
  for (const r of stab) if (!best[r.marker_tissue] || r.effect_size > best[r.marker_tissue].effect_size) best[r.marker_tissue] = r;
  const hard = ["BAT", "SKM-GN", "SKM-VL", "VENACV"].map((t) => best[t]).filter(Boolean);
  const liver = stab.filter((r) => r.marker_tissue === "LIVER").sort((a, b) => b.ovr_score - a.ovr_score);
  document.getElementById("p-hard").replaceChildren(
    "Every hard tissue is hard for the same reason: its best single marker sits only a little above its anatomical neighbour. ",
    hard.map((r) => `${r.gene_symbol} for ${tissueLabel(r.marker_tissue)} is ${fmt(r.effect_size, 2)} log2 CPM above ${r.next_highest_tissue}`).join("; "),
    `. Compare the easy tissues, whose best marker sits ${fmt(Math.min(...["KIDNEY", "ADRNL", "TESTES"].map((t) => best[t].effect_size)), 1)}–${fmt(Math.max(...["KIDNEY", "ADRNL", "TESTES"].map((t) => best[t].effect_size)), 1)} log2 CPM above everything else (${["KIDNEY", "ADRNL", "TESTES"].map((t) => best[t].gene_symbol).join(", ")}). `,
    "The muscles share a fibre programme, vena cava carries perivascular brown fat, and the panel has one gene per tissue to tell them apart.",
  );
  const tissuesAll = Object.keys(best).sort((a, b) => best[a].effect_size - best[b].effect_size);
  await figure(document.getElementById("fig-markers"), {
    title: "Marker strength per tissue: the hard tissues have weak single markers",
    subtitle: "Effect size of each tissue's strongest marker among the 51 stability-selected genes (not necessarily the one the k = 20 panel picked in a given fold): mean log2 CPM in the tissue minus the highest mean of any other tissue; hover for the gene and its runner-up tissue.",
    build: () => {
      const t = tokens(); const p = palette();
      const x = tissuesAll.map((tt) => best[tt].effect_size);
      return { traces: [{ ...bar(x, tissuesAll.map((tt) => tissueLabel(tt)), { horizontal: true, name: "best marker effect", slot: 1, text: tissuesAll.map((tt) => best[tt].gene_symbol), textposition: "outside", hover: "%{customdata}<extra></extra>" }),
                          customdata: tissuesAll.map((tt) => `${tissueLabel(tt)}: ${best[tt].gene_symbol}<br>${fmt(best[tt].effect_size, 2)} log2 CPM above ${best[tt].next_highest_tissue}<br>selected in ${fmt(best[tt].selection_frequency, 2)} of bootstraps`),
                          marker: { color: tissuesAll.map((tt) => (["BAT", "SKM-GN", "SKM-VL", "VENACV"].includes(tt) ? p[1] : p[0])), line: { color: t.surface, width: 2 }, cornerradius: 4 } }],
               layout: { xaxis: { title: { text: "log2 CPM above the next-highest tissue" }, range: [0, Math.max(...x) * 1.25] }, yaxis: { automargin: true, tickfont: { size: 11 } }, margin: { t: 20, l: 10 }, showlegend: false, bargap: 0.3 },
               table: { columns: ["tissue", "gene", "effect_size", "next_highest_tissue", "selection_frequency", "ovr_score"], rows: tissuesAll.map((tt) => ({ tissue: tt, gene: best[tt].gene_symbol, effect_size: best[tt].effect_size, next_highest_tissue: best[tt].next_highest_tissue, selection_frequency: best[tt].selection_frequency, ovr_score: best[tt].ovr_score })) } };
    },
    source: "results/05_panels/TRNSCRPT/stability_k20_annotated.csv (51 genes ever selected in the bootstrap run; orange = the hard tissues)",
    notShow: "tissues whose markers were never among the 51 (none: every tissue has at least one); the effect is on the log scale, so 2 log2 CPM is a four-fold difference.", height: "tall",
  });
  document.getElementById("p-liver").replaceChildren(callout("note", "Liver's instability is redundancy, not weakness",
    `Liver has ${liver.length} near-equivalent candidates (${liver.map((r) => `${r.gene_symbol} ${fmt(r.ovr_score, 1)}`).join(", ")}; one-vs-rest scores from results/05_panels/TRNSCRPT/stability_k20_annotated.csv), so bootstraps split the vote: `
    + `${liver[0].gene_symbol} is selected in ${fmt(liver[0].selection_frequency, 2)} of resamples, ${liver[1]?.gene_symbol} in ${fmt(liver[1]?.selection_frequency, 2)}. Any of them identifies liver; the stability rule simply cannot pick one.`));
}

main().catch((e) => { console.error(e); document.getElementById("lede").textContent = "Failed to load: " + e.message; });
