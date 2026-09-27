import { mountChrome, loadJSON, el, fmt, segmented, control, callout, tableFrom } from "../site.js";
import { figure, bar, line, heatmap, strip, tokens, palette, hexAlpha, template, CONFIG, tissueLabel, mean } from "../charts.js";
import { mountLadder, pick, RUNGS } from "../ladder.js";

const MODELS = ["k20", "k50", "full"];
const MODEL_LABEL = { k20: "20 genes", k50: "50 genes", full: "all genes" };

async function main() {
  await mountChrome("transfer.html");
  const [H, S, BM, GT, REP, GENES] = await Promise.all([
    loadJSON("data/headline.json"), loadJSON("data/shift.json"), loadJSON("data/bodymap.json"), loadJSON("data/gtex.json"),
    loadJSON("data/representations.json"), loadJSON("data/genes.json"),
  ]);
  const ex = H.extras;
  const bm20 = pick(H.ladder, "different_lab", "k20", "marginal"), gt20 = pick(H.ladder, "different_species", "k20", "marginal"), gtfull = pick(H.ladder, "different_species", "full", "marginal");
  const tr20 = pick(H.ladder, "train_control_test_trained", "k20", "marginal"), trFull = pick(H.ladder, "train_control_test_trained", "full", "marginal");
  document.getElementById("lede").replaceChildren(
    `Fit on the ${tr20.n_source_animals} sedentary control animals alone, the 20-gene panel names the tissue of the ${tr20.n_individuals} trained animals at ${fmt(tr20.accuracy)} with coverage ${fmt(tr20.coverage)}: training state is not a shift the panel notices (the all-gene model reaches ${fmt(trFull.accuracy)} but its coverage slips to ${fmt(trFull.coverage)}). `,
    `Fit on MoTrPAC, it names ${fmt(bm20.accuracy)} of adult BodyMap organs and ${fmt(gt20.accuracy)} of human GTEx samples (${fmt(gtfull.accuracy)} with all genes). `,
    `The 90 % guarantee delivers ${fmt(bm20.coverage)} and ${fmt(gt20.coverage)} with the source calibration. Three target animals repair it within species; three donors do not repair it across species.`,
  );

  // ---- ladder --------------------------------------------------------------------------------------
  await mountLadder(document.getElementById("fig-ladder"), H, { full: true });

  // ---- empty-set stack -------------------------------------------------------------------------------
  const emptyState = { model: "k20" };
  const rowsFor = (model) => RUNGS.map((r) => [r.short, pick(H.ladder, r.id, model, "marginal")]).filter(([, r]) => r && !r.pending);
  const e20 = rowsFor("k20");
  document.getElementById("p-empty").replaceChildren(
    "A conformal set fails in two ways: it can be non-empty and wrong, or it can be empty. Across every shift the loss is almost entirely empty sets. ",
    `With the 20-gene panel and the source calibration, ${e20.map(([l, r]) => `${l.toLowerCase()} ${fmt(r.empty, 2)}`).join(", ")} of samples receive no tissue at all, while non-empty wrong sets stay rare. `,
    "That is the desired failure mode of a calibrated model under shift: it says “I do not know” rather than naming the wrong tissue with confidence.",
  );
  const emptyCtl = el("div", { class: "controls" }, [control("Model", segmented(MODELS.map((m) => [m, MODEL_LABEL[m]]), emptyState.model, (v) => { emptyState.model = v; emptyFig.rerender(); }, "model"))]);
  const emptyBuild = () => {
    const t = tokens(); const p = palette();
    const rows = rowsFor(emptyState.model);
    const y = rows.map(([l]) => l);
    const covered = rows.map(([, r]) => r.coverage);
    const empty = rows.map(([, r]) => (r.empty === null || r.empty === undefined ? null : r.empty));
    const wrong = rows.map(([, r], i) => (r.wrong_non_empty !== null && r.wrong_non_empty !== undefined ? r.wrong_non_empty : (empty[i] === null ? null : Math.max(0, 1 - r.coverage - empty[i]))));
    const mk = (name, vals, color, pattern) => ({ type: "bar", orientation: "h", y, x: vals, name, marker: { color, line: { color: t.surface, width: 2 }, pattern: pattern ? { shape: "/", fgcolor: t.ink2, bgcolor: t.surface, size: 6, solidity: 0.35 } : undefined },
                                                  hovertemplate: "%{y}: %{x:.3f}<extra>" + name + "</extra>", text: vals.map((v) => (v !== null && v >= 0.12 ? fmt(v, 2) : "")), textposition: "inside", insidetextanchor: "middle", textfont: { color: pattern ? t.ink : "#fff", size: 12 }, cliponaxis: false });
    return { traces: [mk("true tissue in the set (covered)", covered, p[0]), mk("non-empty but wrong", wrong, p[7]), mk("empty set (abstains)", empty, t.grid, true)],
             layout: { barmode: "stack", xaxis: { range: [0, 1], title: { text: "fraction of test samples" } }, yaxis: { autorange: "reversed", automargin: true }, margin: { t: 40, l: 10 }, bargap: 0.4, legend: { y: 1.18 } },
             table: { columns: ["shift", "model", "covered", "wrong_non_empty", "empty", "n_samples"], rows: rows.map(([l, r], i) => ({ shift: l, model: r.model, covered: covered[i], wrong_non_empty: wrong[i], empty: empty[i], n_samples: r.n_samples })) } };
  };
  const emptyFig = await figure(document.getElementById("fig-empty"), {
    title: "What replaces coverage is abstention, not confident error",
    subtitle: "Fraction of test samples whose α = 0.10 marginal set (source-calibrated) contains the true tissue, is non-empty but wrong, or is empty; per shift.",
    build: emptyBuild, toolbar: emptyCtl,
    source: "results/31_site_regen/06_conformal/TRNSCRPT/scores_*.csv (in-distribution, recomputed), results/08_shift/TRNSCRPT/shift_table.csv, results/12_bodymap/conformal_transfer.csv, results/13_gtex/conformal_transfer.csv",
    notShow: "Mondrian and floored sets, which replace abstention with larger sets (ladder controls above). All three fractions are over the same samples: seen-class vials for the held-out sex, mapped samples for BodyMap and GTEx.",
  });
  emptyFig.rerender = async () => { const b = emptyBuild(); emptyFig.traces = b.traces; emptyFig.table = b.table; await window.Plotly.react(emptyFig.chart, b.traces, { ...template(), ...b.layout }, CONFIG); };

  // ---- recalibration cost ----------------------------------------------------------------------------
  const rb = BM.recalibration, rg = GT.recalibration, rs = S.recalibration;
  const inf = GT.recal_infinite_draws;
  const infK20n3 = inf && inf.k20_n3 ? inf.k20_n3.frac_infinite : null;
  document.getElementById("p-recal").replaceChildren(
    `Recalibrating the threshold on a few target individuals is the standard repair. On BodyMap adults, three animals lift coverage of the 20-gene sets from ${fmt(rb.find((r) => r.model === "k20" && r.n_recal === 3).coverage_source_cal_same_test)} to ${fmt(rb.find((r) => r.model === "k20" && r.n_recal === 3).coverage_recalibrated)} at ${fmt(rb.find((r) => r.model === "k20" && r.n_recal === 3).set_size_recalibrated, 2)} tissues per set: the guarantee is back and the sets are still singletons. `,
    `On GTEx, three donors lift coverage to ${fmt(rg.find((r) => r.model === "k20" && r.n_recal === 3).coverage_recalibrated)}, but the sets hold ${fmt(rg.find((r) => r.model === "k20" && r.n_recal === 3).set_size_recalibrated, 1)} of 19 tissues; five donors bring ${fmt(rg.find((r) => r.model === "k20" && r.n_recal === 5).coverage_recalibrated)} at ${fmt(rg.find((r) => r.model === "k20" && r.n_recal === 5).set_size_recalibrated, 2)} tissues. A set that names most of the body is coverage without information.`,
  );
  document.getElementById("callout-recal").replaceChildren(callout("caveat", "“Three donors restore coverage” is false across species", [
    `The recalibrated GTEx coverage of ${fmt(rg.find((r) => r.model === "k20" && r.n_recal === 3).coverage_recalibrated)} is trivial, not earned. Three donors contribute few mapped samples, so in ${infK20n3 === null ? "pending" : (100 * infK20n3).toFixed(0) + " %"} of the 20 draws the rank ⌈(n + 1)(1 − α)⌉ exceeds n, the threshold is +∞ and every set is the full 19-tissue set (Source: results/31_site_regen/13_gtex/recal_thresholds.csv, per-draw thresholds). `
    + `Read the coverage together with the set size: ${fmt(rg.find((r) => r.model === "k20" && r.n_recal === 3).set_size_recalibrated, 1)} tissues per set. Within species (BodyMap) the same three-animal recalibration gives one tissue per set.`,
  ]));
  await figure(document.getElementById("fig-recal"), {
    title: "Recalibration repairs coverage within species at one tissue per set; across species it repairs it by returning most of the body",
    subtitle: "Coverage vs mean set size after recalibrating the α = 0.10 threshold on n target individuals (● n = 3, ◆ n = 5), one point per target × model × n, mean over 20 draws (BodyMap: all combinations). Shaded: coverage ≥ 0.90 at ≤ 1.5 tissues per set.",
    build: () => {
      const t = tokens(); const p = palette();
      const mk = (rows, name, slot, labelKey) => ({
        type: "scatter", mode: "markers+text", name,
        x: rows.map((r) => r.set_size_recalibrated), y: rows.map((r) => r.coverage_recalibrated),
        text: rows.map((r) => `${MODEL_LABEL[r.model] || r.arm || r.model} · n=${r.n_recal_animals ?? r.n_recal}`), textposition: "top center", textfont: { color: t.ink2, size: 10 },
        marker: { color: p[slot - 1], size: 10, symbol: rows.map((r) => ((r.n_recal_animals ?? r.n_recal) === 5 ? "diamond" : "circle")), line: { color: t.surface, width: 2 } },
        customdata: rows.map((r) => `${name}<br>${r.split ? r.split.replace(/_/g, " ") + " · " : ""}${MODEL_LABEL[r.model] || r.arm}, n = ${r.n_recal_animals ?? r.n_recal}<br>coverage ${fmt(r.coverage_recalibrated ?? r.coverage_target_recalibrated)} (source-calibrated ${fmt(r.coverage_source_cal_same_test ?? r.coverage_target_source_cal_same_test)})<br>mean set size ${fmt(r.set_size_recalibrated, 2)}`),
        hovertemplate: "%{customdata}<extra></extra>",
      });
      const shiftRows = rs.filter((r) => r.arm === "panel_k20" || r.arm === "full").map((r) => ({ ...r, model: r.arm === "panel_k20" ? "k20" : "full", coverage_recalibrated: r.coverage_target_recalibrated, coverage_source_cal_same_test: r.coverage_target_source_cal_same_test }));
      const traces = [mk(rb, "rat BodyMap adults", 1), mk(rg, "human GTEx", 2), mk(shiftRows, "MoTrPAC shifts (sex, time point, training state)", 3)];
      traces.forEach((tr) => { tr.mode = "markers"; tr.text = undefined; });
      const g3 = rg.filter((r) => r.n_recal === 3), g5 = rg.filter((r) => r.n_recal === 5);
      const clusterNotes = [
        { x: Math.max(...rb.map((r) => r.set_size_recalibrated)) + 0.55, xanchor: "left", y: mean(rb.map((r) => r.coverage_recalibrated)), yanchor: "middle", align: "left", text: `BodyMap, 3 or 5 animals:<br>${fmt(Math.min(...rb.map((r) => r.set_size_recalibrated)), 2)}–${fmt(Math.max(...rb.map((r) => r.set_size_recalibrated)), 2)} tissues per set` },
        { x: mean(g3.map((r) => r.set_size_recalibrated)), y: Math.max(...g3.map((r) => r.coverage_recalibrated)) + 0.012, text: `GTEx, 3 donors:<br>${fmt(Math.min(...g3.map((r) => r.set_size_recalibrated)), 1)}–${fmt(Math.max(...g3.map((r) => r.set_size_recalibrated)), 1)} tissues per set`, yanchor: "bottom" },
        { x: mean(g5.map((r) => r.set_size_recalibrated)), y: Math.min(...g5.map((r) => r.coverage_recalibrated)) - 0.012, text: "GTEx, 5 donors", yanchor: "top" },
      ].map((a) => ({ ...a, xref: "x", yref: "y", showarrow: false, font: { color: t.ink2, size: 11 } }));
      return { traces,
               layout: { xaxis: { title: { text: "mean set size (tissues per set)" }, range: [0.5, 13] }, yaxis: { title: { text: "coverage after recalibration" }, range: [0.8, 1.02] },
                         shapes: [{ type: "rect", x0: 0.5, x1: 1.5, y0: 0.9, y1: 1.02, xref: "x", yref: "y", fillcolor: hexAlpha(p[0], 0.08), line: { width: 0 } }],
                         annotations: [{ x: 1.0, y: 1.015, xref: "x", yref: "y", text: "target: guaranteed and informative", showarrow: false, font: { color: t.ink2, size: 11 } }, ...clusterNotes], legend: { y: 1.14 }, margin: { t: 40 } },
               table: { columns: ["target", "model", "n_recal", "coverage_recalibrated", "coverage_source_cal_same_test", "set_size_recalibrated", "draws"],
                        rows: [...rb.map((r) => ({ target: "BodyMap", ...r })), ...rg.map((r) => ({ target: "GTEx", ...r })), ...shiftRows.map((r) => ({ target: r.split, model: r.model, n_recal: r.n_recal_animals, coverage_recalibrated: r.coverage_recalibrated, coverage_source_cal_same_test: r.coverage_source_cal_same_test, set_size_recalibrated: r.set_size_recalibrated, draws: r.repeats }))] } };
    },
    source: "results/12_bodymap/recalibration.csv, results/13_gtex/recalibration.csv, results/08_shift/TRNSCRPT/shift_recalibration.csv",
    notShow: "the spread over draws (the CSVs hold means; the shift file holds an sd); the GTEx points mix finite and infinite thresholds across draws.",
    height: "tall",
  });

  // ---- age ---------------------------------------------------------------------------------------------
  const ages = BM.age_accuracy;
  const a2 = ages.find((r) => r.stage_weeks === 2), a104 = ages.find((r) => r.stage_weeks === 104);
  const eb = await loadJSON("data/expr_bodymap.json");
  document.getElementById("p-age").replaceChildren(
    `The panel was learned on six-month-old rats. In the BodyMap it names ${fmt(a2.k20)} of organs correctly at two weeks of age and ${fmt(a104.k20)} at two years, against ${fmt(ages.find((r) => r.stage_weeks === 21).k20)} in young adults. `,
    `The juvenile failures are explained gene by gene: the panel's only testis marker, Pgk2, is a spermatid glycolytic isozyme (${fmt(ex.pgk2_testes_2wk, 2)} log2 CPM in two-week testes vs ${fmt(ex.pgk2_testes_21wk, 2)} in adults), so a juvenile testis carries no testis evidence and is called by its next-best signal; the juvenile spleen, still erythropoietic, expresses the blood marker Hbq1b and is called blood.`,
  );
  await figure(document.getElementById("fig-age"), {
    title: "Accuracy dips in juveniles and old animals; the full model is more robust than the panel",
    subtitle: "Accuracy on the 9 mapped BodyMap organs by age (sample-weighted; super-class scoring for muscle and brain), per model. Whiskers: 95 % intervals over the animals of each age (cluster bootstrap; exact binomial where every animal is perfect); n in the hover.",
    build: () => ({
      traces: MODELS.map((m, i) => {
        const tr = line(ages.map((r) => String(r.stage_weeks)), ages.map((r) => r[m]), { name: MODEL_LABEL[m], slot: i + 1, hover: "%{customdata}<extra>" + MODEL_LABEL[m] + "</extra>" });
        tr.customdata = ages.map((r) => `${r.stage_weeks} weeks: ${fmt(r[m])}${r[`${m}_ci`] ? ` [${fmt(r[`${m}_ci`][0])}, ${fmt(r[`${m}_ci`][1])}]` : ""}<br>n = ${r.n_samples ?? "—"} samples, ${r.n_animals ?? "—"} animals`);
        if (ages.every((r) => r[`${m}_ci`])) tr.error_y = { type: "data", symmetric: false, array: ages.map((r) => r[`${m}_ci`][1] - r[m]), arrayminus: ages.map((r) => r[m] - r[`${m}_ci`][0]), visible: true, thickness: 1.2, width: 3 };
        return tr;
      }),
      layout: { xaxis: { title: { text: "age (weeks)" }, type: "category" }, yaxis: { range: [0.5, 1.02], title: { text: "accuracy" } }, margin: { t: 40 } },
      table: { columns: ["stage_weeks", "n_samples", "n_animals", "k20", "k20_ci", "k50", "k50_ci", "full", "full_ci"],
               rows: ages.map((r) => ({ stage_weeks: r.stage_weeks, n_samples: r.n_samples, n_animals: r.n_animals, ...Object.fromEntries(MODELS.flatMap((m) => [[m, r[m]], [`${m}_ci`, r[`${m}_ci`] ? r[`${m}_ci`].map((v) => v.toFixed(3)).join(" – ") : ""]])) })) },
    }),
    source: "results/12_bodymap/age_shift_accuracy.csv; intervals and n recomputed from results/31_site_regen/12_bodymap/scores_target_probs.csv",
    notShow: "per-organ detail (results/12_bodymap/accuracy_by_organ.csv): at 2 weeks the panel calls testes SKM-VL and spleen BLOOD in every sample.",
    height: "short",
  });
  // developmental markers: mean log2 CPM vs age in the relevant organ, panel markers solid, confirmatory genes dashed
  const devGenes = [
    { sym: "Pgk2", organ: "Testes", panel: true }, { sym: "Mybph", organ: "Testes", panel: true, note: "the SKM-VL marker, higher in juvenile testis" }, { sym: "Prm1", organ: "Testes", panel: false }, { sym: "Tnp1", organ: "Testes", panel: false }, { sym: "Tnp2", organ: "Testes", panel: false }, { sym: "Acrv1", organ: "Testes", panel: false },
    { sym: "Hbq1b", organ: "Spleen", panel: true }, { sym: "Fcrl5", organ: "Spleen", panel: true }, { sym: "Klf1", organ: "Spleen", panel: false }, { sym: "Alas2", organ: "Spleen", panel: false }, { sym: "Hbb", organ: "Spleen", panel: false },
  ];
  const symIdx = Object.fromEntries(eb.symbols.map((s, i) => [s, i]));
  const AGES = [2, 6, 21, 104];
  const devFigure = async (organ, id, title) => figure(document.getElementById(id), {
    title,
    subtitle: `Mean log2 CPM by age in BodyMap ${organ.toLowerCase()} (8 samples per age; 4 for testes and uterus): panel markers solid with large markers, confirmatory developmental genes dotted.`,
    build: () => {
      const t = tokens();
      const traces = []; const missing = [];
      let slot = 0;
      for (const g of devGenes.filter((g) => g.organ === organ)) {
        const i = symIdx[g.sym];
        if (i === undefined) { missing.push(g.sym); continue; }
        const vals = AGES.map((a) => mean(eb.samples.map((s, k) => (s.organ === organ && s.age_weeks === a ? eb.values[i][k] : NaN))));
        slot += 1;
        traces.push(line(AGES.map(String), vals, { name: g.sym + (g.panel ? " (panel)" : ""), slot: ((slot - 1) % 8) + 1, dash: g.panel ? "solid" : "dot", size: g.panel ? 10 : 6, hover: "%{x} weeks: %{y:.2f}<extra>" + g.sym + "</extra>" }));
      }
      return { traces, layout: { xaxis: { type: "category", title: { text: "age (weeks)" } }, yaxis: { title: { text: "log2 CPM" } }, legend: { y: 1.16, font: { size: 11 } }, margin: { t: 50 },
                                 annotations: missing.length ? [{ xref: "paper", yref: "paper", x: 0, y: -0.25, text: "absent from the export: " + missing.join(", "), showarrow: false, font: { color: t.muted, size: 10 } }] : [] },
               table: { columns: ["gene", "panel_marker", "age_2", "age_6", "age_21", "age_104"], rows: devGenes.filter((g) => g.organ === organ && symIdx[g.sym] !== undefined).map((g) => { const i = symIdx[g.sym]; const m = (a) => mean(eb.samples.map((s, k) => (s.organ === organ && s.age_weeks === a ? eb.values[i][k] : NaN))); return { gene: g.sym, panel_marker: g.panel ? "yes" : "no", age_2: m(2), age_6: m(6), age_21: m(21), age_104: m(104) }; }) } };
    },
    source: "site/data/expr_bodymap.json (log2 CPM from data/external/bodymap_counts.csv); results/12_bodymap/juvenile_marker_check.csv holds the same means for the 20 panel genes",
    notShow: "protein-level confirmation; the genes are read as transcripts only.",
  });
  await devFigure("Testes", "fig-dev-testes", "The juvenile testis has no spermatid transcripts yet: Pgk2 and the other post-meiotic genes appear between 2 and 6 weeks");
  await devFigure("Spleen", "fig-dev-spleen", "The juvenile spleen still makes blood: haemoglobin and erythroid genes are high at 2 weeks, the spleen marker Fcrl5 is not yet");

  // ---- GTEx per tissue ----------------------------------------------------------------------------------
  const gState = { model: "k20" };
  const acc = GT.accuracy_by_tissue;
  const ov = Object.fromEntries(GT.accuracy_overall.map((r) => [r.model, r]));
  document.getElementById("p-gtex").replaceChildren(
    `Across 17 GTEx tissues and ${ex.gtex_donors} donors (${ex.orthologs_in_gtex} one-to-one orthologs of ${ex.motrpac_genes} rat genes), the 20-gene panel is right on ${fmt(ov.k20.accuracy_sample_weighted)} of samples, 50 genes on ${fmt(ov.k50.accuracy_sample_weighted)}, all genes on ${fmt(ov.full.accuracy_sample_weighted)}. `,
    `Two failures are structural rather than noisy: the human heart is called skeletal muscle by the panel (${fmt(ex.gtex_heart_k20)} correct; ${fmt(ex.gtex_heart_k20_to_skm_frac)} called SKM-GN or SKM-VL) and the human ovary is called heart (${fmt(ex.gtex_ovary_k20)} correct; ${fmt(ex.gtex_ovary_k20_top_frac)} called ${ex.gtex_ovary_k20_top}), because single rat markers (Gnb3 for heart, Akr1c3 for ovary) point elsewhere in human. `,
    `A native GTEx 20-gene panel reaches ${fmt(ex.gtex_native_k20)} under donor-grouped CV, so the gap is the transfer, not the task.`,
  );
  const gCtl = el("div", { class: "controls" }, [control("Model", segmented(MODELS.map((m) => [m, MODEL_LABEL[m]]), gState.model, (v) => { gState.model = v; gAcc.rerender(); gConf.rerender(); }, "model"))]);
  const accBuild = () => {
    const rows = acc.filter((r) => r.model === gState.model).sort((a, b) => a.accuracy - b.accuracy);
    return { traces: [{ ...bar(rows.map((r) => r.accuracy), rows.map((r) => r.gtex_tissue), { horizontal: true, slot: 1, name: "accuracy", text: rows.map((r) => fmt(r.accuracy, 2)), textposition: "outside", hover: "%{customdata}<extra></extra>" }),
                        customdata: rows.map((r) => `${r.gtex_tissue} → ${r.rat_classes}<br>accuracy ${fmt(r.accuracy)} (n = ${r.n} samples, ${r.n_donors} donors)<br>most frequent call: ${r.top_prediction} (${fmt(r.top_prediction_frac, 2)})`) }],
             layout: { xaxis: { range: [0, 1.15], title: { text: "accuracy (super-class scoring)" } }, yaxis: { automargin: true, tickfont: { size: 11 } }, margin: { t: 20, l: 10 }, showlegend: false, bargap: 0.3 },
             table: { columns: ["gtex_tissue", "rat_classes", "n", "n_donors", "accuracy", "top_prediction", "top_prediction_frac"], rows } };
  };
  const gAcc = await figure(document.getElementById("fig-gtex-acc"), {
    title: "Most human tissues transfer; heart, ovary, testis and adipose are the panel's blind spots",
    subtitle: "Accuracy per GTEx tissue for the selected model (sample-weighted; ≤ 150 donors per tissue); hover for the most frequent wrong call.",
    build: accBuild, toolbar: gCtl, source: "results/13_gtex/accuracy_by_tissue.csv",
    notShow: "conformal coverage per tissue (results/13_gtex/coverage_by_tissue.csv, in the data table of the ladder).", height: "tall",
  });
  gAcc.rerender = async () => { const b = accBuild(); gAcc.traces = b.traces; gAcc.table = b.table; await window.Plotly.react(gAcc.chart, b.traces, { ...template(), ...b.layout }, CONFIG); };
  const confBuild = () => {
    const c = GT[`confusion_${gState.model}`];
    const rowsum = c.counts.map((r) => r.reduce((a, b) => a + b, 0));
    const z = c.counts.map((r, i) => r.map((v) => (rowsum[i] ? v / rowsum[i] : 0)));
    const keep = c.cols.map((_, j) => c.counts.some((r) => r[j] > 0));
    const cols = c.cols.filter((_, j) => keep[j]);
    const zk = z.map((r) => r.filter((_, j) => keep[j]));
    return { traces: [heatmap(zk, cols, c.rows, { zmin: 0, zmax: 1, ztitle: "row fraction", hover: "%{y} → %{x}: %{z:.2f}<extra></extra>", text: zk.map((r) => r.map((v) => (v >= 0.2 ? v.toFixed(2) : ""))), texttemplate: "%{text}" })],
             layout: { xaxis: { title: { text: "predicted rat tissue" }, tickangle: -45, tickfont: { size: 10 } }, yaxis: { title: { text: "GTEx tissue" }, autorange: "reversed", tickfont: { size: 10 }, automargin: true }, margin: { t: 20, b: 90 } },
             table: { columns: ["gtex_tissue", ...cols], rows: c.rows.map((r, i) => ({ gtex_tissue: r, ...Object.fromEntries(cols.map((cc, j) => [cc, zk[i][j]])) })) } };
  };
  const gConf = await figure(document.getElementById("fig-gtex-conf"), {
    title: "Where the wrong calls go: heart → skeletal muscle, ovary → heart, adipose → brown fat (50 genes)",
    subtitle: "Row-normalised confusion: fraction of each GTEx tissue's samples called each rat tissue (columns with no call dropped); labels shown at ≥ 0.20.",
    build: confBuild, source: "results/13_gtex/confusion_k20.csv, confusion_k50.csv, confusion_full.csv", notShow: "counts (in the data table as fractions; n per tissue in the bar chart above).", height: "tall",
  });
  gConf.rerender = async () => { const b = confBuild(); gConf.traces = b.traces; gConf.table = b.table; await window.Plotly.react(gConf.chart, b.traces, { ...template(), ...b.layout }, CONFIG); };

  // ---- representations ----------------------------------------------------------------------------------
  const ts = REP.target_summary.filter((r) => r.selector === "standard");
  const g20 = (rep) => ts.find((r) => r.target === "GTEx" && r.k === 20 && r.representation === rep)?.accuracy;
  const cpm = REP.cpm_target_summary.filter((r) => r.selector === "standard");
  document.getElementById("p-rep").replaceChildren(
    `Three ways to hand the panel to another dataset: per-gene z-scores within each dataset, within-sample ranks of the panel genes, or which of two panel genes is higher. Within species all three transfer; across species only z-scores do (GTEx, 20 genes: z-score ${fmt(g20("zscore"))}, rank ${fmt(g20("rank"))}, pairs ${fmt(g20("pairs"))}). `,
    "Re-scoring GTEx on log2 CPM instead of TPM, the unit MoTrPAC and BodyMap use, changes the rank and pair accuracies by almost nothing, so the failure is biology (the within-sample ordering of the panel genes differs between rat and human tissues), not units. The ovary is the exception: ranks and pairs get it right where z-scores point at the heart.",
  );
  await figure(document.getElementById("fig-rep"), {
    title: "Only the z-score representation crosses species; the ovary is the one tissue where ranks and pairs win",
    subtitle: "Left: accuracy on mapped tissues by representation (standard selector), BodyMap adults and GTEx, k = 20 and 50. Right: GTEx accuracy per tissue at k = 20 by representation.",
    build: () => {
      const p = palette();
      const reps = ["zscore", "rank", "pairs"];
      const cats = ["BodyMap k20", "BodyMap k50", "GTEx k20", "GTEx k50"];
      const val = (target, k, rep) => ts.find((r) => r.target === (target === "BodyMap" ? "BodyMap adults" : "GTEx") && r.k === k && r.representation === rep)?.accuracy ?? null;
      const left = reps.map((rep, i) => ({ ...bar(cats, [val("BodyMap", 20, rep), val("BodyMap", 50, rep), val("GTEx", 20, rep), val("GTEx", 50, rep)], { name: rep, slot: i + 1, hover: "%{x}: %{y:.3f}<extra>" + rep + "</extra>" }), xaxis: "x", yaxis: "y" }));
      const ta = REP.target_accuracy.filter((r) => r.target === "GTEx" && r.selector === "standard" && r.k === 20);
      const organs = [...new Set(ta.map((r) => r.organ))].sort();
      const right = reps.map((rep, i) => ({ ...bar(organs.map((o) => ta.find((r) => r.organ === o && r.representation === rep)?.accuracy ?? null), organs, { horizontal: true, name: rep, slot: i + 1, showlegend: false, hover: "%{y} (" + rep + "): %{x:.2f}<extra></extra>" }), xaxis: "x2", yaxis: "y2" }));
      return { traces: [...left, ...right],
               layout: { grid: { rows: 1, columns: 2, pattern: "independent" }, xaxis: { domain: [0, 0.42] }, yaxis: { range: [0, 1.05], title: { text: "accuracy" } }, xaxis2: { domain: [0.58, 1], range: [0, 1], title: { text: "GTEx accuracy, k = 20" } }, yaxis2: { automargin: true, tickfont: { size: 9 }, autorange: "reversed" },
                         barmode: "group", legend: { y: 1.14 }, margin: { t: 40, b: 60 } },
               table: { columns: ["target", "representation", "k", "accuracy", "accuracy_macro_over_tissues", "n"], rows: ts } };
    },
    source: "results/14_transfer/target_summary.csv, results/14_transfer/target_accuracy.csv",
    notShow: "the transfer-aware selector (results/14_transfer/target_summary.csv, selector = transfer_aware): it refuses training-regulated and QC-correlated genes at a small in-distribution cost and does not help the species transfer of z-scores.", height: "tall",
  });
  await figure(document.getElementById("fig-cpm"), {
    title: "Matching the unit (CPM instead of TPM) does not rescue ranks or pairs in human",
    subtitle: "GTEx accuracy by representation and panel size, scored on log2 TPM (the run of record) and on log2 CPM from read counts (the MoTrPAC and BodyMap unit).",
    build: () => {
      const reps = ["zscore", "rank", "pairs"];
      const rowsK = [20, 50];
      const cats = rowsK.flatMap((k) => reps.map((r) => `${r} k${k}`));
      const tpm = rowsK.flatMap((k) => reps.map((r) => ts.find((x) => x.target === "GTEx" && x.k === k && x.representation === r)?.accuracy ?? null));
      const cpmv = rowsK.flatMap((k) => reps.map((r) => cpm.find((x) => x.target === "GTEx" && x.k === k && x.representation === r)?.accuracy ?? null));
      return { traces: [bar(cats, tpm, { name: "log2 TPM", slot: 1, hover: "%{x}: %{y:.3f}<extra>TPM</extra>" }), bar(cats, cpmv, { name: "log2 CPM", slot: 4, hover: "%{x}: %{y:.3f}<extra>CPM</extra>" })],
               layout: { barmode: "group", yaxis: { range: [0, 1], title: { text: "GTEx accuracy" } }, margin: { t: 40 }, legend: { y: 1.14 } },
               table: { columns: ["representation", "k", "accuracy_tpm", "accuracy_cpm"], rows: rowsK.flatMap((k) => reps.map((r, i) => ({ representation: r, k, accuracy_tpm: tpm[rowsK.indexOf(k) * 3 + i], accuracy_cpm: cpmv[rowsK.indexOf(k) * 3 + i] }))) } };
    },
    source: "results/14_transfer/target_summary.csv, results/14_transfer_cpm/target_summary.csv", notShow: "coverage under the two units (results/14_transfer_cpm/recalibration.csv): it moves by ≤ 0.01 too.", height: "short",
  });

  // ---- survival tables --------------------------------------------------------------------------------------
  const bsurv = BM.panel_survival, gsurv = GT.panel_survival;
  document.getElementById("p-survival").replaceChildren(
    `Every panel gene is present in the BodyMap (Ensembl ids shared: ${ex.shared_genes_bodymap} of ${ex.motrpac_genes} MoTrPAC genes). `,
    `In human, ${gsurv.find((r) => r.panel.startsWith("k=20")).with_1to1_ortholog_in_gtex} of 20 panel genes and ${gsurv.find((r) => r.panel.startsWith("T5")).with_1to1_ortholog_in_gtex} of the 10 stable-core genes have a one-to-one ortholog in GTEx; the panel is re-selected in ortholog space there. `,
    "Below, the per-gene check on the two targets: the marker's effect (mean z in its tissue minus the highest other tissue) at home and away.",
  );
  document.getElementById("tbl-survival").replaceChildren(el("h3", {}, "Panel survival"), tableFrom({ columns: ["target", "panel", "n_genes", "present", "lost"], rows: [
    ...bsurv.map((r) => ({ target: "BodyMap", panel: r.panel, n_genes: r.n_genes, present: r.present_in_bodymap, lost: "" })),
    ...gsurv.map((r) => ({ target: "GTEx", panel: r.panel, n_genes: r.n_genes, present: r.with_1to1_ortholog_in_gtex, lost: r.lost || "" })),
  ] }), el("p", { class: "small" }, "Source: results/12_bodymap/panel_survival.csv, results/13_gtex/panel_survival.csv"));
  const gcB = BM.gene_check, gcG = GT.gene_check;
  const ffmt = (v) => (v === null || v === undefined ? "" : typeof v === "number" ? v.toFixed(2) : String(v));
  document.getElementById("tbl-genecheck").replaceChildren(el("h3", {}, "Per-gene check of the 20-gene panel on each target"),
    tableFrom({ columns: ["gene", "marker_tissue", "target", "target_organ", "effect_z_source", "effect_z_target", "target_top", "fails", "weakened", "regulated_flag", "qc_flag"], rows: [
      ...gcB.map((r) => ({ gene: r.gene_symbol, marker_tissue: r.marker_tissue, target: "BodyMap", target_organ: r.bodymap_organ, effect_z_source: r.motrpac_effect_z, effect_z_target: r.bodymap_effect_z, target_top: r.bodymap_top_organ, fails: r.fails_in_bodymap, weakened: r.weakened, regulated_flag: r.risk_T7_regulated, qc_flag: r.risk_qc_correlated })),
      ...gcG.map((r) => ({ gene: r.gene_symbol, marker_tissue: r.marker_tissue, target: "GTEx", target_organ: r.target_organ, effect_z_source: r.source_effect_z, effect_z_target: r.target_effect_z, target_top: r.target_top_organ, fails: r.fails_in_target, weakened: r.weakened, regulated_flag: r.risk_T7_regulated, qc_flag: r.risk_qc_correlated })),
    ].map((r) => Object.fromEntries(Object.entries(r).map(([k, v]) => [k, typeof v === "boolean" ? (v ? "yes" : "no") : v]))), format: { effect_z_source: ffmt, effect_z_target: ffmt } }),
    el("p", { class: "small" }, "Source: results/12_bodymap/panel_gene_check.csv, results/13_gtex/panel_gene_check.csv (the GTEx panel is the ortholog-space re-selection, so its 20 genes differ in part)."));
}

main().catch((e) => { console.error(e); document.getElementById("lede").textContent = "Failed to load: " + e.message; });
