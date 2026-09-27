import { mountChrome, loadJSON, el, fmt, segmented, control, select, slider, badge, callout, tableFrom } from "../site.js";
import { figure, bar, strip, line, tokens, palette, organSystem, tissueLabel, TISSUE_NAMES, jitter, median, mean, sd, template, CONFIG } from "../charts.js";
import { setsFor, conformalQuantile, predictSet, animalsNeeded } from "../conformal.js";
import { parseTable, scoreSamples, recalibrate, templateCsv, resultsCsv, exampleRows, minLabelled, MIN_SAMPLES_FOR_WITHIN_Z } from "../score.js";

const AGE_ORDER = [2, 6, 21, 104];
let DATA = {};

async function main() {
  await mountChrome("explore.html");
  const [M, B, G, GENES, PC, CERT] = await Promise.all([
    loadJSON("data/samples_motrpac.json"), loadJSON("data/samples_bodymap.json"), loadJSON("data/samples_gtex.json"),
    loadJSON("data/genes.json"), loadJSON("data/panel_curve.json"), loadJSON("data/certificate.json"),
  ]);
  DATA = { M, B, G, GENES, PC, CERT };
  // expression files load in the background; the strip and gene explorer wait for them
  const exprP = Promise.all([loadJSON("data/expr_motrpac.json"), loadJSON("data/expr_bodymap.json"), loadJSON("data/expr_gtex.json")]).then(([em, eb, eg]) => { DATA.EM = em; DATA.EB = eb; DATA.EG = eg; return DATA; });
  tissueCard(exprP);
  geneExplorer(exprP);
  scoreTool(exprP);
  panelBuilder();
  calculator();
}

// ---- score your own samples ---------------------------------------------------------------------------
function download(name, text, type = "text/csv") {
  const a = el("a", { href: URL.createObjectURL(new Blob([text], { type })), download: name });
  document.body.appendChild(a); a.click(); a.remove();
}

async function scoreTool(exprP) {
  const model = await loadJSON("data/panel_model.json");
  const ctx = model.context || {};
  const intro = document.getElementById("score-intro");
  const ctl = document.getElementById("score-controls");
  const status = document.getElementById("score-status");
  const out = document.getElementById("score-results");
  const ta = document.getElementById("score-text");
  const nGenes = model.genes.length, nShared = model.n_genes_shared_with_all_animal_panel, cal = model.calibration;
  const nVials3 = typeof ctx.recal_n3_mean_cal_samples === "number" ? Math.round(ctx.recal_n3_mean_cal_samples) : null;
  intro.replaceChildren(
    `Score any rat RNA-seq samples with the transfer model that carries the guarantee: the ${nGenes}-gene panel fit on ${cal.n_animals ? 50 - cal.n_animals : "the fit"} MoTrPAC animals, its sets calibrated on the other ${cal.n_animals} (${nShared} of its ${nGenes} genes are in the all-animal panel of the `,
    el("a", { href: "fingerprint.html" }, "Panel page"), "). ",
    `Genes are matched by symbol or Ensembl ID; a missing gene is reported and scored at the MoTrPAC mean. With at least ${MIN_SAMPLES_FOR_WITHIN_Z} samples the genes are z-scored within your set, as the pipeline does for a new dataset; with fewer, the MoTrPAC statistics are used and the calls are less reliable. `,
    nVials3 !== null ? `Samples from another laboratory need their own calibration: in the study "three animals" meant about ${nVials3} labelled samples (${ctx.recal_n3_cal_samples_range[0]}–${ctx.recal_n3_cal_samples_range[1]} over ${ctx.recal_n3_draws} draws), one per mapped organ of each BodyMap rat. ` : "",
    "Everything runs in your browser; nothing is uploaded.",
  );
  const state = { alpha: 0.10, labels: {}, res: null, q: null, rows: [] };
  const need = () => minLabelled(state.alpha);
  const labelledCount = () => state.res ? state.res.results.filter((r) => state.labels[r.id] && model.classes.includes(state.labels[r.id])).length : 0;
  const a = slider(0.05, 0.30, 0.01, state.alpha, (v) => { state.alpha = v; if (state.res) run(); }, (v) => `α = ${v.toFixed(2)}`);
  const file = el("input", { type: "file", accept: ".csv,.tsv,.txt,text/csv,text/tab-separated-values" });
  file.addEventListener("change", () => { const f = file.files[0]; if (f) f.text().then((t) => { ta.value = t; run(); }); });
  const btnScore = el("button", { class: "btn primary", type: "button" }, "Score");
  const btnExample = el("button", { class: "btn", type: "button" }, "Load the example");
  const btnTemplate = el("button", { class: "btn", type: "button" }, "Download the template");
  const btnCard = el("button", { class: "btn", type: "button" }, "Download the panel card (CSV)");
  const linkCardJson = el("span", { class: "note" }, ["or as ", el("a", { href: "data/panel_card.json", download: "panel_card.json" }, "JSON")]);
  btnScore.addEventListener("click", run);
  // the example: every 21-week BodyMap sample, true_tissue pre-filled on 12 of the single-tissue organs
  const example = async () => {
    const [D, SB] = await Promise.all([exprP, loadJSON("data/samples_bodymap.json")]);
    return templateCsv(model, exampleRows(model, D.EB, SB.organ_map));
  };
  btnExample.addEventListener("click", async () => { ta.value = await example(); status.replaceChildren(callout("note", "Example loaded", `${ta.value.trim().split("\n").length - 1} rat BodyMap samples (21 weeks), with true_tissue filled in for 12 of them. Press Score.`)); out.replaceChildren(); });
  btnTemplate.addEventListener("click", async () => download("panel_template.csv", await example()));
  btnCard.addEventListener("click", () => { const l = el("a", { href: "data/panel_card.csv", download: "panel_card.csv" }); document.body.appendChild(l); l.click(); l.remove(); });
  ctl.append(control("Error level α", a.input, a.out), control("Or upload a file", file), btnScore, btnExample, btnTemplate, btnCard, linkCardJson);

  function run() {
    state.q = null;
    const parsed = parseTable(ta.value);
    if (parsed.error) { status.replaceChildren(callout("caveat", "Nothing to score", parsed.error)); out.replaceChildren(); return; }
    state.rows = parsed.rows;
    for (const r of parsed.rows) if (r.label) state.labels[r.id] = r.label;   // a true_tissue column pre-fills the selects
    const res = scoreSamples(model, parsed.rows, { alpha: state.alpha });
    state.res = res;
    render(res, null);
  }
  function render(res, recal) {
    const N = need();
    const bits = [`${res.results.length} samples scored; ${res.matched} of ${nGenes} panel genes matched; ${res.note}.`];
    if (res.missing.length) bits.push(` Missing (scored at the MoTrPAC mean): ${res.missing.map((g) => g.symbol).join(", ")}.`);
    const abstained = res.results.filter((r) => r.setSize === 0).length;
    if (!recal && abstained > 0) {
      bits.push(` ${abstained} of ${res.results.length} samples abstained with the MoTrPAC calibration — expected for another laboratory's samples`,
                typeof ctx.bodymap_coverage_k20 === "number" ? ` (the same effect gives ${fmt(ctx.bodymap_coverage_k20)} coverage on the rat BodyMap)` : "",
                `. Label ≥ ${N} samples and recalibrate`, nVials3 !== null ? `; in the study three animals (≈ ${nVials3} samples) were enough.` : ".");
    }
    if (recal) bits.push(` Recalibrated on ${recal.n} labelled samples: threshold ${recal.q === Infinity ? `+∞ (every set holds all tissues; at α = ${res.alpha.toFixed(2)} a finite threshold needs at least ${N} labelled samples)` : fmt(recal.q)}, coverage on the labelled samples ${fmt(recal.coverage)}.`);
    status.replaceChildren(callout(res.missing.length || res.mode === "source" ? "caveat" : "note", recal ? "Recalibrated" : "Scoring", bits.join("")));
    const rows = res.results.map((r) => ({ sample: r.id, call: r.call, p_call: r.callProb, set: r.set.join(", ") || "∅", set_size: r.setSize, kind: r.kind }));
    const table = tableFrom({ columns: ["sample", "call", "p_call", "set", "set_size", "kind", "true tissue (optional)"], rows: rows.map((r) => ({ ...r, "true tissue (optional)": "" })), format: { p_call: (v) => fmt(v, 3) } });
    // the recalibration control: enabled once N samples are labelled, with a live counter
    const btnRecal = el("button", { class: "btn", type: "button" }, "");
    const counter = el("span", { class: "note" }, "");
    const refresh = () => {
      const k = labelledCount();
      btnRecal.textContent = `Recalibrate on the labelled samples (needs ≥ ${N} at α = ${state.alpha.toFixed(2)})`;
      btnRecal.disabled = k < N;
      counter.textContent = `${k} labelled (minimum ${N})`;
    };
    // a tissue picker per row (pre-filled from a true_tissue column or earlier choices)
    table.querySelectorAll("tbody tr").forEach((tr, i) => {
      const cell = tr.lastElementChild;
      const id = res.results[i].id;
      const sel = select([["", "—"], ...model.classes.map((c) => [c, c])], state.labels[id] || "", (v) => { if (v) state.labels[id] = v; else delete state.labels[id]; refresh(); });
      cell.replaceChildren(sel);
    });
    const cards = el("div", { class: "score-cards" }, res.results.slice(0, 24).map((r) => el("div", { class: "sample-card" }, [
      el("div", { class: "who" }, [el("b", {}, r.id)]),
      badge(r.kind, r.kind === "confident" ? "●" : r.kind === "ambiguous" ? "◐" : "∅", r.kind[0].toUpperCase() + r.kind.slice(1)),
      el("div", { class: "answer" }, r.set.length ? r.set.map((t) => tissueLabel(t)).join("; ") : "no tissue clears the threshold"),
      el("div", { class: "explain" }, `top call ${tissueLabel(r.call)} (p = ${fmt(r.callProb)})`),
    ])));
    const btnCsv = el("button", { class: "btn", type: "button" }, "Download the results (CSV)");
    btnCsv.addEventListener("click", () => download("scored_samples.csv", resultsCsv(res, model)));
    btnRecal.addEventListener("click", () => {
      const rc = recalibrate(model, state.res.results, state.labels, state.alpha);
      if (rc.error) { status.replaceChildren(callout("caveat", "Recalibration", rc.error)); return; }
      const res2 = scoreSamples(model, state.rows, { alpha: state.alpha, q: rc.q });
      state.res = res2; state.q = rc.q;
      render(res2, rc);
    });
    refresh();
    out.replaceChildren(
      res.results.length > 24 ? el("p", { class: "small" }, `Cards for the first 24 samples; the table and the CSV hold all ${res.results.length}.`) : el("span"),
      cards, el("div", { class: "controls" }, [btnCsv, btnRecal, counter]), table,
      el("details", { class: "fig-notes" }, [el("summary", {}, "Source and caveats"),
        el("p", {}, [el("b", {}, "Source: "), "site/data/panel_model.json (results/34_panel_model: the phase-12 transfer model re-fit with the same seed; validated against results/31_site_regen/12_bodymap/scores_target_probs.csv), calibration scores from results/31_site_regen/12_bodymap/scores_calibration.csv; the BodyMap coverage and the size of a three-animal recalibration from results/12_bodymap/conformal_transfer.csv and results/31_site_regen/12_bodymap/recal_thresholds.csv; the example is every 21-week rat BodyMap sample (site/data/expr_bodymap.json)."]),
        el("p", {}, [el("b", {}, "What it does not show: "), "the guarantee is nominal for samples from a new laboratory or species until you recalibrate on labelled samples of your own; z-scoring within your set assumes it spans several tissues (a small upload compresses the marker z-scores and raises the recalibrated threshold); a gene missing from your table is scored at the MoTrPAC mean; abstention on tissues outside the 19 is a tendency, not a guarantee."]),
      ]),
    );
  }
}

// ---- helpers ----------------------------------------------------------------------------------------
function zStats(expr) {
  // per-gene mean and sd over all samples of the dataset (within-dataset z, as the pipeline does)
  if (!expr._z) {
    expr._z = expr.values.map((row) => { const m = mean(row); const s = sd(row) || 1; return { m, s }; });
    expr._idx = Object.fromEntries(expr.genes.map((g, i) => [g, i]));
    expr._sidx = Object.fromEntries(expr.samples.map((s, i) => [s.id, i]));
  }
  return expr;
}

function refProfiles(expr, groupKey) {
  // mean z per group (tissue / organ) for every gene
  const key = "_ref_" + groupKey;
  if (!expr[key]) {
    zStats(expr);
    const groups = {};
    expr.samples.forEach((s, i) => { (groups[s[groupKey]] ||= []).push(i); });
    const out = {};
    for (const [g, idx] of Object.entries(groups)) {
      out[g] = expr.values.map((row, gi) => mean(idx.map((i) => (row[i] - expr._z[gi].m) / expr._z[gi].s)));
    }
    expr[key] = out;
  }
  return expr[key];
}

// ---- (a) the tissue card ------------------------------------------------------------------------------
function tissueCard(exprP) {
  const { M, B, G, GENES } = DATA;
  const classes = M.classes;
  const state = { source: "motrpac", alpha: 0.10, model: "k20", variant: "marginal", calib: "source", mcal: "pooled", recalN: "3", hide: false, sample: null, sub1: null, sub2: null };
  const picker = document.getElementById("picker");
  const ctl = document.getElementById("card-controls");
  const card = document.getElementById("card");
  let figTop, figStrip;

  const sourceSel = segmented([["motrpac", "MoTrPAC held-out vial"], ["bodymap", "rat BodyMap"], ["gtex", "human GTEx"]], state.source, (v) => { state.source = v; buildPicker(); }, "sample source");
  const sub1Wrap = el("span"), sub2Wrap = el("span"), sampleWrap = el("span");
  const randomBtn = el("button", { class: "btn", type: "button" }, "Random sample");
  const hide = el("input", { type: "checkbox" });
  hide.addEventListener("change", () => { state.hide = hide.checked; applyHide(); render(); });
  const hideHint = el("span", { class: "small muted" });
  function applyHide() {
    for (const w of [sub1Wrap, sub2Wrap, sampleWrap]) {
      w.classList.toggle("answer-hidden", state.hide);
      w.querySelectorAll("select").forEach((sel) => { sel.disabled = state.hide; });
    }
    hideHint.textContent = state.hide ? "picker hidden: use Random sample" : "";
  }
  picker.append(control("Source", sourceSel), sub1Wrap, sub2Wrap, sampleWrap, el("span", { class: "control" }, [el("span", {}, " "), randomBtn]),
                el("label", { class: "check" }, [hide, "hide the answer (for demos)"]), hideHint);

  function samplesOf(source) { return source === "motrpac" ? M.samples : source === "bodymap" ? B.samples : G.samples; }
  function buildPicker() {
    const src = state.source;
    sub1Wrap.replaceChildren(); sub2Wrap.replaceChildren(); sampleWrap.replaceChildren();
    if (src === "motrpac") {
      const tissues = [...new Set(M.samples.map((s) => s.tissue))].sort();
      state.sub1 = state.sub1 && tissues.includes(state.sub1) ? state.sub1 : "TESTES";
      sub1Wrap.replaceChildren(control("Tissue", select(tissues.map((t) => [t, tissueLabel(t)]), state.sub1, (v) => { state.sub1 = v; buildSamples(); })));
      buildSamples();
    } else if (src === "bodymap") {
      state.sub1 = AGE_ORDER.includes(Number(state.sub1)) ? Number(state.sub1) : 21;
      const organs = [...new Set(B.samples.map((s) => s.organ))].sort();
      state.sub2 = state.sub2 && organs.includes(state.sub2) ? state.sub2 : "Thymus";
      sub1Wrap.replaceChildren(control("Age", select(AGE_ORDER.map((a) => [String(a), `${a} weeks`]), String(state.sub1), (v) => { state.sub1 = Number(v); buildSamples(); })));
      sub2Wrap.replaceChildren(control("Organ", select(organs.map((o) => [o, B.organ_map[o] ? o : `${o} (not in MoTrPAC)`]), state.sub2, (v) => { state.sub2 = v; buildSamples(); })));
      buildSamples();
    } else {
      const tissues = [...new Set(G.samples.map((s) => s.tissue))].sort();
      state.sub1 = state.sub1 && tissues.includes(state.sub1) ? state.sub1 : "Heart - Left Ventricle";
      if (!state.sample) {  // first visit: a heart sample the all-animal 20-gene panel calls skeletal muscle (the Transfer page's story)
        state.sample = G.samples.find((x) => x.tissue === state.sub1 && /^SKM/.test(x.pred_all_animals.k20)) || null;
      }
      sub1Wrap.replaceChildren(control("GTEx tissue", select(tissues.map((t) => [t, `${t} → ${(G.organ_map[t] || []).join("/")}`]), state.sub1, (v) => { state.sub1 = v; buildSamples(); })));
      buildSamples();
    }
    buildControls();
    applyHide();
  }
  function candidates() {
    const src = state.source;
    if (src === "motrpac") return M.samples.filter((s) => s.tissue === state.sub1);
    if (src === "bodymap") return B.samples.filter((s) => s.age_weeks === state.sub1 && s.organ === state.sub2);
    return G.samples.filter((s) => s.tissue === state.sub1);
  }
  function label(s) {
    if (state.source === "motrpac") return `${s.id} · ${s.sex} · ${s.group} · fold ${s.fold}`;
    if (state.source === "bodymap") return `${s.id} · ${s.sex} · animal ${s.animal}`;
    return `${s.id} · donor ${s.donor}`;
  }
  function buildSamples() {
    const c = candidates();
    if (!c.length) { sampleWrap.replaceChildren(); state.sample = null; render(); return; }
    if (!state.sample || !c.includes(state.sample)) state.sample = c[0];
    sampleWrap.replaceChildren(control("Sample", select(c.map((s) => [s.id, label(s)]), state.sample.id, (v) => { state.sample = c.find((s) => s.id === v); render(); })));
    render();
  }
  randomBtn.addEventListener("click", () => {
    const all = samplesOf(state.source);
    const s = all[Math.floor(Math.random() * all.length)];
    if (state.source === "motrpac") state.sub1 = s.tissue; else if (state.source === "bodymap") { state.sub1 = s.age_weeks; state.sub2 = s.organ; } else state.sub1 = s.tissue;
    state.sample = s;
    buildPicker();
  });

  function buildControls() {
    ctl.replaceChildren();
    // in recalibrated mode the set uses the draw's own threshold (α = 0.10, marginal), so the α and variant controls are switched off
    const off = state.source !== "motrpac" && state.calib === "recal";
    const a = slider(0.05, 0.30, 0.01, off ? 0.10 : state.alpha, (v) => { state.alpha = v; render(); }, (v) => `α = ${v.toFixed(2)} (target coverage ${((1 - v) * 100).toFixed(0)} %)`);
    a.input.disabled = off;
    const alphaCtl = control(off ? "Error level α (fixed at 0.10 by the recalibration draw)" : "Error level α", a.input, a.out);
    alphaCtl.classList.toggle("is-off", off);
    ctl.append(alphaCtl);
    ctl.append(control("Model", segmented([["k20", "20 genes"], ["k50", "50 genes"], ["full", "all genes"]], state.model, (v) => { state.model = v; render(); }, "model")));
    const variantCtl = control(off ? "Prediction set (marginal in recalibrated mode)" : "Prediction set", segmented([["marginal", "marginal"], ["mondrian", "Mondrian"], ["floored", "floored Mondrian"]], state.variant, (v) => { state.variant = v; render(); }, "conformal variant"));
    variantCtl.classList.toggle("is-off", off);
    if (off) variantCtl.querySelectorAll("button").forEach((b) => { b.disabled = true; });
    ctl.append(variantCtl);
    if (state.source === "motrpac") {
      ctl.append(control("Calibration (22 held-out animals)", segmented([["pooled", "pooled vials"], ["one_per_animal", "one vial per animal"]], state.mcal, (v) => { state.mcal = v; render(); }, "calibration")));
    } else {
      ctl.append(control("Calibration", segmented([["source", "15 MoTrPAC animals (source)"], ["recal", "recalibrated on target individuals"]], state.calib, (v) => { state.calib = v; buildControls(); render(); }, "calibration")));
      ctl.append(control("Recalibration n (draw 0 of 20)", segmented([["3", "3"], ["5", "5"]], state.recalN, (v) => { state.recalN = v; render(); }, "recalibration size")));
    }
  }

  function computeSet() {
    const s = state.sample;
    const p = s.p[state.model];
    let set, q, qs = null, note = "";
    if (state.source === "motrpac") {
      const cal = M.calibration[String(s.fold)][state.mcal];
      const r = setsFor({ scores: cal.scores[state.model], y_idx: cal.y_idx, n_classes: classes.length }, p, state.variant, state.alpha);
      set = r.set; q = r.q; qs = r.qs;
      note = `calibration: fold ${s.fold}, ${state.mcal === "pooled" ? cal.scores[state.model].length + " vials" : cal.scores[state.model].length + " vials, one per animal"} of the 22 calibration animals`;
    } else {
      const D = state.source === "bodymap" ? B : G;
      if (state.calib === "recal") {
        const rt = D.recal_thresholds[state.model][state.recalN];
        q = rt.draw0.q === "inf" ? Infinity : rt.draw0.q;
        set = predictSet(p, q);
        note = `recalibrated threshold from draw 0 (α = 0.10 fixed; the variant and α controls do not apply): ${rt.draw0.n_cal_scores} calibration scores from ${state.recalN} target individuals (${rt.draw0.chosen.join(", ")}); ${(100 * rt.frac_infinite).toFixed(0)} % of the 20 draws have an infinite threshold`;
      } else {
        const cal = D.calibration[state.model];
        const r = setsFor({ scores: cal.scores, y_idx: cal.y_idx, n_classes: classes.length }, p, state.variant, state.alpha);
        set = r.set; q = r.q; qs = r.qs;
        note = `calibration: ${cal.n_vials} vials of ${cal.n_animals} held-out MoTrPAC animals (pooled)`;
      }
    }
    return { p, set, q, qs, note };
  }

  function truth() {
    const s = state.sample;
    if (state.source === "motrpac") return { label: tissueLabel(s.tissue), mapped: [s.tissue] };
    if (state.source === "bodymap") { const m = B.organ_map[s.organ]; return { label: `${s.organ}, ${s.age_weeks} weeks` + (m ? ` → ${m.join(" / ")}` : " (no MoTrPAC tissue)"), mapped: m }; }
    const m = G.organ_map[s.tissue];
    return { label: `${s.tissue} → ${(m || []).join(" / ")}`, mapped: m };
  }

  async function render() {
    if (!state.sample) { card.replaceChildren(el("p", {}, "No sample in this selection.")); return; }
    const s = state.sample;
    const { p, set, q, qs, note } = computeSet();
    const chosen = classes.map((c, i) => [c, i]).filter(([, i]) => set[i]);
    const tr = truth();
    const covered = tr.mapped ? tr.mapped.some((t) => set[classes.indexOf(t)]) : null;
    const top = classes.map((c, i) => [c, p[i]]).sort((a, b) => b[1] - a[1]);
    const chips = el("div", { class: "chips" });
    if (!chosen.length) chips.appendChild(el("span", { class: "chip empty" }, "∅ empty set — the fingerprint abstains"));
    for (const [c] of chosen) {
      const os = organSystem(c);
      const ch = el("span", { class: "chip" + (tr.mapped && tr.mapped.includes(c) && !state.hide ? " true" : "") }, [el("span", { class: "dot", style: `background:${os.color}` }), tissueLabel(c)]);
      ch.title = `${os.name} · p = ${fmt(p[classes.indexOf(c)])}`;
      chips.appendChild(ch);
    }
    let kind, icon, sentence;
    if (chosen.length === 0) { kind = "abstains"; icon = "∅"; sentence = `Abstains. At α = ${state.alpha.toFixed(2)} no tissue clears the calibrated threshold, so the set is empty: under the guarantee the model names no tissue rather than guessing.`; }
    else if (chosen.length === 1) { kind = "confident"; icon = "●"; sentence = `Confident. Exactly one tissue clears the threshold: ${tissueLabel(chosen[0][0])}.`; }
    else { kind = "ambiguous"; icon = "◐"; sentence = `Ambiguous. ${chosen.length} tissues clear the threshold, so the set returns all of them.`; }
    const thr = Number.isFinite(q) ? `threshold q̂ = ${fmt(q, 3)}: a tissue enters the set when p ≥ ${fmt(1 - q, 3)}` : "threshold q̂ = +∞ (too few calibration scores for this α): every tissue enters the set";
    const who = el("div", { class: "who" }, [el("b", {}, s.id), " · ", state.source === "motrpac" ? `MoTrPAC vial, ${s.sex}, ${s.group}, out-of-fold (fold ${s.fold})` : state.source === "bodymap" ? `rat BodyMap, ${s.sex}, animal ${s.animal}` : `GTEx v8, donor ${s.donor}`]);
    const topCall = top[0][0];
    const calls = state.source === "motrpac" ? null : el("p", { class: "explain" }, [
      el("b", {}, "Two models, two calls. "),
      `The sets come from the model refit on 35 MoTrPAC animals, calibrated on the other 15 (top call ${tissueLabel(topCall)}); the model fit on all 50 animals, the one the accuracy tables use, calls this sample `, el("b", {}, tissueLabel(s.pred_all_animals[state.model])), ".",
      s.pred_all_animals[state.model] !== topCall ? " The two disagree here: a panel re-selected on 35 animals is not the same panel." : "",
    ]);
    const answer = el("div", { class: "answer" + (state.hide ? " hidden-answer" : "") }, [el("b", {}, "True tissue: "), tr.label]);
    const verdict = covered === null ? badge("unmapped", "◌", "no MoTrPAC tissue to be right about — abstention is the desired behaviour")
      : (state.hide ? badge("unmapped", "?", "answer hidden") : covered ? badge("correct", "✓", "the true tissue is in the set") : badge("wrong", "✗", chosen.length ? "the true tissue is not in the set" : "missed by abstaining"));
    card.replaceChildren(who, el("div", {}, [badge(kind, icon, kind[0].toUpperCase() + kind.slice(1))]), chips, el("p", { class: "explain" }, sentence), el("p", { class: "explain" }, thr + (state.variant !== "marginal" && qs && state.calib !== "recal" ? "; Mondrian per-class thresholds shown in the bars' hover" : "")),
                        answer, el("div", {}, [verdict]), ...(calls ? [calls] : []), el("p", { class: "explain small" }, note));

    // top-5 probabilities
    const t = tokens();
    const top5 = top.slice(0, 5);
    const specTop = {
      title: "Top-5 class probabilities and the calibrated threshold",
      subtitle: `${state.model === "full" ? "all-gene" : state.model.replace("k", "") + "-gene"} model${state.source === "motrpac" ? " (fit on the fold's 18 fit animals)" : " (refit on 35 MoTrPAC animals for calibration)"}; a tissue enters the set when its probability is ≥ 1 − q̂.`,
      build: () => {
        const y = top5.map(([c]) => tissueLabel(c)).reverse(), x = top5.map(([, v]) => v).reverse();
        const inset = top5.map(([c]) => set[classes.indexOf(c)]).reverse();
        const qsArr = qs ? top5.map(([c]) => qs[classes.indexOf(c)]).reverse() : null;
        const traces = [{ ...bar(x, y, { horizontal: true, name: "probability", slot: 1, text: x.map((v) => fmt(v, 3)), textposition: "outside",
                                        hover: "%{y}: p = %{x:.3f}%{customdata}<extra></extra>" }),
                          customdata: y.map((_, i) => (inset[i] ? " · in the set" : "") + (qsArr ? `<br>class threshold 1 − q̂ = ${Number.isFinite(qsArr[i]) ? fmt(1 - qsArr[i], 3) : "−∞"}` : "")),
                          marker: { color: inset.map((b) => (b ? palette()[0] : t.axis)), line: { color: t.surface, width: 2 }, cornerradius: 4 } }];
        const shapes = Number.isFinite(q) && !qs ? [{ type: "line", x0: 1 - q, x1: 1 - q, xref: "x", y0: -0.5, y1: 4.5, yref: "y", line: { color: t.ink2, width: 1, dash: "dash" } }] : [];
        const annotations = Number.isFinite(q) && !qs ? [{ x: 1 - q, y: 4.6, xref: "x", yref: "y", text: `1 − q̂ = ${fmt(1 - q, 2)}`, showarrow: false, font: { color: t.ink2, size: 11 }, xanchor: "left" }] : [];
        return { traces, layout: { xaxis: { range: [0, 1.15], title: { text: "class probability" } }, yaxis: { automargin: true }, shapes, annotations, margin: { t: 20, l: 10 }, showlegend: false, bargap: 0.35 },
                 table: { columns: ["tissue", "probability", "in_set"], rows: top.map(([c, v]) => ({ tissue: c, probability: v, in_set: set[classes.indexOf(c)] ? "yes" : "no" })) } };
      },
      source: state.source === "motrpac" ? "results/31_site_regen/06_conformal/TRNSCRPT/scores_test_probs.csv, scores_calibration.csv" : `results/31_site_regen/${state.source === "bodymap" ? "12_bodymap" : "13_gtex"}/scores_target_probs.csv, scores_calibration.csv, recal_thresholds.csv`,
      notShow: "the 14 lower-ranked tissues (in the data table); under Mondrian variants each class has its own threshold (hover).",
      height: "short",
    };
    if (!figTop) figTop = await figure(document.getElementById("fig-top5"), specTop);
    else { figTop.spec = specTop; const b = specTop.build(); figTop.traces = b.traces; figTop.table = b.table; await window.Plotly.react(figTop.chart, b.traces, { ...template(), ...b.layout }, CONFIG); figTop.root.querySelector(".fig-sub").textContent = specTop.subtitle; figTop.root.querySelector(".fig-notes p").lastChild.textContent = specTop.source; }

    // panel-gene strip (needs the expression exports)
    await exprP;
    const expr = state.source === "motrpac" ? DATA.EM : state.source === "bodymap" ? DATA.EB : DATA.EG;
    zStats(expr);
    const groupKey = state.source === "motrpac" ? "tissue" : state.source === "bodymap" ? "organ" : "tissue";
    const refs = refProfiles(expr, groupKey);
    const panelIds = (state.model === "k50" ? GENES.sets.k50 : GENES.sets.k20).filter((g) => expr._idx[g] !== undefined);
    const gi = Object.fromEntries(GENES.genes.map((g) => [g.id, g]));
    const si = expr._sidx[s.id];
    const sampleZ = panelIds.map((g) => { const i = expr._idx[g]; return (expr.values[i][si] - expr._z[i].m) / expr._z[i].s; });
    const topTissue = top[0][0];
    const refGroupForTissue = (tt) => state.source === "motrpac" ? tt : state.source === "bodymap" ? Object.entries(B.organ_map).find(([, v]) => v && v.includes(tt))?.[0] : Object.entries(G.organ_map).find(([, v]) => v && v.includes(tt))?.[0];
    const topGroup = refGroupForTissue(topTissue);
    const trueGroup = state.source === "motrpac" ? s.tissue : state.source === "bodymap" ? s.organ : s.tissue;
    const specStrip = {
      title: `Where this sample sits on each panel gene, against the tissue reference profiles`,
      subtitle: `Sample z (bars) for the ${panelIds.length} panel genes, with the mean z of the top-scoring tissue (◆) and of the sample's own group (●); marker tissue in the hover.`,
      build: () => {
        const x = panelIds.map((g) => gi[g]?.symbol || g);
        const traces = [{ ...bar(x, sampleZ, { name: "this sample (z)", slot: 1, hover: "%{x}: z = %{y:.2f}<br>%{customdata}<extra>this sample</extra>" }), customdata: panelIds.map((g) => `marker of ${gi[g]?.marker_tissue || "?"}`) }];
        if (topGroup && refs[topGroup]) traces.push(line(x, panelIds.map((g) => refs[topGroup][expr._idx[g]]), { name: `reference: ${topGroup} (top-scoring)`, slot: 2, mode: "markers", symbol: "diamond", size: 9, hover: "%{x}: %{y:.2f}<extra>" + topGroup + " reference</extra>" }));
        if (!state.hide && refs[trueGroup] && trueGroup !== topGroup) traces.push(line(x, panelIds.map((g) => refs[trueGroup][expr._idx[g]]), { name: `reference: ${trueGroup} (sample's group)`, slot: 3, mode: "markers", symbol: "circle", size: 9, hover: "%{x}: %{y:.2f}<extra>" + trueGroup + " reference</extra>" }));
        return { traces, layout: { yaxis: { title: { text: "z within dataset" }, zeroline: true }, xaxis: { tickangle: -45, tickfont: { size: 10 } }, margin: { t: 40, b: 70 }, legend: { y: 1.16 } },
                 table: { columns: ["gene", "marker_tissue", "sample_z", "ref_top", "ref_group"], rows: panelIds.map((g, i) => ({ gene: gi[g]?.symbol || g, marker_tissue: gi[g]?.marker_tissue, sample_z: sampleZ[i], ref_top: topGroup && refs[topGroup] ? refs[topGroup][expr._idx[g]] : null, ref_group: refs[trueGroup] ? refs[trueGroup][expr._idx[g]] : null })) } };
      },
      source: `${state.source === "motrpac" ? "expr_motrpac.json (log2 CPM, io.stack_tissues on data/raw/counts)" : state.source === "bodymap" ? "expr_bodymap.json (log2 CPM from data/external/bodymap_counts.csv)" : "expr_gtex.json (log2 TPM, data/external/gtex_tpm_subset.csv)"}; panel membership from results/12_bodymap/panel_survival.csv`,
      notShow: "the classifier's weights; the strip is descriptive (z of raw expression), the model uses the same z-scores through logistic regression.",
    };
    if (!figStrip) figStrip = await figure(document.getElementById("fig-strip"), specStrip);
    else { figStrip.spec = specStrip; const b = specStrip.build(); figStrip.traces = b.traces; figStrip.table = b.table; await window.Plotly.react(figStrip.chart, b.traces, { ...template(), ...b.layout }, CONFIG); figStrip.root.querySelector(".fig-sub").textContent = specStrip.subtitle; figStrip.root.querySelector("h3").textContent = specStrip.title; }
  }
  buildPicker();
}

// ---- (b) the gene explorer ----------------------------------------------------------------------------
function geneExplorer(exprP) {
  const { GENES, B } = DATA;
  const genes = [...GENES.genes].sort((a, b) => (a.symbol || a.id).localeCompare(b.symbol || b.id));
  const byId = Object.fromEntries(genes.map((g) => [g.id, g]));
  const bySym = Object.fromEntries(genes.map((g) => [g.symbol, g.id]));
  const state = { gene: bySym["Pgk2"] || genes[0].id };
  const ctl = document.getElementById("gene-controls");
  const sel = select(genes.map((g) => [g.id, `${g.symbol || g.id}${g.marker_tissue ? " · " + g.marker_tissue : ""}${g.in_core ? " · core" : g.in_k20 ? " · k20" : ""}`]), state.gene, (v) => { state.gene = v; render(); });
  const feat = (sym) => el("button", { class: "btn", type: "button", onclick: () => { state.gene = bySym[sym]; sel.value = state.gene; render(); } }, sym);
  ctl.append(control("Gene (panel and candidate genes)", sel), el("span", { class: "control" }, [el("span", {}, "Featured"), el("span", {}, [feat("Pgk2"), " ", feat("Gnb3")])]));
  const figs = {};
  async function render() {
    await exprP;
    const { EM, EB, EG } = DATA;
    const g = byId[state.gene];
    const card = document.getElementById("gene-card");
    const dl = el("dl", { class: "kv" });
    const row = (k, v) => { dl.append(el("dt", {}, k), el("dd", {}, v === null || v === undefined ? "—" : String(v))); };
    row("gene", `${g.symbol} (${g.id})`);
    row("marker tissue", g.marker_tissue ? tissueLabel(g.marker_tissue) : "—");
    row("selection frequency (k = 20 bootstraps)", g.freq === null || g.freq === undefined ? "not in the stability run" : fmt(g.freq, 2));
    row("effect (log2 CPM above the next tissue)", g.effect_size === null || g.effect_size === undefined ? "—" : fmt(g.effect_size, 2) + (g.next_highest_tissue ? ` (next: ${g.next_highest_tissue})` : ""));
    row("one-vs-rest score", g.ovr_score === null || g.ovr_score === undefined ? "—" : fmt(g.ovr_score, 2));
    row("r with mRNA fraction (marker tissue)", g.r_pct_mrna === null || g.r_pct_mrna === undefined ? "—" : fmt(g.r_pct_mrna, 2) + (g.qc_flag ? " · QC-correlated flag" : ""));
    row("training-regulated (marker tissue)", g.regulated === null || g.regulated === undefined ? "—" : g.regulated ? "yes (risk flag)" : "no");
    row("in panels", [g.in_core ? "stable core" : null, g.in_k20 ? "k = 20" : null, g.in_k50 ? "k = 50" : null, g.in_developmental ? "developmental marker" : null].filter(Boolean).join(", ") || "none (selected in " + g.n_folds_selected_k20 + " of 5 folds at k = 20)");
    row("BodyMap adults", g.fails_bodymap === null || g.fails_bodymap === undefined ? "not testable (no organ)" : g.fails_bodymap ? `lost (top organ ${g.bodymap_top_organ})` : g.weakened_bodymap ? "weakened" : "holds");
    row("GTEx", g.fails_gtex === null || g.fails_gtex === undefined ? "not testable / not in the ortholog-space panel" : g.fails_gtex ? `lost (top tissue ${g.gtex_top_tissue})` : g.weakened_gtex ? "weakened" : "holds");
    row("human ortholog", g.human_gene || "—");
    card.replaceChildren(el("div", { class: "who" }, [el("b", {}, g.symbol), " annotation card"]), dl,
                         el("details", { class: "fig-notes" }, [el("summary", {}, "Source"), el("p", {}, [el("b", {}, "Source: "), "results/05_panels/TRNSCRPT/stability_k20_annotated.csv, results/12_bodymap/panel_gene_check.csv, results/13_gtex/panel_gene_check.csv"])]));
    const systemOrder = ["CORTEX", "HIPPOC", "HYPOTH", "SKM-GN", "SKM-VL", "HEART", "WAT-SC", "BAT", "COLON", "SMLINT", "OVARY", "TESTES", "BLOOD", "SPLEEN", "VENACV", "LIVER", "KIDNEY", "LUNG", "ADRNL"];
    const t = tokens();
    // MoTrPAC panel
    const specM = {
      title: `${g.symbol} in MoTrPAC: log2 CPM by tissue`,
      subtitle: "Every study vial (899); sex as marker shape (● female, ◆ male), colour = organ system, bars = tissue median.",
      build: () => {
        const i = EM._idx?.[g.id] ?? EM.genes.indexOf(g.id);
        if (i < 0) return { traces: [], layout: { annotations: [{ text: "gene absent from the stacked MoTrPAC matrix", showarrow: false }] } };
        const vals = EM.values[i];
        const tissues = systemOrder.filter((tt) => EM.samples.some((s) => s.tissue === tt));
        const xi = Object.fromEntries(tissues.map((tt, k) => [tt, k]));
        const xs = [], ys = [], cols = [], syms = [], hov = [];
        EM.samples.forEach((s, k) => { xs.push(xi[s.tissue] + jitter(k, 0, 0.28)); ys.push(vals[k]); cols.push(organSystem(s.tissue).color); syms.push(s.sex === "male" ? "diamond" : "circle"); hov.push(`${s.tissue} · ${s.sex} · ${s.group}<br>log2 CPM ${vals[k].toFixed(2)}`); });
        const meds = tissues.map((tt) => median(EM.samples.map((s, k) => (s.tissue === tt ? vals[k] : NaN))));
        return { traces: [
          { ...strip(xs, ys, { colors: cols, symbols: syms, size: 7, name: "vials", hover: "%{customdata}<extra></extra>" }), customdata: hov, showlegend: false },
          { type: "scatter", mode: "markers", x: tissues.map((_, k) => k), y: meds, name: "tissue median", marker: { symbol: "line-ew", size: 22, color: t.ink, line: { width: 2, color: t.ink } }, hovertemplate: "%{text}: median %{y:.2f}<extra></extra>", text: tissues },
        ], layout: { xaxis: { tickvals: tissues.map((_, k) => k), ticktext: tissues, tickangle: -60, tickfont: { size: 10 } }, yaxis: { title: { text: "log2 CPM" } }, margin: { t: 30, b: 80 }, showlegend: false },
          table: { columns: ["tissue", "median_log2_cpm", "n"], rows: tissues.map((tt, k) => ({ tissue: tt, median_log2_cpm: meds[k], n: EM.samples.filter((s) => s.tissue === tt).length })) } };
      },
      source: "site/data/expr_motrpac.json (log2 CPM via io.stack_tissues on data/raw/counts/TRNSCRPT__*.csv)",
      notShow: "the marker's selection frequency or effect size (annotation card); vials are jittered within tissue.",
    };
    const specB = {
      title: `${g.symbol} in the rat BodyMap: by organ and age`,
      subtitle: "Every sample (316), colour = age (2 → 104 weeks, light → dark); thymus and uterus have no MoTrPAC tissue.",
      build: () => {
        const i = EB.genes.indexOf(g.id);
        if (i < 0) return { traces: [], layout: { annotations: [{ text: "gene absent from the BodyMap matrix", showarrow: false }] } };
        const vals = EB.values[i];
        const organs = [...new Set(EB.samples.map((s) => s.organ))].sort();
        const xi = Object.fromEntries(organs.map((o, k) => [o, k]));
        const seq = ["--s2", "--s3", "--s5", "--s7"].map((v) => getComputedStyle(document.documentElement).getPropertyValue(v).trim());
        const traces = AGE_ORDER.map((a, ai) => {
          const idx = EB.samples.map((s, k) => (s.age_weeks === a ? k : -1)).filter((k) => k >= 0);
          return { ...strip(idx.map((k) => xi[EB.samples[k].organ] + (ai - 1.5) * 0.17 + jitter(k, 0, 0.06)), idx.map((k) => vals[k]), { color: seq[ai], size: 7, name: `${a} weeks`, hover: "%{customdata}<extra>" + a + " weeks</extra>" }), customdata: idx.map((k) => `${EB.samples[k].organ} · ${EB.samples[k].sex} · ${a} w<br>log2 CPM ${vals[k].toFixed(2)}`) };
        });
        return { traces, layout: { xaxis: { tickvals: organs.map((_, k) => k), ticktext: organs, tickangle: -45, tickfont: { size: 10 } }, yaxis: { title: { text: "log2 CPM" } }, margin: { t: 40, b: 70 }, legend: { y: 1.16 } },
                 table: { columns: ["organ", "age_weeks", "mean_log2_cpm", "n"], rows: organs.flatMap((o) => AGE_ORDER.map((a) => { const v = EB.samples.map((s, k) => (s.organ === o && s.age_weeks === a ? vals[k] : NaN)).filter(Number.isFinite); return { organ: o, age_weeks: a, mean_log2_cpm: mean(v), n: v.length }; })) } };
      },
      source: "site/data/expr_bodymap.json (log2 CPM from data/external/bodymap_counts.csv)",
      notShow: "animal identity (assumed from the replicate index) or technical replicates (summed).",
    };
    const specG = {
      title: `${g.symbol} in GTEx: log2 TPM by tissue`,
      subtitle: `Every sample (2,485, ≤ 150 donors per tissue) through the 1:1 ortholog${g.human_gene ? " " + g.human_gene : ""}; bars = tissue median.`,
      build: () => {
        const i = EG.genes.indexOf(g.id);
        if (i < 0) return { traces: [], layout: { annotations: [{ text: "no 1:1 human ortholog in GTEx for this gene", showarrow: false, font: { color: t.ink2 } }], xaxis: { visible: false }, yaxis: { visible: false } } };
        const vals = EG.values[i];
        const tissues = [...new Set(EG.samples.map((s) => s.tissue))].sort();
        const xi = Object.fromEntries(tissues.map((tt, k) => [tt, k]));
        const xs = [], ys = [], hov = [];
        EG.samples.forEach((s, k) => { xs.push(xi[s.tissue] + jitter(k, 0, 0.3)); ys.push(vals[k]); hov.push(`${s.tissue}<br>log2 TPM ${vals[k].toFixed(2)}`); });
        const meds = tissues.map((tt) => median(EG.samples.map((s, k) => (s.tissue === tt ? vals[k] : NaN))));
        return { traces: [
          { ...strip(xs, ys, { slot: 1, size: 5, opacity: 0.45, name: "samples", hover: "%{customdata}<extra></extra>" }), customdata: hov, showlegend: false },
          { type: "scatter", mode: "markers", x: tissues.map((_, k) => k), y: meds, name: "tissue median", marker: { symbol: "line-ew", size: 22, color: t.ink, line: { width: 2, color: t.ink } }, hovertemplate: "%{text}: median %{y:.2f}<extra></extra>", text: tissues },
        ], layout: { xaxis: { tickvals: tissues.map((_, k) => k), ticktext: tissues.map((tt) => tt.replace(" - ", "<br>")), tickangle: -60, tickfont: { size: 9 } }, yaxis: { title: { text: "log2(TPM + 1)" } }, margin: { t: 30, b: 110 }, showlegend: false },
          table: { columns: ["tissue", "median_log2_tpm", "n"], rows: tissues.map((tt, k) => ({ tissue: tt, median_log2_tpm: meds[k], n: EG.samples.filter((s) => s.tissue === tt).length })) } };
      },
      source: "site/data/expr_gtex.json (log2 TPM + 1 from data/external/gtex_tpm_subset.csv; orthologs from data/raw/rat_to_human_gene.csv, 1:1 only)",
      notShow: "donor-level structure; GTEx samples of one tissue come from different donors.",
    };
    for (const [key, id, spec] of [["m", "fig-gene-motrpac", specM], ["b", "fig-gene-bodymap", specB], ["g", "fig-gene-gtex", specG]]) {
      if (!figs[key]) figs[key] = await figure(document.getElementById(id), spec);
      else { figs[key].spec = spec; const b = spec.build(); figs[key].traces = b.traces; figs[key].table = b.table; await window.Plotly.react(figs[key].chart, b.traces, { ...template(), ...b.layout }, CONFIG); figs[key].root.querySelector("h3").textContent = spec.title; figs[key].root.querySelector(".fig-sub").textContent = spec.subtitle; }
    }
  }
  render();
}

// ---- (c) the panel builder -----------------------------------------------------------------------------
function panelBuilder() {
  const { PC, GENES } = DATA;
  const grid = PC.curve.map((r) => r.k);
  const order = PC.roundrobin_class_order;
  const gi = Object.fromEntries(GENES.genes.map((g) => [g.id, g]));
  const state = { ki: grid.indexOf(20), fold: 0 };
  const ctl = document.getElementById("builder-controls");
  const s = slider(0, grid.length - 1, 1, state.ki, (v) => { state.ki = v; render(); }, (v) => `k = ${grid[v]} genes`);
  ctl.append(control("Panel size (the evaluated grid)", s.input, s.out), control("Fold whose selection is shown", select([0, 1, 2, 3, 4].map((f) => [String(f), `fold ${f}`]), "0", (v) => { state.fold = Number(v); render(); })));
  let fig;
  async function render() {
    const k = grid[state.ki];
    const row = PC.curve[state.ki];
    const genes = PC.selected_by_fold[String(k)][String(state.fold)];
    const covered = order.slice(0, Math.min(k, order.length));
    document.getElementById("builder-summary").replaceChildren(
      el("div", { class: "who" }, [el("b", {}, `k = ${k}`), ` · balanced accuracy ${fmt(row.roundrobin_mean)} ± ${fmt(row.roundrobin_sd)} (5 folds, ${row.n_train_animals} train / ${row.n_test_animals} test animals) · F-test selector at the same k: ${fmt(row.fclassif_mean)}`]),
      el("p", { class: "explain" }, `${covered.length} of 19 tissues have a marker after ${k} picks (one tissue per pick, classes in sorted order).`),
    );
    const stripEl = document.getElementById("builder-strip");
    stripEl.replaceChildren(...order.map((t) => el("div", { class: "cell" + (covered.includes(t) ? " on" : ""), title: tissueLabel(t) }, t)));
    const chips = document.getElementById("builder-genes");
    chips.replaceChildren(...genes.map((g) => { const a = gi[g]; const sym = PC.gene_symbols[g] || a?.symbol || g; const mt = a?.marker_tissue || PC.marker_tissue_known[g]; const os = mt ? organSystem(mt) : null;
      return el("span", { class: "chip", title: mt ? `marker of ${mt}` : "marker tissue not annotated for this gene" }, [os ? el("span", { class: "dot", style: `background:${os.color}` }) : null, sym + (mt ? ` · ${mt}` : "")]); }));
    const spec = {
      title: "Accuracy vs panel size, with the current k marked",
      subtitle: "Mean ± sd balanced accuracy over 5 animal-grouped folds, round-robin selector vs F-test; log x.",
      build: () => {
        const t = tokens();
        const ks = PC.curve.map((r) => r.k);
        return { traces: [
          line(ks, PC.curve.map((r) => r.roundrobin_mean), { name: "round-robin", slot: 1, sd: PC.curve.map((r) => r.roundrobin_sd || 0) }),
          line(ks, PC.curve.map((r) => r.fclassif_mean), { name: "F-test", slot: 2, sd: PC.curve.map((r) => r.fclassif_sd || 0) }),
        ], layout: { xaxis: { type: "log", tickvals: ks, ticktext: ks.map(String), title: { text: "k" } }, yaxis: { range: [0, 1.05], title: { text: "balanced accuracy" } },
                     shapes: [{ type: "line", x0: k, x1: k, xref: "x", y0: 0, y1: 1, yref: "paper", line: { color: t.ink2, width: 1, dash: "dot" } }], margin: { t: 40 } },
          table: { columns: ["k", "roundrobin_mean", "roundrobin_sd", "fclassif_mean", "fclassif_sd"], rows: PC.curve } };
      },
      source: "results/05_panels/TRNSCRPT/panel_curve.csv, panel_curve_fclassif.csv, selected_by_fold.csv",
      notShow: "which gene was picked for which tissue: the pipeline stores the selected ids per fold, not the pick order; the tissue strip follows the selector's construction (verified on the annotated genes).",
      height: "short",
    };
    if (!fig) fig = await figure(document.getElementById("fig-builder"), spec);
    else { fig.spec = spec; const b = spec.build(); fig.traces = b.traces; fig.table = b.table; await window.Plotly.react(fig.chart, b.traces, { ...template(), ...b.layout }, CONFIG); }
  }
  render();
}

// ---- (d) the animals-needed calculator -----------------------------------------------------------------
function calculator() {
  const { CERT } = DATA;
  const state = { alpha: 0.10, delta: 0.10 };
  const ctl = document.getElementById("calc-controls");
  const a = slider(0.01, 0.30, 0.01, state.alpha, (v) => { state.alpha = v; render(); }, (v) => `α = ${v.toFixed(2)}`);
  const d = slider(0.01, 0.30, 0.01, state.delta, (v) => { state.delta = v; render(); }, (v) => `δ = ${v.toFixed(2)}`);
  ctl.append(control("Error level α", a.input, a.out), control("Confidence level δ (certificate holds with probability 1 − δ)", d.input, d.out));
  function render() {
    const n = animalsNeeded(state.alpha, state.delta);
    const out = document.getElementById("calc-out");
    const refs = CERT.sizing.map((r) => ({ ...r, agrees: animalsNeeded(r.alpha, r.delta) === r.n_zero_error_needed }));
    out.replaceChildren(
      el("div", { class: "who" }, [el("b", {}, `n ≥ ${n} calibration animals`), ` with zero calibration errors certify a panel at (α, δ) = (${state.alpha.toFixed(2)}, ${state.delta.toFixed(2)}).`]),
      el("p", { class: "explain" }, `This study calibrates on ${CERT.sizing[0].n_cal_this_run} animals per fold (one vial per animal).`),
      tableFrom({ columns: ["alpha", "delta", "n_zero_error_needed", "n_cal_this_run", "formula_agrees"], rows: refs.map((r) => ({ alpha: r.alpha, delta: r.delta, n_zero_error_needed: r.n_zero_error_needed, n_cal_this_run: r.n_cal_this_run, formula_agrees: r.agrees ? "yes" : "no" })) }),
      el("details", { class: "fig-notes" }, [el("summary", {}, "Source"), el("p", {}, [el("b", {}, "Source: "), "results/06_conformal/TRNSCRPT/sizing_table.csv (the four reference points); the slider applies n = ⌈ln δ / ln(1 − α)⌉, the same formula the pipeline uses."])]),
    );
  }
  render();
}

main().catch((e) => { console.error(e); document.querySelector(".lede").textContent = "Failed to load: " + e.message; });
