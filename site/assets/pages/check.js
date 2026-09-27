// "Check samples" (site/index.html). The DOM side of the product: input, banners, results table, reference map, sample
// drawer, calibration to a lab, downloads, and the two evidence notes. The logic lives in ../check-core.js (pure, tested
// by tests/test_check_core.js); every number shown comes from site/data (panel_model.json, headline.json,
// product.json with its pv_* provenance entries, expr_motrpac.json, expr_bodymap.json).
import { mountChrome, loadJSON, el, fmt, pct, callout, tableFrom, control, slider, select, badge } from "../site.js";
import { figure, tokens, palette, organSystem, template, CONFIG } from "../charts.js";
import { exampleRows, recalibrate, minLabelled, templateCsv } from "../score.js";
import { compositionNote, modelNote } from "../notes.js";
import { parseUpload, runCheck, resultsToCsv, exampleTable, acceptedList, CLASS_NAMES, contrast, explainSentence, project, referenceStats, orderByPriority, labelledForRecal, setLabel, sdLabel, exampleInputTables,
         WITHIN_MIN_SAMPLES, WITHIN_MIN_TISSUES } from "../check-core.js";

const nm = (c) => CLASS_NAMES[c] || c;
const pctA = () => `${Math.round(100 * (1 - state.alpha))} %`;   // 1 − α as a percentage, e.g. "90 %"
const state = { model: null, samples: null, parsed: null, check: null, alpha: 0.1, scaling: "auto", recal: null, filterFlagged: false, sort: { key: null, dir: 1 }, selected: null, source: null, swapped: [] };
let D = {};

async function main() {
  await mountChrome("index.html");
  const [model, H, PR, PC] = await Promise.all([loadJSON("data/panel_model.json"), loadJSON("data/headline.json"), loadJSON("data/product.json"), loadJSON("data/panel_curve.json")]);
  state.model = model; state.alpha = model.alpha_default ?? 0.1;
  D = { H, PR, PC };
  document.getElementById("how-min-labelled").textContent = String(minLabelled(model.alpha_default ?? 0.1));
  trustStrip();
  document.getElementById("model-note").textContent = modelNote(model, H.extras.n_animals);
  acceptedNames();
  expectNote();
  flagNote();
  document.getElementById("btn-example").addEventListener("click", runExample);
  document.getElementById("btn-upload").addEventListener("click", () => { document.getElementById("input").scrollIntoView({ behavior: "smooth" }); document.getElementById("dropzone").focus({ preventScroll: true }); });
  document.getElementById("btn-load-example").addEventListener("click", runExample);
  dropZone();
  exampleTables();
  acceptedLink();
  document.getElementById("btn-check").addEventListener("click", () => {
    const text = document.getElementById("in-text").value;
    if (!text.trim()) { setProgress("Paste a table first, or choose a file."); return; }
    state.swapped = []; document.getElementById("example-note").replaceChildren();
    parseInWorker({ text }, "pasted table");
  });
  document.getElementById("in-file").addEventListener("change", (e) => { const f = e.target.files && e.target.files[0]; if (f) readFile(f); });
  document.getElementById("btn-template").addEventListener("click", () => download("tissue_check_template.csv", templateCsv(model).replace(",true_tissue", ",claimed_tissue"), "text/csv"));
  document.getElementById("drawer-close").addEventListener("click", closeDrawer);
  document.addEventListener("keydown", (e) => { if (e.key === "Escape" && !document.getElementById("drawer").hidden) closeDrawer(); });
  if (new URLSearchParams(location.search).get("example") === "1") runExample();
}

// ---- input: one file at a time, from the file picker or the drop zone (same parse path) ----------------------------------
function readFile(f, extra = 0) {
  state.swapped = []; document.getElementById("example-note").replaceChildren();
  parseInWorker({ file: f }, f.name);
  if (extra) setTimeout(() => setProgress(`${document.getElementById("progress").textContent} (${extra} more file${extra === 1 ? " was" : "s were"} dropped and ignored: one file at a time.)`), 0);
}

function dropZone() {
  const dz = document.getElementById("dropzone"), input = document.getElementById("in-file"), choose = document.getElementById("choose-file");
  dz.setAttribute("tabindex", "0"); dz.setAttribute("role", "button");
  dz.setAttribute("aria-label", "Drop a CSV or TSV file here, or press Enter to choose a file");
  const open = (e) => { if (e.key === "Enter" || e.key === " ") { e.preventDefault(); input.click(); } };
  dz.addEventListener("keydown", (e) => { if (e.target === dz) open(e); });
  choose.setAttribute("tabindex", "0"); choose.setAttribute("role", "button");
  choose.addEventListener("keydown", open);
  const on = (e) => { e.preventDefault(); dz.classList.add("dragover"); };
  const off = () => dz.classList.remove("dragover");
  dz.addEventListener("dragenter", on); dz.addEventListener("dragover", on);
  dz.addEventListener("dragleave", (e) => { if (!dz.contains(e.relatedTarget)) off(); });
  dz.addEventListener("drop", (e) => {
    e.preventDefault(); off();
    const files = [...((e.dataTransfer && e.dataTransfer.files) || [])];
    if (files.length) readFile(files[0], files.length - 1);
  });
}

// ---- "Show example input tables": built from the example rows and the exported real counts --------------------------------
function exampleTables() {
  const det = document.getElementById("example-tables");
  let built = false;
  det.addEventListener("toggle", async () => {
    if (!det.open || built) return;
    built = true;
    const [expr, bm] = await Promise.all([loadJSON("data/expr_bodymap.json"), loadJSON("data/bodymap.json")]);
    const rows = exampleRows(state.model, expr, bm.organ_map, { perOrgan: 8, nLabelled: 0 });
    const T = exampleInputTables(state.model, rows, D.PR.example_counts);
    document.getElementById("ex-table-a").replaceChildren(tableFrom(T.a));
    document.getElementById("ex-table-a-src").textContent = `Two samples of the example (rat BodyMap adults), the first ${T.a.columns.length - 3} panel genes, values rounded to 2 decimals; a real table has all ${state.model.genes.length}.`;
    document.getElementById("ex-table-b").replaceChildren(tableFrom(T.b));
    document.getElementById("ex-table-b-src").textContent = `Real counts: ${D.PR.example_counts.source}. A real matrix has one row (or column) per gene, every gene.`;
  });
}

// ---- the "list of accepted names" link opens and shows the list ---------------------------------------------------------
function acceptedLink() {
  const show = () => { const d = document.getElementById("accepted"); d.open = true; d.scrollIntoView({ behavior: "smooth", block: "start" }); };
  document.getElementById("accepted-link").addEventListener("click", (e) => { e.preventDefault(); show(); history.replaceState(null, "", "#accepted"); });
  if (location.hash === "#accepted") setTimeout(show, 0);
  window.addEventListener("hashchange", () => { if (location.hash === "#accepted") show(); });
}

// ---- trust strip --------------------------------------------------------------------------------------------------
// percent with one decimal where needed: 97.6 %, 100 %
const pct1 = (v) => `${(100 * v).toFixed(1).replace(/\.0$/, "")} %`;

function trustStrip() {
  const t = Object.fromEntries(D.H.tiles.map((x) => [x.id, x]));
  const ood = D.PR.ood;
  const nOod = ood.Thymus.n + ood.Uterus.n;
  const nEmpty = Math.round(ood.Thymus.frac_empty * ood.Thymus.n) + Math.round(ood.Uterus.frac_empty * ood.Uterus.n);
  const bl = D.PC.baselines.map((b) => b.balanced_accuracy_mean);
  document.getElementById("trust").replaceChildren(
    el("b", {}, `${pct1(t.tile_acc_k20.value)} balanced accuracy`), ` across ${D.H.extras.n_tissues} rat tissues on held-out MoTrPAC animals, with the panel selected inside each fold (simple all-gene baselines on the same folds: ${pct1(Math.min(...bl))}–${pct1(Math.max(...bl))}); `,
    "the panel names ", el("b", {}, pct1(t.tile_bodymap_k20.value)), " of mapped adult organs in another laboratory's rats (rat BodyMap); and on organs it never saw it abstains (",
    el("b", {}, `${nEmpty} of ${nOod}`), " samples in a mixed upload). ",
    el("a", { href: "science.html#sec-accuracy" }, "How it was measured"), " · ", el("a", { href: "transfer.html" }, "Transfer to another lab and species"), ".",
  );
}

function acceptedNames() {
  document.getElementById("accepted-list").replaceChildren(
    (() => { const t = tableFrom({ columns: ["code", "tissue", "also accepted"], rows: acceptedList().map((a) => ({ code: a.code, tissue: a.name, "also accepted": a.also.join(", ") })) }); t.classList.add("wrap-last"); return t; })(),
    el("p", { class: "small" }, "Also accepted: brain (cortex, hippocampus or hypothalamus), skeletal muscle or muscle (either muscle), adipose or fat (brown or white), intestine (colon or small intestine). Any other name is reported as “tissue not in reference”."));
}

// ---- input ----------------------------------------------------------------------------------------------------------
function setProgress(msg) { document.getElementById("progress").textContent = msg; }

async function runExample() {
  const [expr, bm] = await Promise.all([loadJSON("data/expr_bodymap.json"), loadJSON("data/bodymap.json")]);
  const rows = exampleRows(state.model, expr, bm.organ_map, { perOrgan: 8, nLabelled: 0 });
  const exTbl = exampleTable(state.model, rows);
  state.swapped = exTbl.swapped;
  document.getElementById("in-text").value = exTbl.csv;
  document.getElementById("in-format").value = "auto";
  const organs = [...new Set(rows.map((r) => r.organ))];
  const unseen = organs.filter((o) => !bm.organ_map[o]);
  document.getElementById("example-note").replaceChildren(callout("note", "The example", [
    `${rows.length} rat BodyMap samples from 21-week-old animals (Yu et al. 2014, GEO GSE53960: another laboratory, processed independently), ${organs.length} organs, as 20-gene log2 CPM with the organ as the claimed tissue. `
    + `To demonstrate mismatch detection we SWAPPED the claimed labels of two samples: ${exTbl.swapped[0]} (really ${exTbl.organs[0].toLowerCase()}, labelled ${exTbl.organs[1].toLowerCase()}) and ${exTbl.swapped[1]} (really ${exTbl.organs[1].toLowerCase()}, labelled ${exTbl.organs[0].toLowerCase()}). `
    + `${unseen.map((o) => o.toLowerCase()).join(" and ").replace(/^./, (c) => c.toUpperCase())} are organs the model never saw.`,
  ]));
  parseInWorker({ text: exTbl.csv }, "the example");
}

function clearResults() {
  for (const id of ["results", "map-sec", "calibrate", "drawer"]) document.getElementById(id).hidden = true;
  document.getElementById("banners").replaceChildren();
  state.check = null; state.samples = null;
}

function parseInWorker(payload, sourceName) {
  state.source = sourceName;
  setProgress(`Reading ${sourceName}…`);
  const format = document.getElementById("in-format").value;
  const opts = { format };
  const done = (r) => {
    setProgress(`Read ${sourceName}: ${r.samples.length} sample${r.samples.length === 1 ? "" : "s"}, ${r.format === "counts" ? `${r.nFeatures} gene features (raw counts, converted to log2 CPM)` : "20-gene log2 CPM"}, ${r.orientation === "genes_x_samples" ? "genes in rows" : "samples in rows"}.`);
    state.parsed = r; state.samples = r.samples; state.recal = null; state.selected = null;
    closeDrawer();
    analyse();
  };
  let w = null;
  try { w = new Worker(new URL("../check-worker.js", import.meta.url), { type: "module" }); } catch (e) { w = null; }
  if (!w) {   // no module workers: parse on the main thread
    const go = async () => { const text = payload.text ?? (await payload.file.text()); const r = parseUpload(text, state.model, opts); if (r.error) { clearResults(); setProgress(`Could not read ${sourceName}: ${r.error}`); } else done(r); };
    go();
    return;
  }
  w.onmessage = (e) => {
    const m = e.data;
    if (m.type === "progress") setProgress(m.stage === "reading" ? `Reading ${sourceName}…` : `Parsing ${sourceName}: ${Math.round(100 * m.f)} %`);
    else if (m.type === "error") { clearResults(); setProgress(`Could not read ${sourceName}: ${m.message}`); w.terminate(); }
    else if (m.type === "done") { w.terminate(); done(m.result); }
  };
  w.onerror = (e) => { clearResults(); setProgress(`Could not read ${sourceName}: ${e.message || "worker error"}`); w.terminate(); };
  w.postMessage({ ...payload, model: state.model, opts });
}

// ---- analysis -------------------------------------------------------------------------------------------------------
function analyse() {
  if (!state.samples || !state.samples.length) return;
  state.check = runCheck(state.model, state.samples, { alpha: state.alpha, scaling: state.scaling, q: state.recal && Number.isFinite(state.recal.q) ? state.recal.q : null });
  banners();
  document.getElementById("results").hidden = false;
  document.getElementById("map-sec").hidden = false;
  document.getElementById("calibrate").hidden = false;
  resultControls();
  guide();
  summary();
  table();
  map();
  calibrateSection();
  if (state.selected) openDrawer(state.selected, false);
}

function banners() {
  const P = state.parsed, C = state.check, K = D.PR.scaling_key;
  const out = [];
  if (P.species === "human") out.push(callout("caveat", "These look like human genes", [el("span", {}, ["This is a rat model. Human samples can be scored through the matching gene symbols, but expect many abstentions and weaker calls; see the ", el("a", { href: "transfer.html" }, "GTEx results"), " for how the fingerprint does in human tissue."])]));
  else if (P.species === "mouse") out.push(callout("caveat", "These look like mouse genes", "This is a rat model; mouse genes are matched by symbol only. Expect weaker calls and more abstentions."));
  if (P.missing.length) out.push(callout("caveat", `${P.missing.length} of ${state.model.genes.length} panel genes are missing`, `Not found in your table: ${P.missing.map((g) => g.symbol).join(", ")}. A missing gene is scored at the reference mean, which weakens the calls.`));
  for (const w of P.warnings || []) out.push(callout("caveat", "Check the input", w));
  const s = C.scaling;
  if (s.mode === "reference") {
    out.push(callout("caveat", "Reference scaling", [
      `${s.reason.replace(/^./, (c) => c.toUpperCase())}${s.forced ? " (chosen by you; the automatic rule would scale within your set)" : ""}. `,
      `Measured on rat BodyMap adults from another laboratory, reference scaling named ${pct(K["adult_21wk.reference.accuracy"])} of mapped organs correctly with ${pct(K["adult_21wk.reference.coverage"])} coverage, but gave sets to unseen thymus that within-set scaling left empty. It assumes your log2 CPM is on the pipeline's scale (log2 CPM on the total library of bulk RNA-seq counts, as the pipeline computes it). For the most reliable calls, upload a set that mixes at least ${WITHIN_MIN_TISSUES} tissues and ${WITHIN_MIN_SAMPLES} samples.`,
    ]));
  } else {
    out.push(callout(s.caveat ? "caveat" : "note", "Within-set scaling", `${s.reason.replace(/^./, (c) => c.toUpperCase())}.${s.forced ? " Chosen by you; the automatic rule would use reference scaling." : ""}${s.caveat ? ` Caution: ${s.caveat}.` : ""}`));
  }
  if (state.recal) out.push(callout("note", "Calibrated to your lab", state.recal.error ? state.recal.error : `The prediction sets use a threshold recalibrated on ${state.recal.n} of your labelled samples (see “Calibrate to my lab”).`));
  document.getElementById("banners").replaceChildren(...out);
}

function resultControls() {
  const { input, out } = slider(0.05, 0.3, 0.01, state.alpha, (v) => { state.alpha = v; if (state.recal) state.recal = recalFor(); analyse(); }, (v) => `α = ${v.toFixed(2)} (${Math.round(100 * (1 - v))} % sets)`);
  input.id = "alpha"; input.setAttribute("aria-label", "error rate alpha");
  const sc = select([["auto", "automatic"], ["within", "within your set"], ["reference", "MoTrPAC reference"]], state.scaling, (v) => { state.scaling = v; analyse(); });
  sc.id = "scaling"; sc.setAttribute("aria-label", "gene scaling");
  document.getElementById("result-controls").replaceChildren(control("Error rate", input, out), control("Gene scaling", sc));
}

const STATUS_BADGE = { Confident: ["confident", "●"], Ambiguous: ["ambiguous", "◐"], Unknown: ["abstains", "∅"] };
const CLAIM_BADGE = { Consistent: ["correct", "✓"], Mismatch: ["wrong", "✗"], "Can't confirm": ["unmapped", "?"], "Not in reference": ["unmapped", "–"] };
const sBadge = (s) => badge(...STATUS_BADGE[s], s);
const cBadge = (s) => (s ? badge(...CLAIM_BADGE[s], s) : el("span", { class: "small" }, "no claim"));

// What to look at first: the swapped pair of the example, and a one-click recalibration when most sets are empty.
function guide() {
  const C = state.check, c = C.counts, out = [];
  const mism = C.results.filter((r) => r.claimStatus === "Mismatch");
  if (state.swapped.length) {
    const hit = state.swapped.filter((id) => mism.some((r) => r.id === id));
    out.push(el("p", {}, [el("b", {}, `The two swapped samples ${hit.length === state.swapped.length ? "are both flagged Mismatch" : `: ${hit.length} of ${state.swapped.length} flagged Mismatch`}`),
      ` and sit at the top of the table (${state.swapped.join(" and ")}); select one to see why.`]));
  } else if (mism.length) {
    out.push(el("p", {}, [el("b", {}, `${mism.length} claimed label${mism.length === 1 ? " does" : "s do"} not match the call`), "; they sit at the top of the table."]));
  }
  const nLab = Object.keys(labelledMap()).length, need = minLabelled(state.alpha);
  if (!state.recal && c.Unknown > c.total / 2) {
    const b = el("button", { class: "btn primary", type: "button" }, "Recalibrate on my labelled samples");
    if (nLab < need) b.setAttribute("disabled", "");
    b.addEventListener("click", () => { state.recal = recalFor(); analyse(); document.getElementById("results").scrollIntoView({ behavior: "smooth" }); });
    out.push(el("p", {}, [el("b", {}, `${c.Unknown} of ${c.total} prediction sets are empty.`),
      " That is the expected first result for samples from another laboratory: the threshold was calibrated on MoTrPAC animals, so the model abstains rather than guess. ",
      nLab >= need ? `Recalibrating on your ${nLab} labelled samples adapts the threshold to your lab (see “Calibrate to my lab” below)${state.swapped.length ? `; the ${state.swapped.length} deliberately swapped samples are left out, since their labels are wrong on purpose` : ""}. ` : `Label at least ${need} samples to recalibrate to your lab. `]), b);
  } else if (state.recal && !state.recal.error) {
    out.push(el("p", {}, `Recalibrated on ${state.recal.n} of your labelled samples: the flags below use your lab's threshold (observed coverage ${fmt(state.recal.coverage, 3)} on those samples, not a guarantee for new ones).`));
  }
  document.getElementById("results-guide").replaceChildren(...(out.length ? [callout("note", "Start here", out)] : []));
}

function summary() {
  const c = state.check.counts;
  const tile = (value, label, sub) => el("div", { class: "tile" }, [el("div", { class: "value" }, String(value)), el("div", { class: "label" }, label), sub ? el("div", { class: "sub" }, sub) : null]);
  const tiles = [tile(c.total, "samples checked", `${state.check.scaling.mode === "within" ? "within-set" : "reference"} scaling · α = ${state.alpha.toFixed(2)}`),
                 tile(c.flagged, "flagged", "mismatch, can't confirm, unknown tissue or not in reference"),
                 tile(`${c.Confident} · ${c.Ambiguous} · ${c.Unknown}`, "confident · ambiguous · unknown", `one tissue · several · none in the ${setLabel(state.alpha)}`)];
  if (c.claimed) tiles.push(tile(`${c.Consistent} · ${c.Mismatch} · ${c["Can't confirm"]}`, "consistent · mismatch · can't confirm", `of ${c.claimed} claimed labels${c["Not in reference"] ? `; ${c["Not in reference"]} not in the reference` : ""}`));
  document.getElementById("summary").replaceChildren(...tiles);
}

function table() {
  const C = state.check;
  const onlyFlag = el("input", { type: "checkbox", id: "flagged-only" });
  onlyFlag.checked = state.filterFlagged;
  onlyFlag.addEventListener("change", () => { state.filterFlagged = onlyFlag.checked; table(); });
  const dl = el("button", { class: "btn", type: "button" }, "Download CSV");
  dl.addEventListener("click", () => download("tissue_check.csv", resultsToCsv(C, state.model), "text/csv"));
  const rp = el("button", { class: "btn", type: "button" }, "Download report");
  rp.addEventListener("click", () => download("tissue_check_report.html", reportHtml(), "text/html"));
  document.getElementById("table-controls").replaceChildren(el("label", { class: "control inline", for: "flagged-only" }, [onlyFlag, el("span", {}, "Flagged only")]), dl, rp);
  const cols = [["id", "sample"], ["label", "claimed"], ["call", "call"], ["set", setLabel(state.alpha)], ["status", "status"], ["claimStatus", "claim"], ["callProb", "top p"], ["runnerUp", "runner-up"], ...(C.results.some((r) => r.missingGenes.length) ? [["missing", "missing"]] : [])];
  const hasMissing = cols.some(([k]) => k === "missing");
  let rows = orderByPriority(C.results.filter((r) => !state.filterFlagged || r.flagged));
  if (state.sort.key) {
    const k = state.sort.key, d = state.sort.dir;
    const val = (r) => (k === "set" ? r.set.length : k === "runnerUp" ? r.runnerUpProb : k === "missing" ? r.missingGenes.length : r[k] ?? "");
    rows = [...rows].sort((a, b) => (val(a) < val(b) ? -d : val(a) > val(b) ? d : 0));
  }
  const thead = el("thead", {}, el("tr", {}, cols.map(([k, lab]) => {
    const b = el("button", { type: "button", class: "sort" }, lab + (state.sort.key === k ? (state.sort.dir > 0 ? " ▲" : " ▼") : ""));
    b.addEventListener("click", () => { state.sort = { key: k, dir: state.sort.key === k ? -state.sort.dir : 1 }; table(); });
    const th = el("th", { scope: "col", class: k === "callProb" ? "num" : "" }, b);
    if (state.sort.key === k) th.setAttribute("aria-sort", state.sort.dir > 0 ? "ascending" : "descending");
    return th;
  })));
  const tbody = el("tbody");
  for (const r of rows) {
    const tr = el("tr", { tabindex: "0", class: "clickable" + (r.flagged ? " flagged" : "") + (state.selected === r.id ? " selected" : ""), "aria-label": `${r.id}: ${r.status}${r.claimStatus ? ", " + r.claimStatus : ""}. Open details.` }, [
      el("td", {}, r.id), el("td", { class: "wrap" }, r.label ?? ""), el("td", { class: "wrap" }, nm(r.call)),
      el("td", { class: "wrap" }, r.set.length ? r.set.map(nm).join("; ") : "none"),
      el("td", {}, sBadge(r.status)), el("td", {}, cBadge(r.claimStatus)),
      el("td", { class: "num" }, fmt(r.callProb, 2)), el("td", { class: "wrap" }, `${nm(r.runnerUp)} (${fmt(r.runnerUpProb, 2)})`),
      hasMissing ? el("td", {}, r.missingGenes.length ? el("span", { class: "missing-mark", title: `missing: ${r.missingGenes.join(", ")}` }, `⚠ ${r.missingGenes.length} gene${r.missingGenes.length === 1 ? "" : "s"}`) : "") : null,
    ]);
    const open = () => openDrawer(r.id, true);
    tr.addEventListener("click", open);
    tr.addEventListener("keydown", (e) => { if (e.key === "Enter" || e.key === " ") { e.preventDefault(); open(); } });
    tbody.appendChild(tr);
  }
  const wrap = el("div", { class: "table-wrap" }, el("table", { class: "data results" }, [thead, tbody]));
  document.getElementById("tbl-results").replaceChildren(rows.length ? wrap : el("p", {}, "No flagged samples."));
  document.getElementById("status-key").textContent = `Rows are ordered with the most actionable first (Mismatch, Can't confirm, not in reference, Unknown); click a column header to sort. The ${setLabel(state.alpha)} is a conformal prediction set: every tissue whose probability clears a threshold calibrated on held-out MoTrPAC animals (see the glossary on The science page). Status: Confident = one tissue in the set; Ambiguous = several; Unknown = none (the model abstains). Claim: Consistent = your label is in the set; Mismatch = the set is not empty and your label is not in it; Can't confirm = empty set. Select a row (click, or Enter) for the explanation.`;
}

// ---- reference map ------------------------------------------------------------------------------------------------------
async function map() {
  const pca = D.PR.pca, C = state.check;
  await figure(document.getElementById("fig-map"), {
    title: "Where your samples sit among the MoTrPAC reference tissues",
    subtitle: `The first two principal components of the 20-gene z-space fitted on ${pca.reference.x.length} MoTrPAC reference vials (${pct(pca.explained[0])} and ${pct(pca.explained[1])} of the variance); reference vials coloured by organ system, your samples as outlined markers (diamonds are flagged).`,
    build: () => {
      const t = tokens();
      const byT = {};
      pca.reference.tissue.forEach((tt, i) => { (byT[tt] ||= []).push(i); });
      const traces = Object.entries(byT).sort().map(([tt, idx]) => ({ type: "scatter", mode: "markers", name: nm(tt), x: idx.map((i) => pca.reference.x[i]), y: idx.map((i) => pca.reference.y[i]),
        marker: { size: 5, color: organSystem(tt).color, opacity: 0.45 }, hovertemplate: `${nm(tt)} (reference)<extra></extra>`, legendgroup: organSystem(tt).name }));
      const pts = C.results.map((r) => ({ r, xy: project(r.z, pca) }));
      for (const [flag, sym, name] of [[false, "circle-open", "your samples"], [true, "diamond-open", "your samples, flagged"]]) {
        const sel = pts.filter((p) => p.r.flagged === flag);
        if (!sel.length) continue;
        traces.push({ type: "scatter", mode: "markers", name, x: sel.map((p) => p.xy[0]), y: sel.map((p) => p.xy[1]), customdata: sel.map((p) => p.r.id),
          marker: { size: 11, symbol: sym, color: t.ink, line: { width: 2, color: t.ink } },
          hovertemplate: sel.map((p) => `${p.r.id}<br>call ${nm(p.r.call)} · ${p.r.status}${p.r.claimStatus ? " · " + p.r.claimStatus : ""}<extra></extra>`) });
      }
      return { traces, layout: { xaxis: { title: { text: "PC1" }, zeroline: false }, yaxis: { title: { text: "PC2" }, zeroline: false }, legend: { font: { size: 10 } }, margin: { t: 20 } },
               table: { columns: ["sample", "PC1", "PC2", "call", "status"], rows: pts.map((p) => ({ sample: p.r.id, PC1: p.xy[0], PC2: p.xy[1], call: p.r.call, status: p.r.status })) } };
    },
    source: "results_product/40_product/pca_loadings.csv and pca_reference_coords.csv (scripts/40_product_validation.py, from site/data/expr_motrpac.json and the model's MoTrPAC z-scoring); your samples are projected in the browser (tests/test_check_core.js checks the projection against Python).",
    notShow: `how the model decides: two components carry only part of the variance of 20 genes, so a sample can sit between clouds and still be called confidently. Under within-set scaling your samples' z-scores are centred on your set, not on MoTrPAC.`,
    height: "tall",
  });
  const chart = document.querySelector("#fig-map .chart");
  if (chart && chart.on) chart.on("plotly_click", (ev) => { const id = ev.points && ev.points[0] && ev.points[0].customdata; if (id) openDrawer(id, true); });
}

// ---- drawer -----------------------------------------------------------------------------------------------------------
let REF = null;
async function openDrawer(id, focus) {
  const r = state.check.results.find((x) => x.id === id);
  if (!r) return;
  state.selected = id;
  document.querySelectorAll("#tbl-results tr.selected").forEach((t) => t.classList.remove("selected"));
  const dr = document.getElementById("drawer");
  dr.hidden = false;
  document.getElementById("drawer-h").textContent = `Sample ${r.id}`;
  const claimCls = r.claim && r.claim.classes.length ? r.claim.classes : [];
  const other = claimCls.length && !claimCls.includes(r.call) ? claimCls[0] : r.runnerUp;
  document.getElementById("drawer-body").replaceChildren(...[
    el("div", { class: "drawer-badges" }, [sBadge(r.status), cBadge(r.claimStatus)]),
    el("p", {}, [el("b", {}, "Claimed: "), r.label ? `${r.label}${claimCls.length ? ` (${claimCls.map(nm).join(" or ")})` : ` (not one of the ${state.model.classes.length} reference tissues)`}` : "none", " · ", el("b", {}, "Call: "), `${nm(r.call)} (p = ${fmt(r.callProb, 2)})`, " · ", el("b", {}, `${setLabel(state.alpha).replace(/^./, (c) => c.toUpperCase())}: `), r.set.length ? r.set.map(nm).join("; ") : "empty"]),
    el("p", { class: "explain" }, explainSentence(r, state.model)),
    r.missingGenes.length ? el("p", { class: "small" }, `Missing genes (scored at the reference mean): ${r.missingGenes.join(", ")}.`) : null,
  ].filter(Boolean));
  if (focus) dr.focus();
  // probabilities vs the calibrated threshold
  const classes = state.model.classes;
  const order = r.p.map((v, k) => [v, k]).sort((a, b) => b[0] - a[0]).slice(0, 5).map((x) => x[1]);
  for (const c of claimCls) { const k = classes.indexOf(c); if (!order.includes(k)) order.push(k); }
  await figure(document.getElementById("fig-probs"), {
    title: "Probabilities and the calibrated threshold",
    subtitle: `Top five tissues (and the claimed one); a tissue enters the ${setLabel(state.alpha)} when its bar crosses the line.`,
    build: () => {
      const t = tokens(), p = palette();
      const ys = order.map((k) => nm(classes[k])).reverse(), xs = order.map((k) => r.p[k]).reverse();
      const inSet = order.map((k) => r.set.includes(classes[k])).reverse();
      const thr = Number.isFinite(r.threshold) ? r.threshold : 0;
      return { traces: [{ type: "bar", orientation: "h", x: xs, y: ys, marker: { color: inSet.map((b) => (b ? p[0] : t.muted)) }, text: xs.map((v, i) => `${fmt(v, 3)}${inSet[i] ? " · in set" : ""}`), textposition: "outside", hovertemplate: "%{y}: %{x:.3f}<extra></extra>" }],
               layout: { xaxis: { range: [0, 1.25], title: { text: "probability" } }, yaxis: { automargin: true }, showlegend: false, margin: { t: 10, l: 10 }, bargap: 0.35,
                         shapes: [{ type: "line", x0: thr, x1: thr, y0: -0.5, y1: ys.length - 0.5, line: { color: t.ink, width: 2, dash: "dash" } }],
                         annotations: [{ x: thr, y: ys.length - 0.5, text: `threshold ${fmt(thr, 3)}`, showarrow: false, yanchor: "bottom", font: { size: 11, color: t.ink2 } }] },
               table: { columns: ["tissue", "probability", "in_set"], rows: order.map((k) => ({ tissue: classes[k], probability: r.p[k], in_set: r.set.includes(classes[k]) ? "yes" : "no" })) } };
    },
    source: "site/data/panel_model.json (the 20-gene transfer model and its 15-animal calibration scores), scored in the browser; the threshold is 1 − q̂ at your α" + (state.recal ? ", recalibrated on your labelled samples" : ""),
    notShow: `classes outside the top five unless claimed; the set is a ${pctA()} guarantee over samples like the calibration animals, not a per-sample probability.`,
    height: "short",
  });
  // why X, not Y
  const ct = contrast(r.z, state.model, r.call, other);
  const terms = [...ct.terms].sort((a, b) => Math.abs(b.diff) - Math.abs(a.diff));
  await figure(document.getElementById("fig-contrib"), {
    title: `Why ${nm(r.call)}, not ${nm(other)}`,
    subtitle: `Per gene: coefficient × z-score for ${nm(r.call)} minus the same for ${nm(other)}; bars to the right favour ${nm(r.call)}. With the intercepts they sum to the difference in log-odds (${fmt(ct.logitDiff, 2)}).`,
    build: () => {
      const t = tokens(), p = palette();
      const ys = terms.map((x) => `${x.symbol} (${x.z >= 0 ? "+" : ""}${fmt(x.z, 1)} sd)`).reverse(), xs = terms.map((x) => x.diff).reverse();
      return { traces: [{ type: "bar", orientation: "h", x: xs, y: ys, marker: { color: xs.map((v) => (v >= 0 ? p[0] : p[3])) }, hovertemplate: "%{y}: %{x:.2f}<extra></extra>" }],
               layout: { xaxis: { title: { text: `favours ${nm(other)} ←  → favours ${nm(r.call)}` }, zeroline: true, zerolinecolor: t.ink2 }, yaxis: { automargin: true, tickfont: { size: 10 } }, showlegend: false, margin: { t: 10, l: 10 }, bargap: 0.25 },
               table: { columns: ["gene", "z", `contribution_${r.call}`, `contribution_${other}`, "difference"], rows: terms.map((x) => ({ gene: x.symbol, z: x.z, [`contribution_${r.call}`]: x.value, [`contribution_${other}`]: x.value - x.diff, difference: x.diff })) } };
    },
    source: "the model's coefficients (site/data/panel_model.json) times this sample's z-scores; contributions plus intercept equal the logit (tests/test_check_core.js)",
    notShow: "causation: a coefficient says how the classifier uses a gene, not why the tissue expresses it; genes are correlated, so contributions are shared among them.",
    height: "tall",
  });
  // gene values vs the reference
  if (!REF) REF = referenceStats(await loadJSON("data/expr_motrpac.json"), state.model);
  const sample = state.samples.find((s) => s.id === r.id);
  const pair = [r.call, ...(other !== r.call ? [other] : [])];
  await figure(document.getElementById("fig-genes"), {
    title: "This sample's genes beside the reference tissues",
    subtitle: `log2 CPM of each panel gene: this sample (dots) against the MoTrPAC median and interquartile range of ${pair.map(nm).join(" and ")}.`,
    build: () => {
      const t = tokens(), p = palette();
      const genes = state.model.genes;
      const ys = genes.map((g) => g.symbol);
      const traces = pair.map((c, k) => ({ type: "scatter", mode: "markers", name: `${nm(c)} (MoTrPAC median, IQR)`, y: ys, x: genes.map((g) => REF[c][g.id]?.median ?? null),
        error_x: { type: "data", symmetric: false, array: genes.map((g) => (REF[c][g.id] ? REF[c][g.id].q3 - REF[c][g.id].median : 0)), arrayminus: genes.map((g) => (REF[c][g.id] ? REF[c][g.id].median - REF[c][g.id].q1 : 0)), color: k ? p[3] : p[0], thickness: 3, width: 0 },
        marker: { symbol: k ? "square" : "diamond", size: 8, color: k ? p[3] : p[0] }, hovertemplate: `%{y}: median %{x:.2f}<extra>${nm(c)}</extra>` }));
      traces.push({ type: "scatter", mode: "markers", name: "this sample", y: ys, x: genes.map((g) => sample.values[g.id]), marker: { size: 10, symbol: "circle-open", color: t.ink, line: { width: 2 } }, hovertemplate: "%{y}: %{x:.2f}<extra>this sample</extra>" });
      return { traces, layout: { xaxis: { title: { text: "log2 CPM" } }, yaxis: { automargin: true, autorange: "reversed", tickfont: { size: 10 } }, legend: { orientation: "h", y: 1.08 }, margin: { t: 30, l: 10 } },
               table: { columns: ["gene", "sample", ...pair.flatMap((c) => [`${c}_median`, `${c}_q1`, `${c}_q3`])], rows: genes.map((g) => Object.fromEntries([["gene", g.symbol], ["sample", sample.values[g.id]], ...pair.flatMap((c) => [[`${c}_median`, REF[c][g.id]?.median], [`${c}_q1`, REF[c][g.id]?.q1], [`${c}_q3`, REF[c][g.id]?.q3]])])) } };
    },
    source: `site/data/expr_motrpac.json (MoTrPAC log2 CPM, total library, ${Object.values(REF).reduce((a, t) => a + (Object.values(t)[0]?.n || 0), 0)} vials); the sample's own log2 CPM as uploaded or converted from counts`,
    notShow: "batch or depth differences between your lab and MoTrPAC, which shift whole columns; the classifier uses z-scores, not these raw levels.",
    height: "tall",
  });
  table();
}

function closeDrawer() {
  const dr = document.getElementById("drawer");
  if (dr.hidden) return;
  dr.hidden = true;
  const id = state.selected;
  state.selected = null;
  if (id && state.check) table();
}

// ---- calibrate to my lab -------------------------------------------------------------------------------------------------
function labelledMap() {
  return labelledForRecal(state.check.results, state.swapped);   // the example's two swapped labels are wrong on purpose
}
function recalFor() {
  // thresholds come from the samples scored with the MoTrPAC calibration (q override off)
  const base = runCheck(state.model, state.samples, { alpha: state.alpha, scaling: state.scaling });
  const rc = recalibrate(state.model, base.results, labelledMap(), state.alpha);
  if (rc.error) return { error: rc.error };
  return { ...rc, error: Number.isFinite(rc.q) ? null : `At α = ${state.alpha.toFixed(2)} a finite threshold needs at least ${minLabelled(state.alpha)} labelled samples; with ${rc.n} every set would hold all tissues.` };
}
function calibrateSection() {
  const nLab = Object.keys(labelledMap()).length, need = minLabelled(state.alpha);
  document.getElementById("cal-intro").textContent = `The ${pctA()} guarantee holds for samples like the MoTrPAC calibration animals. For a new lab, recalibrate the threshold on your own samples with verified labels: here, the samples whose claimed label names exactly one of the ${state.model.classes.length} tissues (${nLab} in this upload). At α = ${state.alpha.toFixed(2)} this needs at least ${need} labelled samples. The recalibrated coverage is observed on your labelled samples, not guaranteed for new samples. When the labelled samples come from a few animals (several tissues each), the coverage is pooled within those animals and is not a guarantee for new ones. Recalibrate only on labels you trust: a swapped label enlarges every set.${state.swapped.length ? ` In the example the ${state.swapped.length} deliberately swapped samples (${state.swapped.join(", ")}) are left out of the recalibration.` : ""}`;
  const b = el("button", { class: "btn", type: "button" }, state.recal ? "Recalibrate again" : "Recalibrate on my labelled samples");
  if (nLab < need) b.setAttribute("disabled", "");
  b.addEventListener("click", () => { state.recal = recalFor(); analyse(); });
  const reset = el("button", { class: "btn", type: "button" }, "Use the MoTrPAC calibration");
  if (!state.recal) reset.setAttribute("disabled", "");
  reset.addEventListener("click", () => { state.recal = null; analyse(); });
  document.getElementById("cal-controls").replaceChildren(b, reset);
  const c = state.check.counts;
  if (!state.recal && c.Unknown > c.total / 2) document.getElementById("cal-intro").prepend(el("b", {}, `${c.Unknown} of ${c.total} sets are empty. `), "That is typical for samples from another laboratory: the threshold was calibrated on MoTrPAC animals, and the model abstains rather than guess. ");
  const st = document.getElementById("cal-status");
  if (!state.recal) st.textContent = nLab < need ? `${nLab} labelled samples: add labels to at least ${need} to recalibrate.` : `Using the MoTrPAC calibration (${state.model.calibration.n_animals} animals).`;
  else if (state.recal.error) st.textContent = state.recal.error;
  else st.textContent = `Recalibrated on ${state.recal.n} labelled samples: threshold ${fmt(1 - state.recal.q, 3)} (a tissue enters the set when p ≥ that); coverage on the labelled samples ${fmt(state.recal.coverage, 3)}. The labelled samples are used both to set and to check the threshold, so this coverage is optimistic.`;
}

// ---- evidence notes -------------------------------------------------------------------------------------------------------
function expectNote() {
  const K = D.PR.scaling_key;
  const mappedOrgans = D.PR.scaling.filter((r) => r.subset === "adult_21wk" && r.mode === "within_all" && r.mapped && r.organ !== "ALL_MAPPED").length;
  const row = (label, mode) => ({ scaling: label, "organs named correctly": pct(K[`adult_21wk.${mode}.accuracy`]), "90 % sets that hold the organ": pct(K[`adult_21wk.${mode}.coverage`]), "empty sets": pct(K[`adult_21wk.${mode}.frac_empty`]) });
  document.getElementById("expect-body").replaceChildren(
    el("p", {}, `The model z-scores each gene before scoring. Within a mixed upload that is what the pipeline does for a new dataset; within a single-tissue upload it erases the very differences the model reads. We measured the three options on ${K["adult_21wk.within_all.n"]} rat BodyMap samples from 21-week-old animals of another laboratory (${mappedOrgans} organs that map to a reference tissue):`),
    tableFrom({ columns: ["scaling", "organs named correctly", "90 % sets that hold the organ", "empty sets"],
                rows: [row("the whole mixed set together (the pipeline)", "within_all"), row("MoTrPAC reference means and scales", "reference"), row("each organ alone (a single-tissue upload)", "within_organ_alone")] }),
    el("p", {}, `So the page never scales a small or single-tissue upload within itself: with fewer than ${WITHIN_MIN_SAMPLES} samples or fewer than ${WITHIN_MIN_TISSUES} tissues it switches to reference scaling and says so. Reference scaling held up in this test, better than within-set scaling on coverage, but it does not fail safe on organs the model never saw: thymus got an empty set in ${pct(K["adult_21wk.reference.Thymus.frac_empty"])} of samples under reference scaling and ${pct(K["adult_21wk.within_all.Thymus.frac_empty"])} under within-set scaling. Across all four BodyMap ages the pattern is the same (reference ${pct(K["all_ages.reference.accuracy"])} correct, within-set ${pct(K["all_ages.within_all.accuracy"])}, each organ alone ${pct(K["all_ages.within_organ_alone.accuracy"])}).`),
    el("p", {}, `Within-set scaling also depends on what else you upload. ${compositionNote(D.PR)}`),
    el("p", { class: "small" }, "Source: results_product/40_product/scaling.csv and composition.csv (scripts/40_product_validation.py; the whole-set scaling reproduces the pipeline's BodyMap probabilities to machine precision). Provenance: pv_scaling_* in site/data/provenance.json."),
  );
}

function flagNote() {
  const F = D.PR.flag_key, V = D.PR.venacv;
  const settings = [["motrpac_heldout", "MoTrPAC, held-out animals"], ["bodymap_adult_21wk", "rat BodyMap, 21-week adults"], ["bodymap_all_ages", "rat BodyMap, all ages"]];
  const rows = [];
  for (const [s, lab] of settings) for (const a of [0.05, 0.1, 0.2]) rows.push({ data: lab, "α": a.toFixed(2), samples: F[`${s}.${a}.n_samples`], "false Mismatch (correct label)": pct(F[`${s}.${a}.false_mismatch_rate`], 1), "Can't confirm (correct label)": pct(F[`${s}.${a}.cant_confirm_rate`], 1), "swapped labels flagged Mismatch": pct(F[`${s}.${a}.swap_vial_detection_rate`], 1), "swaps with ≥ 1 of the pair flagged": pct(F[`${s}.${a}.swap_pair_detection_rate`], 1) });
  const maxMiss = Math.max(...D.PR.flag.map((f) => f.swap_vial_missed_rate));
  const nV = V.length, nFl = V.filter((v) => v.consortium_flagged).length, nMis = V.filter((v) => v.status_claimed_venacv === "mismatch").length, nCc = V.filter((v) => v.status_claimed_venacv === "cant_confirm").length;
  document.getElementById("flag-body").replaceChildren(
    el("p", {}, `A correctly labelled sample should almost never be flagged Mismatch, and a swapped label should be. We measured both on existing held-out scores: every correctly labelled sample, and ${F["motrpac_heldout.0.1.n_swaps"]} simulated swaps between two samples of different tissues (fixed seed), at three error rates.`),
    (() => { const t = tableFrom({ columns: ["data", "α", "samples", "false Mismatch (correct label)", "Can't confirm (correct label)", "swapped labels flagged Mismatch", "swaps with ≥ 1 of the pair flagged"], rows }); t.classList.add("wrap-head"); return t; })(),
    el("p", {}, `${maxMiss === 0 ? "A swapped label is never called Consistent in these data" : `At most ${pct(maxMiss, 2)} of swapped labels are called Consistent in any row`}; the swaps that are not flagged Mismatch mostly get an empty set (Can't confirm), so they are still not confirmed. In another laboratory the price is more Can't confirm on correct labels; recalibrating on a few of your own labelled samples raises coverage (see the Transfer page).`),
    el("p", {}, `A real case: of the ${D.H.extras.venacv_vials} held-out MoTrPAC vena cava vials, the model calls ${nV} brown fat, and the consortium had flagged ${nFl} of those ${nV} as contaminated with brown fat. Checked with the claimed label “vena cava” at α = 0.10, ${nMis} are flagged Mismatch and ${nCc} Can't confirm${nMis + nCc === nV ? ": none is passed as consistent" : ""}.`),
    (() => { const VA = D.PR.venacv_all, f = VA.flagged, u = VA.unflagged;
      return el("p", {}, `The other direction: of all ${f.n_vials + u.n_vials} held-out vena cava vials, the consortium flagged ${f.n_vials}; the model calls ${f.n_called_bat} of those ${f.n_vials} brown fat and ${u.n_called_bat} of the ${u.n_vials} unflagged. Claimed as vena cava, ${f.n_consistent} of the ${f.n_vials} flagged vials passes as Consistent (${f.n_mismatch} Mismatch, ${f.n_cant_confirm} Can't confirm), and ${u.n_mismatch} of the ${u.n_vials} unflagged is flagged Mismatch (${u.n_consistent} Consistent, ${u.n_cant_confirm} Can't confirm).`); })(),
    el("p", { class: "small" }, "Source: results_product/40_product/flag_rates.csv, venacv_cases.csv and venacv_summary.csv (scripts/40_product_validation.py). MoTrPAC rows use the pipeline's held-out 20-gene models (panel re-selected in each animal-grouped fold, each fold's own calibration scores; results_frozen/31_site_regen/06_conformal/TRNSCRPT/scores_*.csv): the same method as the model on this page, not the identical fit. BodyMap rows use this page's model exactly. Provenance: pv_flag_*, pv_venacv_cases."),
  );
}

// ---- downloads ---------------------------------------------------------------------------------------------------------
function download(name, text, type) {
  const a = el("a", { href: URL.createObjectURL(new Blob([text], { type })), download: name });
  document.body.appendChild(a); a.click(); a.remove();
}
const esc = (s) => String(s ?? "").replace(/[&<>"']/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c]));

function reportHtml() {
  const C = state.check, c = C.counts, m = state.model;
  const rows = C.results.map((r) => `<tr class="${r.flagged ? "flag" : ""}"><td>${esc(r.id)}</td><td>${esc(r.label ?? "")}</td><td>${esc(nm(r.call))}</td><td>${esc(r.set.map(nm).join("; ") || "none")}</td><td>${esc(r.status)}</td><td>${esc(r.claimStatus ?? "")}</td><td class="n">${r.callProb.toFixed(2)}</td><td>${esc(nm(r.runnerUp))} (${r.runnerUpProb.toFixed(2)})</td></tr>`).join("\n");
  const flagged = C.results.filter((r) => r.flagged).map((r) => `<li><b>${esc(r.id)}</b> — ${esc(r.status)}${r.claimStatus ? ", " + esc(r.claimStatus) : ""}. ${esc(explainSentence(r, m))}</li>`).join("\n");
  const when = new Date().toISOString().slice(0, 16).replace("T", " ");
  return `<!doctype html><html lang="en"><head><meta charset="utf-8"><title>Tissue check report</title>
<style>body{font:14px/1.45 system-ui,sans-serif;color:#111;margin:2rem;max-width:60rem}h1{font-size:1.5rem}table{border-collapse:collapse;width:100%;font-size:12px}th,td{border-bottom:1px solid #ccc;padding:.25rem .4rem;text-align:left}td.n{text-align:right}tr.flag td{background:#fff3e0}.meta{color:#555;font-size:12px}@media print{body{margin:1cm}tr{break-inside:avoid}}</style></head><body>
<h1>Is this sample the tissue you think it is? — report</h1>
<p class="meta">Generated ${esc(when)} UTC from ${esc(state.source)} by https://corygardner.github.io/motrpac_challenge/ (Tissue Fingerprints, The Rat PAC). Model: the 20-gene transfer model (site/data/panel_model.json, generated ${esc(m._meta.generated)}, git ${esc(m._meta.git_hash)}). α = ${C.alpha.toFixed(2)}; ${esc(C.scaling.mode === "within" ? "within-set" : "MoTrPAC reference")} scaling (${esc(C.scaling.reason)})${state.recal && !state.recal.error ? `; threshold recalibrated on ${state.recal.n} labelled samples` : `; MoTrPAC calibration (${m.calibration.n_animals} animals)`}.</p>
<p><b>${c.total}</b> samples; <b>${c.flagged}</b> flagged. Confident ${c.Confident}, ambiguous ${c.Ambiguous}, unknown ${c.Unknown}.${c.claimed ? ` Of ${c.claimed} claimed labels: consistent ${c.Consistent}, mismatch ${c.Mismatch}, can't confirm ${c["Can't confirm"]}, not in reference ${c["Not in reference"]}.` : ""}</p>
${flagged ? `<h2>Flagged samples</h2><ul>${flagged}</ul>` : ""}
<h2>All samples</h2><table><thead><tr><th>sample</th><th>claimed</th><th>call</th><th>${esc(setLabel(C.alpha))}</th><th>status</th><th>claim</th><th>top p</th><th>runner-up</th></tr></thead><tbody>
${rows}
</tbody></table>
<p class="meta">Status: Confident = one tissue in the set; Ambiguous = several; Unknown = none. Claim: Consistent = the claimed tissue is in the set; Mismatch = the set is not empty and does not hold it; Can't confirm = empty set. The 90 % guarantee holds for samples like the MoTrPAC calibration animals; in another laboratory recalibrate on your own labelled samples. Rat model; ${m.classes.length} reference tissues.</p>
</body></html>`;
}

main().catch((e) => { console.error(e); setProgress("The page could not load its data: " + e.message); });
