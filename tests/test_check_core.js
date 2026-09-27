// The "Check samples" core (site/assets/check-core.js): status rules, count → log2 CPM parity with the pipeline on real
// BodyMap counts, reference-map projection parity with Python, contributions + intercept = logit, synonym mapping, the
// scaling and species guards, CSV safety, and the example (a known label swap is flagged).
// Run: node tests/test_check_core.js
import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import { gunzipSync } from "node:zlib";
import { fileURLToPath } from "node:url";
import { dirname, join } from "node:path";

const here = dirname(fileURLToPath(import.meta.url));
const K = await import(join(here, "..", "site", "assets", "check-core.js"));
const S = await import(join(here, "..", "site", "assets", "score.js"));
const read = (n) => JSON.parse(readFileSync(join(here, "..", "site", "data", n), "utf8"));
const model = read("panel_model.json"), product = read("product.json"), exprM = read("expr_motrpac.json"), exprB = read("expr_bodymap.json");
const organMap = read("bodymap.json").organ_map;
let n = 0;
const check = (c, m) => { n += 1; assert.ok(c, m); };
const close = (a, b, tol, m) => check(Math.abs(a - b) <= tol, `${m}: ${a} vs ${b}`);

// ---- status rules -------------------------------------------------------------------------------------------------
check(K.sampleStatus(1) === "Confident" && K.sampleStatus(2) === "Ambiguous" && K.sampleStatus(5) === "Ambiguous" && K.sampleStatus(0) === "Unknown", "sample status");
const liver = K.mapClaim("liver");
check(K.claimStatus(["LIVER"], liver) === "Consistent", "claimed in set → Consistent");
check(K.claimStatus(["LIVER", "KIDNEY"], liver) === "Consistent", "claimed in an ambiguous set → Consistent");
check(K.claimStatus(["KIDNEY"], liver) === "Mismatch", "set non-empty without the claim → Mismatch");
check(K.claimStatus([], liver) === "Can't confirm", "empty set → Can't confirm");
check(K.claimStatus(["LIVER"], K.mapClaim("")) === null, "no claim → no claim status");
check(K.claimStatus(["LIVER"], K.mapClaim("thymus")) === "Not in reference", "claim outside the 19 → Not in reference");
check(K.claimStatus(["HIPPOC"], K.mapClaim("brain")) === "Consistent" && K.claimStatus(["LIVER"], K.mapClaim("Brain")) === "Mismatch", "super-class claims");
check(K.isFlagged({ status: "Confident", claimStatus: "Mismatch" }) && K.isFlagged({ status: "Unknown", claimStatus: null }) && !K.isFlagged({ status: "Ambiguous", claimStatus: "Consistent" }), "flagged rule");

// ---- synonym mapping ------------------------------------------------------------------------------------------------
const expect = { "Liver": ["LIVER"], "adrenal gland": ["ADRNL"], "Brown fat": ["BAT"], "BAT": ["BAT"], "white_adipose": ["WAT-SC"], "vena cava": ["VENACV"], "SKM-GN": ["SKM-GN"],
                 "gastrocnemius": ["SKM-GN"], "Testis": ["TESTES"], "Testes": ["TESTES"], "small intestine": ["SMLINT"], "Hippocampus": ["HIPPOC"], "left ventricle": ["HEART"],
                 "Muscle": ["SKM-GN", "SKM-VL"], "brain": ["CORTEX", "HIPPOC", "HYPOTH"], "Kidneys": ["KIDNEY"], "lungs": ["LUNG"], "whole blood": ["BLOOD"] };
for (const [k, v] of Object.entries(expect)) check(K.mapClaim(k).classes.join() === v.join(), `synonym ${k} → ${v}`);
for (const k of ["Thymus", "Uterus", "pancreas", "skin"]) check(K.mapClaim(k).kind === "not_in_reference", `${k} is not in the reference`);
for (const c of model.classes) check(K.mapClaim(c).classes.join() === c, `class code ${c} maps to itself`);
check(K.acceptedList().length === model.classes.length && K.acceptedList().every((a) => model.classes.includes(a.code)), "accepted list covers exactly the 19 classes");

// ---- count → log2 CPM parity with src/tfp/io.py log_cpm on real counts (both orientations) ----------------------------------
const fx = JSON.parse(readFileSync(join(here, "fixtures", "product_parity.json"), "utf8"));
const gxs = gunzipSync(readFileSync(join(here, "fixtures", "bodymap_counts_subset.csv.gz"))).toString("utf8");
const parsed = K.parseUpload(gxs, model);
check(parsed.format === "counts" && parsed.orientation === "genes_x_samples" && parsed.samples.length === fx.samples.length, `count matrix detected (${parsed.format}, ${parsed.orientation})`);
let maxD = 0;
fx.samples.forEach((sid, s) => {
  const got = parsed.samples.find((x) => x.id === sid);
  check(got !== undefined, `sample ${sid}`);
  close(got.libSize, fx.library_size[s], 1e-6, `library size ${sid}`);
  fx.genes.forEach((g, j) => { maxD = Math.max(maxD, Math.abs(got.values[g] - fx.log2cpm[s][j])); });
});
check(maxD < 1e-9, `log2 CPM parity genes × samples, max |Δ| ${maxD}`);
// transpose to samples × genes and parse again
const L = gxs.trim().split("\n").map((l) => l.split(","));
const T = L[0].map((_, c) => L.map((r) => r[c]).join(",")).join("\n");
const pT = K.parseUpload(T, model);
check(pT.orientation === "samples_x_genes" && pT.format === "counts", "transposed matrix detected");
let maxT = 0;
fx.samples.forEach((sid, s) => { const got = pT.samples.find((x) => x.id.replace(/"/g, "") === sid); fx.genes.forEach((g, j) => { maxT = Math.max(maxT, Math.abs(got.values[g] - fx.log2cpm[s][j])); }); });
check(maxT < 1e-9, `log2 CPM parity samples × genes, max |Δ| ${maxT}`);
check(parsed.species === "rat", `species of the BodyMap counts: ${parsed.species}`);

// ---- species guard ------------------------------------------------------------------------------------------------------
check(K.guessSpecies(["ENSG00000141510", "ENSG00000012048", "ENSG00000139618", "ENSG00000157764", "ENSG00000146648", "ENSG00000105974"]) === "human", "ENSG ids → human");
check(K.guessSpecies(["TP53", "BRCA1", "ALB", "UMOD", "MYH7", "SFTPC"]) === "human", "upper-case symbols → human");
check(K.guessSpecies(["Tp53", "Brca1", "Alb", "Umod", "Myh7", "Sftpc"]) === "rat or mouse", "capitalised symbols → rat or mouse");
check(K.guessSpecies(model.genes.map((g) => g.id)) === "rat", "ENSRNOG ids → rat");

// ---- scaling guard --------------------------------------------------------------------------------------------------------
check(K.chooseScaling(8, 3, 0).mode === "within" && K.chooseScaling(80, 0, 9).mode === "within", "≥ 8 samples and ≥ 3 tissues → within");
check(K.chooseScaling(7, 5, 5).mode === "reference" && K.chooseScaling(20, 2, 2).mode === "reference" && K.chooseScaling(1, 1, 1).mode === "reference", "otherwise reference");

// ---- contributions + intercept = logit, and the softmax of the logits is the probability ---------------------------------------
const rows = S.exampleRows(model, exprB, organMap, { perOrgan: 8, nLabelled: 0 });
const ex = K.exampleTable(model, rows);
check(ex.swapped.length === 2, "the example swaps two named samples");
const pe = K.parseUpload(ex.csv, model);
check(pe.format === "log2cpm" && pe.orientation === "samples_x_genes" && pe.missing.length === 0 && pe.samples.length === rows.length, "the example parses as a 20-gene table");
const res = K.runCheck(model, pe.samples, { alpha: 0.1 });
check(res.scaling.mode === "within", "the example is scaled within the set");
for (const r of res.results) {
  const lg = K.logits(r.z, model);
  const sm = S.softmax(lg);
  for (let c = 0; c < model.classes.length; c++) {
    const ct = K.contributions(r.z, model, c);
    close(ct.logit, lg[c], 1e-9, `contributions + intercept = logit (${r.id}, ${model.classes[c]})`);
    close(sm[c], r.p[c], 1e-12, `softmax(logit) = p (${r.id}, ${model.classes[c]})`);
  }
  const d = K.contrast(r.z, model, r.call, r.runnerUp);
  close(d.terms.reduce((a, t) => a + t.diff, d.interceptDiff), d.logitDiff, 1e-9, `contrast sums (${r.id})`);
}
// the example reproduces score.js and the pipeline's own sets for the same rows
const sc = S.scoreSamples(model, pe.samples.map((s) => ({ id: s.id, values: s.values })), { alpha: 0.1 });
res.results.forEach((r, i) => { check(r.call === sc.results[i].call && r.set.join() === sc.results[i].set.join(), `same call and set as score.js (${r.id})`); });
// the swapped samples are flagged; the unseen organs never get Consistent
for (const id of ex.swapped) { const r = res.results.find((x) => x.id === id); check(r.claimStatus === "Mismatch", `swapped sample ${id} is flagged Mismatch (${r.claimStatus})`); }
for (const r of res.results.filter((x) => /^(Thymus|Uterus)/.test(x.id))) check(r.claimStatus === "Not in reference" && r.flagged, `${r.id}: not in reference and flagged`);
check(res.results.filter((r) => r.claimStatus === "Mismatch").every((r) => ex.swapped.includes(r.id)), "only the swapped samples are Mismatch in the example");
// α drives the sets: a larger α never gives larger sets
const r20 = K.runCheck(model, pe.samples, { alpha: 0.2 }), r05 = K.runCheck(model, pe.samples, { alpha: 0.05 });
res.results.forEach((r, i) => check(r20.results[i].setSize <= r.setSize && r.setSize <= r05.results[i].setSize, `set size monotone in α (${r.id})`));

// ---- reference-map projection parity with Python (scripts/40_product_validation.py) -----------------------------------------------
const gi = Object.fromEntries(exprM.genes.map((g, i) => [g, i]));
check(product.pca.genes.join() === model.genes.map((g) => g.id).join(), "PCA genes are the model genes in order");
let maxP = 0;
exprM.samples.forEach((s, j) => {
  const X = [model.genes.map((g) => exprM.values[gi[g.id]][j])];
  const z = S.zscoreSource(X, model).Z[0];
  const [x, y] = K.project(z, product.pca);
  maxP = Math.max(maxP, Math.abs(x - product.pca.reference.x[j]), Math.abs(y - product.pca.reference.y[j]));
  if (j < 5) check(product.pca.reference.id[j] === s.id, "reference ids in order");
});
check(maxP < 1e-5, `projection parity over ${exprM.samples.length} reference vials, max |Δ| ${maxP}`);

// ---- reference statistics and CSV safety ------------------------------------------------------------------------------------------------
const st = K.referenceStats(exprM, model);
check(Object.keys(st).length === model.classes.length && st.LIVER[model.genes[0].id].q1 <= st.LIVER[model.genes[0].id].median, "reference median/IQR per tissue");
check(K.csvCell("=HYPERLINK(1)") === "'=HYPERLINK(1)" && K.csvCell("+1") === "'+1" && K.csvCell("@x") === "'@x" && K.csvCell(-1) === "-1", "formula injection neutralised, numbers untouched");
check(K.csvCell('a,"b"') === '"a,""b"""', "quoting");
const csv = K.resultsToCsv(res, model);
check(csv.split("\n")[0].startsWith("sample,claimed,") && csv.trim().split("\n").length === res.results.length + 1, "results CSV shape");
check(typeof K.explainSentence(res.results[0], model) === "string" && K.explainSentence(res.results[0], model).startsWith("Called "), "explanation sentence");

// ---- missing values: blank, whitespace, NA and a short row are all missing, never zero -------------------------------
check(K.parseCell("") === null && K.parseCell("  ") === null && K.parseCell("NA") === null && K.parseCell("NaN") === null && K.parseCell("null") === null
      && K.parseCell("-") === null && K.parseCell("abc") === null && K.parseCell(undefined) === null && K.parseCell("0") === 0 && K.parseCell(" 2.5 ") === 2.5, "parseCell");
{
  const lines = ex.csv.trim().split("\n");
  const head = lines[0].split(",");
  let same = 0, cases = 0, named = 0;
  const base = lines.slice(1).map((l) => l.split(","));
  for (let i = 0; i < base.length; i++) {
    for (let j = 1; j <= model.genes.length; j++) {
      const mk = (val) => [lines[0], ...base.map((c, k) => (k === i ? c.map((x, jj) => (jj === j ? val : x)) : c).join(","))].join("\n");
      const rb = K.runCheck(model, K.parseUpload(mk(""), model).samples, { alpha: 0.1 });
      const rn = K.runCheck(model, K.parseUpload(mk("NA"), model).samples, { alpha: 0.1 });
      cases += 1;
      if (rb.results.every((r, k) => r.call === rn.results[k].call && r.set.join() === rn.results[k].set.join() && r.p.every((v, c) => Math.abs(v - rn.results[k].p[c]) < 1e-12))) same += 1;
      if (rb.results[i].missingGenes.includes(head[j])) named += 1;
    }
  }
  check(cases === base.length * model.genes.length && same === cases, `blank and NA give identical probabilities, calls and sets in ${same} of ${cases} single-cell cases`);
  check(named === cases, `the blanked gene is named in the sample's missing-gene list in ${named} of ${cases} cases`);
  // a truncated row: its last cells are missing (null), with a warning naming the row
  const trunc = [lines[0], ...base.map((c, k) => (k === 3 ? c.slice(0, 10) : c).join(","))].join("\n");
  const pt = K.parseUpload(trunc, model);
  const s3 = pt.samples[3];
  check(model.genes.slice(9).every((g) => s3.values[g.id] === null) && model.genes.slice(0, 9).every((g) => s3.values[g.id] !== null), "a truncated row gives nulls for its missing cells");
  check(pt.warnings.some((w) => w.includes(base[3][0]) && w.includes("shorter than the header")), "a truncated row is named in a warning");
  // a linear-scale table (2^x − 1 of the example) warns; the example itself does not
  const lin = [lines[0], ...base.map((c) => c.map((x, jj) => (jj >= 1 && jj <= model.genes.length && x !== "" ? String(2 ** Number(x) - 1) : x)).join(","))].join("\n");
  check(K.parseUpload(lin, model).warnings.some((w) => w.includes("linear or TPM")), "a linear-scale table warns");
  check(!K.parseUpload(ex.csv, model).warnings.some((w) => w.includes("linear or TPM")), "the log2 CPM example does not warn");
  // a count matrix with an empty cell leaves it out of the library size and says so
  const gl = gxs.trim().split("\n");
  const holed = [gl[0], gl[1].split(",").map((x, jj) => (jj === 1 ? "" : x)).join(","), ...gl.slice(2)].join("\n");
  const ph = K.parseUpload(holed, model);
  check(ph.nMissing === 1 && ph.warnings.some((w) => w.includes("left out of the library size")), "an empty count cell is left out of the library size, with a warning");
}

// ---- input checks and the display order ----------------------------------------------------------------------------
check(K.parseUpload("a,b\n1,2", model).error && K.parseUpload("a,b\n1,2", model).error.includes("none of the"), "a table without any panel gene is refused");
{
  const lines = ex.csv.trim().split("\n");
  const few = lines.map((l) => l.split(",").filter((_, j) => j === 0 || j > 12).join(",")).join("\n");     // keeps 8 of the 20 genes + the claim
  check(K.parseUpload(few, model).warnings.some((w) => w.includes("of the 20 panel genes were found")), "fewer than half the panel genes warns");
  const holed = [lines[0], ...lines.slice(1).map((l, k) => (k === 0 ? l.split(",").map((x, j) => (j === 3 ? "" : x)).join(",") : l))].join("\n");
  const ph = K.parseUpload(holed, model);
  check(ph.partialMissing.length === 1 && ph.warnings.some((w) => w.includes("empty or non-numeric panel-gene value") && w.includes(ph.partialMissing[0])), "a blank panel-gene cell is named in a warning");
  check(K.parseUpload(ex.csv, model).partialMissing.length === 0, "the example has no partial missing values");
}
{
  const ord = K.orderByPriority(res.results);
  const pr = ord.map(K.priority);
  check(pr.every((v, i) => i === 0 || pr[i - 1] <= v), "priority order is non-decreasing");
  check(ex.swapped.every((id) => ord.slice(0, ex.swapped.length).some((r) => r.id === id)), "the swapped samples come first in the example");
  const same = ord.filter((r) => K.priority(r) === 5).map((r) => res.results.indexOf(r));
  check(same.every((v, i) => i === 0 || same[i - 1] < v), "upload order kept within a priority group");
}

// ---- the explanation names only genes that matter; the example never recalibrates on its swapped labels ----------------
{
  const kid = res.results.find((r) => r.id === ex.swapped[1]);            // really kidney, claimed liver
  check(kid.call === "KIDNEY" && kid.claimStatus === "Mismatch", "the swapped kidney is called kidney and flagged");
  const other = kid.claim.classes[0];
  const d = K.drivers(K.contrast(kid.z, model, kid.call, other).terms).map((t) => t.symbol);
  check(d.join() === "Umod,Cfhr1", `the swapped kidney's drivers are exactly Umod and Cfhr1 (got ${d.join(", ")})`);
  const sent = K.explainSentence(kid, model);
  check(sent.includes("Umod (+") && sent.includes("Cfhr1 (") && !sent.includes("Trim29") && sent.includes(" sd)"), `the sentence names Umod and Cfhr1 in sd units only: ${sent}`);
  for (const r of res.results) {
    const ds = K.drivers(K.contrast(r.z, model, r.call, r.runnerUp).terms);
    if (ds.length) check(ds.every((t) => t.diff >= Math.max(K.DRIVER_MIN, K.DRIVER_REL * ds[0].diff) - 1e-12), `driver floor (${r.id})`);
  }
  check(K.explainSentence({ ...kid, z: kid.z.map(() => 0) }, model).includes("no single gene dominates"), "no qualifying driver → no single gene dominates");
  check(K.setLabel(0.1) === "90 % prediction set" && K.setLabel(0.2) === "80 % prediction set" && K.setLabel(0.05) === "95 % prediction set", "set label follows α");
  for (const a of [0.05, 0.1, 0.2]) {
    const base = K.runCheck(model, pe.samples, { alpha: a });
    const lab = K.labelledForRecal(base.results, ex.swapped);
    check(ex.swapped.every((id) => !(id in lab)), `α = ${a}: the swapped samples are left out of the recalibration`);
    const rc = S.recalibrate(model, base.results, lab, a);
    const after = K.runCheck(model, pe.samples, { alpha: a, q: rc.q });
    for (const id of ex.swapped) check(after.results.find((r) => r.id === id).claimStatus === "Mismatch", `α = ${a}: ${id} stays Mismatch after recalibration (q = ${rc.q})`);
  }
}

// ---- the scaling report describes the mode actually used -------------------------------------------------------------
{
  const fr = K.runCheck(model, pe.samples, { alpha: 0.1, scaling: "reference" }).scaling;
  check(fr.mode === "reference" && fr.auto === "within" && fr.forced && fr.reason.includes("MoTrPAC reference") && !fr.reason.includes("within your") && fr.caveat === null, "forced reference: its own reason, the automatic choice reported");
  const fw = K.runCheck(model, pe.samples.slice(0, 5), { alpha: 0.1, scaling: "within" }).scaling;
  check(fw.mode === "within" && fw.auto === "reference" && fw.forced && fw.caveat && fw.reason.includes("within your"), "forced within on 5 samples: a caveat");
  const au = K.runCheck(model, pe.samples, { alpha: 0.1 }).scaling;
  check(!au.forced && au.caveat === null && au.mode === "within" && au.reason.includes("within your"), "automatic within: no caveat");
}

// ---- the landing page's example input tables equal their source data ---------------------------------------------------
{
  const T = K.exampleInputTables(model, rows, product.example_counts);
  check(T.a.columns.length === 6 && T.a.rows.length === 2 && T.a.rows[0].claimed_tissue === "adrenal" && T.a.rows[1].claimed_tissue === "liver", "table A: an adrenal and a liver sample");
  for (const r of T.a.rows) {
    const src = rows.find((x) => x.id === r.sample);
    model.genes.slice(0, 3).forEach((g) => check(r[g.symbol] === Number(src.values[g.id]).toFixed(2), `table A value ${r.sample} ${g.symbol}`));
    const j = exprB.samples.findIndex((sm) => sm.id === r.sample), gi2 = exprB.genes.indexOf(model.genes[0].id);
    check(Math.abs(Number(r[model.genes[0].symbol]) - exprB.values[gi2][j]) <= 0.005 + 1e-9, `table A equals expr_bodymap.json (${r.sample})`);
  }
  const fxL = gunzipSync(readFileSync(join(here, "fixtures", "bodymap_counts_subset.csv.gz"))).toString("utf8").trim().split("\n").map((l) => l.split(","));
  const col = Object.fromEntries(fxL[0].map((h, i) => [h, i]));
  for (const r of T.b.rows.slice(0, -1)) {
    const line = fxL.find((l) => l[0] === r.gene_id);
    check(/^ENSRNOG\d{11}$/.test(r.gene_id), `table B uses 11-digit rat Ensembl IDs (${r.gene_id})`);
    for (const smp of product.example_counts.samples) check(r[smp] === line[col[smp]], `table B count ${r.gene_id} ${smp} equals the fixture`);
  }
}

console.log(`ok: ${n} assertions (count → CPM max |Δ| ${maxD.toExponential(1)}, projection max |Δ| ${maxP.toExponential(1)}, ${res.results.length} example samples)`);
