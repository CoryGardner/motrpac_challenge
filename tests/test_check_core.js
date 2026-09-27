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

console.log(`ok: ${n} assertions (count → CPM max |Δ| ${maxD.toExponential(1)}, projection max |Δ| ${maxP.toExponential(1)}, ${res.results.length} example samples)`);
