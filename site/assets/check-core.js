// "Check samples": pure functions (no DOM) behind site/index.html, tested by tests/test_check_core.js.
// Upload parsing (20-gene log2 CPM table or a full count matrix in either orientation), the pipeline's count → log2 CPM
// transform, claimed-label mapping, the species and scaling guards, the status rules, per-gene contributions, the
// reference-map projection and CSV output. The arithmetic of the call itself is score.js + conformal.js, unchanged.
import { zscoreWithin, zscoreSource, probabilities } from "./score.js";
import { setsFor } from "./conformal.js";

// ---- claimed labels ------------------------------------------------------------------------------------------------
export const CLASS_NAMES = {
  ADRNL: "adrenal gland", BAT: "brown adipose", BLOOD: "whole blood", COLON: "colon", CORTEX: "cerebral cortex", HEART: "heart",
  HIPPOC: "hippocampus", HYPOTH: "hypothalamus", KIDNEY: "kidney", LIVER: "liver", LUNG: "lung", OVARY: "ovary",
  "SKM-GN": "gastrocnemius", "SKM-VL": "vastus lateralis", SMLINT: "small intestine", SPLEEN: "spleen", TESTES: "testes",
  VENACV: "vena cava", "WAT-SC": "white adipose (subcutaneous)",
};
// accepted names per class (lower case, compared after normalisation); the class code itself is always accepted
export const SYNONYMS = {
  ADRNL: ["adrenal", "adrenal gland", "adrenals", "adrenal glands"],
  BAT: ["brown adipose", "brown adipose tissue", "brown fat", "bat", "interscapular brown adipose", "ibat"],
  BLOOD: ["blood", "whole blood", "peripheral blood"],
  COLON: ["colon", "large intestine"],
  CORTEX: ["cortex", "cerebral cortex", "frontal cortex", "brain cortex", "neocortex"],
  HEART: ["heart", "cardiac", "left ventricle", "ventricle", "myocardium"],
  HIPPOC: ["hippocampus", "hippocampal"],
  HYPOTH: ["hypothalamus"],
  KIDNEY: ["kidney", "kidneys", "renal"],
  LIVER: ["liver", "hepatic"],
  LUNG: ["lung", "lungs", "pulmonary"],
  OVARY: ["ovary", "ovaries"],
  "SKM-GN": ["gastrocnemius", "gastrocnemius muscle", "gastroc", "skm gn"],
  "SKM-VL": ["vastus lateralis", "vastus", "skm vl"],
  SMLINT: ["small intestine", "jejunum", "ileum", "duodenum", "small bowel"],
  SPLEEN: ["spleen", "splenic"],
  TESTES: ["testes", "testis", "testicle"],
  VENACV: ["vena cava", "venacava", "inferior vena cava", "ivc"],
  "WAT-SC": ["white adipose", "white adipose tissue", "subcutaneous white adipose", "subcutaneous adipose", "white fat", "wat", "inguinal fat", "inguinal adipose"],
};
// names that cover several classes (the model's super-classes in the BodyMap transfer)
export const SUPERCLASSES = {
  brain: ["CORTEX", "HIPPOC", "HYPOTH"], "whole brain": ["CORTEX", "HIPPOC", "HYPOTH"],
  muscle: ["SKM-GN", "SKM-VL"], "skeletal muscle": ["SKM-GN", "SKM-VL"], "muscle skeletal": ["SKM-GN", "SKM-VL"], quadriceps: ["SKM-GN", "SKM-VL"],
  adipose: ["BAT", "WAT-SC"], fat: ["BAT", "WAT-SC"], "adipose tissue": ["BAT", "WAT-SC"],
  intestine: ["COLON", "SMLINT"], gut: ["COLON", "SMLINT"],
};
const normLabel = (s) => String(s ?? "").trim().toLowerCase().replace(/[_\-.]+/g, " ").replace(/\s+/g, " ").replace(/^"|"$/g, "");
const LOOKUP = (() => {
  const m = new Map();
  for (const [c, names] of Object.entries(SYNONYMS)) { m.set(normLabel(c), [c]); for (const n of names) m.set(normLabel(n), [c]); }
  for (const [n, cs] of Object.entries(SUPERCLASSES)) m.set(normLabel(n), cs);
  return m;
})();

/** Map a claimed label to model classes: { raw, classes, kind: "class" | "superclass" | "not_in_reference" | "none" }. */
export function mapClaim(raw) {
  const n = normLabel(raw);
  if (!n || n === "na" || n === "nan" || n === "null") return { raw: raw ?? null, classes: [], kind: "none" };
  const hit = LOOKUP.get(n) || LOOKUP.get(n.replace(/ tissue$/, "")) || LOOKUP.get(n.replace(/s$/, ""));
  if (!hit) return { raw, classes: [], kind: "not_in_reference" };
  return { raw, classes: hit.slice(), kind: hit.length === 1 ? "class" : "superclass" };
}

export function acceptedList() {
  return Object.keys(SYNONYMS).map((c) => ({ code: c, name: CLASS_NAMES[c], also: SYNONYMS[c] }));
}

// ---- parsing ----------------------------------------------------------------------------------------------------------
export const LABEL_KEY = /^(true[ _]?tissue|claimed[ _]?tissue|tissue|claimed|label)$/i;
const ENS = /^ENS([A-Z]{0,4})G\d{6,}/i;
const splitterFor = (line) => (line.includes("\t") ? "\t" : line.includes(";") && !line.includes(",") ? ";" : ",");
const clean = (c) => c.trim().replace(/^"|"$/g, "");
const MISSING = /^(|na|n\/a|nan|null|none|-|\.)$/i;
/** One numeric cell: "", whitespace, NA, NaN, null, "-" and anything non-numeric (or a cell past the end of a short row)
 *  → null; never 0. (Number("") is 0, which would read an empty cell as "not expressed".) */
export function parseCell(raw) {
  if (raw === undefined || raw === null) return null;
  const c = clean(String(raw));
  if (MISSING.test(c)) return null;
  const v = Number(c);
  return Number.isFinite(v) ? v : null;
}
const geneKey = (s) => String(s).trim().toLowerCase().replace(/\.\d+$/, "");

/** Species from gene identifiers: Ensembl prefixes first, then symbol case (rat/mouse: Capitalised; human: UPPER). */
export function guessSpecies(ids) {
  let ensHuman = 0, ensRat = 0, ensMouse = 0, ens = 0, upper = 0, cap = 0, sym = 0;
  for (const raw of ids) {
    const s = String(raw).trim();
    const m = s.match(ENS);
    if (m) { ens += 1; const p = m[1].toUpperCase(); if (p === "") ensHuman += 1; else if (p === "RNO") ensRat += 1; else if (p === "MUS") ensMouse += 1; continue; }
    if (!/^[A-Za-z][A-Za-z0-9-]{1,}$/.test(s)) continue;
    sym += 1;
    if (s === s.toUpperCase() && /[A-Z]/.test(s)) upper += 1; else if (/^[A-Z][a-z0-9]/.test(s)) cap += 1;
  }
  if (ens >= Math.max(5, 0.5 * (ens + sym))) {
    if (ensHuman > 0.5 * ens) return "human";
    if (ensMouse > 0.5 * ens) return "mouse";
    if (ensRat > 0.5 * ens) return "rat";
    return "unknown";
  }
  if (sym >= 5 && upper > 0.8 * sym) return "human";
  if (sym >= 5 && cap > 0.5 * sym) return "rat or mouse";
  return "unknown";
}

/**
 * Parse an upload into the panel genes' log2 CPM per sample, streaming over lines (a full count matrix is never held
 * as a matrix: only library sizes and the 20 panel rows/columns are kept).
 *   format "auto": a table with ≥ minFeaturesForCounts gene features is a count matrix, otherwise 20-gene log2 CPM.
 *   orientation "auto": the panel genes are looked for in the header and in the first column; failing that, the long
 *   dimension is genes.
 * Returns { format, orientation, samples: [{ id, label, values: {geneId: log2cpm|null}, libSize }], nFeatures, species,
 *           missing: [gene], warnings: [..] }.
 */
export function parseUpload(text, model, { format = "auto", orientation = "auto", minFeaturesForCounts = 200, onProgress = null } = {}) {
  const lines = String(text).replace(/\r/g, "").split("\n");
  let h = 0;
  while (h < lines.length && (!lines[h].trim() || lines[h].startsWith("#"))) h += 1;
  if (h >= lines.length) return { error: "the table is empty" };
  const sep = splitterFor(lines[h]);
  const header = lines[h].split(sep).map(clean);
  const body = [];
  for (let i = h + 1; i < lines.length; i++) if (lines[i].trim() && !lines[i].startsWith("#")) body.push(i);
  if (!body.length) return { error: "the table needs a header row and at least one data row" };
  const panelKeys = new Map();
  model.genes.forEach((g) => { panelKeys.set(geneKey(g.id), g.id); panelKeys.set(geneKey(g.symbol), g.id); });
  const firstCol = body.slice(0, 5000).map((i) => clean(lines[i].split(sep, 1)[0]));
  const inHeader = header.slice(1).filter((c) => panelKeys.has(geneKey(c))).length;
  const inFirst = firstCol.filter((c) => panelKeys.has(geneKey(c))).length;
  let orient = orientation;
  if (orient === "auto") {
    if (inHeader > inFirst) orient = "samples_x_genes";
    else if (inFirst > inHeader) orient = "genes_x_samples";
    else orient = body.length > header.length - 1 ? "genes_x_samples" : "samples_x_genes";
  }
  const warnings = [];
  const out = { orientation: orient, warnings };
  if (orient === "samples_x_genes") {
    const labelCol = header.findIndex((c, i) => i > 0 && LABEL_KEY.test(c));
    const geneCols = header.map((c, i) => i).filter((i) => i > 0 && i !== labelCol);
    const nFeatures = geneCols.length;
    const fmt = format === "auto" ? (nFeatures >= minFeaturesForCounts ? "counts" : "log2cpm") : format;
    const panelCol = {};
    for (const i of geneCols) { const g = panelKeys.get(geneKey(header[i])); if (g && panelCol[g] === undefined) panelCol[g] = i; }
    const samples = [];
    let nonInt = 0, neg = 0, htseq = 0, nMissing = 0;
    const shortRows = [];
    const htseqCol = new Set(geneCols.filter((i) => header[i].startsWith("__")));
    body.forEach((li, k) => {
      const cells = lines[li].split(sep);
      const id = clean(cells[0] ?? "");
      if (!id) return;
      let lib = 0;
      const raw = {};
      if (cells.length < header.length) shortRows.push(id);
      if (fmt === "counts") {
        for (const i of geneCols) {
          if (htseqCol.has(i)) continue;
          const v = parseCell(cells[i]);
          if (v === null) { nMissing += 1; continue; }
          lib += v; if (v < 0) neg += 1; else if (v !== Math.floor(v)) nonInt += 1;
        }
      }
      for (const g of model.genes) {
        const i = panelCol[g.id];
        raw[g.id] = i === undefined ? null : parseCell(cells[i]);
      }
      const label = labelCol > 0 ? clean(cells[labelCol] ?? "") || null : null;
      samples.push({ id, label, raw, libSize: fmt === "counts" ? lib : null });
      if (onProgress && k % 200 === 0) onProgress(k / body.length);
    });
    htseq = htseqCol.size;
    Object.assign(out, { format: fmt, nFeatures, samples, geneIds: geneCols.map((i) => header[i]), nonInt, neg, htseq, nMissing, shortRows });
  } else {
    const sampleCols = header.map((c, i) => i).filter((i) => i > 0);
    const nFeatures = body.length;
    const fmt = format === "auto" ? (nFeatures >= minFeaturesForCounts ? "counts" : "log2cpm") : format;
    const lib = new Array(sampleCols.length).fill(0);
    const raw = sampleCols.map(() => ({}));
    let labels = null, nonInt = 0, neg = 0, htseq = 0, nMissing = 0;
    const shortRows = [];
    const geneIds = [];
    body.forEach((li, k) => {
      const cells = lines[li].split(sep);
      const key = clean(cells[0] ?? "");
      if (LABEL_KEY.test(key)) { labels = sampleCols.map((i) => clean(cells[i] ?? "") || null); return; }
      if (key.startsWith("__")) { htseq += 1; return; }
      geneIds.push(key);
      if (cells.length < header.length) shortRows.push(key);
      const g = panelKeys.get(geneKey(key));
      for (let s = 0; s < sampleCols.length; s++) {
        const v = parseCell(cells[sampleCols[s]]);
        if (v === null) { if (fmt === "counts") nMissing += 1; }
        else if (fmt === "counts") { lib[s] += v; if (v < 0) neg += 1; else if (v !== Math.floor(v)) nonInt += 1; }
        if (g && raw[s][g] === undefined) raw[s][g] = v;
      }
      if (onProgress && k % 500 === 0) onProgress(k / body.length);
    });
    const samples = sampleCols.map((i, s) => ({ id: header[i], label: labels ? labels[s] : null,
      raw: Object.fromEntries(model.genes.map((g) => [g.id, raw[s][g.id] ?? null])), libSize: fmt === "counts" ? lib[s] : null }));
    Object.assign(out, { format: fmt, nFeatures, samples, geneIds, nonInt, neg, htseq, nMissing, shortRows });
  }
  if (onProgress) onProgress(1);
  // values: log2 CPM of the panel genes
  for (const s of out.samples) {
    s.values = out.format === "counts" ? countsToLog2Cpm(s.raw, s.libSize) : { ...s.raw };
  }
  out.missing = model.genes.filter((g) => out.samples.every((s) => s.values[g.id] === null));
  if (!out.samples.length) return { error: "no sample rows or columns were found" };
  if (out.missing.length === model.genes.length) {
    return { error: `none of the ${model.genes.length} panel genes (${model.genes.slice(0, 3).map((g) => g.symbol).join(", ")}, …) was found in the header or the first column, as a gene symbol or an Ensembl ID. Check that the table has one column (or row) per gene and a header row.` };
  }
  if (out.missing.length > model.genes.length / 2) {
    warnings.push(`only ${model.genes.length - out.missing.length} of the ${model.genes.length} panel genes were found; the calls rest on those alone and are unreliable.`);
  }
  // samples with a missing value in a gene the table does have (a blank or NA cell), not a gene absent from the table
  const present = model.genes.filter((g) => !out.missing.includes(g));
  out.partialMissing = out.samples.filter((s) => present.some((g) => s.values[g.id] === null)).map((s) => s.id);
  if (out.partialMissing.length) {
    const names = out.partialMissing.slice(0, 5).join(", ") + (out.partialMissing.length > 5 ? ` and ${out.partialMissing.length - 5} more` : "");
    warnings.push(`${out.partialMissing.length} sample${out.partialMissing.length === 1 ? " has" : "s have"} an empty or non-numeric panel-gene value (${names}); each such gene is scored at the reference mean, which weakens that sample's call. They are marked in the table's "missing" column.`);
  }
  out.species = guessSpecies(out.geneIds);
  if (out.shortRows.length) {
    const names = out.shortRows.slice(0, 5).join(", ") + (out.shortRows.length > 5 ? ` and ${out.shortRows.length - 5} more` : "");
    warnings.push(`${out.shortRows.length} row${out.shortRows.length === 1 ? " is" : "s are"} shorter than the header (${names}); the missing cells are treated as missing values, not zeros.`);
  }
  if (out.format === "log2cpm") {
    const vals = out.samples.flatMap((s) => Object.values(s.values)).filter((v) => v !== null);
    const lo = vals.length ? Math.min(...vals) : 0, hi = vals.length ? Math.max(...vals) : 0;
    if (lo < LOG2CPM_MIN || hi > LOG2CPM_MAX) warnings.push(`values from ${lo.toFixed(2)} to ${hi.toFixed(2)} are impossible for log2(CPM + 1) (expected ${LOG2CPM_MIN} to ${LOG2CPM_MAX}): these look like linear or TPM values; take log2(x + 1) of CPM, or upload raw counts.`);
  }
  if (out.format === "counts") {
    if (out.nMissing) warnings.push(`${out.nMissing} empty or non-numeric cells in the count matrix were left out of the library size.`);
    if (out.nonInt) warnings.push(`${out.nonInt} non-integer values in what looks like a count matrix; the page treats them as counts. If these are TPM or normalised values the calls are not comparable.`);
    if (out.neg) warnings.push(`${out.neg} negative values: a count matrix has none. Check the file.`);
    if (out.htseq) warnings.push(`${out.htseq} HTSeq summary rows or columns (names starting with "__") were left out of the library size.`);
    if (out.nFeatures < 10000) warnings.push(`only ${out.nFeatures} gene features: the library size is the sum over the genes given, so a filtered matrix shifts every CPM.`);
  }
  return out;
}

/** The pipeline's transform (src/tfp/io.py log_cpm): CPM = count / library × 1e6, library = the sum over every gene
 *  given, then log2(CPM + 1); a library of 0 gives null. */
export function countsToLog2Cpm(raw, libSize) {
  const out = {};
  for (const [g, v] of Object.entries(raw)) out[g] = v === null || !Number.isFinite(v) || !(libSize > 0) ? null : Math.log2((v / libSize) * 1e6 + 1);
  return out;
}

// ---- guards and status rules -------------------------------------------------------------------------------------------
/** Plausible range of log2(CPM + 1): ≥ 0 by construction (−0.5 allows rounding), ≤ 20 (2^20 CPM exceeds a whole library). */
export const LOG2CPM_MIN = -0.5, LOG2CPM_MAX = 20;
export const WITHIN_MIN_SAMPLES = 8;
export const WITHIN_MIN_TISSUES = 3;

/** Within-set z-scoring only for ≥ 8 samples spanning ≥ 3 tissues (claimed, or called under reference scaling). */
export function chooseScaling(nSamples, nClaimedTissues, nPredictedTissues) {
  const tissues = Math.max(nClaimedTissues, nPredictedTissues);
  const ok = nSamples >= WITHIN_MIN_SAMPLES && tissues >= WITHIN_MIN_TISSUES;
  return { mode: ok ? "within" : "reference", nSamples, tissues,
           reason: ok ? `${nSamples} samples spanning ${tissues} tissues: genes z-scored within your set, as the pipeline does for a new dataset`
                      : `${nSamples} sample${nSamples === 1 ? "" : "s"} spanning ${tissues} tissue${tissues === 1 ? "" : "s"} (within-set scaling needs ≥ ${WITHIN_MIN_SAMPLES} samples and ≥ ${WITHIN_MIN_TISSUES} tissues): genes z-scored with the MoTrPAC reference means and scales` };
}

/** The scaling actually used, with a reason that describes it, and what the automatic rule would have chosen.
 *  { mode, auto, forced, reason, autoReason, caveat } — caveat is set when within-set scaling is forced on a set too
 *  small or too uniform for it. */
export function scalingReport(auto, mode) {
  const n = auto.nSamples, t = auto.tissues;
  const nS = `${n} sample${n === 1 ? "" : "s"} spanning ${t} tissue${t === 1 ? "" : "s"}`;
  const reason = mode === "within"
    ? `genes z-scored within your ${nS}, as the pipeline does for a new mixed dataset`
    : `genes z-scored with the MoTrPAC reference means and scales (${nS})`;
  const forced = mode !== auto.mode;
  const caveat = forced && mode === "within"
    ? `within-set scaling needs at least ${WITHIN_MIN_SAMPLES} samples spanning at least ${WITHIN_MIN_TISSUES} tissues; on a smaller or more uniform set it erases the differences the model reads (on BodyMap, scaling each organ alone named almost no organ correctly)`
    : null;
  return { ...auto, mode, auto: auto.mode, forced, reason, autoReason: auto.reason, caveat };
}

/** Sample status from the prediction set: Confident (1 class), Ambiguous (≥ 2), Unknown (empty). */
export function sampleStatus(setSize) {
  return setSize === 0 ? "Unknown" : setSize === 1 ? "Confident" : "Ambiguous";
}

/** Claim status: Consistent (a claimed class is in the set), Mismatch (set non-empty, no claimed class in it),
 *  Can't confirm (empty set); no claim → null; a claim outside the 19 tissues → "Not in reference". */
export function claimStatus(set, claim) {
  if (!claim || claim.kind === "none") return null;
  if (claim.kind === "not_in_reference") return "Not in reference";
  if (!set.length) return "Can't confirm";
  return claim.classes.some((c) => set.includes(c)) ? "Consistent" : "Mismatch";
}

/** Display order of the results table: the most actionable first (Mismatch, then Can't confirm, Not in reference, Unknown
 *  without a claim, Ambiguous, Consistent / Confident), keeping the upload order within each group. */
export function priority(r) {
  if (r.claimStatus === "Mismatch") return 0;
  if (r.claimStatus === "Can't confirm") return 1;
  if (r.claimStatus === "Not in reference") return 2;
  if (r.status === "Unknown") return 3;
  if (r.status === "Ambiguous") return 4;
  return 5;
}
export function orderByPriority(results) {
  return results.map((r, i) => [r, i]).sort((a, b) => priority(a[0]) - priority(b[0]) || a[1] - b[1]).map(([r]) => r);
}

export function isFlagged(r) {
  return r.claimStatus === "Mismatch" || r.claimStatus === "Can't confirm" || r.claimStatus === "Not in reference" || r.status === "Unknown";
}

// ---- explanation -------------------------------------------------------------------------------------------------------
export function logits(z, model) {
  return model.coef.map((w, c) => w.reduce((a, wj, j) => a + wj * z[j], 0) + model.intercept[c]);
}

/** Per-gene contributions to one class's logit: coef[c][j] · z[j]; their sum plus the intercept is the logit. */
export function contributions(z, model, cls) {
  const c = typeof cls === "number" ? cls : model.classes.indexOf(cls);
  const terms = model.genes.map((g, j) => ({ id: g.id, symbol: g.symbol, z: z[j], coef: model.coef[c][j], value: model.coef[c][j] * z[j] }));
  const intercept = model.intercept[c];
  return { cls: model.classes[c], terms, intercept, logit: terms.reduce((a, t) => a + t.value, intercept) };
}

/** "Why X, not Y": per gene, contribution to X minus contribution to Y (positive favours X). */
export function contrast(z, model, a, b) {
  const A = contributions(z, model, a), B = contributions(z, model, b);
  return { a: A.cls, b: B.cls, terms: A.terms.map((t, j) => ({ ...t, diff: t.value - B.terms[j].value })), interceptDiff: A.intercept - B.intercept, logitDiff: A.logit - B.logit };
}

/** One plain-language sentence naming the genes that drive the call (and, with a claim, why not the claim). */
/** "90 % prediction set" at α = 0.10: every label of the sets derives from α. */
export function setLabel(alpha) {
  return `${Math.round(100 * (1 - alpha))} % prediction set`;
}

/** A gene is named as a driver only if its contribution to "why X, not Y" is at least max(DRIVER_MIN, DRIVER_REL × the
 *  largest); written as the chart writes it, "Umod (+2.9 sd)". */
export const DRIVER_MIN = 0.1, DRIVER_REL = 0.1;
export function drivers(terms, key = "diff", n = 3) {
  const pos = [...terms].filter((t) => t[key] > 0).sort((x, y) => y[key] - x[key]);
  if (!pos.length) return [];
  const floor = Math.max(DRIVER_MIN, DRIVER_REL * pos[0][key]);
  return pos.filter((t) => t[key] >= floor).slice(0, n);
}
export const sdLabel = (t) => `${t.symbol} (${t.z >= 0 ? "+" : ""}${t.z.toFixed(1)} sd)`;

/** One plain-language sentence naming the genes that drive the call (and, with a claim, why not the claim). */
export function explainSentence(r, model, names = CLASS_NAMES) {
  const nm = (c) => names[c] || c;
  const list = (xs) => (xs.length > 1 ? `${xs.slice(0, -1).join(", ")} and ${xs[xs.length - 1]}` : xs.join(""));
  const other = r.claim && r.claim.classes.length && !r.claim.classes.includes(r.call) ? r.claim.classes[0] : r.runnerUp;
  const ct = contrast(r.z, model, r.call, other);
  const d = drivers(ct.terms);
  const base = `Called ${nm(r.call)} (p = ${r.callProb.toFixed(2)})`;
  const why = d.length ? ` mainly because of ${list(d.map(sdLabel))}, which favour${d.length === 1 ? "s" : ""} ${nm(r.call)} over ${nm(other)}` : `; no single gene dominates the call over ${nm(other)}`;
  const lab = setLabel(r.alpha ?? model.alpha_default ?? 0.1);
  const setTxt = r.set.length === 0 ? `; no tissue reaches the calibrated threshold, so the ${lab} is empty` : r.set.length > 1 ? `; the ${lab} keeps ${list(r.set.map(nm))}` : "";
  return `${base}${why}${setTxt}.`;
}

/** The samples a recalibration may use: a claim naming exactly one reference tissue, minus `exclude` (the example's
 *  deliberately swapped samples, whose labels are wrong on purpose). Returns { id: class }. */
export function labelledForRecal(results, exclude = []) {
  const skip = new Set(exclude), m = {};
  for (const r of results) if (!skip.has(r.id) && r.claim && r.claim.kind === "class") m[r.id] = r.claim.classes[0];
  return m;
}

// ---- reference map -----------------------------------------------------------------------------------------------------
export function project(z, pca) {
  return pca.loadings.map((L) => L.reduce((a, l, j) => a + l * (z[j] - pca.centre[j]), 0));
}

/** Median, first and third quartile of each panel gene per tissue in the MoTrPAC reference (expr_motrpac.json). */
export function referenceStats(expr, model) {
  const gi = Object.fromEntries(expr.genes.map((g, i) => [g, i]));
  const byT = {};
  expr.samples.forEach((s, j) => { (byT[s.tissue] ||= []).push(j); });
  const q = (arr, p) => { const a = [...arr].sort((x, y) => x - y); const h = (a.length - 1) * p, lo = Math.floor(h); return a[lo] + (a[Math.min(lo + 1, a.length - 1)] - a[lo]) * (h - lo); };
  const out = {};
  for (const [t, idx] of Object.entries(byT)) {
    out[t] = {};
    for (const g of model.genes) {
      if (gi[g.id] === undefined) continue;
      const v = idx.map((j) => expr.values[gi[g.id]][j]).filter((x) => x !== null && Number.isFinite(x));
      out[t][g.id] = { median: q(v, 0.5), q1: q(v, 0.25), q3: q(v, 0.75), n: v.length };
    }
  }
  return out;
}

// ---- the check ---------------------------------------------------------------------------------------------------------
/**
 * Score parsed samples. opts: alpha, scaling ("auto" | "within" | "reference"), q (a recalibrated threshold that
 * replaces the calibration). Returns { scaling, alpha, q, results: [...], counts }.
 */
export function runCheck(model, samples, opts = {}) {
  const alpha = opts.alpha ?? model.alpha_default ?? 0.1;
  const X = samples.map((s) => model.genes.map((g) => (s.values[g.id] === null || s.values[g.id] === undefined ? null : s.values[g.id])));
  const claims = samples.map((s) => mapClaim(s.label));
  const Pref = probabilities(zscoreSource(X, model).Z, model);
  const predRef = new Set(Pref.map((p) => model.classes[p.indexOf(Math.max(...p))]));
  const claimed = new Set(claims.filter((c) => c.kind === "class" || c.kind === "superclass").map((c) => c.classes.join("|")));
  const auto = chooseScaling(samples.length, claimed.size, predRef.size);
  const mode = opts.scaling && opts.scaling !== "auto" ? opts.scaling : auto.mode;
  const Z = mode === "within" ? zscoreWithin(X).Z : zscoreSource(X, model).Z;
  const P = mode === "within" ? probabilities(Z, model) : Pref;
  const cal = { scores: model.calibration.scores, y_idx: model.calibration.y_idx, n_classes: model.classes.length };
  const results = P.map((p, i) => {
    const r = opts.q !== undefined && opts.q !== null ? { set: p.map((v) => v >= 1 - opts.q), q: opts.q } : setsFor(cal, p, "marginal", alpha);
    const set = model.classes.filter((_, k) => r.set[k]);
    const order = p.map((v, k) => [v, k]).sort((a, b) => b[0] - a[0]);
    const res = { id: samples[i].id, label: samples[i].label, claim: claims[i], z: Z[i], p, q: r.q, threshold: 1 - r.q, alpha,
                  call: model.classes[order[0][1]], callProb: order[0][0], runnerUp: model.classes[order[1][1]], runnerUpProb: order[1][0],
                  set, setSize: set.length, missingGenes: model.genes.filter((g, j) => X[i][j] === null).map((g) => g.symbol) };
    res.status = sampleStatus(set.length);
    res.claimStatus = claimStatus(set, claims[i]);
    res.flagged = isFlagged(res);
    return res;
  });
  const counts = { total: results.length, Confident: 0, Ambiguous: 0, Unknown: 0, Consistent: 0, Mismatch: 0, "Can't confirm": 0, "Not in reference": 0, flagged: 0, claimed: 0 };
  for (const r of results) { counts[r.status] += 1; if (r.claimStatus) { counts[r.claimStatus] += 1; counts.claimed += 1; } if (r.flagged) counts.flagged += 1; }
  return { scaling: scalingReport(auto, mode), alpha, q: results[0]?.q, results, counts };
}

// ---- CSV ---------------------------------------------------------------------------------------------------------------
/** One CSV cell: quoted when needed, and neutralised against spreadsheet formula injection (a leading = + - @ tab CR). */
export function csvCell(v) {
  if (v === null || v === undefined) return "";
  let s = String(v);
  if (typeof v !== "number" && /^[=+\-@\t\r]/.test(s)) s = "'" + s;
  return /[",\n\r;]/.test(s) || s !== s.trim() ? `"${s.replace(/"/g, '""')}"` : s;
}

export function resultsToCsv(check, model) {
  const head = ["sample", "claimed", "claimed_classes", "call", "p_call", "runner_up", "p_runner_up", "set", "set_size", "status", "claim_status", "flagged", "alpha", "threshold", "scaling", ...model.classes.map((c) => `p_${c}`)];
  const lines = check.results.map((r) => [r.id, r.label ?? "", r.claim.classes.join(";"), r.call, r.callProb.toFixed(4), r.runnerUp, r.runnerUpProb.toFixed(4), r.set.join(";"), r.setSize, r.status, r.claimStatus ?? "", r.flagged ? "yes" : "no", check.alpha, Number.isFinite(r.threshold) ? r.threshold.toFixed(4) : "-inf", check.scaling.mode, ...r.p.map((v) => v.toFixed(5))].map(csvCell).join(","));
  return [head.join(","), ...lines].join("\n") + "\n";
}

/** The example: BodyMap 21-week samples as a 20-gene log2 CPM table with a claimed_tissue column (the organ), with the
 *  claims of two named samples swapped. rows: from score.js exampleRows. Returns { csv, swapped: [idA, idB] }. */
export function exampleTable(model, rows, swapOrgans = ["Liver", "Kidney"]) {
  const claims = Object.fromEntries(rows.map((r) => [r.id, r.organ]));
  const a = rows.find((r) => r.organ === swapOrgans[0]), b = rows.find((r) => r.organ === swapOrgans[1]);
  if (a && b) { claims[a.id] = b.organ; claims[b.id] = a.organ; }
  const head = ["sample", ...model.genes.map((g) => g.symbol), "claimed_tissue"].join(",");
  const lines = rows.map((r) => [r.id, ...model.genes.map((g) => (r.values[g.id] === null || r.values[g.id] === undefined ? "" : Number(r.values[g.id]).toFixed(4))), claims[r.id]].join(","));
  return { csv: [head, ...lines].join("\n") + "\n", swapped: a && b ? [a.id, b.id] : [], organs: swapOrgans };
}

/** The two "example input" tables of the landing page, built from real data (no typed values).
 *  A: two samples of the example (the first of each of `organs`), the first `nGenes` panel genes in log2 CPM rounded
 *     to 2 decimals, "…", and the claimed tissue (the organ). rows: from score.js exampleRows (the "Try the example" rows).
 *  B: the exported real counts (product.json example_counts): genes × samples, Ensembl IDs, then a "…" row. */
export function exampleInputTables(model, rows, exampleCounts, { organs = ["Adrenal", "Liver"], nGenes = 3 } = {}) {
  const genes = model.genes.slice(0, nGenes);
  const picked = organs.map((o) => rows.find((r) => r.organ === o)).filter(Boolean);
  const a = { columns: ["sample", ...genes.map((g) => g.symbol), "…", "claimed_tissue"],
              rows: picked.map((r) => Object.fromEntries([["sample", r.id], ...genes.map((g) => [g.symbol, r.values[g.id] === null ? "" : Number(r.values[g.id]).toFixed(2)]), ["…", "…"], ["claimed_tissue", r.organ.toLowerCase()]])) };
  const ec = exampleCounts;
  const b = { columns: ["gene_id", ...ec.samples, "…"],
              rows: [...ec.genes.map((g, i) => Object.fromEntries([["gene_id", g], ...ec.samples.map((smp, j) => [smp, String(ec.counts[i][j])]), ["…", "…"]])),
                     Object.fromEntries([["gene_id", "…"], ...ec.samples.map((smp) => [smp, "…"]), ["…", ""]])] };
  return { a, b };
}
