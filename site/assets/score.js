// "Score your own samples": pure functions (no DOM) that turn a pasted table of log2 CPM values into tissue calls and
// conformal prediction sets with the exported 20-gene transfer model (site/data/panel_model.json). Used by the Explorer
// page and by tests/test_score_tool.js, which checks that the rat BodyMap samples come out exactly as the pipeline
// scored them. The arithmetic is the pipeline's: z-score each gene within the uploaded set (population sd, as
// scikit-learn's StandardScaler does), then softmax(z · Wᵀ + b); sets use conformal.js with the 15-animal calibration.
import { conformalQuantile, lacScores, setsFor } from "./conformal.js";

export const MIN_SAMPLES_FOR_WITHIN_Z = 8;

/** Parse CSV/TSV text: first row = header, first column = sample id, other columns numeric (blank/NA → null). */
export function parseTable(text) {
  const lines = String(text).replace(/\r/g, "").split("\n").map((l) => l.trim()).filter((l) => l && !l.startsWith("#"));
  if (lines.length < 2) return { header: [], rows: [], error: "the table needs a header row and at least one sample row" };
  const sep = lines[0].includes("\t") ? "\t" : lines[0].includes(";") && !lines[0].includes(",") ? ";" : ",";
  const split = (l) => l.split(sep).map((c) => c.trim().replace(/^"|"$/g, ""));
  const header = split(lines[0]);
  const rows = [];
  for (const l of lines.slice(1)) {
    const cells = split(l);
    if (!cells[0]) continue;
    const values = {};
    header.slice(1).forEach((h, i) => {
      const v = cells[i + 1];
      const num = v === undefined || v === "" || /^(na|nan|null)$/i.test(v) ? null : Number(v);
      values[h] = Number.isFinite(num) ? num : null;
    });
    rows.push({ id: cells[0], values });
  }
  return { header, rows, error: rows.length ? null : "no sample rows" };
}

const norm = (s) => String(s).trim().toLowerCase().replace(/\.\d+$/, "");

/** Match header columns to the model's genes by Ensembl id or symbol (case-insensitive; a version suffix is ignored). */
export function matchGenes(header, model) {
  const byKey = new Map();
  header.slice(1).forEach((h) => byKey.set(norm(h), h));
  const columns = {};
  const missing = [];
  for (const g of model.genes) {
    const h = byKey.get(norm(g.id)) ?? byKey.get(norm(g.symbol)) ?? null;
    columns[g.id] = h;
    if (h === null) missing.push(g);
  }
  return { columns, missing, matched: model.genes.length - missing.length };
}

/** rows × genes matrix in the model's feature order; null where the value is absent. */
export function buildMatrix(rows, columns, model) {
  return rows.map((r) => model.genes.map((g) => (columns[g.id] === null || columns[g.id] === undefined ? null : r.values[columns[g.id]] ?? null)));
}

/** z-score each column within the matrix (population sd); columns with no spread or no values give 0. */
export function zscoreWithin(X) {
  const n = X.length, m = X[0] ? X[0].length : 0;
  const mean = [], sd = [];
  for (let j = 0; j < m; j++) {
    const v = X.map((r) => r[j]).filter((x) => x !== null && Number.isFinite(x));
    const mu = v.length ? v.reduce((a, b) => a + b, 0) / v.length : 0;
    const s = v.length ? Math.sqrt(v.reduce((a, b) => a + (b - mu) ** 2, 0) / v.length) : 0;
    mean.push(mu); sd.push(s);
  }
  const Z = X.map((r) => r.map((x, j) => (x === null || !Number.isFinite(x) || sd[j] === 0 ? 0 : (x - mean[j]) / sd[j])));
  return { Z, mean, sd, n };
}

/** z-score with the MoTrPAC statistics the model was trained on (the fallback for small uploads). */
export function zscoreSource(X, model) {
  const { mean, scale } = model.source_z;
  return { Z: X.map((r) => r.map((x, j) => (x === null || !Number.isFinite(x) || !scale[j] ? 0 : (x - mean[j]) / scale[j]))), mean, sd: scale, n: X.length };
}

export function softmax(logits) {
  const mx = Math.max(...logits);
  const e = logits.map((v) => Math.exp(v - mx));
  const s = e.reduce((a, b) => a + b, 0);
  return e.map((v) => v / s);
}

export function probabilities(Z, model) {
  return Z.map((z) => softmax(model.coef.map((w, c) => w.reduce((acc, wj, j) => acc + wj * z[j], 0) + model.intercept[c])));
}

/**
 * Score parsed rows. Returns { mode, note, missing, results: [{ id, p, call, callProb, set, setSize, kind, q }] }.
 * opts: alpha (0.10), variant ("marginal" | "mondrian" | "floored"), minN (8), calibration override ({scores, y_idx}).
 */
export function scoreSamples(model, rows, opts = {}) {
  const alpha = opts.alpha ?? model.alpha_default ?? 0.1;
  const variant = opts.variant || "marginal";
  const minN = opts.minN ?? MIN_SAMPLES_FOR_WITHIN_Z;
  const header = ["id", ...Object.keys(rows[0]?.values || {})];
  const { columns, missing, matched } = matchGenes(header, model);
  const X = buildMatrix(rows, columns, model);
  const within = rows.length >= minN;
  const { Z } = within ? zscoreWithin(X) : zscoreSource(X, model);
  const P = probabilities(Z, model);
  const cal = opts.calibration || model.calibration;
  const calibration = { scores: cal.scores, y_idx: cal.y_idx, n_classes: model.classes.length };
  const results = P.map((p, i) => {
    const r = opts.q !== undefined && opts.q !== null ? { set: p.map((v) => v >= 1 - opts.q), q: opts.q } : setsFor(calibration, p, variant, alpha);
    const chosen = model.classes.filter((_, k) => r.set[k]);
    const top = p.indexOf(Math.max(...p));
    return { id: rows[i].id, p, call: model.classes[top], callProb: p[top], set: chosen, setSize: chosen.length,
             kind: chosen.length === 0 ? "abstains" : chosen.length === 1 ? "confident" : "ambiguous", q: r.q };
  });
  const note = within
    ? `genes z-scored within your ${rows.length} samples (as the pipeline does for a new dataset)`
    : `fewer than ${minN} samples, so the genes were z-scored with the MoTrPAC means and scales instead of within your set; calls are less reliable`;
  return { mode: within ? "within" : "source", note, missing, matched, alpha, variant, results };
}

/** Recalibrate on user-labelled samples: LAC threshold at α from those samples' scores; returns coverage on them. */
export function recalibrate(model, results, labels, alpha) {
  const idx = [], probs = [];
  for (const r of results) {
    const t = labels[r.id];
    if (t && model.classes.includes(t)) { idx.push(model.classes.indexOf(t)); probs.push(r.p); }
  }
  if (idx.length < 3) return { error: `label at least 3 samples (${idx.length} labelled)` };
  const scores = lacScores(probs, idx);
  const q = conformalQuantile(scores, alpha);
  const covered = probs.filter((p, i) => p[idx[i]] >= 1 - q).length;
  return { q, n: idx.length, coverage: covered / idx.length, scores };
}

export function templateCsv(model, examples = []) {
  const head = ["sample", ...model.genes.map((g) => g.symbol)].join(",");
  const lines = examples.map((e) => [e.id, ...model.genes.map((g) => (e.values[g.id] === null || e.values[g.id] === undefined ? "" : Number(e.values[g.id]).toFixed(4)))].join(","));
  return [head, ...lines].join("\n") + "\n";
}

export function resultsCsv(res, model) {
  const head = ["sample", "call", "p_call", "set", "set_size", "kind", ...model.classes.map((c) => `p_${c}`)].join(",");
  const lines = res.results.map((r) => [r.id, r.call, r.callProb.toFixed(4), `"${r.set.join(";")}"`, r.setSize, r.kind, ...r.p.map((v) => v.toFixed(5))].join(","));
  return [head, ...lines].join("\n") + "\n";
}
