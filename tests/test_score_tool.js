// The browser scoring tool must reproduce the pipeline: the rat BodyMap samples, run through score.js with the exported
// 20-gene transfer model, get the same probabilities, calls and prediction sets as site/data/samples_bodymap.json
// (which holds the pipeline's own scores from the --save-scores regeneration).
// Run: node tests/test_score_tool.js
import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import { fileURLToPath } from "node:url";
import { dirname, join } from "node:path";

const here = dirname(fileURLToPath(import.meta.url));
const S = await import(join(here, "..", "site", "assets", "score.js"));
const cf = await import(join(here, "..", "site", "assets", "conformal.js"));
const read = (n) => JSON.parse(readFileSync(join(here, "..", "site", "data", n), "utf8"));
const model = read("panel_model.json"), expr = read("expr_bodymap.json"), bm = read("samples_bodymap.json");

let n = 0;
const check = (c, m) => { n += 1; assert.ok(c, m); };
check(model.genes.length === 20 && model.coef.length === model.classes.length && model.coef[0].length === 20, "model shape");
check(model.classes.join() === bm.classes.join(), "class order matches the BodyMap export");
const gi = Object.fromEntries(expr.genes.map((g, i) => [g, i]));
for (const g of model.genes) check(g.id in gi, `expression exported for ${g.symbol} (${g.id})`);

// every BodyMap sample as an upload row (log2 CPM of the 20 genes), z-scored within the upload as the pipeline does
const rows = expr.samples.map((s, j) => ({ id: s.id, values: Object.fromEntries(model.genes.map((g) => [g.id, expr.values[gi[g.id]][j]])) }));
const res = S.scoreSamples(model, rows, { alpha: 0.1 });
check(res.mode === "within" && res.missing.length === 0, "all 20 genes matched, z-scored within the upload");
const byId = Object.fromEntries(res.results.map((r) => [r.id, r]));
let compared = 0, maxDiff = 0;
for (const s of bm.samples) {
  const r = byId[s.id];
  check(r !== undefined, `scored ${s.id}`);
  const pRef = s.p.k20;
  for (let k = 0; k < pRef.length; k++) maxDiff = Math.max(maxDiff, Math.abs(pRef[k] - r.p[k]));
  const refCall = model.classes[pRef.indexOf(Math.max(...pRef))];
  check(r.call === refCall, `${s.id}: call ${r.call} vs ${refCall}`);
  const refSet = cf.setsFor({ scores: bm.calibration.k20.scores, y_idx: bm.calibration.k20.y_idx, n_classes: model.classes.length }, pRef, "marginal", 0.1).set;
  check(refSet.map((b) => (b ? 1 : 0)).join() === model.classes.map((c) => (r.set.includes(c) ? 1 : 0)).join(), `${s.id}: prediction set`);
  compared += 1;
}
check(compared === bm.samples.length && compared >= 300, `compared ${compared} samples`);
check(maxDiff < 1e-3, `max |Δp| ${maxDiff} (the exported expression is rounded to 4 decimals)`);
// the same calibration the pipeline used
check(model.calibration.scores.length === bm.calibration.k20.scores.length, "calibration scores match the BodyMap export");
// small-upload fallback and recalibration
const small = S.scoreSamples(model, rows.slice(0, 3), { alpha: 0.1 });
check(small.mode === "source" && small.results.length === 3, "fewer than 8 samples → MoTrPAC z-statistics");
const labels = Object.fromEntries(bm.samples.slice(0, 6).map((s) => [s.id, model.classes[s.p.k20.indexOf(Math.max(...s.p.k20))]]));
const rc = S.recalibrate(model, res.results, labels, 0.1);
check(rc.n === 6 && rc.q === Infinity, "6 labelled samples at α = 0.10 → rank 7 > 6 → +∞ threshold (full sets), as the pipeline's rule says");
console.log(`ok: ${n} assertions over ${compared} BodyMap samples; max |Δp| = ${maxDiff.toExponential(2)}`);
