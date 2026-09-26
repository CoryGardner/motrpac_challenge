// The JS conformal port must reproduce the Python prediction sets exactly, including infinite thresholds.
// Run: node tests/test_site_conformal.js
import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import { fileURLToPath } from "node:url";
import { dirname, join } from "node:path";

const here = dirname(fileURLToPath(import.meta.url));
const cf = await import(join(here, "..", "site", "assets", "conformal.js"));

let n = 0;
function check(cond, msg) { n += 1; assert.ok(cond, msg); }

// ---- unit cases: the rank rule ---------------------------------------------------------------
const s22 = Array.from({ length: 22 }, (_, i) => (i + 1) / 100);        // 0.01 .. 0.22
check(cf.conformalQuantile(s22, 0.10) === 0.21, "n=22, α=0.10: rank ceil(23·0.9)=21 → 21st smallest = 0.21");
check(cf.conformalQuantile(s22, 0.05) === 0.22, "n=22, α=0.05: rank ceil(23·0.95)=22 ≤ 22 → the largest score (regime 19–38 at α=0.05)");
// (n+1)(1-α) exactly an integer: n=19, α=0.05 → 20·0.95 = 19.000000000000004 in floating point; the −1e-9 tolerance keeps rank 19
const s19 = Array.from({ length: 19 }, (_, i) => (i + 1) / 100);
check(cf.conformalQuantile(s19, 0.05) === 0.19, "n=19, α=0.05: rank 19 (tolerance), not +∞");
check(cf.conformalQuantile(Array.from({ length: 18 }, (_, i) => (i + 1) / 100), 0.05) === Infinity, "n=18, α=0.05: rank 19 > 18 → +∞");
check(cf.conformalQuantile([], 0.1) === Infinity, "no calibration scores → +∞");
check(cf.conformalQuantile([0.5, 0.1, 0.9], 0.5) === 0.5, "unsorted input is sorted first (rank 2 of 3)");
// LAC set with an infinite threshold contains every class
check(cf.predictSet([0.0, 0.2, 0.8], Infinity).every(Boolean), "+∞ threshold → full set");
check(cf.predictSet([0.0, 0.2, 0.8], 0.25).join() === "false,false,true", "p ≥ 1−q rule");
// Mondrian fallback and floor
const qs = cf.perClassQuantiles([0.1, 0.2, 0.3, 0.4], [0, 0, 1, 1], 0.10, 3, 0.35, null);
check(qs[2] === 0.35, "class with no calibration score falls back to the marginal threshold");
check(qs[0] === Infinity && qs[1] === Infinity, "two scores per class at α=0.10: rank 3 > 2 → +∞");
const qf = cf.perClassQuantiles([0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 0.13, 0.14, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9],
  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1], 0.10, 2, 0.5, 0.5);
check(qf[0] === 0.5 && qf[1] === 0.9, "floored: per-class threshold never below the marginal one");

// ---- fixtures: exact reproduction of Python sets ------------------------------------------------
const fx = JSON.parse(readFileSync(join(here, "..", "site", "data", "conformal_fixtures.json"), "utf8"));
let cases = 0;
for (const c of fx.cases) {
  const cal = fx.calibrations[c.calibration];            // {scores: [...], y_idx: [...], n_classes}
  const r = cf.setsFor(cal, c.p, c.variant, c.alpha);
  const got = r.set.map((b) => (b ? 1 : 0));
  assert.deepEqual(got, c.expected_set, `case ${c.id}: ${c.calibration} ${c.variant} α=${c.alpha}`);
  if (c.expected_q === "inf") check(r.q === Infinity, `case ${c.id}: q is +∞`);
  else check(Math.abs(r.q - c.expected_q) < 1e-12, `case ${c.id}: q matches (${r.q} vs ${c.expected_q})`);
  cases += 1;
}
check(cases >= 100, `fixture has ${cases} cases (need ≥ 100)`);
const inf = fx.cases.filter((c) => c.expected_q === "inf").length;
check(inf > 0, "fixture includes at least one infinite-threshold case");
console.log(`ok: ${n} assertions, ${cases} fixture cases (${inf} with infinite threshold)`);
