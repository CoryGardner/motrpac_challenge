// Split-conformal prediction sets for classification — a faithful port of src/tfp/conformal.py
// (conformal_quantile, conformal_quantile_per_class, predict_sets, predict_sets_conditional; LAC scores).
//
// The threshold rule is the textbook one (fixed in the pipeline on 2026-09-25): the ceil((n+1)(1−α))-th
// smallest of the n calibration scores, and +∞ (every class in every set) when that rank exceeds n.
// The −1e-9 tolerance mirrors the Python code: (n+1)(1−α) can be an integer plus floating-point noise.

export function conformalQuantile(scores, alpha) {
  const n = scores.length;
  if (n === 0) return Infinity;
  const k = Math.ceil((n + 1) * (1 - alpha) - 1e-9);
  if (k > n) return Infinity;
  const sorted = Array.from(scores, Number).sort((a, b) => a - b);
  return sorted[Math.max(k, 1) - 1];
}

// Nonconformity of the true class under LAC: 1 − p̂(true class).
export function lacScores(probs, yIdx) {
  return probs.map((p, i) => 1 - p[yIdx[i]]);
}

// One threshold per class from that class's calibration scores (Mondrian). A class with no calibration
// score gets `fallback` (the marginal threshold); with `floor` every class threshold is raised to at
// least the floor (the marginal-floor Mondrian variant).
export function perClassQuantiles(scores, yIdx, alpha, nClasses, fallback = null, floor = null) {
  const q = new Array(nClasses).fill(NaN);
  for (let c = 0; c < nClasses; c++) {
    const sc = [];
    for (let i = 0; i < scores.length; i++) if (yIdx[i] === c) sc.push(scores[i]);
    if (sc.length) q[c] = conformalQuantile(sc, alpha);
  }
  for (let c = 0; c < nClasses; c++) {
    if (Number.isNaN(q[c]) && fallback !== null) q[c] = fallback;
    if (floor !== null && !Number.isNaN(q[c])) q[c] = Math.max(q[c], floor);
  }
  return q;
}

// LAC set with one threshold: class j is in the set when p[j] ≥ 1 − q (q = +∞ → every class).
export function predictSet(p, q) {
  const thr = 1 - q;
  return p.map((v) => v >= thr);
}

// LAC set with a class-specific threshold.
export function predictSetConditional(p, qs) {
  return p.map((v, j) => v >= 1 - qs[j]);
}

// Convenience used by the explorer and the fixture test.
// calibration: {scores: number[], y_idx: number[], n_classes: number}; variant: marginal | mondrian | floored.
export function setsFor(calibration, p, variant, alpha) {
  const q = conformalQuantile(calibration.scores, alpha);
  if (variant === "marginal") return { set: predictSet(p, q), q, qs: null };
  const nClasses = calibration.n_classes ?? p.length;
  const qs = perClassQuantiles(calibration.scores, calibration.y_idx, alpha, nClasses, q, variant === "floored" ? q : null);
  return { set: predictSetConditional(p, qs), q, qs };
}

// Zero-error sizing: calibration units needed so that a panel with no calibration error is certified at
// (α, δ): n ≥ ln δ / ln(1 − α).
export function animalsNeeded(alpha, delta) {
  return Math.ceil(Math.log(delta) / Math.log(1 - alpha));
}
