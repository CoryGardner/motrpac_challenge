import { mountChrome, loadJSON, el, fmt } from "../site.js";

async function main() {
  await mountChrome("methods.html");
  const [M, H, CERT] = await Promise.all([loadJSON("data/manifest.json"), loadJSON("data/headline.json"), loadJSON("data/certificate.json")]);
  const opa = H.ladder.find((r) => r.rung_id === "in_distribution" && r.model === "full" && r.variant === "marginal" && r.calibration === "one_per_animal");
  const pooled = H.ladder.find((r) => r.rung_id === "in_distribution" && r.model === "full" && r.variant === "marginal" && r.calibration === "pooled");
  document.getElementById("p-quantile-fix").replaceChildren(
    `The rule was corrected on 2026-09-25: the previous code took one rank too high for n ≥ 19 at α = 0.10 and returned the largest score instead of +∞ when the rank exceeded n. Phases 06, 08 and 12–14 were rerun; the in-distribution coverage of the full model at α = 0.10 is ${fmt(opa?.coverage)} with one vial per animal and ${fmt(pooled?.coverage)} with pooled vials (results/06_conformal/TRNSCRPT/coverage.csv). Every pre-fix number that moved is listed in docs/NUMBERS_RECONCILIATION.md.`,
  );
  const dl = document.getElementById("versions");
  const row = (k, v) => dl.append(el("dt", {}, k), el("dd", {}, v));
  for (const [k, v] of Object.entries(M.versions || {})) row(k, v);
  for (const [k, v] of Object.entries(M.data || {})) row(k, v);
  row("results git hash", M.git_hash);
  row("site data generated", M.generated);
  row("results phases used", Object.keys(M.phases).join(", "));
  row("zero-error sizing (α, δ → n)", CERT.sizing.map((r) => `(${r.alpha}, ${r.delta}) → ${r.n_zero_error_needed}`).join("; "));
}

main().catch((e) => { console.error(e); });
