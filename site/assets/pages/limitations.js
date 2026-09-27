import { mountChrome, loadJSON } from "../site.js";
import { minLabelled } from "../score.js";
import { recalNote, compositionNote } from "../notes.js";

async function main() {
  await mountChrome("limitations.html");
  try {
    const H = await loadJSON("data/headline.json");
    const ex = H.extras;
    document.getElementById("n-orth").textContent = `${ex.orthologs_in_gtex.toLocaleString()} of ${ex.motrpac_genes.toLocaleString()}`;
    document.getElementById("n-inf").textContent = typeof ex.gtex_recal_k20_n3_frac_inf === "number" ? `${(100 * ex.gtex_recal_k20_n3_frac_inf).toFixed(0)} %` : "pending";
  } catch (e) {
    console.error(e);
  }
  try {
    const [PR, H] = await Promise.all([loadJSON("data/product.json"), loadJSON("data/headline.json")]);
    document.getElementById("lim-recal").textContent = recalNote(PR, H.design, minLabelled(H.design.alpha));
    document.getElementById("lim-comp").textContent = compositionNote(PR);
  } catch (e) {
    console.error(e);
  }
  // the multiomic follow-up paragraph (branch multiomic-overnight); every value from site/data/multiomic.json
  try {
    const M = await loadJSON("data/multiomic.json");
    const pc1 = M.scales.rows.find((r) => r.PC === "PC1");
    const jr = M.ladder_species.find((r) => r.target.startsWith("protein → Jiang 2020 (cleaned") && r.model === "k20");
    const r5 = M.recalibration.find((r) => r.source.startsWith("Phase 3 protein → Jiang 2020 (cleaned") && r.model === "k20" && r.n_recal === 5);
    const set = (id, v) => { const n = document.getElementById(id); if (n) n.textContent = v; };
    set("mo-r2", pc1.r2_rii.toFixed(3)); set("mo-r2-ratio", pc1.r2_ratio.toFixed(4));
    set("mo-acc", jr.accuracy.toFixed(3)); set("mo-n", String(jr.n_samples)); set("mo-donors", String(jr.n_individuals)); set("mo-ci", `${jr.accuracy_ci[0].toFixed(2)}–${jr.accuracy_ci[1].toFixed(2)}`);
    set("mo-cov", jr.coverage.toFixed(3)); set("mo-cov5", r5.coverage_recalibrated.toFixed(3)); set("mo-size5", r5.set_size_recalibrated.toFixed(2));
    set("mo-pairs", `${M.design_jiang.n_pairs_estimable} of ${M.design_jiang.n_pairs_total}`);
  } catch (e) {
    console.error(e);
  }
}

main().catch((e) => console.error(e));
