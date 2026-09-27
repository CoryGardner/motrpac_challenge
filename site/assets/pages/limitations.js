import { mountChrome, loadJSON } from "../site.js";

async function main() {
  await mountChrome("limitations.html");
  try {
    const [P, H] = await Promise.all([loadJSON("data/provenance.json"), loadJSON("data/headline.json")]);
    const rec = P._meta.reconciliation;
    document.getElementById("n-moved").textContent = rec ? `${rec.changed} of ${rec.comparable} comparable numbers changed (its §2)` : "see docs/NUMBERS_RECONCILIATION.md for the count";
    const ex = H.extras;
    document.getElementById("n-orth").textContent = `${ex.orthologs_in_gtex.toLocaleString()} of ${ex.motrpac_genes.toLocaleString()}`;
    document.getElementById("n-inf").textContent = typeof ex.gtex_recal_k20_n3_frac_inf === "number" ? `${(100 * ex.gtex_recal_k20_n3_frac_inf).toFixed(0)} %` : "pending";
  } catch (e) {
    document.getElementById("n-moved").textContent = "see docs/NUMBERS_RECONCILIATION.md for the count";
  }
}

main().catch((e) => console.error(e));
