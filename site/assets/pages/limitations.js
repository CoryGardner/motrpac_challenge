import { mountChrome, loadJSON } from "../site.js";

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
}

main().catch((e) => console.error(e));
