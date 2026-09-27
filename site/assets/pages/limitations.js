import { mountChrome, loadJSON } from "../site.js";

async function main() {
  await mountChrome("limitations.html");
  try {
    const P = await loadJSON("data/provenance.json");
    // count of provenance entries whose pre-fix value differs is written by the reconciliation; here we count entries
    // that cite a phase rerun after the fix (06, 08, 12, 13, 14) as the population that could have moved
    const rerun = P.entries.filter((e) => e.file && /results\/(06|08|12|13|14)_/.test(e.file)).length;
    document.getElementById("n-moved").textContent = `102 of 182 comparable (docs/NUMBERS_RECONCILIATION.md §2; ${rerun} exported entries cite a rerun phase)`;
  } catch (e) {
    document.getElementById("n-moved").textContent = "see docs/NUMBERS_RECONCILIATION.md for the count";
  }
}

main().catch((e) => console.error(e));
