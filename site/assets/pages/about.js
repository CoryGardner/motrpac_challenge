import { mountChrome, loadJSON, el } from "../site.js";

async function main() {
  await mountChrome("about.html");
  try {
    const M = await loadJSON("data/manifest.json");
    const dl = document.getElementById("build-info");
    const row = (k, v) => dl.append(el("dt", {}, k), el("dd", {}, v));
    row("results git hash", M.git_hash);
    row("site data generated", M.generated);
    row("site data files", String((M.site_data_files || []).length || "see manifest"));
  } catch (e) { /* manifest missing: nothing to show */ }
}

main().catch((e) => console.error(e));
