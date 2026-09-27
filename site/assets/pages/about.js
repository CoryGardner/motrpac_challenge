import { mountChrome, loadJSON, el } from "../site.js";

// Labels for the manifest's data_access keys; anything unlisted falls back to the key itself.
const DATA_LABELS = { motrpac: "MoTrPAC data", bodymap: "Rat BodyMap", gtex: "GTEx" };

function row(dl, label, value) {
  dl.append(el("dt", {}, label), el("dd", {}, value));
}

// what · how · accessed <date> · documented in <code>file</code>
function accessValue(a) {
  const text = [a.what, a.how, a.date ? `accessed ${a.date}` : null].filter(Boolean).join(" · ");
  return a.documented_in ? [`${text} · documented in `, el("code", {}, a.documented_in)] : text;
}

async function main() {
  await mountChrome("about.html");
  const details = document.getElementById("versions-details");
  if (details && location.hash === "#versions") details.open = true;   // a deep link shows the list
  const dl = document.getElementById("versions-kv");
  if (!dl) return;
  try {
    const M = await loadJSON("data/manifest.json");
    for (const [name, version] of Object.entries(M.versions || {})) row(dl, name, String(version));
    for (const [key, access] of Object.entries(M.data_access || {})) row(dl, DATA_LABELS[key] || key, accessValue(access || {}));
    row(dl, "results snapshot", M.results_snapshot || `none (read from ${M.results_dir || "results"}/)`);
    row(dl, "export git hash", M.git_hash || "pending");
    row(dl, "generated", M.generated || "pending");
    row(dl, "site data files", String((M.site_data_files || []).length || "see manifest"));
  } catch (e) {
    row(dl, "manifest", "pending");
  }
}

main().catch((e) => console.error(e));
