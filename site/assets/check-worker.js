// Parses an upload off the main thread (site/index.html): reads the File, streams it through parseUpload and reports
// progress. Only the panel genes' log2 CPM, library sizes and gene identifiers come back, never the full matrix.
import { parseUpload } from "./check-core.js";

self.onmessage = async (e) => {
  const { file, text, model, opts } = e.data;
  try {
    self.postMessage({ type: "progress", stage: "reading", f: 0 });
    const body = text ?? (await file.text());
    let last = -1;
    const r = parseUpload(body, model, { ...opts, onProgress: (f) => {
      const pct = Math.floor(f * 100);
      if (pct !== last) { last = pct; self.postMessage({ type: "progress", stage: "parsing", f }); }
    } });
    if (r.error) { self.postMessage({ type: "error", message: r.error }); return; }
    // the gene identifier list can be long; the page only needs the species guess and the count
    const { geneIds, ...rest } = r;
    self.postMessage({ type: "done", result: { ...rest, nGeneIds: geneIds ? geneIds.length : 0 } });
  } catch (err) {
    self.postMessage({ type: "error", message: String((err && err.message) || err) });
  }
};
