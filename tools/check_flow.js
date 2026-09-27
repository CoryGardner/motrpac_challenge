// The Check samples demo flow, run by `make screenshots` after the page renders (needs Chrome + playwright):
// "Try the example" → the "Start here" note names the swapped pair → both sit at the top of the table as Mismatch →
// a row opens the sample drawer → "Recalibrate on my labelled samples" runs; plus drag-and-drop of the real count
// fixture, the "list of accepted names" link, the example input tables, and a refused table. Exit 1 on any failure.
// Usage: node tools/check_flow.js --base http://localhost:8000
import { chromium } from "playwright";
import fs from "node:fs";
import path from "node:path";
import zlib from "node:zlib";
import { fileURLToPath } from "node:url";
const __dirname = path.dirname(fileURLToPath(import.meta.url));

const args = {};
for (let i = 2; i < process.argv.length; i += 2) args[process.argv[i].replace(/^--/, "")] = process.argv[i + 1];
const base = args.base || "http://localhost:8000";
const fails = [];
const check = (ok, msg) => { console.log(`${ok ? "ok  " : "FAIL"} ${msg}`); if (!ok) fails.push(msg); };

(async () => {
  const browser = await chromium.launch({ channel: "chrome" });
  for (const theme of ["light", "dark"]) {
    const page = await browser.newPage({ viewport: { width: 1440, height: 900 }, colorScheme: theme });
    const errors = [];
    page.on("pageerror", (e) => errors.push(String(e)));
    page.on("console", (m) => { if (m.type() === "error") errors.push(m.text()); });
    await page.goto(`${base}/index.html`);
    await page.waitForSelector("#trust b");
    // the demo flow
    await page.click("#btn-example");
    await page.waitForSelector("#results-guide .callout", { timeout: 15000 });
    const guide = await page.textContent("#results-guide");
    check(/both flagged Mismatch/.test(guide), `${theme}: "Start here" names the swapped pair`);
    const top = await page.$$eval("#tbl-results tbody tr", (rs) => rs.slice(0, 2).map((r) => r.innerText));
    check(top.length === 2 && top.every((t) => t.includes("Mismatch")), `${theme}: the swapped pair tops the table`);
    await page.click("#tbl-results tbody tr");
    await page.waitForSelector("#drawer:not([hidden]) #fig-contrib .chart");
    check(/mainly because of|no single gene dominates/.test(await page.textContent("#drawer-body")), `${theme}: the drawer explains the call`);
    await page.click("#results-guide button");
    await page.waitForFunction(() => /Recalibrated on/.test(document.getElementById("results-guide").textContent));
    const after = await page.$$eval("#tbl-results tbody tr", (rs) => rs.slice(0, 2).map((r) => r.innerText));
    check(after.every((t) => t.includes("Mismatch")), `${theme}: both swaps stay Mismatch after recalibrating`);
    // the accepted-names link
    await page.click("#accepted-link");
    await page.waitForTimeout(400);
    check(await page.$eval("#accepted", (d) => d.open), `${theme}: "list of accepted names" opens the list`);
    // the example input tables
    await page.click("#example-tables summary");
    await page.waitForSelector("#ex-table-b table");
    check((await page.$$("#ex-table-a tbody tr")).length === 2 && (await page.$$("#ex-table-b tbody tr")).length >= 3, `${theme}: example input tables built`);
    // drag and drop the real count fixture (unzipped)
    const csv = zlib.gunzipSync(fs.readFileSync(path.join(__dirname, "..", "tests", "fixtures", "bodymap_counts_subset.csv.gz"))).toString("utf8");
    const dt = await page.evaluateHandle((text) => { const d = new DataTransfer(); d.items.add(new File([text], "bodymap_counts_subset.csv", { type: "text/csv" })); return d; }, csv);
    await page.dispatchEvent("#dropzone", "dragover", { dataTransfer: dt });
    check(await page.$eval("#dropzone", (d) => d.classList.contains("dragover")), `${theme}: the drop zone highlights on dragover`);
    await page.dispatchEvent("#dropzone", "drop", { dataTransfer: dt });
    await page.waitForFunction(() => /Read bodymap_counts_subset\.csv/.test(document.getElementById("progress").textContent), null, { timeout: 20000 });
    check(/raw counts/.test(await page.textContent("#progress")) && (await page.$$("#tbl-results tbody tr")).length === 4, `${theme}: the dropped count matrix is read and scored`);
    // a table without panel genes is refused
    await page.fill("#in-text", "a,b\n1,2");
    await page.click("#btn-check");
    await page.waitForFunction(() => /none of the/.test(document.getElementById("progress").textContent));
    check(await page.$eval("#results", (r) => r.hidden), `${theme}: a table without panel genes is refused`);
    check(errors.length === 0, `${theme}: no console errors${errors.length ? ": " + errors.join(" | ") : ""}`);
    await page.close();
  }
  // the #accepted hash on load
  const p2 = await browser.newPage();
  await p2.goto(`${base}/index.html#accepted`);
  await p2.waitForSelector("#trust b");
  await p2.waitForTimeout(300);
  check(await p2.$eval("#accepted", (d) => d.open), "#accepted in the URL opens the list on load");
  await browser.close();
  console.log(fails.length ? `\n${fails.length} flow check(s) failed` : "\nall flow checks passed");
  process.exit(fails.length ? 1 : 0);
})();
