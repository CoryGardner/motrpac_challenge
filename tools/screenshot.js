// Render pass: screenshot every page at 1440×900 and 390×844 in both themes, collect console errors,
// measure load time and the visible words in <main>, and optionally block the CDN to prove the vendored Plotly fallback works.
// Usage: node tools/screenshot.js [--base http://localhost:8765] [--pages index,explore] [--block-cdn] [--strict-words] [--out site/_screenshots]
// --strict-words fails a render whose visible word count exceeds the page's BUDGET.
// Needs playwright (npm install in tools/) and a Chrome/Chromium: it uses the system Chrome channel.
import { mkdirSync, writeFileSync } from "node:fs";
import { join } from "node:path";
import { createRequire } from "node:module";

const require = createRequire(import.meta.url);
const { chromium } = require("playwright");

// --key=value, --key value, or a bare --flag
const args = {};
for (let i = 2, argv = process.argv; i < argv.length; i += 1) {
  const a = argv[i];
  if (!a.startsWith("--")) continue;
  const eq = a.indexOf("=");
  if (eq > 0) args[a.slice(2, eq)] = a.slice(eq + 1);
  else if (i + 1 < argv.length && !argv[i + 1].startsWith("--")) { args[a.slice(2)] = argv[i + 1]; i += 1; }
  else args[a.slice(2)] = true;
}
const base = args.base || "http://localhost:8765";
const out = args.out || "site/_screenshots";
const pages = (args.pages ? String(args.pages).split(",") : ["index", "explore", "transfer", "exercise", "fingerprint", "identifiability", "methods", "limitations", "about"]);
// visible-word budgets per page (main.innerText), enforced with --strict-words
const BUDGET = { index: 700, explore: 1000, transfer: 1300, exercise: 900, fingerprint: 700, identifiability: 900, methods: 1200, limitations: 700, about: 650 };
const strictWords = !!args["strict-words"];
const sizes = [{ name: "desktop", width: 1440, height: 900 }, { name: "phone", width: 390, height: 844 }];
const themes = ["light", "dark"];
mkdirSync(out, { recursive: true });

const browser = await chromium.launch({ channel: "chrome", headless: true });
const report = [];
let failures = 0;
for (const page of pages) {
  for (const theme of themes) {
    for (const size of sizes) {
      const ctx = await browser.newContext({ viewport: { width: size.width, height: size.height }, colorScheme: theme, deviceScaleFactor: 1 });
      if (args["block-cdn"]) await ctx.route(/cdn\.jsdelivr\.net|cdnjs\.cloudflare\.com/, (route) => route.abort());
      const p = await ctx.newPage();
      const errors = [];
      p.on("console", (m) => {
        if (m.type() !== "error") return;
        // in offline mode the blocked CDN request is the expected failure that triggers the vendored fallback
        if (args["block-cdn"] && /Failed to load resource: net::ERR_FAILED/.test(m.text())) return;
        errors.push(m.text());
      });
      p.on("pageerror", (e) => errors.push("pageerror: " + e.message));
      p.on("requestfailed", (r) => { if (!args["block-cdn"] || !/jsdelivr|cdnjs/.test(r.url())) errors.push("requestfailed: " + r.url()); });
      const t0 = Date.now();
      const url = `${base}/${page}.html?theme=${theme}`;
      try {
        await p.goto(url, { waitUntil: "networkidle", timeout: 60000 });
        await p.waitForFunction(() => document.querySelectorAll(".chart .plot-container").length > 0 || document.querySelectorAll(".chart").length === 0, null, { timeout: 30000 }).catch(() => {});
        await p.waitForTimeout(600);
      } catch (e) {
        errors.push("navigation: " + e.message);
      }
      const ms = Date.now() - t0;
      const timing = await p.evaluate(() => { const n = performance.getEntriesByType("navigation")[0]; return n ? { domContentLoaded: Math.round(n.domContentLoadedEventEnd), load: Math.round(n.loadEventEnd) } : null; }).catch(() => null);
      const charts = await p.evaluate(() => ({ charts: document.querySelectorAll(".chart").length, rendered: document.querySelectorAll(".chart .plot-container").length,
        emptyCharts: [...document.querySelectorAll(".chart")].filter((c) => !c.querySelector(".plot-container")).length,
        overflowX: document.documentElement.scrollWidth > document.documentElement.clientWidth + 1,
        imgsNoAlt: [...document.images].filter((i) => !i.hasAttribute("alt")).length,
        words: (() => { const m = document.querySelector("main"); return m ? m.innerText.trim().split(/\s+/).filter(Boolean).length : 0; })(),
        plotly: !!window.Plotly, plotlySrc: [...document.scripts].filter((s) => /plotly/.test(s.src)).map((s) => s.src) })).catch(() => ({}));
      const file = join(out, `${page}-${theme}-${size.name}${args["block-cdn"] ? "-offline" : ""}.png`);
      await p.screenshot({ path: file, fullPage: true }).catch((e) => errors.push("screenshot: " + e.message));
      const budget = BUDGET[page] ?? null;
      const overBudget = budget !== null && (charts.words ?? 0) > budget;
      const row = { page, theme, size: size.name, ms, timing, ...charts, budget, overBudget, errors };
      if (errors.length || charts.emptyCharts || charts.overflowX || (strictWords && overBudget)) failures += 1;
      report.push(row);
      const wordsNote = `words ${String(charts.words ?? "?").padStart(4)}${budget !== null ? `/${budget}` : ""}${overBudget ? (strictWords ? " OVER BUDGET" : " (over budget)") : ""}`;
      console.log(`${page.padEnd(16)} ${theme.padEnd(5)} ${size.name.padEnd(7)} ${String(ms).padStart(5)} ms  charts ${charts.rendered}/${charts.charts}  overflowX ${charts.overflowX}  ${wordsNote}  errors ${errors.length}${errors.length ? "\n    " + errors.join("\n    ") : ""}`);
      await ctx.close();
    }
  }
}
await browser.close();
writeFileSync(join(out, "report.json"), JSON.stringify(report, null, 1));
console.log(failures ? `\n${failures} page renders with problems` : "\nall page renders clean");
process.exit(failures ? 1 : 0);
