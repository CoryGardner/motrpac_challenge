// Shared chrome and components: nav, footer, theme toggle, data loader, stat tile, callouts, tables.
// All text goes through textContent (never innerHTML with data).

export const REPO_URL = "https://github.com/CoryGardner/motrpac_challenge";
const PAGES = [
  ["index.html", "Home"], ["explore.html", "Explorer"], ["fingerprint.html", "Panel"], ["transfer.html", "Transfer"],
  ["identifiability.html", "Identifiability"], ["exercise.html", "Exercise"], ["methods.html", "Methods"], ["limitations.html", "Limitations"],
  ["about.html", "About"], ["multiomic.html", "Multiomic · follow-up (branch)"],
];
const cache = new Map();

export function el(tag, attrs = {}, children = []) {
  const node = document.createElement(tag);
  for (const [k, v] of Object.entries(attrs || {})) {
    if (v === undefined || v === null) continue;
    if (k === "class") node.className = v;
    else if (k === "dataset") Object.assign(node.dataset, v);
    else if (k.startsWith("on") && typeof v === "function") node.addEventListener(k.slice(2), v);
    else node.setAttribute(k, v);
  }
  const kids = Array.isArray(children) ? children : [children];
  for (const c of kids) {
    if (c === null || c === undefined) continue;
    node.appendChild(typeof c === "string" || typeof c === "number" ? document.createTextNode(String(c)) : c);
  }
  return node;
}

export async function loadJSON(path) {
  if (!cache.has(path)) {
    cache.set(path, fetch(path).then((r) => {
      if (!r.ok) throw new Error(`${path}: ${r.status}`);
      return r.json();
    }));
  }
  return cache.get(path);
}

// ---- theme ---------------------------------------------------------------------------------------
export function initTheme() {
  const params = new URLSearchParams(location.search);
  let theme = params.get("theme");
  if (!theme) {
    try { theme = localStorage.getItem("tfp-theme"); } catch (e) { theme = null; }
  }
  if (theme === "dark" || theme === "light") document.documentElement.setAttribute("data-theme", theme);
}
function currentTheme() {
  const set = document.documentElement.getAttribute("data-theme");
  if (set) return set;
  return window.matchMedia && window.matchMedia("(prefers-color-scheme: dark)").matches ? "dark" : "light";
}
export function toggleTheme() {
  const next = currentTheme() === "dark" ? "light" : "dark";
  document.documentElement.setAttribute("data-theme", next);
  try { localStorage.setItem("tfp-theme", next); } catch (e) { /* private mode */ }
  document.dispatchEvent(new CustomEvent("themechange", { detail: next }));
  updateToggle();
  syncThemeColor();
}
// Keep <meta name="theme-color"> (the browser chrome colour on phones) equal to the page background token.
function syncThemeColor() {
  const page = getComputedStyle(document.documentElement).getPropertyValue("--page").trim();
  if (!page) return;
  document.querySelectorAll('meta[name="theme-color"]').forEach((m) => m.setAttribute("content", page));
}
function updateToggle() {
  const b = document.querySelector(".theme-toggle");
  if (b) b.textContent = currentTheme() === "dark" ? "☀ Light" : "☾ Dark";
}

// ---- chrome ---------------------------------------------------------------------------------------
export async function mountChrome(page) {
  initTheme();
  syncThemeColor();
  window.addEventListener("load", syncThemeColor, { once: true });
  const header = document.querySelector("header.site-header") || document.body.insertBefore(el("header", { class: "site-header" }), document.body.firstChild);
  const logo = el("img", { class: "mark", src: "assets/brand/badge-192.png", alt: "", width: 44, height: 44, decoding: "async" });
  const nav = el("nav", { class: "primary", "aria-label": "Site" });
  for (const [href, label] of PAGES) {
    const a = el("a", { href }, label);
    if (href === page) a.setAttribute("aria-current", "page");
    nav.appendChild(a);
  }
  const toggle = el("button", { class: "theme-toggle", type: "button", "aria-label": "Toggle colour theme", onclick: toggleTheme }, "");
  header.replaceChildren(el("div", { class: "wrap" }, [el("a", { class: "brand", href: "index.html" }, [logo, "Tissue Fingerprints"]), nav, toggle]));
  updateToggle();
  document.body.insertBefore(el("a", { class: "skip", href: "#main" }, "Skip to content"), document.body.firstChild);
  const footer = document.querySelector("footer.site-footer") || document.body.appendChild(el("footer", { class: "site-footer" }));
  footer.replaceChildren(el("div", { class: "wrap" }, [
    el("span", {}, ["Code and data exports: ", el("a", { href: REPO_URL }, "github.com/CoryGardner/motrpac_challenge"), " (MIT)."]),
    el("span", {}, ["Every number on this site is read from the result tables; ", el("code", {}, "site/data/provenance.json"), " records where."]),
  ]));
  window.addEventListener("resize", () => { /* Plotly handles responsive */ });
}

// ---- components -----------------------------------------------------------------------------------
export function fmt(v, digits = 3) {
  if (v === null || v === undefined || Number.isNaN(v)) return "pending";
  if (v === "inf") return "+∞";
  if (typeof v === "string") return v;
  if (Number.isInteger(v) && digits === 0) return String(v);
  return Number(v).toFixed(digits);
}
export function pct(v, digits = 0) {
  if (v === null || v === undefined) return "pending";
  return (100 * v).toFixed(digits) + " %";
}

export function statTile(t) {
  const tile = el("div", { class: "tile", id: t.id });
  const value = el("div", { class: "value" });
  if (t.value === null || t.value === undefined) {
    value.textContent = "pending";
    tile.classList.add("pending");
  } else if (t.format === "of") {
    value.append(document.createTextNode(String(t.value)), el("span", { class: "of" }, `of ${t.total}`));
  } else {
    value.textContent = t.format === "pct1" ? pct(t.value, 1) : fmt(t.value, Number(t.format || 3));
  }
  tile.append(value);
  if (t.value !== null && t.value !== undefined && t.format !== "of") {
    if (t.sd !== undefined && t.sd !== null) tile.append(el("div", { class: "ci" }, `± ${fmt(t.sd, 3)} sd over folds`));
    else if (t.ci) tile.append(el("div", { class: "ci" }, `95 % interval [${fmt(t.ci[0], 2)}, ${fmt(t.ci[1], 2)}]`));
  }
  tile.append(el("div", { class: "label" }, t.label), el("div", { class: "sub" }, t.sub || ""),
    el("details", { class: "tile-src" }, [el("summary", {}, "Source"), el("div", {}, t.source || "pending")]));
  return tile;
}

export function callout(kind, title, body) {
  const box = el("div", { class: `callout ${kind}` });
  box.appendChild(el("div", { class: "kind" }, title || kind));
  if (typeof body === "string") box.appendChild(el("p", {}, body));
  else if (Array.isArray(body)) body.forEach((b) => box.appendChild(typeof b === "string" ? el("p", {}, b) : b));
  else if (body) box.appendChild(body);
  return box;
}

export function pendingBlock(what, reason) {
  return el("div", { class: "pending-block" }, [el("div", { class: "kind" }, "Pending"), el("div", {}, [el("b", {}, what + ". "), reason])]);
}

export function tableFrom({ columns, rows, format = {} }) {
  const wrap = el("div", { class: "table-wrap" });
  const table = el("table", { class: "data" });
  const thead = el("thead");
  const tr = el("tr");
  const numeric = columns.map((c) => rows.some((r) => typeof r[c] === "number"));
  columns.forEach((c, i) => tr.appendChild(el("th", { class: numeric[i] ? "num" : "", scope: "col" }, c)));
  thead.appendChild(tr);
  table.appendChild(thead);
  const tbody = el("tbody");
  for (const r of rows) {
    const row = el("tr");
    columns.forEach((c, i) => {
      let v = r[c];
      if (typeof v === "number") v = format[c] ? format[c](v) : (Number.isInteger(v) ? String(v) : v.toFixed(3));
      else if (v === null || v === undefined) v = "";
      else if (v === "inf") v = "+∞";
      row.appendChild(el("td", { class: numeric[i] ? "num" : "" }, String(v)));
    });
    tbody.appendChild(row);
  }
  table.appendChild(tbody);
  wrap.appendChild(table);
  return wrap;
}

export function control(label, input, valueEl) {
  const c = el("label", { class: "control" }, [el("span", {}, label), input]);
  if (valueEl) c.appendChild(valueEl);
  return c;
}

export function select(options, value, onchange) {
  const s = el("select", { onchange: (e) => onchange(e.target.value) });
  for (const o of options) {
    const [v, lab] = Array.isArray(o) ? o : [o, o];
    const opt = el("option", { value: v }, lab);
    if (v === value) opt.selected = true;
    s.appendChild(opt);
  }
  return s;
}

export function segmented(options, value, onchange, ariaLabel) {
  const seg = el("div", { class: "seg", role: "group", "aria-label": ariaLabel || "" });
  const buttons = options.map(([v, lab]) => {
    const b = el("button", { type: "button", "aria-pressed": String(v === value) }, lab);
    b.addEventListener("click", () => {
      buttons.forEach((bb) => bb.setAttribute("aria-pressed", "false"));
      b.setAttribute("aria-pressed", "true");
      onchange(v);
    });
    seg.appendChild(b);
    return b;
  });
  return seg;
}

export function slider(min, max, step, value, onchange, format = (v) => v) {
  const out = el("span", { class: "val" }, format(value));
  const input = el("input", { type: "range", min, max, step, value });
  input.addEventListener("input", () => { out.textContent = format(Number(input.value)); onchange(Number(input.value)); });
  return { input, out };
}

export function sourceList(files) {
  return (files || []).join("; ");
}

export function badge(kind, icon, text) {
  return el("span", { class: `badge ${kind}` }, [el("span", { class: "icon", "aria-hidden": "true" }, icon), text]);
}

// Print: open the collapsed figure notes and tile sources so the paper copy carries them.
window.addEventListener("beforeprint", () => document.querySelectorAll("details.fig-notes, details.tile-src").forEach((d) => { d.open = true; }));
