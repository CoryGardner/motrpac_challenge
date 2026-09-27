// Chart template and helpers — every chart on the site goes through `figure()`.
// Rules enforced here: one y-axis, tokens from theme.css, thin marks (2px lines, ≥ 8px markers, ≤ 24px bars
// with a 2px surface gap), hairline solid gridlines, hover on every mark, text in ink tokens, a legend for
// ≥ 2 series, a data-table toggle on every figure, and a caption with "Source:" and "What it does not show:".

import { el, tableFrom } from "./site.js";

export const PLOTLY_VERSION = "2.35.2";
const registry = new Set();   // live figures, re-rendered on theme change

function cssVar(name) {
  return getComputedStyle(document.documentElement).getPropertyValue(name).trim();
}

export function tokens() {
  return {
    page: cssVar("--page"), surface: cssVar("--surface"), ink: cssVar("--ink"), ink2: cssVar("--ink-2"), muted: cssVar("--muted"),
    grid: cssVar("--grid"), axis: cssVar("--axis"), font: cssVar("--font"),
    good: cssVar("--good"), warning: cssVar("--warning"), serious: cssVar("--serious"), critical: cssVar("--critical"),
    divNeg: cssVar("--div-neg"), divMid: cssVar("--div-mid"), divPos: cssVar("--div-pos"),
  };
}

export function palette() {
  return [1, 2, 3, 4, 5, 6, 7, 8].map((i) => cssVar(`--c${i}`));
}

export function sequential() {
  return [1, 2, 3, 4, 5, 6, 7].map((i) => cssVar(`--s${i}`));
}

export function sequentialScale() {
  const s = sequential();
  const t = tokens();
  const n = s.length;
  // the lightest step sits just above the surface so zero cells read as "nothing"
  return [[0, t.surface], ...s.map((c, i) => [(i + 1) / n, c])];
}

export function divergingScale() {
  const t = tokens();
  return [[0, t.divNeg], [0.5, t.divMid], [1, t.divPos]];
}

// organ systems, palette order fixed by the design spec
export const ORGAN_SYSTEMS = [
  { name: "brain", tissues: ["CORTEX", "HIPPOC", "HYPOTH"] },
  { name: "muscle", tissues: ["SKM-GN", "SKM-VL", "HEART"] },
  { name: "adipose", tissues: ["WAT-SC", "BAT"] },
  { name: "gut", tissues: ["COLON", "SMLINT"] },
  { name: "gonad", tissues: ["OVARY", "TESTES"] },
  { name: "circulation & immune", tissues: ["BLOOD", "SPLEEN", "VENACV"] },
  { name: "visceral", tissues: ["LIVER", "KIDNEY", "LUNG"] },
  { name: "endocrine", tissues: ["ADRNL"] },
];
export const TISSUE_NAMES = {
  ADRNL: "adrenal gland", BAT: "brown adipose", BLOOD: "whole blood", COLON: "colon", CORTEX: "cerebral cortex", HEART: "heart",
  HIPPOC: "hippocampus", HYPOTH: "hypothalamus", KIDNEY: "kidney", LIVER: "liver", LUNG: "lung", OVARY: "ovary", PLASMA: "plasma",
  "SKM-GN": "gastrocnemius", "SKM-VL": "vastus lateralis", SMLINT: "small intestine", SPLEEN: "spleen", TESTES: "testes",
  VENACV: "vena cava", "WAT-SC": "white adipose",
};

export function organSystem(tissue) {
  const p = palette();
  for (let i = 0; i < ORGAN_SYSTEMS.length; i++) {
    if (ORGAN_SYSTEMS[i].tissues.includes(tissue)) return { index: i + 1, name: ORGAN_SYSTEMS[i].name, color: p[i] };
  }
  return { index: 0, name: "other", color: tokens().muted };
}

export function tissueLabel(t) {
  return TISSUE_NAMES[t] ? `${t} · ${TISSUE_NAMES[t]}` : t;
}

// Plotly layout template built from the live tokens
export function template(extra = {}) {
  const t = tokens();
  const axis = {
    gridcolor: t.grid, gridwidth: 1, zerolinecolor: t.axis, zerolinewidth: 1, linecolor: t.axis, linewidth: 1,
    tickfont: { color: t.ink2, size: 12 }, titlefont: { color: t.ink2, size: 12 }, tickcolor: t.axis, automargin: true,
    showline: false, ticks: "", fixedrange: false,
  };
  return {
    paper_bgcolor: t.surface, plot_bgcolor: t.surface,
    font: { family: t.font, color: t.ink2, size: 12 },
    colorway: palette(),
    margin: { l: 56, r: 20, t: 16, b: 48 },
    xaxis: { ...axis }, yaxis: { ...axis },
    legend: { orientation: "h", x: 0, y: 1.12, yanchor: "bottom", font: { color: t.ink2, size: 12 }, bgcolor: "rgba(0,0,0,0)" },
    hoverlabel: { bgcolor: t.page, bordercolor: t.axis, font: { family: t.font, color: t.ink, size: 12 } },
    hovermode: "closest",
    bargap: 0.35, bargroupgap: 0.08,
    ...extra,
  };
}

export const CONFIG = { responsive: true, displaylogo: false, displayModeBar: false, scrollZoom: false };

function plotlyReady(timeoutMs = 20000) {
  return new Promise((resolve, reject) => {
    const t0 = Date.now();
    (function poll() {
      if (window.Plotly) return resolve(window.Plotly);
      if (Date.now() - t0 > timeoutMs) return reject(new Error("Plotly did not load (CDN and vendored fallback)"));
      setTimeout(poll, 40);
    })();
  });
}

/**
 * figure(container, spec)
 *  spec = { title, subtitle, traces | build(), layout, source, notShow, table: {columns, rows} | null, height, id }
 *  `build()` (optional) returns {traces, layout} and is re-run on theme change so colours follow the tokens.
 */
export async function figure(container, spec) {
  const root = typeof container === "string" ? document.getElementById(container) : container;
  root.classList.add("figure");
  root.replaceChildren();
  const head = el("div", { class: "fig-head" }, [el("h3", {}, spec.title || "")]);
  if (spec.subtitle) head.appendChild(el("p", { class: "fig-sub" }, spec.subtitle));
  root.appendChild(head);
  const toolbar = el("div", { class: "fig-toolbar" });
  root.appendChild(toolbar);
  if (spec.toolbar) toolbar.appendChild(spec.toolbar);
  const chart = el("div", { class: "chart" + (spec.height === "short" ? " short" : spec.height === "tall" ? " tall" : ""), role: "img", "aria-label": spec.alt || spec.title || "chart" });
  root.appendChild(chart);
  const tableWrap = el("div", { class: "fig-table", hidden: "" });
  const btn = el("button", { class: "btn", type: "button", "aria-pressed": "false", "aria-controls": tableWrap.id = `tbl-${Math.random().toString(36).slice(2, 8)}` }, "Show data table");
  toolbar.appendChild(btn);
  root.appendChild(tableWrap);
  const cap = el("p", { class: "fig-caption" });
  cap.appendChild(el("b", {}, "Source: "));
  cap.appendChild(document.createTextNode(spec.source || "pending"));
  if (spec.notShow) {
    cap.appendChild(el("br"));
    cap.appendChild(el("b", {}, "What it does not show: "));
    cap.appendChild(document.createTextNode(spec.notShow));
  }
  root.appendChild(cap);

  const state = { root, chart, spec, tableWrap };
  btn.addEventListener("click", () => {
    const open = btn.getAttribute("aria-pressed") !== "true";
    btn.setAttribute("aria-pressed", String(open));
    btn.textContent = open ? "Hide data table" : "Show data table";
    if (open) {
      tableWrap.hidden = false;
      tableWrap.replaceChildren(tableFrom(state.table || tableFromTraces(state.traces)));
    } else {
      tableWrap.hidden = true;
    }
  });
  registry.add(state);
  await render(state);
  return state;
}

function tableFromTraces(traces) {
  const rows = [];
  for (const tr of traces || []) {
    const name = tr.name || tr.type;
    if (tr.type === "heatmap") {
      (tr.y || []).forEach((yy, i) => (tr.x || []).forEach((xx, j) => rows.push({ series: name, row: yy, column: xx, value: tr.z[i][j] })));
    } else {
      const xs = tr.x || [], ys = tr.y || [];
      const n = Math.max(xs.length, ys.length);
      for (let i = 0; i < n; i++) rows.push({ series: name, x: xs[i], y: ys[i], ...(tr.error_y && tr.error_y.array ? { sd: tr.error_y.array[i] } : {}) });
    }
  }
  const columns = rows.length ? Object.keys(rows[0]) : ["series"];
  return { columns, rows };
}

async function render(state) {
  const Plotly = await plotlyReady();
  const built = state.spec.build ? state.spec.build() : { traces: state.spec.traces, layout: state.spec.layout || {} };
  state.traces = built.traces;
  state.table = built.table || state.spec.table || null;
  const layout = mergeLayout(template(), built.layout || {});
  await Plotly.react(state.chart, built.traces, layout, CONFIG);
  if (!state.tableWrap.hidden) state.tableWrap.replaceChildren(tableFrom(state.table || tableFromTraces(state.traces)));
}

function mergeLayout(base, extra) {
  const out = { ...base, ...extra };
  for (const k of Object.keys(extra)) {
    if (extra[k] && typeof extra[k] === "object" && !Array.isArray(extra[k]) && base[k] && typeof base[k] === "object") out[k] = { ...base[k], ...extra[k] };
  }
  // any secondary axis (xaxis2/yaxis2, used only for small multiples) inherits the axis styling
  for (const k of Object.keys(extra)) {
    if (/^[xy]axis\d+$/.test(k)) out[k] = { ...base[k[0] + "axis"], ...extra[k] };
  }
  return out;
}

export async function rerenderAll() {
  for (const s of registry) {
    try { await render(s); } catch (e) { console.error(e); }
  }
}
document.addEventListener("themechange", () => { rerenderAll(); });

// ---- mark helpers ------------------------------------------------------------------------------
export function bar(x, y, opts = {}) {
  const t = tokens();
  const p = palette();
  return {
    type: "bar", x, y, name: opts.name, orientation: opts.horizontal ? "h" : "v",
    marker: { color: opts.color || p[(opts.slot || 1) - 1], line: { color: t.surface, width: 2 }, cornerradius: 4, ...(opts.marker || {}) },
    width: opts.width, text: opts.text, textposition: opts.text ? (opts.textposition || "outside") : undefined, textfont: { color: t.ink2 },
    cliponaxis: false,
    hovertemplate: opts.hover || (opts.horizontal ? "%{y}: %{x}<extra>%{fullData.name}</extra>" : "%{x}: %{y}<extra>%{fullData.name}</extra>"),
    error_y: opts.sd ? { type: "data", array: opts.sd, visible: true, color: t.ink2, thickness: 1.5, width: 4 } : undefined,
    error_x: opts.sdx ? { type: "data", array: opts.sdx, visible: true, color: t.ink2, thickness: 1.5, width: 4 } : undefined,
    customdata: opts.customdata, showlegend: opts.showlegend,
  };
}

export function line(x, y, opts = {}) {
  const t = tokens();
  const p = palette();
  const color = opts.color || p[(opts.slot || 1) - 1];
  return {
    type: "scatter", mode: opts.mode || "lines+markers", x, y, name: opts.name,
    line: { color, width: 2, dash: opts.dash || "solid", shape: opts.shape || "linear" },
    marker: { color, size: opts.size || 8, line: { color: t.surface, width: 2 }, symbol: opts.symbol || "circle" },
    hovertemplate: opts.hover || "%{x}: %{y:.3f}<extra>%{fullData.name}</extra>",
    text: opts.text, textposition: opts.textposition || "top center", textfont: { color: t.ink2, size: 11 },
    customdata: opts.customdata, showlegend: opts.showlegend, legendgroup: opts.legendgroup,
    error_y: opts.sd ? { type: "data", array: opts.sd, visible: true, color, thickness: 1.5, width: 4 } : undefined,
    fill: opts.fill, fillcolor: opts.fillcolor,
  };
}

export function band(x, lo, hi, opts = {}) {
  // a translucent band between lo and hi (mean ± sd), in the series hue at 12 % opacity
  const p = palette();
  const color = opts.color || p[(opts.slot || 1) - 1];
  return {
    type: "scatter", x: [...x, ...[...x].reverse()], y: [...hi, ...[...lo].reverse()], fill: "toself", fillcolor: hexAlpha(color, 0.12),
    line: { color: "rgba(0,0,0,0)", width: 0 }, hoverinfo: "skip", showlegend: false, name: (opts.name || "") + " band", legendgroup: opts.legendgroup,
  };
}

export function heatmap(z, x, y, opts = {}) {
  const t = tokens();
  return {
    type: "heatmap", z, x, y, colorscale: opts.diverging ? divergingScale() : sequentialScale(), zmin: opts.zmin, zmax: opts.zmax, zmid: opts.zmid,
    xgap: 2, ygap: 2, hoverongaps: false, showscale: opts.showscale !== false,
    colorbar: { thickness: 10, len: 0.8, outlinewidth: 0, tickfont: { color: t.ink2, size: 11 }, title: { text: opts.ztitle || "", font: { color: t.ink2, size: 11 } } },
    text: opts.text, texttemplate: opts.texttemplate, textfont: { color: t.ink, size: 10 },
    hovertemplate: opts.hover || "%{y} → %{x}: %{z}<extra></extra>",
  };
}

export function strip(xs, ys, opts = {}) {
  const t = tokens();
  const p = palette();
  return {
    type: "scatter", mode: "markers", x: xs, y: ys, name: opts.name,
    marker: { color: opts.colors || opts.color || p[(opts.slot || 1) - 1], size: opts.size || 8, symbol: opts.symbols || opts.symbol || "circle", opacity: opts.opacity || 0.75,
              line: { color: t.surface, width: 1.5 } },
    hovertemplate: opts.hover || "%{x}: %{y:.2f}<extra></extra>", customdata: opts.customdata, showlegend: opts.showlegend !== false, legendgroup: opts.legendgroup,
  };
}

export function hexAlpha(hex, a) {
  const h = hex.replace("#", "");
  const r = parseInt(h.slice(0, 2), 16), g = parseInt(h.slice(2, 4), 16), b = parseInt(h.slice(4, 6), 16);
  return `rgba(${r},${g},${b},${a})`;
}

// a dashed reference line (e.g. 1 − α = 0.90) as a layout shape + annotation
export function refLine(y, label, opts = {}) {
  const t = tokens();
  return {
    shapes: [{ type: "line", xref: "paper", x0: 0, x1: 1, yref: "y", y0: y, y1: y, line: { color: t.ink2, width: 1, dash: "dash" } }],
    annotations: label ? [{ xref: "paper", x: opts.x ?? 1, xanchor: "right", yref: "y", y, yanchor: "bottom", text: label, showarrow: false, font: { color: t.ink2, size: 11 } }] : [],
  };
}

// jitter for strip plots, deterministic
export function jitter(i, n, width = 0.32) {
  const g = ((i * 9301 + 49297) % 233280) / 233280;
  return (g - 0.5) * 2 * width;
}

export function sd(arr) {
  const v = arr.filter((x) => Number.isFinite(x));
  if (v.length < 2) return 0;
  const m = v.reduce((a, b) => a + b, 0) / v.length;
  return Math.sqrt(v.reduce((a, b) => a + (b - m) ** 2, 0) / (v.length - 1));
}
export function mean(arr) {
  const v = arr.filter((x) => Number.isFinite(x));
  return v.length ? v.reduce((a, b) => a + b, 0) / v.length : NaN;
}
export function median(arr) {
  const v = arr.filter((x) => Number.isFinite(x)).sort((a, b) => a - b);
  if (!v.length) return NaN;
  const m = Math.floor(v.length / 2);
  return v.length % 2 ? v[m] : (v[m - 1] + v[m]) / 2;
}
