// Data-driven SVG diagrams: "the study in one picture" (home) and the in-fold pipeline (methods).
// Every label and number arrives in `slots` from the page script (read from site/data/*.json); the builders
// never contain a result number. Theme comes from CSS tokens (theme.css .svg-diagram rules); the layout switches
// to a stacked column below 640 px and re-renders when the viewport crosses that width.
const NS = "http://www.w3.org/2000/svg";
const NARROW = "(max-width: 640px)";

function svg(tag, attrs = {}, children = []) {
  const n = document.createElementNS(NS, tag);
  for (const [k, v] of Object.entries(attrs)) if (v !== undefined && v !== null) n.setAttribute(k, String(v));
  for (const c of children) n.appendChild(typeof c === "string" ? document.createTextNode(c) : c);
  return n;
}

function text(x, y, str, attrs = {}) {
  return svg("text", { x, y, ...attrs }, [str]);
}

function arrowDefs(prefix) {
  return svg("defs", {}, [svg("marker", { id: `${prefix}-arr`, viewBox: "0 0 10 10", refX: 9, refY: 5, markerWidth: 7, markerHeight: 7, orient: "auto-start-reverse" },
                                [svg("path", { d: "M0 0 L10 5 L0 10 z", fill: "currentColor" })])]);
}

function arrow(prefix, x1, y1, x2, y2) {
  return svg("path", { class: "arrow", d: `M${x1} ${y1} L${x2} ${y2}`, "stroke-width": 2, fill: "none", "marker-end": `url(#${prefix}-arr)`, style: "color: var(--ink-2); stroke: currentColor" });
}

/** Wrap a string into lines of at most `max` characters (word boundaries). */
export function wrap(str, max) {
  const out = [];
  let cur = "";
  for (const w of String(str).split(/\s+/)) {
    if ((cur + " " + w).trim().length > max && cur) { out.push(cur); cur = w; } else cur = (cur + " " + w).trim();
  }
  if (cur) out.push(cur);
  return out;
}

/**
 * The study in one picture. `nodes`: [{ label: string (≤ 6 words), num: string | null, unit: string | null, verdict?: boolean }].
 * Returns the <svg>; re-renders on viewport changes. `alt` becomes the <title>; `desc` the <desc>.
 */
export function mountOverview(container, nodes, { alt, desc, prefix = "ov" } = {}) {
  const build = () => {
    const narrow = window.matchMedia(NARROW).matches;
    const n = nodes.length;
    const root = svg("svg", { class: "svg-diagram wide overview-svg", role: "img", "aria-labelledby": `${prefix}-title`, "aria-describedby": `${prefix}-desc` });
    root.appendChild(svg("title", { id: `${prefix}-title` }, [alt || ""]));
    root.appendChild(svg("desc", { id: `${prefix}-desc` }, [desc || ""]));
    root.appendChild(arrowDefs(prefix));
    if (!narrow) {
      // every box reads in full: the label wraps at 17 characters (up to 4 lines), the unit at 19 (up to 2), and all
      // boxes take the height of the tallest content
      const W = 960, w = 118, gap = (W - 2 * 12 - n * w) / (n - 1), y0 = 12;
      const parts = nodes.map((nd) => ({ lines: wrap(nd.label, 17).slice(0, 4), units: nd.num !== null && nd.num !== undefined && nd.unit ? wrap(nd.unit, 19).slice(0, 2) : [] }));
      const labelH = Math.max(...parts.map((p) => p.lines.length)) * 15;
      const unitH = Math.max(...parts.map((p) => p.units.length)) * 12;
      const hasNum = nodes.some((nd) => nd.num !== null && nd.num !== undefined);
      const h = 18 + labelH + (hasNum ? 30 : 0) + unitH + 10;
      root.setAttribute("viewBox", `0 0 ${W} ${h + 2 * y0}`);
      nodes.forEach((nd, i) => {
        const x = 12 + i * (w + gap);
        root.appendChild(svg("rect", { class: nd.verdict ? "node verdict" : "node", x, y: y0, width: w, height: h, rx: 8 }));
        parts[i].lines.forEach((l, j) => root.appendChild(text(x + w / 2, y0 + 20 + j * 15, l, { "text-anchor": "middle", "font-size": 12 })));
        if (nd.num !== null && nd.num !== undefined) {
          const yNum = y0 + 18 + labelH + 22;
          root.appendChild(text(x + w / 2, yNum, nd.num, { class: "num", "text-anchor": "middle", "font-size": 18 }));
          parts[i].units.forEach((u, j) => root.appendChild(text(x + w / 2, yNum + 15 + j * 12, u, { class: "note", "text-anchor": "middle", "font-size": 11 })));
        }
        if (i < n - 1) root.appendChild(arrow(prefix, x + w + 2, y0 + h / 2, x + w + gap - 3, y0 + h / 2));
      });
    } else {
      const W = 360, w = 320, h = 64, step = 82, H = n * step;
      root.setAttribute("viewBox", `0 0 ${W} ${H}`);
      nodes.forEach((nd, i) => {
        const y = 4 + i * step, x = 20;
        root.appendChild(svg("rect", { class: nd.verdict ? "node verdict" : "node", x, y, width: w, height: h, rx: 8 }));
        const lines = wrap(nd.label, nd.num ? 22 : 34).slice(0, 3);
        lines.forEach((l, j) => root.appendChild(text(x + 12, y + (lines.length > 2 ? 18 : 22) + j * (lines.length > 2 ? 14 : 15), l, { "font-size": lines.length > 2 ? 12 : 13 })));
        if (nd.num !== null && nd.num !== undefined) {
          root.appendChild(text(x + w - 12, y + 26, nd.num, { class: "num", "text-anchor": "end", "font-size": 17 }));
          if (nd.unit) wrap(nd.unit, 24).slice(0, 2).forEach((u, j) => root.appendChild(text(x + w - 12, y + 42 + j * 12, u, { class: "note", "text-anchor": "end", "font-size": 10.5 })));
        }
        if (i < n - 1) root.appendChild(arrow(prefix, x + w / 2, y + h + 2, x + w / 2, y + step - 4));
      });
    }
    return root;
  };
  let current = build();
  container.replaceChildren(current);
  const mq = window.matchMedia(NARROW);
  mq.addEventListener("change", () => { const next = build(); current.replaceWith(next); current = next; });
  return current;
}

/**
 * The in-fold pipeline (Methods). `s`: { nTrain, nTest, nTotal, prefilter, cGrid: [..], nFit, nCal, nFolds }.
 */
export function mountPipeline(container, s, { prefix = "pl" } = {}) {
  const steps = [
    { label: ["median", "imputation"], note: ["fold medians"] },
    { label: ["variance", `prefilter (${s.prefilter.toLocaleString("en-US")})`], note: [`top ${s.prefilter.toLocaleString("en-US")} genes by variance`] },
    { label: ["z-score", "per gene"], note: ["fold mean and sd"] },
    { label: ["round-robin", "selector (k genes)"], note: ["one tissue per pick,", "classes in turn"] },
    { label: ["logistic", "regression (L2)"], note: [`C from {${s.cGrid.join(", ")}} by an`, "inner animal-grouped grid search"] },
  ];
  const alt = `Pipeline diagram: median imputation, variance prefilter (${s.prefilter}), z-score, round-robin selector and an L2 logistic regression, all fit inside the training animals of the fold (${s.nTrain} of ${s.nTotal}); the ${s.nTest} held-out animals are scored only. For the conformal guarantee the training animals are split again into ${s.nFit} fit and ${s.nCal} calibration animals.`;
  const build = () => {
    const narrow = window.matchMedia(NARROW).matches;
    const root = svg("svg", { class: "svg-diagram wide pipeline-svg", role: "img", "aria-labelledby": `${prefix}-title` });
    root.appendChild(svg("title", { id: `${prefix}-title` }, [alt]));
    root.appendChild(arrowDefs(prefix));
    if (!narrow) {
      root.setAttribute("viewBox", "0 0 960 300");
      root.appendChild(text(24, 24, `training animals of the fold (${s.nTrain} of ${s.nTotal}): every step below is fit here only`, { "font-size": 14, "font-weight": 600 }));
      root.appendChild(svg("rect", { class: "fold", x: 12, y: 36, width: 772, height: 196, rx: 10, fill: "none" }));
      const xs = [28, 174, 320, 466, 612];
      steps.forEach((st, i) => {
        const x = xs[i];
        root.appendChild(svg("rect", { class: "box", x, y: 70, width: 128, height: 64, rx: 8 }));
        st.label.forEach((l, j) => root.appendChild(text(x + 64, 96 + j * 18, l, { "text-anchor": "middle", "font-size": 13 })));
        st.note.forEach((l, j) => root.appendChild(text(x + 64, 156 + j * 16, l, { class: "note", "text-anchor": "middle", "font-size": 12 })));
        if (i < steps.length - 1) root.appendChild(arrow(prefix, x + 130, 102, xs[i + 1] - 4, 102));
      });
      root.appendChild(text(28, 214, `inner split for the guarantee: ${s.nFit} fit animals, ${s.nCal} calibration animals (one vial each)`, { class: "note", "font-size": 12 }));
      root.appendChild(svg("rect", { class: "box", x: 812, y: 70, width: 136, height: 64, rx: 8 }));
      root.appendChild(text(880, 96, "held-out animals", { "text-anchor": "middle", "font-size": 13 }));
      root.appendChild(text(880, 114, `(${s.nTest}): scored only`, { "text-anchor": "middle", "font-size": 13 }));
      root.appendChild(arrow(prefix, 742, 102, 808, 102));
      root.appendChild(text(880, 156, "→ balanced accuracy,", { class: "note", "text-anchor": "middle", "font-size": 12 }));
      root.appendChild(text(880, 172, "prediction sets, coverage", { class: "note", "text-anchor": "middle", "font-size": 12 }));
      root.appendChild(text(12, 268, `The ${s.nFolds} outer folds partition the ${s.nTotal} animals; the panel is re-selected inside each fold.`, { class: "note", "font-size": 12 }));
    } else {
      const items = [...steps.map((st) => ({ label: st.label.join(" "), note: st.note.join(" "), inFold: true })),
                     { label: `held-out animals (${s.nTest}): scored only`, note: "→ balanced accuracy, prediction sets, coverage", inFold: false }];
      const step = 78, H = 30 + items.length * step + 40;
      root.setAttribute("viewBox", `0 0 360 ${H}`);
      root.appendChild(text(12, 18, `training animals of the fold (${s.nTrain} of ${s.nTotal})`, { "font-size": 13, "font-weight": 600 }));
      root.appendChild(svg("rect", { class: "fold", x: 8, y: 26, width: 344, height: steps.length * step + 4, rx: 10, fill: "none" }));
      items.forEach((it, i) => {
        const y = 34 + i * step;
        root.appendChild(svg("rect", { class: "box", x: 20, y, width: 320, height: 40, rx: 8 }));
        root.appendChild(text(180, y + 25, it.label, { "text-anchor": "middle", "font-size": 13 }));
        wrap(it.note, 52).slice(0, 1).forEach((l) => root.appendChild(text(180, y + 56, l, { class: "note", "text-anchor": "middle", "font-size": 11 })));
        if (i < items.length - 1) root.appendChild(arrow(prefix, 180, y + 62, 180, y + step - 2));
      });
      root.appendChild(text(12, H - 8, `inner split for the guarantee: ${s.nFit} fit / ${s.nCal} calibration animals`, { class: "note", "font-size": 11 }));
    }
    return root;
  };
  let current = build();
  container.replaceChildren(current);
  const mq = window.matchMedia(NARROW);
  mq.addEventListener("change", () => { const next = build(); current.replaceWith(next); current = next; });
  return current;
}
