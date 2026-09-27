import { mountChrome, loadJSON, el, statTile, fmt, pct } from "../site.js";
import { mountOverview } from "../overview.js";
import { figure, line, band, tokens, palette } from "../charts.js";
import { mountLadder, pick } from "../ladder.js";

// Every number on this page is read from site/data/*.json; nothing below is typed by hand.
async function main() {
  await mountChrome("index.html");
  const [H, PC, N, G] = await Promise.all([loadJSON("data/headline.json"), loadJSON("data/panel_curve.json"), loadJSON("data/nesting.json"), loadJSON("data/genes.json")]);
  const ex = H.extras, d = H.design, acc = H.accuracy;
  const tile = Object.fromEntries(H.tiles.map((t) => [t.id, t]));
  const idk = pick(H.ladder, "in_distribution", "k20", "marginal", "pooled");
  const opa = pick(H.ladder, "in_distribution", "k20", "marginal", "one_per_animal");
  const tr = pick(H.ladder, "train_control_test_trained", "k20", "marginal");
  const mf = pick(H.ladder, "train_male_test_female", "k20", "marginal");
  const bm = pick(H.ladder, "different_lab", "k20", "marginal");
  const gt = pick(H.ladder, "different_species", "k20", "marginal");
  const guarantee = pct(1 - d.alpha);                       // "90 %"
  const bridge = ex.bridge_sum_ratio_all_genes_pool99;      // Σ V_batch / Σ V_tissue over all genes, bridging pool 99
  const bmAcc = tile.tile_bodymap_k20.value;
  const kOf = (key) => Number(String(key).replace(/^k/, ""));   // the panel size encoded in an accuracy key ("k50" → 50)
  const perSet = (v) => `${fmt(v, 2)} tissue${fmt(v, 2) === "1.00" ? "" : "s"} per set`;

  // ---- lede ---------------------------------------------------------------------------------------------------
  document.getElementById("lede").replaceChildren(
    el("b", {}, "Yes."),
    ` ${d.k_panel === 20 ? "Twenty" : String(d.k_panel)} genes identify ${ex.n_tissues} rat tissues at ${fmt(tile.tile_acc_k20.value)} balanced accuracy and name ${bmAcc === 1 ? "every mapped organ" : fmt(bmAcc) + " of the mapped organs"} in another laboratory's rats. `,
    `Its ${guarantee} guarantee holds within the study, abstains rather than guesses beyond it, and three animals from the new laboratory restore it.`,
  );

  // ---- tiles --------------------------------------------------------------------------------------------------
  const tiles = document.getElementById("tiles");
  for (const t of H.tiles) tiles.appendChild(statTile(t));

  // ---- the study in one picture (labels ≤ 4 lines of 17 characters, units ≤ 2 lines of 19, so every box reads in full)
  {
    const nodes = [
      { label: `MoTrPAC RNA-seq, ${ex.n_tissues} tissues`, num: String(ex.n_vials), unit: "vials" },
      { label: "Animal-grouped folds, whole animals held out", num: String(d.n_outer_folds), unit: "folds" },
      { label: "Panel selected inside each fold", num: fmt(acc.k20.mean), unit: `balanced accuracy, ${d.k_panel} genes` },
      { label: `${guarantee} sets calibrated on held-out animals`, num: String(d.n_cal_animals), unit: "calibration animals" },
      { label: "Shift ladder: state, sex, lab, species", num: fmt(bm.accuracy), unit: "mapped organs named, other lab" },
      { label: "Batch measured directly on bridging pools", num: pct(bridge, 1), unit: "of the variance separating tissues" },
      { label: "Sufficient, transferable, calibratable", num: null, unit: null, verdict: true },
    ];
    const alt = `The study in one picture: MoTrPAC RNA-seq of ${ex.n_tissues} tissues (${ex.n_vials} vials) → ${d.n_outer_folds} animal-grouped folds with whole animals held out → a ${d.k_panel}-gene panel selected inside each fold (${fmt(acc.k20.mean)} balanced accuracy) → ${guarantee} conformal sets calibrated on ${d.n_cal_animals} held-out animals → a shift ladder (training state, sex, laboratory, species; ${fmt(bm.accuracy)} of mapped organs named in another laboratory) → batch measured directly on bridging pools (${pct(bridge, 1)} of the variance separating tissues) → verdict: sufficient, transferable, calibratable.`;
    mountOverview(document.getElementById("overview-svg"), nodes, { alt, desc: nodes.map((n) => n.label + (n.num ? ` ${n.num} ${n.unit || ""}` : "")).join("; ") });
    document.getElementById("overview-caption").replaceChildren(
      `MoTrPAC's rat endurance-training study makes this possible: ${ex.n_tissues} tissues from the same ${ex.n_animals} animals, ${ex.n_vials} RNA-seq vials, and reference RNA pools on every extraction plate. `,
      `Inside animal-grouped folds we selected the panel, calibrated its ${guarantee} sets on held-out animals, carried both up a ladder of shifts, and measured batch directly on those pools.`,
    );
  }

  // ---- the transfer ladder (hero; shared with the Transfer page) ----------------------------------------------
  await mountLadder(document.getElementById("fig-ladder"), H, { full: false, shortTitle: true });

  // ---- three one-line points under the ladder (full content width, not the 72ch prose measure) -----------------
  const point = (text) => el("li", { style: "max-width: none" }, text);
  document.getElementById("key-points").replaceChildren(
    point(`A class-aware round-robin selector makes ${d.k_panel} genes sufficient (${fmt(acc.k20.mean)} ± ${fmt(acc.k20.sd)} on ${ex.n_tissues} tissues); a univariate F-test reaches ${fmt(ex.acc_fclassif_k20)}.`),
    point(`The ${guarantee} guarantee holds in the study, abstains rather than guesses beyond it; three same-species animals restore it (${fmt(bm.recal_n3)}).`),
    point(`MoTrPAC's design makes the check possible: reference RNA pools on every extraction plate let batch be measured directly, and the answer (${pct(bridge, 1)}) supports the biology reading.`),
  );

  // ---- 1. twenty genes are enough: paragraph + panel curve ----------------------------------------------------
  document.getElementById("p-accuracy").replaceChildren(
    `Under ${d.n_outer_folds} animal-grouped folds the round-robin selector reaches ${fmt(acc.k20.mean)} ± ${fmt(acc.k20.sd)} balanced accuracy at ${d.k_panel} genes, `,
    `${fmt(acc.k50.mean)} at ${kOf("k50")} and ${fmt(acc.full.mean)} with all genes. `,
    `A univariate F-test at the same size reaches ${fmt(ex.acc_fclassif_k20)}: the selector, not the classifier, is the result.`,
  );
  // the smallest panel whose mean is within one fold-sd of the 20-gene value (where the curve levels off)
  const curveAt = (k) => PC.curve.find((r) => r.k === k);
  const kLevel = PC.curve.filter((r) => r.roundrobin_mean >= acc.k20.mean - acc.k20.sd).map((r) => r.k).sort((a, b) => a - b)[0];
  await figure(document.getElementById("fig-curve"), {
    title: `A class-aware selector reaches ${fmt(curveAt(kLevel).roundrobin_mean)} by ${kLevel} genes and ${fmt(acc.k20.mean)} by ${d.k_panel}; a univariate F-test at ${d.k_panel} genes reaches ${fmt(ex.acc_fclassif_k20)}`,
    subtitle: `Mean ± sd balanced accuracy over ${d.n_outer_folds} animal-grouped folds vs panel size (log scale); logreg_l2 on the selected genes.`,
    build: () => {
      const c = PC.curve;
      const ks = c.map((r) => r.k);
      const rr = c.map((r) => r.roundrobin_mean), rs = c.map((r) => r.roundrobin_sd || 0);
      const fc = c.map((r) => r.fclassif_mean), fs = c.map((r) => r.fclassif_sd || 0);
      return {
        traces: [
          band(ks, rr.map((v, i) => v - rs[i]), rr.map((v, i) => v + rs[i]), { slot: 1, legendgroup: "rr" }),
          line(ks, rr, { name: "round-robin selector", slot: 1, legendgroup: "rr", hover: "k = %{x}: %{y:.3f}<extra>round-robin</extra>" }),
          band(ks, fc.map((v, i) => v - fs[i]), fc.map((v, i) => v + fs[i]), { slot: 2, legendgroup: "fc" }),
          line(ks, fc, { name: "F-test selector", slot: 2, legendgroup: "fc", hover: "k = %{x}: %{y:.3f}<extra>F-test</extra>" }),
        ],
        layout: { xaxis: { type: "log", title: { text: "panel size k (genes)" }, tickvals: ks, ticktext: ks.map(String) }, yaxis: { range: [0, 1.05], title: { text: "balanced accuracy" } }, margin: { t: 40 } },
        table: { columns: ["k", "roundrobin_mean", "roundrobin_sd", "fclassif_mean", "fclassif_sd", "n_folds", "n_train_animals", "n_test_animals"], rows: c },
      };
    },
    source: "results/05_panels/TRNSCRPT/panel_curve.csv, results/05_panels/TRNSCRPT/panel_curve_fclassif.csv",
    notShow: "which genes: the selection changes per fold (see the Panel page for the stable core and the Explorer's panel builder).",
  });

  // ---- 2. the guarantee travels honestly: paragraph + empty-set stack ------------------------------------------
  document.getElementById("p-guarantee").replaceChildren(
    `In-distribution the ${d.k_panel}-gene sets cover ${fmt(idk?.coverage)} of held-out vials (${fmt(opa?.coverage)} with one vial per animal), and ${fmt(tr?.coverage)} on `,
    el("a", { href: "exercise.html" }, "trained animals"),
    " with the panel fit on sedentary controls only. ",
    `Beyond the study they abstain rather than guess: coverage ${fmt(bm?.coverage)} in another laboratory and ${fmt(gt?.coverage)} in human, with wrong non-empty sets at ${fmt(bm?.wrong_non_empty)} and ${fmt(gt?.wrong_non_empty)}. `,
    `Three target animals restore ${fmt(bm?.recal_n3)} within species at ${perSet(bm?.recal_n3_size)}.`,
  );
  await figure(document.getElementById("fig-empty"), {
    title: "Under shift the sets abstain rather than guess: the shortfall is empty sets; wrong confident sets stay rare",
    subtitle: `Fraction of test samples whose α = ${fmt(d.alpha, 2)} set (${d.k_panel}-gene panel, marginal, source-calibrated) holds the true tissue, is wrong, or is empty.`,
    build: () => {
      const t = tokens();
      const p = palette();
      const rows = [["In-distribution", idk], ["Held-out sex (M → F)", mf], ["Other laboratory (BodyMap)", bm], ["Other species (GTEx)", gt]].filter(([, r]) => r && !r.pending);
      const y = rows.map(([l]) => l);
      const covered = rows.map(([, r]) => r.coverage);
      const empty = rows.map(([, r]) => (r.empty === null || r.empty === undefined ? null : r.empty));
      const wrong = rows.map(([, r], i) => (r.wrong_non_empty !== null && r.wrong_non_empty !== undefined ? r.wrong_non_empty : (empty[i] === null ? null : Math.max(0, 1 - r.coverage - empty[i]))));
      const mk = (name, vals, color, pattern) => ({ type: "bar", orientation: "h", y, x: vals, name, marker: { color, line: { color: t.surface, width: 2 }, pattern: pattern ? { shape: "/", fgcolor: t.ink2, bgcolor: t.surface, size: 6, solidity: 0.35 } : undefined },
                                                    hovertemplate: "%{y}: %{x:.3f}<extra>" + name + "</extra>", text: vals.map((v) => (v !== null && v >= 0.12 ? fmt(v, 2) : "")), textposition: "inside", insidetextanchor: "middle", textfont: { color: pattern ? t.ink : "#fff" }, cliponaxis: false });
      return {
        traces: [mk("true tissue in the set (covered)", covered, p[0]), mk("non-empty but wrong", wrong, p[7]), mk("empty set (abstains)", empty, t.grid, true)],
        layout: { barmode: "stack", xaxis: { range: [0, 1], title: { text: "fraction of test samples" } }, yaxis: { autorange: "reversed", automargin: true }, margin: { t: 40, l: 10 }, bargap: 0.4, legend: { y: 1.18 } },
        table: { columns: ["shift", "covered", "wrong_non_empty", "empty"], rows: rows.map(([l], i) => ({ shift: l, covered: covered[i], wrong_non_empty: wrong[i], empty: empty[i] })) },
      };
    },
    source: "results/31_site_regen/06_conformal/TRNSCRPT/scores_*.csv (in-distribution k20, recomputed), results/08_shift/TRNSCRPT/shift_table.csv, results/12_bodymap/conformal_transfer.csv, results/13_gtex/conformal_transfer.csv",
    notShow: "Mondrian and floored sets, which trade abstention for larger sets (Transfer page). All three fractions are over the same samples: seen-class vials for the held-out sex, mapped samples for BodyMap and GTEx.",
  });

  // ---- 3. biology, not processing: paragraph + batch measured per panel gene on the bridging pools --------------
  const est = N.estimable_pairs.find((r) => r.assay === "TRNSCRPT");
  document.getElementById("p-identifiability").replaceChildren(
    `As in any multi-tissue design, each tissue was processed as a unit: ${ex.n_plates} extraction plates, ${ex.n_lib_batches} library batches and ${ex.n_flowcells} flowcells hold whole tissues, `,
    `${est.n_pairs_estimable} of ${est.n_pairs_total} tissue pairs ${est.n_pairs_estimable === 1 ? "is" : "are"} contrastable inside one batch, so within-study accuracy needs an outside check. `,
    `Two external facts make the fingerprint credible as biology: ${fmt(bmAcc)} of mapped adult organs named in a laboratory with none of these batches, `,
    `and batch measured directly on MoTrPAC's reference RNA pools: ${pct(bridge, 1)} of the tissue-separating variance on the pool run on ${ex.bridge_n_plates_pool99} plates at both sites, and ${fmt(100 * ex.bridge_pools_min_sum_ratio_all_genes, 1)}–${pct(ex.bridge_pools_max_sum_ratio_all_genes, 1)} across all ${ex.bridge_n_pools} bridging pools.`,
  );
  // one row per bridging reference pool (all genes): Σ V_batch / Σ V_tissue in %, from nesting.json bridge.pools
  const pools = [...(N.bridge.pools || [])].sort((a, b) => a.sum_ratio_batch_over_tissue - b.sum_ratio_batch_over_tissue);
  const poolLabel = (r) => `${r.pool_tissue.replace(" Powder", "").toLowerCase()} pool ${r.pool_type}`;
  await figure(document.getElementById("fig-bridge"), {
    title: "On every bridging pool, batch is a few percent of the tissue signal",
    subtitle: "Between-plate variance of a reference RNA pool over the variance of the 19 tissue means, all genes; the two gastrocnemius pools ran on both sequencing sites.",
    build: () => {
      const t = tokens();
      const p = palette();
      const y = pools.map(poolLabel);
      const x = pools.map((r) => 100 * r.sum_ratio_batch_over_tissue);
      const both = pools.map((r) => r.sites.includes(";"));
      return {
        traces: [
          { type: "bar", orientation: "h", y, x, name: "Σ V_batch / Σ V_tissue", marker: { color: t.grid, line: { color: t.surface, width: 0 } }, width: 0.08, hoverinfo: "skip", showlegend: false },
          { type: "scatter", mode: "markers+text", y, x, name: "Σ V_batch / Σ V_tissue", marker: { color: both.map((b) => (b ? p[0] : p[1])), size: 12, line: { color: t.surface, width: 2 } },
            text: pools.map((r) => `${fmt(100 * r.sum_ratio_batch_over_tissue, 1)} % · ${r.n_plates} plates${r.sites.includes(";") ? ", both sites" : ""}`), textposition: "middle right", textfont: { color: t.ink2, size: 11 }, cliponaxis: false,
            customdata: pools.map((r) => [r.n_plates, r.sites, r.n_genes]), hovertemplate: "%{y}: %{x:.2f} % of the tissue-separating variance<br>%{customdata[0]} plates (%{customdata[1]}), %{customdata[2]} genes<extra></extra>", showlegend: false },
          // legend-only traces: the marker colour says whether the pool crossed both sequencing sites
          { type: "scatter", mode: "markers", x: [null], y: [null], name: "both sequencing sites", marker: { color: p[0], size: 12 }, hoverinfo: "skip" },
          { type: "scatter", mode: "markers", x: [null], y: [null], name: "one site", marker: { color: p[1], size: 12 }, hoverinfo: "skip" },
        ],
        layout: {
          xaxis: { title: { text: "batch variance as % of the tissue-separating variance (all genes)" }, range: [0, 10], ticksuffix: " %" },
          yaxis: { automargin: true, categoryorder: "array", categoryarray: y, autorange: "reversed" },
          annotations: [{ xref: "x", x: 10, xanchor: "right", yref: "paper", y: 0, yanchor: "bottom", text: "tissue signal = 100 %", showarrow: false, font: { color: t.muted, size: 11 } }],
          legend: { orientation: "h", y: 1.02, yanchor: "bottom", x: 0 }, margin: { t: 40, l: 10, r: 20, b: 56 },
        },
        table: { columns: ["pool_tissue", "pool_type", "n_plates", "sites", "n_genes", "sum_ratio_batch_over_tissue"], rows: pools,
                 format: { sum_ratio_batch_over_tissue: (v) => v.toFixed(4) } },
      };
    },
    source: "results/16_identifiability/bridge_variance.csv (gene_set all_genes; one row per reference pool run on more than one extraction plate)",
    notShow: "the per-gene picture, on the Identifiability page: a muscle-derived pool measures batch only on the genes it expresses, so on markers of other tissues the ratio reads zero. The 100 % reference is the variance that separates the 19 tissue means, which every ratio is taken against.",
  });
}

main().catch((e) => { console.error(e); document.getElementById("lede").textContent = "Failed to load site data: " + e.message; });
