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
  await mountLadder(document.getElementById("fig-ladder"), H, { full: false });

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
    `and batch measured at ${pct(bridge, 1)} of tissue-separating variance on bridging pools run on ${ex.bridge_n_plates_pool99} plates at both sites.`,
  );
  // the two bridging pools that crossed every plate: pool_bid 80001 is reference pool 99 (the headline pool), 80000 is pool 88
  const k20 = new Set(G.sets.k20);
  const geneInfo = Object.fromEntries(G.genes.map((g) => [g.id, g]));
  const poolRows = (bid) => N.bridge.per_gene.filter((r) => r.pool_bid === bid && k20.has(r.feature_ID));
  const pools = [poolRows(80001), poolRows(80000)].filter((rows) => rows.length);
  const poolName = (rows) => `pool ${rows[0].pool_type}`;
  const geneLabel = (r) => `${r.gene_symbol || geneInfo[r.feature_ID]?.symbol || r.feature_ID} · ${geneInfo[r.feature_ID]?.marker_tissue || "—"}`;
  // rows top to bottom: by the headline pool's ratio, then the other pool's, then symbol
  const ratioIn = (rows, id) => rows.find((r) => r.feature_ID === id)?.ratio_batch_over_tissue ?? 0;
  const order = [...k20].filter((id) => pools[0].some((r) => r.feature_ID === id))
    .sort((a, b) => ratioIn(pools[0], b) - ratioIn(pools[0], a) || (pools[1] ? ratioIn(pools[1], b) - ratioIn(pools[1], a) : 0) || geneLabel(pools[0].find((r) => r.feature_ID === a)).localeCompare(geneLabel(pools[0].find((r) => r.feature_ID === b))));
  const labels = order.map((id) => geneLabel(pools[0].find((r) => r.feature_ID === id)));
  await figure(document.getElementById("fig-bridge"), {
    title: `Batch measured on MoTrPAC's bridging pools: ${pct(bridge, 1)} of tissue-separating variance`,
    subtitle: "V_batch / V_tissue per panel gene; the all-gene ratio is the headline number.",
    build: () => {
      const t = tokens();
      const p = palette();
      const traces = [];
      pools.forEach((rows, i) => {
        const name = poolName(rows), color = p[i];
        const top = Math.max(...rows.map((r) => r.ratio_batch_over_tissue));
        for (const expressed of [true, false]) {
          const sel = order.map((id) => rows.find((r) => r.feature_ID === id)).filter((r) => r && Boolean(r.expressed_in_pool) === expressed);
          if (!sel.length) continue;
          traces.push({
            type: "bar", orientation: "h", name: `${name}${expressed ? "" : ", not expressed (hollow)"}`,
            y: sel.map(geneLabel), x: sel.map((r) => r.ratio_batch_over_tissue),
            offsetgroup: name, alignmentgroup: "pools", legendgroup: name,
            marker: expressed ? { color, line: { color: t.surface, width: 2 }, cornerradius: 4 } : { color: "rgba(0,0,0,0)", line: { color, width: 1.5 }, cornerradius: 4 },
            text: sel.map((r) => (expressed && r.ratio_batch_over_tissue === top ? fmt(r.ratio_batch_over_tissue, 3) : "")), textposition: "outside", textfont: { color: t.ink2, size: 11 }, cliponaxis: false,
            customdata: sel.map((r) => [r.v_batch, r.v_tissue, r.mean_log2cpm_in_pool, r.n_plates]),
            hovertemplate: "%{y}: V_batch / V_tissue = %{x:.4f}<br>V_batch (across plates) %{customdata[0]:.4f}, V_tissue (across tissue means) %{customdata[1]:.3f}<br>mean log2 CPM in the pool %{customdata[2]:.2f}, %{customdata[3]} plates<extra>" + name + "</extra>",
          });
        }
      });
      return {
        traces,
        layout: {
          barmode: "group", height: 26 * order.length + 130,
          xaxis: { title: { text: "V_batch / V_tissue per gene" }, rangemode: "tozero" },
          yaxis: { autorange: "reversed", automargin: true, categoryorder: "array", categoryarray: labels, tickfont: { size: 11 } },
          shapes: [{ type: "line", xref: "x", x0: bridge, x1: bridge, yref: "paper", y0: 0, y1: 1, line: { color: t.ink2, width: 1, dash: "dash" } }],
          annotations: [{ xref: "x", x: bridge, xanchor: "left", yref: "paper", y: 0.02, yanchor: "bottom", xshift: 6, text: `all genes, ${poolName(pools[0])}: ${pct(bridge, 1)}`, showarrow: false, font: { color: t.ink2, size: 11 } }],
          // the four legend entries wrap to four rows on a phone, so the top margin grows there
          legend: { y: 1.02, yanchor: "bottom" }, margin: { t: window.matchMedia("(max-width: 640px)").matches ? 120 : 64, l: 10, b: 56 },
        },
        table: { columns: ["gene", "marker_tissue", "feature_ID", "pool", "ratio_batch_over_tissue", "v_batch", "v_tissue", "mean_log2cpm_in_pool", "expressed_in_pool", "n_plates"],
                 rows: pools.flatMap((rows) => order.map((id) => rows.find((r) => r.feature_ID === id)).filter(Boolean).map((r) => ({ gene: r.gene_symbol, marker_tissue: geneInfo[r.feature_ID]?.marker_tissue || "", feature_ID: r.feature_ID, pool: `${poolName(rows)} (${r.pool_bid})`, ratio_batch_over_tissue: r.ratio_batch_over_tissue, v_batch: r.v_batch, v_tissue: r.v_tissue, mean_log2cpm_in_pool: r.mean_log2cpm_in_pool, expressed_in_pool: r.expressed_in_pool ? "yes" : "no", n_plates: r.n_plates }))),
                 format: { ratio_batch_over_tissue: (v) => v.toFixed(4), v_batch: (v) => v.toFixed(4), v_tissue: (v) => v.toFixed(3), mean_log2cpm_in_pool: (v) => v.toFixed(2) } },
      };
    },
    source: "results/16_identifiability/bridge_variance_per_gene.csv, bridge_variance.csv",
    notShow: "genes the pool does not express (mean log2 CPM < 1) read zero or near it; the pools are muscle-derived, so batch on other tissues' markers is measurable only where they are expressed.",
    height: "tall",
  });
}

main().catch((e) => { console.error(e); document.getElementById("lede").textContent = "Failed to load site data: " + e.message; });
