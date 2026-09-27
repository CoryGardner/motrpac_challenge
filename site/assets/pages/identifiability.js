import { mountChrome, loadJSON, el, fmt, tableFrom, callout, pendingBlock, statTile } from "../site.js";
import { figure, bar, heatmap, strip, tokens, palette, tissueLabel, hexAlpha } from "../charts.js";
import { pick } from "../ladder.js";

const LAYER_LABEL = { TRNSCRPT: "RNA-seq", METHYL: "RRBS methylation", ATAC: "ATAC-seq", PROT: "proteomics (TMT)", PHOSPHO: "phosphoproteomics (TMT)", ACETYL: "acetylproteomics (TMT)", UBIQ: "ubiquitylome (TMT)", IMMUNO: "immunoassays (Luminex)", METAB: "metabolomics" };
const VAR_LABEL = { GET_site: "site", RNA_extr_plate_ID: "extraction plate", RNA_extr_date: "extraction date", DNA_extr_plate_ID: "extraction plate", DNA_extr_date: "extraction date", Lib_prep_date: "library prep date", Lib_batch_ID: "library batch", Seq_date: "sequencing date", Seq_flowcell_ID: "flowcell", Seq_flowcell_lane: "lane", Seq_batch: "sequencing batch", Sample_batch: "sample batch", Nuclei_extr_date: "nuclei extraction date", Tagmentation_date: "tagmentation date", PCR_date: "PCR date", plex_id: "TMT plex (tissue × label)", tmt11_channel: "TMT channel", plate_id: "assay plate", panel_name: "assay panel" };
const pairsText = (s) => s.split(";").map((p) => p.split("|").map((t) => t.toLowerCase()).join(" and ")).join("; ");

async function main() {
  await mountChrome("identifiability.html");
  const [H, N, Q, E, PC, PR] = await Promise.all([loadJSON("data/headline.json"), loadJSON("data/nesting.json"), loadJSON("data/qc_baseline.json"), loadJSON("data/eda.json"), loadJSON("data/panel_curve.json"), loadJSON("data/product.json")]);
  const RB3 = PR.recal_draws["bodymap.3"];
  const ex = H.extras;
  document.getElementById("tiles").replaceChildren(...(H.tiles_identifiability || []).map((t) => statTile(t)));
  const est = Object.fromEntries(N.estimable_pairs.map((r) => [r.assay, r]));
  const rna = est.TRNSCRPT;
  const immunoMax = Math.max(...N.nesting.IMMUNO.filter((r) => r.variable === "plate_id").map((r) => r.max_tissues_per_level));
  document.getElementById("p-nesting").replaceChildren(
    `In RNA-seq, ${ex.n_plates} extraction plates, ${ex.n_lib_batches} library batches and ${ex.n_flowcells} flowcells each hold whole tissues, and a pair can be contrasted inside a batch only when both tissues share a level of all three: `,
    `${rna.n_pairs_estimable} of ${rna.n_pairs_total} RNA-seq pairs ${rna.n_pairs_estimable === 1 ? "does" : "do"} (${pairsText(rna.estimable_pairs)}, also the sex contrast, since each is single-sex), ${est.METHYL.n_pairs_estimable} of ${est.METHYL.n_pairs_total} in methylation, ${est.ATAC.n_pairs_estimable} of ${est.ATAC.n_pairs_total} in ATAC-seq, `,
    `${est.PROT.n_pairs_estimable} of ${est.PROT.n_pairs_total} in TMT proteomics (a plex is ten samples plus one tissue's reference pool: nested by construction) and ${est.IMMUNO.n_pairs_estimable} of ${est.IMMUNO.n_pairs_total} in the immunoassays, whose Luminex plates hold up to ${immunoMax} tissues. `,
    "The metabolomics tables carry no batch variable.",
  );
  // nesting heatmap: rows layer · variable, columns: Cramér's V, fraction of pairs sharing a level, levels per tissue, tissues per level
  const rows = [];
  for (const [assay, tab] of Object.entries(N.nesting)) for (const r of tab) rows.push({ assay, ...r });
  // every processing variable (one heatmap row each) sits behind a toggle; the estimable-pairs chart below carries the section
  const nestingFig = document.getElementById("fig-nesting");
  const nestingSpec = {
    title: "Each processing variable tracks tissue almost one-to-one, as in any multi-tissue design; the immunoassay plates mix tissues",
    subtitle: "Per layer and processing variable: Cramér's V with tissue (1 = determined), fraction of tissue pairs sharing a level, median levels per tissue, most tissues in one level.",
    build: () => {
      const y = rows.map((r) => `${LAYER_LABEL[r.assay] || r.assay} · ${VAR_LABEL[r.variable] || r.variable}`);
      const cols = ["Cramér's V", "pairs sharing a level", "levels per tissue (median, scaled)", "tissues per level (max, scaled)"];
      const z = rows.map((r) => [r.cramers_v, r.n_pairs_sharing_level / r.n_pairs_total, Math.min(1, (r.median_levels_per_tissue - 1) / 5), Math.min(1, (r.max_tissues_per_level - 1) / (r.n_tissues - 1 || 1))]);
      const text = rows.map((r) => [r.cramers_v.toFixed(2), `${r.n_pairs_sharing_level}/${r.n_pairs_total}`, String(r.median_levels_per_tissue), String(r.max_tissues_per_level)]);
      return { traces: [heatmap(z, cols, y, { zmin: 0, zmax: 1, showscale: false, text, texttemplate: "%{text}", hover: "%{y}<br>%{x}: %{text}<extra></extra>" })],
               layout: { height: 22 * rows.length + 110, yaxis: { autorange: "reversed", automargin: true, tickfont: { size: 10 } }, xaxis: { side: "top", tickfont: { size: 11 } }, margin: { t: 60, l: 10, b: 20 } },
               table: { columns: ["assay", "variable", "n_levels", "n_tissues", "median_levels_per_tissue", "max_levels_per_tissue", "max_tissues_per_level", "n_levels_shared", "tissues_in_one_level", "cramers_v", "n_pairs_sharing_level", "n_pairs_total"], rows } };
    },
    source: "results/16_identifiability/nesting_<ASSAY>.csv (recomputed from data/raw/meta/*.csv, study vials only)",
    notShow: "metabolomics (its tables carry no batch variable) or the reference-standard vials; colour is a value scale per column, the printed numbers are the data.", height: "tall",
  };
  let nestingState = null;
  const nestingBtn = el("button", { class: "btn", type: "button", "aria-pressed": "false", "aria-controls": "fig-nesting" }, "Show every processing variable");
  nestingBtn.addEventListener("click", async () => {
    const show = nestingBtn.getAttribute("aria-pressed") !== "true";
    nestingBtn.setAttribute("aria-pressed", String(show));
    nestingBtn.textContent = show ? "Hide the processing variables" : "Show every processing variable";
    nestingFig.hidden = !show;
    if (!show) return;
    if (!nestingState) nestingState = await figure(nestingFig, nestingSpec);   // first click: render, now that the container has a width
    if (window.Plotly) await window.Plotly.Plots.resize(nestingState.chart);   // a chart laid out while hidden (theme re-render) has none
  });
  document.getElementById("nesting-toggle").replaceChildren(nestingBtn);
  await figure(document.getElementById("fig-estimable"), {
    title: `Tissue pairs contrastable inside a batch, per layer: ${est.TRNSCRPT.n_pairs_estimable} in RNA-seq (${pairsText(est.TRNSCRPT.estimable_pairs)}), ${est.METHYL.n_pairs_estimable + est.ATAC.n_pairs_estimable + est.PROT.n_pairs_estimable + est.PHOSPHO.n_pairs_estimable} in the epigenome and TMT layers, ${est.IMMUNO.n_pairs_estimable} in the immunoassays`,
    subtitle: "Tissue pairs sharing a level of every processing variable of the layer (variables in the hover), out of all pairs.",
    build: () => {
      const order = ["TRNSCRPT", "METHYL", "ATAC", "PROT", "PHOSPHO", "ACETYL", "UBIQ", "IMMUNO"].filter((a) => est[a]);
      const y = order.map((a) => LAYER_LABEL[a]);
      const frac = order.map((a) => est[a].n_pairs_estimable / est[a].n_pairs_total);
      return { traces: [{ ...bar(frac, y, { horizontal: true, name: "estimable fraction", slot: 1, text: order.map((a) => `${est[a].n_pairs_estimable} of ${est[a].n_pairs_total}`), textposition: "outside", hover: "%{customdata}<extra></extra>" }),
                          customdata: order.map((a) => `${LAYER_LABEL[a]}: ${est[a].n_pairs_estimable} of ${est[a].n_pairs_total} pairs (${est[a].n_tissues} tissues, ${est[a].n_samples} samples)<br>variables: ${est[a].variables_used.replace(/;/g, ", ")}${est[a].estimable_pairs ? "<br>pairs: " + est[a].estimable_pairs.replace(/\|/g, "–").replace(/;/g, ", ") : ""}`) }],
               layout: { xaxis: { range: [0, 0.25], title: { text: "fraction of tissue pairs estimable" } }, yaxis: { autorange: "reversed", automargin: true }, margin: { t: 20, l: 10 }, showlegend: false },
               table: { columns: ["assay", "n_tissues", "n_samples", "n_pairs_total", "n_pairs_estimable", "estimable_pairs", "variables_used"], rows: N.estimable_pairs } };
    },
    source: "results/16_identifiability/estimable_pairs.csv (layers: " + Object.entries(N.layers).filter(([, v]) => v.status === "recomputed").map(([k]) => LAYER_LABEL[k]).join(", ") + ")",
    notShow: "whether an estimable pair is confounded with something else: the one RNA-seq pair is ovary vs testes, which is also female vs male. ATAC-seq is not an exception: its nuclei-extraction, tagmentation and PCR dates cross tissues, but each of its flowcells holds one tissue."
      + (Object.values(N.layers).some((v) => v.status !== "recomputed") ? " Not measured: " + Object.entries(N.layers).filter(([, v]) => v.status !== "recomputed").map(([k, v]) => `${LAYER_LABEL[k] || k} (${v.reason})`).join("; ") + "." : ""), height: "short",
  });

  // ---- QC-only ------------------------------------------------------------------------------------------
  const qs = Object.fromEntries(Q.summary.map((r) => [r.features, r]));
  const pf = Q.per_fold;
  const geneBase = PC.baselines.find((r) => r.model === "logreg_l2");
  document.getElementById("p-qc").replaceChildren(
    `A multinomial logistic regression that sees only the per-library QC numbers, never a gene, identifies the tissue on the fingerprint's animal-grouped folds (balanced accuracy): ${fmt(qs.technical.bal_acc_mean)} ± ${fmt(qs.technical.bal_acc_sd)} from ${qs.technical.n_features} technical numbers (RIN, adapter and duplication rates, GC, depth), ${fmt(qs.composition.bal_acc_mean)} ± ${fmt(qs.composition.bal_acc_sd)} from ${qs.composition.n_features} composition fractions and ${fmt(qs.all.bal_acc_mean)} ± ${fmt(qs.all.bal_acc_sd)} from both, against ${fmt(Q.info.chance)} by chance and ${fmt(geneBase.balanced_accuracy_mean)} for the all-gene model. `,
    `Library QC numbers alone identify the tissue at ${fmt(ex.qc_all)} — partly because tissues genuinely differ in composition (see the note below).`,
  );
  await figure(document.getElementById("fig-qc"), {
    title: `Library QC numbers alone identify the tissue at ${fmt(qs.all.bal_acc_mean)}: technical numbers, composition fractions, and both together`,
    subtitle: "Balanced accuracy per fold (dots) and mean ± sd (bars) by feature set; chance and the all-gene model for reference.",
    build: () => {
      const t = tokens(); const p = palette();
      const sets = ["technical", "composition", "all"];
      const x = [...sets.map((s) => `${s} (${qs[s].n_features})`), "gene-based model"];
      const y = [...sets.map((s) => qs[s].acc_mean), geneBase.balanced_accuracy_mean];
      const sdv = [...sets.map((s) => qs[s].acc_sd), geneBase.balanced_accuracy_std];
      const dots = strip(pf.map((r) => sets.indexOf(r.features) + (r.fold - 2) * 0.05), pf.map((r) => r.balanced_accuracy), { color: t.ink2, size: 7, name: "per fold", hover: "%{customdata}<extra></extra>" });
      dots.customdata = pf.map((r) => `${r.features}, fold ${r.fold}: ${fmt(r.balanced_accuracy)} (${r.n_test_animals} test animals, ${r.n_test_vials} vials)`);
      dots.x = pf.map((r) => r.features === "technical" ? 0 + (r.fold - 2) * 0.06 : r.features === "composition" ? 1 + (r.fold - 2) * 0.06 : 2 + (r.fold - 2) * 0.06);
      return { traces: [{ ...bar(x, y, { name: "mean balanced accuracy", slot: 1, sd: sdv, text: y.map((v) => fmt(v)), hover: "%{x}: %{y:.3f}<extra>mean</extra>" }), marker: { color: [p[0], p[0], p[0], t.axis], line: { color: t.surface, width: 2 }, cornerradius: 4 } }, { ...dots, xaxis: "x", x: dots.x.map((v) => x[Math.round(v)]) }],
               layout: { yaxis: { range: [0, 1.08], title: { text: "balanced accuracy" } }, shapes: [{ type: "line", xref: "paper", x0: 0, x1: 1, yref: "y", y0: Q.info.chance, y1: Q.info.chance, line: { color: t.ink2, width: 1, dash: "dash" } }],
                         annotations: [{ xref: "paper", x: 1, xanchor: "right", yref: "y", y: Q.info.chance, yanchor: "bottom", text: "chance (1/19)", showarrow: false, font: { color: t.ink2, size: 11 } }], margin: { t: 40 }, legend: { y: 1.12 } },
               table: { columns: ["features", "n_features", "acc_mean", "acc_sd", "bal_acc_mean", "bal_acc_sd", "n_folds", "n_test_animals_mean"], rows: Q.summary } };
    },
    source: "results/16_identifiability/qc_only_summary.csv, qc_only_per_fold.csv (QC table: " + Q.info.qc_source + "); results/04_baselines/TRNSCRPT/summary.csv",
    notShow: "tuning: the QC model is untuned (C = 1), so it is a floor, not the best a QC-only model could do.",
  });
  document.getElementById("qc-caveat").replaceChildren(callout("caveat", "The composition fractions are not pure processing",
    "The composition set (mitochondrial, globin, rRNA, intronic and intergenic read fractions, chrX/chrY) is read biologically by MoTrPAC itself: a higher mitochondrial fraction in trained muscle is mitochondrial biogenesis. The technical set is closer to pure processing, though RIN and duplication also depend on the tissue's RNA; the two are shown separately and neither number is a verdict on its own."));

  // ---- bridge ---------------------------------------------------------------------------------------------
  const br = N.bridge;
  const bridgeEl = document.getElementById("bridge-block");
  if (!br || br.status !== "recomputed") {
    bridgeEl.replaceChildren(pendingBlock("Bridge-sample variance measurement", (br && br.reason) || "not computed (scripts/16_identifiability.py --bridge)"));
  } else {
    const sum = br.summary;
    const row = (bid, gs) => sum.find((r) => r.pool_bid === bid && r.gene_set === gs);
    const g99 = row(80001, "all_genes"), g88 = row(80000, "all_genes"), e99 = row(80001, "all_genes_expressed_in_pool"), p99 = row(80001, "panel_k20_expressed_in_pool");
    const others = sum.filter((r) => r.gene_set === "all_genes" && ![80000, 80001].includes(r.pool_bid));
    const pc = (v) => `${(100 * v).toFixed(1)} %`;
    bridgeEl.replaceChildren(
      el("p", {}, [
        `The consortium ran ${br.info.n_ref_vials} reference-standard vials from ${br.info.n_pools} RNA pools; ${br.info.n_bridging_pools} bridging pools were run on more than one extraction plate, and two gastrocnemius-derived pools (types 99 and 88) were run on ${g99.n_plates} plates at both sequencing sites. `,
        `For those two pools the variance of a gene's log2 CPM across plates, summed over all ${g99.n_genes.toLocaleString()} genes, is `, el("b", {}, pc(g99.sum_ratio_batch_over_tissue)), ` and ${pc(g88.sum_ratio_batch_over_tissue)} of the summed variance of the 19 tissue means`,
        ` (${pc(e99.sum_ratio_batch_over_tissue)} over the ${e99.n_genes.toLocaleString()} genes the pool expresses; median per-gene ratio ${fmt(g99.median_ratio_batch_over_tissue, 3)}), ${p99 ? pc(p99.sum_ratio_batch_over_tissue) : "—"} on the ${p99 ? p99.n_genes : 0} panel genes the pool expresses`,
        others.length ? `, and the within-site pools agree in order of magnitude (${others.map((r) => `${r.pool_tissue.replace(" Powder", "").toLowerCase()} pool ${r.pool_type}, ${r.n_plates} plates: ${pc(r.sum_ratio_batch_over_tissue)}`).join("; ")}).` : ".",
      ]),
    );
    const fig = el("div");
    bridgeEl.appendChild(fig);
    const pg = br.per_gene.filter((r) => r.pool_bid === 80001 || r.pool_bid === 80000);
    await figure(fig, {
      title: "On the panel genes a bridging pool expresses, batch is a few percent of the tissue-separating variance",
      subtitle: "Per panel gene: between-plate variance of the two gastrocnemius pools (6 plates, both sites) over the variance of the 19 tissue means; genes the pool does not express (mean log2 CPM < 1) are hollow and read zero.",
      build: () => {
        const t = tokens(); const p = palette();
        const genes = [...new Set(pg.map((r) => r.gene_symbol))].sort((a, b) => (pg.find((r) => r.gene_symbol === a && r.pool_bid === 80001)?.ratio_batch_over_tissue ?? 0) - (pg.find((r) => r.gene_symbol === b && r.pool_bid === 80001)?.ratio_batch_over_tissue ?? 0));
        // one filled trace (expressed in the pool) and one hollow trace (not expressed) per pool, so the legend shows both
        const traces = [];
        [[80001, "pool 99", 1], [80000, "pool 88", 2]].forEach(([bid, name, slot]) => {
          const color = p[slot - 1];
          for (const expressed of [true, false]) {
            const rows = genes.map((g) => pg.find((r) => r.gene_symbol === g && r.pool_bid === bid)).filter((r) => r && Boolean(r.expressed_in_pool) === expressed);
            if (!rows.length) continue;
            traces.push({
              type: "bar", orientation: "h", name: `${name}${expressed ? "" : ", not expressed (hollow)"}`,
              y: rows.map((r) => r.gene_symbol), x: rows.map((r) => r.ratio_batch_over_tissue),
              offsetgroup: name, alignmentgroup: "pools", legendgroup: name,
              marker: expressed ? { color, line: { color: t.surface, width: 2 }, cornerradius: 4 } : { color: "rgba(0,0,0,0)", line: { color, width: 1.5 }, cornerradius: 4 },
              customdata: rows.map((r) => `${r.gene_symbol} (${name}): batch/tissue ${fmt(r.ratio_batch_over_tissue, 4)}<br>V_batch ${fmt(r.v_batch, 4)}, V_tissue ${fmt(r.v_tissue, 2)}, within-tissue ${fmt(r.v_within_tissue, 3)}<br>mean log2 CPM in pool ${fmt(r.mean_log2cpm_in_pool, 2)}${expressed ? "" : " (not expressed)"}`),
              hovertemplate: "%{customdata}<extra></extra>",
            });
          }
        });
        return { traces, layout: { barmode: "group", xaxis: { title: { text: "V_batch / V_tissue" }, rangemode: "tozero" }, yaxis: { automargin: true, tickfont: { size: 10 }, categoryorder: "array", categoryarray: genes }, margin: { t: 40, l: 10 }, legend: { y: 1.1 }, bargap: 0.25 },
                 table: { columns: ["gene_symbol", "pool_bid", "pool_type", "expressed_in_pool", "mean_log2cpm_in_pool", "v_batch", "v_tissue", "v_within_tissue", "ratio_batch_over_tissue", "n_plates"], rows: pg } };
      },
      source: "results/16_identifiability/bridge_variance_per_gene.csv (reference vials from the portal per-tissue RSEM count files; study-vial tissue means from the pipeline's stacked matrix). Definition: " + br.info.definition + "; study vials and reference vials use the same unit, log2(CPM + 1) on the total library (scripts/16_identifiability.py --bridge).",
      notShow: "genes the pool does not express (hollow bars): a muscle-derived pool measures batch only on the genes it expresses, markers of other tissues (Umod, Pgk2, Hbq1b, …) read zero on every plate and contribute no batch variance, which is why the expressed-in-pool sets are the fair comparison. The plate-to-plate variance of one pool also contains ordinary technical replicate noise, so it is an upper bound on the systematic batch effect for those genes. The liver and hippocampus pools (2–3 plates, one site) are in the data table of the summary.", height: "tall",
    });
  }

  // ---- two ways processing enters the data -------------------------------------------------------------------
  const vp = Object.fromEntries(["TRNSCRPT", "PROT", "METAB"].map((a) => [a, E[`variance_${a}`]]));
  const pc1 = (a) => vp[a].find((r) => r.PC === "PC1");
  const diag = E.prot_diagnostic[0];
  document.getElementById("p-modes").replaceChildren(
    `In RNA-seq and metabolomics the tissue signal is present and large (tissue explains ${fmt(pc1("TRNSCRPT").R2_tissue, 3)} of the first transcript PC and ${fmt(pc1("METAB").R2_tissue, 3)} of the first metabolite PC) and shares its axis with the batch that processed each tissue. `,
    `In TMT proteomics the quantification removes it: values are ratios to a per-tissue reference pool, so tissue explains ${fmt(pc1("PROT").R2_tissue, 3)} of the first protein PC and the tissue is identified by each plex's missingness pattern (a missingness-only classifier reaches ${fmt(diag.missingness_indicators_only, 3)} on fold 0, ${fmt(diag.per_tissue_means_removed, 3)} once per-tissue means are removed, against ${fmt(diag.chance_balanced, 3)} by chance).`,
  );
  await figure(document.getElementById("fig-modes"), {
    title: "Tissue signal is present in RNA and metabolites and, by design, absent from TMT ratios (per-tissue reference pools)",
    subtitle: "Variance explained (R²) by tissue, sex and training group for the first three principal components of each stacked matrix.",
    build: () => {
      const assays = ["TRNSCRPT", "PROT", "METAB"];
      const labels = { TRNSCRPT: "RNA-seq", PROT: "proteomics (TMT)", METAB: "metabolomics" };
      const rows = assays.flatMap((a) => vp[a].filter((r) => ["PC1", "PC2", "PC3"].includes(r.PC)).map((r) => ({ assay: a, ...r })));
      const x = rows.map((r) => `${labels[r.assay]} ${r.PC}`);
      return { traces: [bar(x, rows.map((r) => r.R2_tissue), { name: "tissue", slot: 1, hover: "%{x}: R² tissue %{y:.3f}<extra></extra>" }), bar(x, rows.map((r) => r.R2_sex), { name: "sex", slot: 5, hover: "%{x}: R² sex %{y:.3f}<extra></extra>" }), bar(x, rows.map((r) => r.R2_group), { name: "training group", slot: 4, hover: "%{x}: R² group %{y:.3f}<extra></extra>" })],
               layout: { barmode: "group", yaxis: { range: [0, 1.05], title: { text: "R² of the PC" } }, xaxis: { tickangle: -30, tickfont: { size: 10 } }, margin: { t: 40, b: 80 }, legend: { y: 1.12 } },
               table: { columns: ["assay", "PC", "explained", "R2_tissue", "R2_sex", "R2_group"], rows } };
    },
    source: "results/03_eda/variance_partition_{TRNSCRPT,PROT,METAB}.csv; the missingness-only diagnostic in the paragraph: results/04_baselines/PROT/diagnostic_accuracy.csv", notShow: "the fraction of variance each PC carries (in the data table) and the batch covariates (results/03_eda/batch_partition_*.csv): in RNA-seq the plate, library batch and flowcell explain the same PCs as tissue, because they are the same partition.",
  });

  // ---- resolution ------------------------------------------------------------------------------------------
  const bm = H.tiles.find((t) => t.id === "tile_bodymap_k20"), cov = H.tiles.find((t) => t.id === "tile_bodymap_cov_k20");
  const bm20 = pick(H.ladder, "different_lab", "k20", "marginal");
  document.getElementById("p-resolution").replaceChildren(
    "Within MoTrPAC a classifier at the fingerprint's accuracy could in principle be reading the batch; the rat BodyMap was collected, extracted and sequenced by another laboratory, with none of the MoTrPAC plates, batches or flowcells, and the MoTrPAC-trained 20-gene panel still names ",
    el("b", {}, fmt(bm.value)), " of its mapped adult organs. Only tissue biology can do that; what stays study-specific, as a batch-influenced quantity would, is the confidence scale: with MoTrPAC thresholds the 90 % sets cover ", el("b", {}, fmt(cov.value)),
    ` of the same adults, the shortfall almost all empty sets, and three target animals restore observed coverage (${fmt(bm20?.recal_n3)}; per draw ${fmt(RB3.min_coverage, 2)}–${fmt(RB3.max_coverage, 2)}). `,
    el("a", { href: "transfer.html" }, "The Transfer page has the full ladder →"),
  );
}

main().catch((e) => { console.error(e); document.querySelector(".lede").textContent = "Failed to load: " + e.message; });
