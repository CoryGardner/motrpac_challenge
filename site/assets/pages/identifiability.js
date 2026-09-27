import { mountChrome, loadJSON, el, fmt, tableFrom, callout, pendingBlock } from "../site.js";
import { figure, bar, heatmap, strip, tokens, palette, tissueLabel, hexAlpha } from "../charts.js";

const LAYER_LABEL = { TRNSCRPT: "RNA-seq", METHYL: "RRBS methylation", ATAC: "ATAC-seq", PROT: "proteomics (TMT)", PHOSPHO: "phosphoproteomics (TMT)", ACETYL: "acetylproteomics (TMT)", UBIQ: "ubiquitylome (TMT)", IMMUNO: "immunoassays (Luminex)", METAB: "metabolomics" };
const VAR_LABEL = { GET_site: "site", RNA_extr_plate_ID: "extraction plate", RNA_extr_date: "extraction date", DNA_extr_plate_ID: "extraction plate", DNA_extr_date: "extraction date", Lib_prep_date: "library prep date", Lib_batch_ID: "library batch", Seq_date: "sequencing date", Seq_flowcell_ID: "flowcell", Seq_flowcell_lane: "lane", Seq_batch: "sequencing batch", Sample_batch: "sample batch", Nuclei_extr_date: "nuclei extraction date", Tagmentation_date: "tagmentation date", PCR_date: "PCR date", plex_id: "TMT plex (tissue × label)", tmt11_channel: "TMT channel", plate_id: "assay plate", panel_name: "assay panel" };

async function main() {
  await mountChrome("identifiability.html");
  const [H, N, Q, BV, E, PC] = await Promise.all([loadJSON("data/headline.json"), loadJSON("data/nesting.json"), loadJSON("data/qc_baseline.json"), loadJSON("data/batch_verdict.json"), loadJSON("data/eda.json"), loadJSON("data/panel_curve.json")]);
  const ex = H.extras;
  const est = Object.fromEntries(N.estimable_pairs.map((r) => [r.assay, r]));
  const rna = est.TRNSCRPT;
  document.getElementById("p-nesting").replaceChildren(
    `In RNA-seq, ${ex.n_plates} extraction plates, ${ex.n_lib_batches} library batches and ${ex.n_flowcells} flowcells each hold whole tissues: the median tissue spans one level of every processing variable. `,
    `A tissue pair is estimable when the two tissues share a level of every one of them, so that a contrast exists inside a batch; ${rna.n_pairs_estimable} of ${rna.n_pairs_total} RNA-seq pairs does (${rna.estimable_pairs.replace("|", " vs ").toLowerCase()}, the one pair that is also the sex contrast, since each is single-sex). `,
    `Methylation and ATAC-seq have none of 28. TMT proteomics is nested by construction: a plex is ten samples plus one tissue's reference pool. The immunoassays are the exception: their Luminex plates hold two to four tissues, so ${est.IMMUNO.n_pairs_estimable} of ${est.IMMUNO.n_pairs_total} pairs are estimable there. Metabolomics carries no batch variable in the export.`,
  );
  // nesting heatmap: rows layer · variable, two columns of colour: Cramér's V and fraction of pairs sharing a level
  const rows = [];
  for (const [assay, tab] of Object.entries(N.nesting)) for (const r of tab) rows.push({ assay, ...r });
  await figure(document.getElementById("fig-nesting"), {
    title: "Every processing variable is a near-perfect proxy for tissue, except the immunoassay plates",
    subtitle: "Per omic layer and processing variable: Cramér's V between the variable and tissue (1 = determined), the fraction of tissue pairs that share at least one level, the median number of levels a tissue spans, and the most tissues any level holds.",
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
    notShow: "metabolomics (no batch variable in the export) or the reference-standard vials; colour is a value scale per column, the printed numbers are the data.", height: "tall",
  });
  await figure(document.getElementById("fig-estimable"), {
    title: "Estimable tissue pairs per layer: one in RNA-seq, none in the epigenome or TMT layers, sixteen in the immunoassays",
    subtitle: "Pairs of tissues that share a level of every processing variable of the layer (the variables used are listed in the hover), out of all pairs the layer has.",
    build: () => {
      const order = ["TRNSCRPT", "METHYL", "ATAC", "PROT", "PHOSPHO", "ACETYL", "UBIQ", "IMMUNO"].filter((a) => est[a]);
      const y = order.map((a) => LAYER_LABEL[a]);
      const frac = order.map((a) => est[a].n_pairs_estimable / est[a].n_pairs_total);
      return { traces: [{ ...bar(frac, y, { horizontal: true, name: "estimable fraction", slot: 1, text: order.map((a) => `${est[a].n_pairs_estimable} of ${est[a].n_pairs_total}`), textposition: "outside", hover: "%{customdata}<extra></extra>" }),
                          customdata: order.map((a) => `${LAYER_LABEL[a]}: ${est[a].n_pairs_estimable} of ${est[a].n_pairs_total} pairs (${est[a].n_tissues} tissues, ${est[a].n_samples} samples)<br>variables: ${est[a].variables_used.replace(/;/g, ", ")}${est[a].estimable_pairs ? "<br>pairs: " + est[a].estimable_pairs.replace(/\|/g, "–").replace(/;/g, ", ") : ""}`) }],
               layout: { xaxis: { range: [0, 0.25], title: { text: "fraction of tissue pairs estimable" } }, yaxis: { autorange: "reversed", automargin: true }, margin: { t: 20, l: 10 }, showlegend: false },
               table: { columns: ["assay", "n_tissues", "n_samples", "n_pairs_total", "n_pairs_estimable", "estimable_pairs", "variables_used"], rows: N.estimable_pairs } };
    },
    source: "results/16_identifiability/estimable_pairs.csv", notShow: "whether an estimable pair is confounded with something else: the one RNA-seq pair is ovary vs testes, which is also female vs male.", height: "short",
  });
  const unavailable = Object.entries(N.layers).filter(([, v]) => v.status !== "recomputed");
  document.getElementById("layers-note").replaceChildren(callout("note", "What was recomputed here and what was not", [
    `Recomputed from the metadata export in this repository (scripts/16_identifiability.py): ${Object.entries(N.layers).filter(([, v]) => v.status === "recomputed").map(([k]) => LAYER_LABEL[k]).join(", ")}. `
    + (unavailable.length ? `Not available: ${unavailable.map(([k, v]) => `${LAYER_LABEL[k] || k} (${v.reason})`).join("; ")}. ` : "")
    + "ATAC-seq is not an exception: its nuclei-extraction, tagmentation and PCR dates cross tissues, but each flowcell holds one tissue, which closes it. The parallel identifiability audit (phase 21) is not present in this copy of the results; its bridge-sample measurement is marked pending below.",
  ]));

  // ---- QC-only ------------------------------------------------------------------------------------------
  const qs = Object.fromEntries(Q.summary.map((r) => [r.features, r]));
  const pf = Q.per_fold;
  const geneBase = PC.baselines.find((r) => r.model === "logreg_l2");
  document.getElementById("p-qc").replaceChildren(
    `A multinomial logistic regression that never sees a gene, only the consortium's per-library QC numbers, identifies the tissue on the same animal-grouped folds as the fingerprint: ${fmt(qs.technical.acc_mean)} ± ${fmt(qs.technical.acc_sd)} from ${qs.technical.n_features} purely technical numbers (RIN, adapter and duplication rates, GC, read depth), ${fmt(qs.composition.acc_mean)} ± ${fmt(qs.composition.acc_sd)} from ${qs.composition.n_features} composition fractions, ${fmt(qs.all.acc_mean)} ± ${fmt(qs.all.acc_sd)} from both, against ${fmt(Q.info.chance)} by chance and ${fmt(geneBase.balanced_accuracy_mean)} for the gene-based model. `,
    "The consortium's QC table is, in effect, a tissue label.",
  );
  await figure(document.getElementById("fig-qc"), {
    title: "With no gene at all, library QC numbers identify the tissue almost as well as the fingerprint",
    subtitle: "Accuracy per fold (dots) and mean ± sd (bars) of a QC-only classifier on the phase-04 animal-grouped folds, by feature set; chance and the all-gene logistic regression for reference.",
    build: () => {
      const t = tokens(); const p = palette();
      const sets = ["technical", "composition", "all"];
      const x = [...sets.map((s) => `${s} (${qs[s].n_features})`), "gene-based model"];
      const y = [...sets.map((s) => qs[s].acc_mean), geneBase.balanced_accuracy_mean];
      const sdv = [...sets.map((s) => qs[s].acc_sd), geneBase.balanced_accuracy_std];
      const dots = strip(pf.map((r) => sets.indexOf(r.features) + (r.fold - 2) * 0.05), pf.map((r) => r.accuracy), { color: t.ink2, size: 7, name: "per fold", hover: "%{customdata}<extra></extra>" });
      dots.customdata = pf.map((r) => `${r.features}, fold ${r.fold}: ${fmt(r.accuracy)} (${r.n_test_animals} test animals, ${r.n_test_vials} vials)`);
      dots.x = pf.map((r) => r.features === "technical" ? 0 + (r.fold - 2) * 0.06 : r.features === "composition" ? 1 + (r.fold - 2) * 0.06 : 2 + (r.fold - 2) * 0.06);
      return { traces: [{ ...bar(x, y, { name: "mean accuracy", slot: 1, sd: sdv, text: y.map((v) => fmt(v)), hover: "%{x}: %{y:.3f}<extra>mean</extra>" }), marker: { color: [p[0], p[0], p[0], t.axis], line: { color: t.surface, width: 2 }, cornerradius: 4 } }, { ...dots, xaxis: "x", x: dots.x.map((v) => x[Math.round(v)]) }],
               layout: { yaxis: { range: [0, 1.08], title: { text: "accuracy" } }, shapes: [{ type: "line", xref: "paper", x0: 0, x1: 1, yref: "y", y0: Q.info.chance, y1: Q.info.chance, line: { color: t.ink2, width: 1, dash: "dash" } }],
                         annotations: [{ xref: "paper", x: 1, xanchor: "right", yref: "y", y: Q.info.chance, yanchor: "bottom", text: "chance (1/19)", showarrow: false, font: { color: t.ink2, size: 11 } }], margin: { t: 40 }, legend: { y: 1.12 } },
               table: { columns: ["features", "n_features", "acc_mean", "acc_sd", "bal_acc_mean", "bal_acc_sd", "n_folds", "n_test_animals_mean"], rows: Q.summary } };
    },
    source: "results/16_identifiability/qc_only_summary.csv, qc_only_per_fold.csv (QC table: " + Q.info.qc_source + "); results/04_baselines/TRNSCRPT/summary.csv",
    notShow: "tuning: the QC model is untuned (C = 1), so it is a floor, not the best a QC-only model could do.",
  });
  document.getElementById("qc-caveat").replaceChildren(callout("caveat", "The composition fractions are not pure processing",
    "The composition set (mitochondrial, globin, rRNA, intronic and intergenic read fractions, chrX/chrY) is read biologically by MoTrPAC itself: a higher mitochondrial fraction in trained muscle is what mitochondrial biogenesis produces. The technical set is closer to pure processing, but RIN and duplication also depend on the tissue's RNA. Which is why the two are shown separately, and why neither number is a verdict on its own."));

  // ---- bridge ---------------------------------------------------------------------------------------------
  document.getElementById("bridge-block").replaceChildren(pendingBlock("Bridge-sample variance measurement",
    "The 36 reference-standard vials (Sample_category = ref in the consortium QC table, two per extraction plate) are not in the pipeline's data export (data/raw/counts holds study vials only), and the parallel identifiability audit that measured them (phase 21) is not present in this copy of the results. The measurement will appear here when results/21_identifiability/ lands; the site regenerates from it (make site)."));

  // ---- two failure modes ------------------------------------------------------------------------------------
  const vp = Object.fromEntries(["TRNSCRPT", "PROT", "METAB"].map((a) => [a, E[`variance_${a}`]]));
  const pc1 = (a) => vp[a].find((r) => r.PC === "PC1");
  const diag = E.prot_diagnostic[0];
  document.getElementById("p-modes").replaceChildren(
    `The confound takes two forms. In RNA-seq and metabolomics the tissue signal is present and enormous (tissue explains ${fmt(pc1("TRNSCRPT").R2_tissue, 3)} of the first transcript PC and ${fmt(pc1("METAB").R2_tissue, 3)} of the first metabolite PC) but it cannot be separated from the batch that processed each tissue. `,
    `In TMT proteomics the tissue signal is removed by the quantification itself: values are ratios to a per-tissue reference pool, so tissue explains ${fmt(pc1("PROT").R2_tissue, 3)} of the first protein PC, and what identifies the tissue is the missingness pattern of each plex (a missingness-only classifier reaches ${fmt(diag.missingness_indicators_only, 3)} on fold 0, ${fmt(diag.per_tissue_means_removed, 3)} once per-tissue means are removed, against ${fmt(diag.chance_balanced, 3)} by chance; Source: results/04_baselines/PROT/diagnostic_accuracy.csv).`,
  );
  await figure(document.getElementById("fig-modes"), {
    title: "Present but confounded (RNA, metabolites) vs removed by quantification (TMT proteomics)",
    subtitle: "Variance explained (R²) by tissue, sex and training group for the first three principal components of each stacked matrix, and the fraction of variance each PC carries.",
    build: () => {
      const assays = ["TRNSCRPT", "PROT", "METAB"];
      const labels = { TRNSCRPT: "RNA-seq", PROT: "proteomics (TMT)", METAB: "metabolomics" };
      const rows = assays.flatMap((a) => vp[a].filter((r) => ["PC1", "PC2", "PC3"].includes(r.PC)).map((r) => ({ assay: a, ...r })));
      const x = rows.map((r) => `${labels[r.assay]} ${r.PC}`);
      return { traces: [bar(x, rows.map((r) => r.R2_tissue), { name: "tissue", slot: 1, hover: "%{x}: R² tissue %{y:.3f}<extra></extra>" }), bar(x, rows.map((r) => r.R2_sex), { name: "sex", slot: 5, hover: "%{x}: R² sex %{y:.3f}<extra></extra>" }), bar(x, rows.map((r) => r.R2_group), { name: "training group", slot: 4, hover: "%{x}: R² group %{y:.3f}<extra></extra>" })],
               layout: { barmode: "group", yaxis: { range: [0, 1.05], title: { text: "R² of the PC" } }, xaxis: { tickangle: -30, tickfont: { size: 10 } }, margin: { t: 40, b: 80 }, legend: { y: 1.12 } },
               table: { columns: ["assay", "PC", "explained", "R2_tissue", "R2_sex", "R2_group"], rows } };
    },
    source: "results/03_eda/variance_partition_{TRNSCRPT,PROT,METAB}.csv", notShow: "batch covariates (results/03_eda/batch_partition_*.csv): in RNA-seq the plate, library batch and flowcell explain the same PCs as tissue, because they are the same partition.",
  });

  // ---- verdict ---------------------------------------------------------------------------------------------
  const conc = BV.conclusion;
  const nUnres = conc.filter((r) => String(r.verdict).startsWith("unresolvable")).length;
  document.getElementById("p-verdict").replaceChildren(
    `Within tissue, control vs 8-week-trained animals separate at AUROC ≈ 1 in every omic. Whether that is training or a collection batch was checked with covariate-only classifiers (collection, library, depth, QC, plex): in ${conc.length - nUnres} of ${conc.length} tissues no recorded covariate separates the arms and the omic separation is attributed to training; in ${nUnres} (${conc.filter((r) => String(r.verdict).startsWith("unresolvable")).map((r) => r.tissue).join(", ")}) a processing variable separates them too and the question is unresolvable from the inside.`,
  );
  document.getElementById("tbl-verdict").replaceChildren(tableFrom({ columns: ["tissue", "max_auroc_collection_library_plex", "max_auroc_qc_metrics", "max_auroc_all_covariates", "null_p95_logreg_max", "batch_variables_separating", "qc_variables_separating", "verdict"],
    rows: conc.map((r) => ({ ...r, batch_variables_separating: r.batch_variables_separating || "", qc_variables_separating: (r.qc_variables_separating || "").replace(/;/g, ", ") })), format: { max_auroc_collection_library_plex: (v) => v.toFixed(2), max_auroc_qc_metrics: (v) => v.toFixed(2), max_auroc_all_covariates: (v) => v.toFixed(2), null_p95_logreg_max: (v) => v.toFixed(2) } }),
    el("p", { class: "small" }, "Source: results/07_fusion/batch_conclusion.csv (AUROC of covariate-only classifiers, animal-grouped CV; null = 95th percentile of a within-sex label permutation)."));
  const rev = BV.time_course_revision;
  document.getElementById("verdict-revision").replaceChildren(callout(rev.status === "interpreted" ? "caveat" : "pending", "A revision is pending: the QC set behind this table contains biologically-read covariates", [
    rev.note + (rev.verdicts ? " Per-tissue readings of the time-course investigation (interpreted, source " + rev.source + "): " + Object.entries(rev.verdicts).map(([t, v]) => `${t}: ${v}`).join("; ") + "." : ""),
  ]));

  // ---- resolution ------------------------------------------------------------------------------------------
  const bm = H.tiles.find((t) => t.id === "tile_bodymap_k20"), cov = H.tiles.find((t) => t.id === "tile_bodymap_cov_k20");
  document.getElementById("p-resolution").replaceChildren(
    "Neither argument settles it alone. Within MoTrPAC, a classifier that reached the fingerprint's accuracy could in principle be reading the batch. But the rat BodyMap was collected, extracted and sequenced by another laboratory, so none of the MoTrPAC plates, batches or flowcells exist there, and the MoTrPAC-trained 20-gene panel still names ",
    el("b", {}, fmt(bm.value)), " of the adult organs correctly. A batch signature cannot do that; tissue biology can. What does not travel is the confidence scale: with MoTrPAC thresholds the 90 % sets cover ", el("b", {}, fmt(cov.value)),
    " of the same adults, almost all of the shortfall through empty sets. The ranking of tissues is biology; the calibration is study-specific, as a batch-influenced quantity would be. ",
    el("a", { href: "transfer.html" }, "The transfer page has the full ladder →"),
  );
}

main().catch((e) => { console.error(e); document.querySelector(".lede").textContent = "Failed to load: " + e.message; });
