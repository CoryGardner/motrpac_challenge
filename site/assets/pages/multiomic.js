// The multiomic follow-up page (branch multiomic-overnight). Every number is read from site/data/multiomic.json, whose values
// trace to results_multiomic/ and results_frozen/ files through the mo_* entries of site/data/provenance.json.
import { mountChrome, loadJSON, el, fmt, statTile, callout, tableFrom, segmented, control, badge } from "../site.js";
import { figure, bar, line, refLine, tokens, palette, hexAlpha, template, CONFIG, TISSUE_NAMES } from "../charts.js";

const MODEL_LABEL = { k20: "20 features", k50: "50 features", full: "all features" };
const perSet = (v, n) => `${fmt(v, 2)} of ${n} classes per set`;
const ci = (c, d = 2) => (c && c[0] !== null && c[0] !== undefined ? ` [${fmt(c[0], d)}, ${fmt(c[1], d)}]` : "");

async function main() {
  await mountChrome("multiomic.html");
  const M = await loadJSON("data/multiomic.json");
  const S = M.scales, X = M.ladder_species, W = M.ladder_within, RC = M.recalibration, R0 = M.recalibration_n0, FU = M.fusion, D = M.design, RP = M.rna_protein, SA = M.sato;
  const pick = (arr, pred) => arr.find(pred);
  const jr20 = pick(X, (r) => r.target.startsWith("protein → Jiang 2020 (cleaned") && r.model === "k20");
  const jw20 = pick(X, (r) => r.target.startsWith("protein → Jiang 2020 (raw") && r.model === "k20");
  const gt20 = pick(X, (r) => r.target.startsWith("RNA → GTEx") && r.model === "k20");
  const rna42 = pick(X, (r) => r.layer === "RNA" && r.target.startsWith("RNA, 7-class") && r.model === "k20");
  const rec = (src, model, n) => RC.find((r) => r.source.startsWith(src) && r.model === model && r.n_recal === n);
  const jrec5 = rec("Phase 3 protein → Jiang 2020 (cleaned", "k20", 5), jrec3 = rec("Phase 3 protein → Jiang 2020 (cleaned", "k20", 3);
  const grec3 = rec("frozen phase 13", "k20", 3), brec3 = rec("frozen phase 12", "k20", 3);
  const pc1 = S.rows.find((r) => r.PC === "PC1");
  const mw = M.metabolites.find((m) => m.leg === "deep_mw");
  const DJ = M.design_jiang;

  document.getElementById("status-line").replaceChildren(badge("ambiguous", "◔", "follow-up on the branch multiomic-overnight — not part of the frozen submission; merging is a human decision"));
  document.getElementById("lede").replaceChildren(
    `The submission left proteomics and metabolomics out of the cross-tissue fingerprint because the distributed proteomics are ratios to per-tissue reference pools. On the portal's reporter-ion intensities tissue explains R² ${fmt(pc1.r2_rii)} of the first principal component (${fmt(pc1.r2_ratio, 4)} on the ratios). `,
    `A 20-protein panel selected on that scale names the tissue of ${fmt(jr20.accuracy)} of ${jr20.n_samples} human TMT samples from ${jr20.n_individuals} GTEx donors (chance ${fmt(1 / 7)}; ${fmt(jw20.accuracy)} of ${jw20.n_samples} when both sides are processed the same way), and its 90 % sets cover ${fmt(jr20.coverage)} with MoTrPAC calibration and ${fmt(jrec5.coverage_recalibrated)} after recalibrating on five donors — at ${perSet(jrec5.set_size_recalibrated, jrec5.n_classes_label_space)}: the GTEx pattern, not the BodyMap one. `,
    `A 20-metabolite panel names the organ of ${fmt(mw.acc_k20)} of ${mw.n_mapped} samples from ${mw.n_individuals} mice in a mouse aging atlas, and in one external TMT design ${DJ.n_pairs_estimable} of ${DJ.n_pairs_total} tissue pairs share a run.`,
  );

  // ---- tiles --------------------------------------------------------------------------------------
  document.getElementById("tiles").replaceChildren(
    statTile({ id: "tile-r2", value: pc1.r2_rii, format: 3, label: "tissue R² of PC1, reporter-ion scale", sub: `${fmt(pc1.r2_ratio, 4)} on the distributed ratios; permutation null 95th pct ${fmt(pc1.r2_rii_null95)}; ${S.n_vials} vials, ${S.n_animals} animals, ${S.n_proteins_inner} proteins`, source: "results_multiomic/01_rii/variance_partition.csv, variance_partition_ratio.csv (mo_r2_rii_PC1, mo_r2_ratio_PC1)" }),
    statTile({ id: "tile-acc", value: jr20.accuracy, ci: jr20.accuracy_ci, format: 3, label: "20-protein panel on human TMT samples", sub: `${jr20.n_samples} Jiang 2020 samples, ${jr20.n_individuals} donors, 5 mapped tissues; chance ${fmt(1 / 7, 2)}; ${fmt(jw20.accuracy)} on the raw-ppm scale`, source: "results_multiomic/03_prot_transfer/accuracy_overall.csv (mo_jiang_relative_k20_accuracy; donor bootstrap)" }),
    statTile({ id: "tile-cov", value: jrec5.coverage_recalibrated, format: 3, label: "coverage after recalibrating on 5 donors", sub: `${fmt(jr20.coverage)} with MoTrPAC calibration → ${fmt(jrec5.coverage_recalibrated)} at ${perSet(jrec5.set_size_recalibrated, jrec5.n_classes_label_space)}`, source: "results_multiomic/03_prot_transfer/conformal_transfer.csv, recalibration.csv (mo_jiang_relative_k20_coverage; recalibration table)" }),
    statTile({ id: "tile-design", value: DJ.n_pairs_estimable, format: "of", total: DJ.n_pairs_total, label: "tissue pairs estimable within a TMT run (Jiang 2020)", sub: `up to ${DJ.max_tissues_per_run} tissues per run, ${DJ.n_runs} runs, Cramér's V ${fmt(DJ.cramers_v)}; MoTrPAC proteomics: 0 of ${DJ.motrpac_prot_pairs_total}, V ${fmt(DJ.motrpac_prot_cramers_v)}`, source: "results_multiomic/06_external_identifiability/design_comparison.csv (mo_jiang_pairs_est)" }),
  );

  // ---- caveats -------------------------------------------------------------------------------------
  document.getElementById("callout-caveats").replaceChildren(callout("caveat", "Read every number with these caveats", [
    `Within-study accuracies are context, not findings: on the reporter-ion scale, as in the audit, one TMT plex holds one tissue plus that tissue's reference pool (plex is nested in tissue), and the missing-value pattern alone classifies the tissue at ${fmt(S.missingness_outer_acc)}. `
    + "The human atlases are adult donors of both sexes, post-mortem, with their own TMT reference and search pipeline, so species, age, death and processing shift together. "
    + "Metabolite matches across platforms are RefMet name matches, not identical analyte measurements. Every recalibrated coverage is shown with its mean set size and the number of classes in its label space; a set that holds most of the label space restores the number, not the information.",
  ]));

  // ---- 1. same proteins, two scales ----------------------------------------------------------------
  document.getElementById("p-scales").replaceChildren(
    `Rolling the peptide reporter-ion intensities up to proteins and dividing each channel by its total gives a value comparable across plexes and tissues in the sense “fraction of this sample's quantified protein mass”. On that matrix tissue explains R² ${fmt(pc1.r2_rii)} of PC1 (label-permutation null 95th percentile ${fmt(pc1.r2_rii_null95)}); on the distributed ratios, ${fmt(pc1.r2_ratio, 4)}. `,
    `The result is the same on the ${S.n_proteins_complete} proteins with no missing value (${fmt(S.rows.find((r) => r.PC === "PC1").r2_complete)}), so it is not an imputation artefact; a 20-protein panel reaches balanced accuracy ${fmt(S.panel_k20_bal_acc)} under animal-grouped folds (permutation null ${fmt(S.panel_k20_null95)}), and ${S.stability_ge_0_8} of its members are selected in ≥ 80 % of animal bootstraps.`,
  );
  await figure(document.getElementById("fig-scales"), {
    title: "The same proteins carry a tissue axis on the reporter-ion scale and none on the ratio scale",
    subtitle: "Tissue R² of the first three principal components: distributed ratios to per-tissue reference pools vs reporter-ion log2 ppm (all proteins quantified in every tissue, and the complete-protein subset); dashes: permutation null 95th percentile.",
    build: () => {
      const t = tokens(); const p = palette();
      const x = S.rows.map((r) => r.PC);
      const mk = (vals, name, slot) => ({ ...bar(x, vals, { name, slot, text: vals.map((v) => fmt(v, 3)), hover: "%{x}: %{y:.4f}<extra>" + name + "</extra>" }) });
      const traces = [mk(S.rows.map((r) => r.r2_ratio), "distributed ratios", 8), mk(S.rows.map((r) => r.r2_rii), "reporter-ion log2 ppm", 1), mk(S.rows.map((r) => r.r2_complete), "reporter-ion, complete proteins only", 3),
                      { type: "scatter", mode: "markers", name: "permutation null (95th pct)", x, y: S.rows.map((r) => r.r2_rii_null95), marker: { symbol: "line-ew", size: 22, color: t.ink2, line: { width: 2, color: t.ink2 } }, hovertemplate: "%{x}: null95 %{y:.3f}<extra></extra>" }];
      return { traces, layout: { barmode: "group", yaxis: { range: [0, 1.08], title: { text: "tissue R² of the component" } }, margin: { t: 40 }, legend: { y: 1.14 } },
               table: { columns: ["PC", "r2_ratio", "r2_rii", "r2_complete", "r2_rii_null95", "explained_ratio", "explained_rii"], rows: S.rows } };
    },
    source: "results_multiomic/01_rii/variance_partition.csv, variance_partition_ratio.csv, variance_partition_complete.csv (phase-03 PCA code on the stacked matrix)",
    notShow: "whether the axis is tissue or plex: one plex holds one tissue, so R2_plex_id equals R2_tissue here (variance_partition.csv), exactly the audit's confound; an external design settles it (section 5).",
    height: "short",
  });

  // ---- 2. ladder protein vs RNA ---------------------------------------------------------------------
  const lad = { model: "k20" };
  const RUNGS = [
    { id: "in_distribution", label: "Held-out animals<br>(MoTrPAC)", w: true }, { id: "train_control_test_trained", label: "Trained animals<br>(fit on controls)", w: true },
    { id: "train_male_test_female", label: "Held-out sex<br>M → F", w: true }, { id: "train_female_test_male", label: "Held-out sex<br>F → M", w: true },
    { id: "species", label: "Other species<br>protein → Jiang, RNA → GTEx" }, { id: "same42", label: "Other species, same 42 samples<br>(7-class RNA vs protein)" },
  ];
  const rungRow = (rung, layer, model) => {
    if (rung.w) return W.find((r) => r.rung_id === rung.id && r.layer === layer && r.model === model);
    if (rung.id === "species") return layer === "protein" ? X.find((r) => r.target.startsWith("protein → Jiang 2020 (cleaned") && r.model === model) : X.find((r) => r.target.startsWith("RNA → GTEx") && r.model === model);
    return X.find((r) => r.layer === layer && r.target.startsWith(`${layer}, 7-class`) && r.model === model);
  };
  const ladBuild = () => {
    const t = tokens();
    const x = RUNGS.map((r) => r.label);
    const series = [["protein", "accuracy", 1], ["protein", "coverage", 3], ["RNA", "accuracy", 2], ["RNA", "coverage", 4]];
    const traces = series.map(([layer, key, slot]) => {
      const rows = RUNGS.map((r) => rungRow(r, layer, lad.model));
      const vals = rows.map((r) => (r && r[key] !== null && r[key] !== undefined ? r[key] : null));
      const hov = rows.map((r) => (r ? `${layer} ${key} ${fmt(r[key])}${r[key + "_ci"] ? " 95 % CI" + ci(r[key + "_ci"]) : r[key + "_sd"] ? " ± " + fmt(r[key + "_sd"]) + " sd" : ""}<br>n = ${r.n_samples} samples, ${r.n_individuals} animals/donors; ${r.n_classes} classes${key === "coverage" ? "<br>mean set size " + fmt(r.set_size, 2) + ", empty " + fmt(r.empty, 2) : ""}` : "not run"));
      const plus = rows.map((r) => (r && r[key + "_ci"] ? r[key + "_ci"][1] - r[key] : r && r[key + "_sd"] ? r[key + "_sd"] : 0));
      const minus = rows.map((r) => (r && r[key + "_ci"] ? r[key] - r[key + "_ci"][0] : r && r[key + "_sd"] ? r[key + "_sd"] : 0));
      return { ...bar(x, vals, { name: `${layer} ${key}`, slot, text: vals.map((v) => (v === null ? "" : fmt(v, 2))), hover: "%{customdata}<extra></extra>" }), customdata: hov,
               error_y: { type: "data", symmetric: false, array: plus, arrayminus: minus, visible: true, color: t.ink2, thickness: 1.2, width: 3 } };
    });
    const ref = refLine(0.9, "1 − α = 0.90");
    const rows = [];
    for (const r of RUNGS) for (const layer of ["protein", "RNA"]) { const q = rungRow(r, layer, lad.model); if (q) rows.push({ rung: r.label.replace(/<br>/g, " "), layer, model: lad.model, accuracy: q.accuracy, coverage: q.coverage, set_size: q.set_size, empty: q.empty, n_samples: q.n_samples, n_individuals: q.n_individuals, n_classes: q.n_classes, source: q.source }); }
    return { traces, layout: { barmode: "group", yaxis: { range: [0, 1.1], title: { text: "fraction" } }, xaxis: { tickfont: { size: 10 } }, shapes: ref.shapes, annotations: ref.annotations, legend: { y: 1.16 }, margin: { t: 48, b: 70 } }, table: { columns: ["rung", "layer", "model", "accuracy", "coverage", "set_size", "empty", "n_samples", "n_individuals", "n_classes", "source"], rows } };
  };
  const idP = W.find((r) => r.rung_id === "in_distribution" && r.layer === "protein" && r.model === "k20"), ctP = W.find((r) => r.rung_id === "train_control_test_trained" && r.layer === "protein" && r.model === "k20"), ctR = W.find((r) => r.rung_id === "train_control_test_trained" && r.layer === "RNA" && r.model === "k20");
  const mfP = W.find((r) => r.rung_id === "train_male_test_female" && r.layer === "protein" && r.model === "k20"), mfR = W.find((r) => r.rung_id === "train_male_test_female" && r.layer === "RNA" && r.model === "k20");
  const p42 = X.find((r) => r.layer === "protein" && r.target.startsWith("protein, 7-class") && r.model === "k20");
  document.getElementById("p-ladder").replaceChildren(
    `Inside MoTrPAC the 20-protein fingerprint (7 tissues) behaves like the 20-gene one (19 tissues): held-out animals ${fmt(idP.accuracy)} accuracy and ${fmt(idP.coverage)} coverage at ${fmt(idP.set_size, 2)} per set; fit on the sedentary controls and tested on the trained animals, ${fmt(ctP.accuracy)} / ${fmt(ctP.coverage)} (RNA ${fmt(ctR.accuracy)} / ${fmt(ctR.coverage)}); held-out sex M → F ${fmt(mfP.accuracy)} / ${fmt(mfP.coverage)} (RNA ${fmt(mfR.accuracy)} / ${fmt(mfR.coverage)}). `,
    `Across species the protein rung is lower than the RNA rung on its own target (${fmt(jr20.accuracy)} on Jiang vs ${fmt(gt20.accuracy)} on GTEx), and on the same 42 human samples the 7-class RNA fingerprint reaches ${fmt(rna42.accuracy)} where protein reaches ${fmt(p42.accuracy)}${ci(p42.accuracy_ci)}; both collapse in coverage under MoTrPAC calibration.`,
  );
  const ladCtl = el("div", { class: "controls" }, [control("Model", segmented([["k20", "20 features"], ["full", "all features"]], lad.model, (v) => { lad.model = v; ladFig.rerender(); }, "model"))]);
  const ladFig = await figure(document.getElementById("fig-ladder"), {
    title: "Protein and RNA fingerprints hold inside MoTrPAC and both lose their guarantee across species; the protein rung is lower",
    subtitle: "Per rung: accuracy and coverage of α = 0.10 sets calibrated on the source, protein (reporter-ion scale, 7 classes) beside RNA (log2 CPM, 19 classes; 7 on the same-42 rung). Whiskers: fold sd within MoTrPAC, 95 % donor-bootstrap intervals across species.",
    build: ladBuild, toolbar: ladCtl,
    source: "results_multiomic/09_extensions/ladder_side_by_side.csv (within MoTrPAC; RNA rows from results_frozen/08_shift/TRNSCRPT/shift_table.csv, results_frozen/05_panels, results_frozen/06_conformal), results_multiomic/03_prot_transfer/accuracy_overall.csv + conformal_transfer.csv + coverage_ci.csv, results_frozen/13_gtex/accuracy_overall.csv + conformal_transfer.csv, results_multiomic/05_fusion_transfer/fusion_transfer_summary.csv",
    notShow: "the raw-ppm processing of the human atlas (accuracy " + fmt(jw20.accuracy) + ci(jw20.accuracy_ci) + " on " + jw20.n_samples + " samples incl. technical replicates; recalibration section) and the two atlases with one sample per tissue (results_multiomic/03_prot_transfer/wang2019, geiger2013). Kidney and adipose have no human protein target; the GTEx rung has 17 tissues, the Jiang rung 5.",
    height: "tall",
  });
  ladFig.rerender = async () => { const b = ladBuild(); ladFig.traces = b.traces; ladFig.table = b.table; await window.Plotly.react(ladFig.chart, b.traces, { ...template(), ...b.layout }, CONFIG); };
  const ptRows = M.per_tissue_k20.filter((r) => ["jiang_relative_k20", "jiang_rawppm_k20", "same42_RNA_k20", "same42_protein_k20"].includes(r.run));
  document.getElementById("tbl-per-tissue").replaceChildren(el("details", { class: "fig-notes" }, [el("summary", {}, "Per-tissue accuracy at k20 with n samples and n donors (verification table)"),
    tableFrom({ columns: ["run", "layer", "rat_class", "jiang_tissue", "n_samples", "n_donors", "accuracy", "top_prediction"], rows: ptRows }),
    el("p", {}, [el("b", {}, "Source: "), "results_multiomic/08_verification/per_tissue_k20.csv"])]));

  // ---- 3. recalibration -----------------------------------------------------------------------------
  const SRC = [["Phase 3 protein → Jiang 2020 (cleaned relative)", "protein → Jiang (cleaned relative), 7 classes", 1], ["Phase 3 protein → Jiang 2020 (raw ppm)", "protein → Jiang (raw ppm), 7 classes", 3],
               ["frozen phase 13 RNA → GTEx", "RNA → GTEx, 19 classes", 2], ["Phase 4 deep_mw", "metabolites → mouse aging atlas, 9 classes", 4], ["Phase 4 deep_sato", "metabolites → Sato 2022 (deep), 9 classes", 5], ["Phase 4 hilic_sato", "metabolites → Sato 2022 (HILIC+), 19 classes", 7]];
  const recSeries = (key) => SRC.map(([src, name, slot]) => {
    const n0 = R0.find((r) => r.source === src && r.model === "k20");
    const pts = [{ n: 0, cov: n0 ? n0.coverage : null, size: n0 ? n0.set_size : null, ncls: n0 ? n0.n_classes_label_space : null }, ...[3, 5].map((n) => { const r = rec(src, "k20", n); return r ? { n, cov: r.coverage_recalibrated, size: r.set_size_recalibrated, ncls: r.n_classes_label_space, tests: r.n_test_individuals, draws: r.draws } : null; }).filter(Boolean)];
    const tr = line(pts.map((p) => String(p.n)), pts.map((p) => (key === "cov" ? p.cov : p.size)), { name, slot, hover: "%{customdata}<extra></extra>" });
    tr.customdata = pts.map((p) => `${name}<br>${p.n === 0 ? "MoTrPAC calibration" : p.n + " target individuals, " + p.draws + " draws, tested on " + p.tests}<br>coverage ${fmt(p.cov)}, mean set size ${fmt(p.size, 2)} of ${p.ncls} classes`);
    return tr;
  });
  document.getElementById("p-recal").replaceChildren(
    `Recalibrating the threshold on a few target individuals restores the coverage number everywhere; the set size says what it costs. On Jiang, three donors bring the protein sets to ${fmt(jrec3.coverage_recalibrated)} at ${perSet(jrec3.set_size_recalibrated, jrec3.n_classes_label_space)} and five to ${fmt(jrec5.coverage_recalibrated)} at ${perSet(jrec5.set_size_recalibrated, jrec5.n_classes_label_space)} — the GTEx pattern (RNA: ${fmt(grec3.coverage_recalibrated)} at ${perSet(grec3.set_size_recalibrated, grec3.n_classes_label_space)} with three donors), not the BodyMap pattern (${fmt(brec3.coverage_recalibrated)} at ${perSet(brec3.set_size_recalibrated, brec3.n_classes_label_space)}). `,
    "The metabolite recalibrations behave the same way; the HILIC+ → Sato sets hold most of the 19 classes.",
  );
  await figure(document.getElementById("fig-recal-cov"), {
    title: "Coverage is restored by three to five target individuals in every transfer",
    subtitle: "Coverage of the α = 0.10 sets (20-feature panels) vs number of target individuals used to recalibrate (0 = MoTrPAC calibration); mean over draws.",
    build: () => { const ref = refLine(0.9, "0.90"); return { traces: recSeries("cov"), layout: { xaxis: { type: "category", title: { text: "target individuals used for recalibration" } }, yaxis: { range: [0, 1.05], title: { text: "coverage" } }, shapes: ref.shapes, annotations: ref.annotations, legend: { y: 1.22, font: { size: 11 } }, margin: { t: 60 } } }; },
    source: "results_multiomic/08_verification/recalibration_set_sizes.csv (every recalibrated coverage in the report, with its set size and label-space size), results_multiomic/03_prot_transfer/conformal_transfer.csv, rawppm/conformal_transfer.csv, results_multiomic/04_metab_transfer/<leg>/conformal_transfer.csv, results_frozen/13_gtex/conformal_transfer.csv (the n = 0 points)",
    notShow: "the spread over draws (means only); the recalibration label of a mapped sample is its best-scoring mapped class.", height: "short",
  });
  await figure(document.getElementById("fig-recal-size"), {
    title: "…at set sizes that range from one tissue to most of the label space",
    subtitle: "Mean set size after the same recalibration; read against the number of classes in the hover (7 protein, 9 or 19 metabolite, 19 RNA).",
    build: () => ({ traces: recSeries("size"), layout: { xaxis: { type: "category", title: { text: "target individuals used for recalibration" } }, yaxis: { title: { text: "mean set size (classes per set)" }, rangemode: "tozero" }, legend: { y: 1.22, font: { size: 11 } }, margin: { t: 60 } } }),
    source: "results_multiomic/08_verification/recalibration_set_sizes.csv and the conformal_transfer.csv files above",
    notShow: "set sizes of the Mondrian and floored variants (full sets on these small calibration sets).", height: "short",
  });

  // ---- 4. fusion ------------------------------------------------------------------------------------
  const fu = { model: "k20" };
  const LAYERS = ["RNA", "protein", "late_mean", "stacked_LR"];
  const LAYER_LABEL = { RNA: "RNA alone", protein: "protein alone", late_mean: "late fusion (mean probability)", stacked_LR: "stacked fusion (LR on MoTrPAC)" };
  const fuBuild = () => {
    const t = tokens();
    const rows = LAYERS.map((l) => FU.find((r) => r.model === fu.model && r.layer === l));
    const x = LAYERS.map((l) => LAYER_LABEL[l].replace(" (", "<br>("));
    const mk = (key, name, slot, ciLo, ciHi) => ({ ...bar(x, rows.map((r) => r[key]), { name, slot, text: rows.map((r) => fmt(r[key], 2)), hover: "%{customdata}<extra></extra>" }),
      customdata: rows.map((r) => `${name} ${fmt(r[key])}${ciLo && r[ciLo] !== null ? " [" + fmt(r[ciLo], 2) + ", " + fmt(r[ciHi], 2) + "]" : ""}<br>n = ${r.n_mapped} samples, ${r.n_donors} donors; mean set size ${fmt(r.avg_set_size, 2)} of 7`),
      error_y: ciLo ? { type: "data", symmetric: false, array: rows.map((r) => (r[ciHi] ?? r[key]) - r[key]), arrayminus: rows.map((r) => r[key] - (r[ciLo] ?? r[key])), visible: true, color: t.ink2, thickness: 1.2, width: 3 } : undefined });
    const ref = refLine(0.9, "0.90");
    return { traces: [mk("accuracy", "accuracy", 1, "acc_ci95_low_donor_boot", "acc_ci95_high_donor_boot"), mk("coverage_motrpac_cal", "coverage (MoTrPAC calibration)", 2, "cov_ci95_low_donor_boot", "cov_ci95_high_donor_boot"), mk("frac_empty", "empty sets", 8)],
             layout: { barmode: "group", yaxis: { range: [0, 1.1], title: { text: "fraction" } }, shapes: ref.shapes, annotations: ref.annotations, legend: { y: 1.14 }, margin: { t: 48, b: 60 } },
             table: { columns: ["model", "layer", "n_mapped", "n_donors", "accuracy", "acc_ci95_low_donor_boot", "acc_ci95_high_donor_boot", "coverage_motrpac_cal", "frac_empty", "avg_set_size", "ood_frac_empty"], rows: FU.filter((r) => r.model === fu.model) } };
  };
  const f20 = Object.fromEntries(LAYERS.map((l) => [l, FU.find((r) => r.model === "k20" && r.layer === l)]));
  const cs20 = M.confusion_structure.find((r) => r.model === "k20" && r.where.startsWith("Jiang"));
  document.getElementById("p-fusion").replaceChildren(
    `Fusion is judged by transfer, not by within-study accuracy (where it was 0 of 7 tissues in the submission). On the ${f20.RNA.n_mapped} Jiang samples with both layers (${f20.RNA.n_donors} donors, five tissues), at k20 the RNA layer alone reaches accuracy ${fmt(f20.RNA.accuracy)} with coverage ${fmt(f20.RNA.coverage_motrpac_cal)}, protein ${fmt(f20.protein.accuracy)}${ci([f20.protein.acc_ci95_low_donor_boot, f20.protein.acc_ci95_high_donor_boot])} with ${fmt(f20.protein.coverage_motrpac_cal)}, and both fusions ${fmt(f20.late_mean.accuracy)} with ${fmt(f20.late_mean.coverage_motrpac_cal)}: the weaker layer pulls the fusion down, so the pre-registered rule (accuracy and coverage at least the better single layer, empty sets no higher) fails. `,
    `The two layers' confusion structures on these samples correlate at ${fmt(cs20.pearson_soft_offdiag)} (soft off-diagonal Pearson).`,
  );
  const fuCtl = el("div", { class: "controls" }, [control("Model", segmented([["k20", "20 features"], ["k50", "50 features"], ["full", "all features"]], fu.model, (v) => { fu.model = v; fuFig.rerender(); }, "model"))]);
  const fuFig = await figure(document.getElementById("fig-fusion"), {
    title: "Under the species shift fusion is not more robust than RNA alone",
    subtitle: "Same 42 Jiang samples, 7-class label space, one MoTrPAC calibration set for every model (α = 0.10): accuracy (whiskers: 95 % donor bootstrap), coverage and empty-set fraction.",
    build: fuBuild, toolbar: fuCtl,
    source: "results_multiomic/05_fusion_transfer/fusion_transfer_summary.csv, confusion_structure.csv",
    notShow: "kidney and adipose (no human target); a precise fusion benefit — two of the five classes have ≤ 5 samples, and every coverage under MoTrPAC calibration is a collapse, so 'coverage ≥' compares collapses. The protein arm uses the authors' cleaned relative scale (conservative; the raw-ppm scale is higher, section 3).",
  });
  fuFig.rerender = async () => { const b = fuBuild(); fuFig.traces = b.traces; fuFig.table = b.table; await window.Plotly.react(fuFig.chart, b.traces, { ...template(), ...b.layout }, CONFIG); };

  // ---- 5. design comparison --------------------------------------------------------------------------
  const dsort = [...D].sort((a, b) => (a.cramers_v ?? 0) - (b.cramers_v ?? 0));
  const labelD = (r) => `${r.dataset} · ${r.batch_variable}`;
  document.getElementById("p-design").replaceChildren(
    `The audit found every MoTrPAC processing batch holding one tissue. Applied to the external datasets' metadata, the same framing finds one TMT design that is not nested: in Jiang 2020 each run holds up to ${DJ.max_tissues_per_run} tissues from several donors (${DJ.n_runs} runs; Cramér's V run × tissue ${fmt(DJ.cramers_v)}; ${DJ.n_pairs_estimable} of ${DJ.n_pairs_total} tissue pairs share a run), against MoTrPAC proteomics (V ${fmt(DJ.motrpac_prot_cramers_v)}, 0 of ${DJ.motrpac_prot_pairs_total}). `,
    "Nesting is a choice of design, not a property of the assay. Wang 2019 (one run per tissue), the mouse aging atlas (one batch per organ) and MoTrPAC's other layers are nested; Sato's ROUND is a per-tissue-table label and is listed, not claimed.",
  );
  await figure(document.getElementById("fig-design-v"), {
    title: "Cramér's V between tissue and the batch variable: 1 means each batch holds one tissue",
    subtitle: "Every MoTrPAC layer (from the audit) and every external dataset obtained; the batch variable is named after the dataset.",
    build: () => ({ traces: [{ ...bar(dsort.map((r) => r.cramers_v), dsort.map(labelD), { horizontal: true, slot: 1, name: "Cramér's V", text: dsort.map((r) => fmt(r.cramers_v, 2)), textposition: "outside", hover: "%{customdata}<extra></extra>" }),
                               customdata: dsort.map((r) => `${labelD(r)}<br>V ${fmt(r.cramers_v)}; ${r.n_levels} levels, up to ${r.max_tissues_per_level} tissues per level<br>${r.n_samples} samples, ${r.n_tissues} tissues`) }],
                    layout: { xaxis: { range: [0, 1.15], title: { text: "Cramér's V (tissue × batch)" } }, yaxis: { automargin: true, tickfont: { size: 10 } }, margin: { t: 20, l: 10 }, showlegend: false, bargap: 0.3 },
                    table: { columns: ["dataset", "batch_variable", "n_samples", "n_tissues", "n_levels", "max_tissues_per_level", "cramers_v", "n_pairs_estimable", "n_pairs_total", "frac_pairs_estimable", "source"], rows: D } }),
    source: "results_multiomic/06_external_identifiability/design_comparison.csv (external datasets recomputed from their sample metadata; MoTrPAC rows from results_frozen/16_identifiability)",
    notShow: "Geiger 2013 (no per-sample batch variable shipped); whether MoTrPAC's fingerprint is or is not batch — only that the confound is MoTrPAC's sample-to-plex allocation, not TMT.", height: "tall",
  });
  await figure(document.getElementById("fig-design-pairs"), {
    title: "Fraction of tissue pairs that share a batch (estimable within one run)",
    subtitle: "Pairs of tissues with at least one common level of the batch variable, over all pairs.",
    build: () => ({ traces: [{ ...bar(dsort.map((r) => r.frac_pairs_estimable ?? 0), dsort.map(labelD), { horizontal: true, slot: 3, name: "estimable fraction", text: dsort.map((r) => `${r.n_pairs_estimable}/${r.n_pairs_total}`), textposition: "outside", hover: "%{customdata}<extra></extra>" }),
                               customdata: dsort.map((r) => `${labelD(r)}<br>${r.n_pairs_estimable} of ${r.n_pairs_total} pairs estimable`) }],
                    layout: { xaxis: { range: [0, 1.15], title: { text: "estimable tissue pairs / all pairs" } }, yaxis: { automargin: true, tickfont: { size: 10 } }, margin: { t: 20, l: 10 }, showlegend: false, bargap: 0.3 } }),
    source: "results_multiomic/06_external_identifiability/design_comparison.csv, estimable_pairs.csv",
    notShow: "which pairs (estimable_pairs.csv lists them); MoTrPAC IMMUNO (Luminex plates) was already the audit's one crossed layer.", height: "tall",
  });

  // ---- 6. RNA markers at the protein level ------------------------------------------------------------
  document.getElementById("p-markers").replaceChildren(
    `On the ${RP.n_genes} genes quantified in both layers, the cross-tissue RNA–protein Spearman correlation (tissue means over the same animals) has median ${fmt(RP.spearman_median)}; a mismatched-pair null (RNA of one gene against the protein of another) has median ${fmt(RP.null_median)} and 95th percentile ${fmt(RP.null_q95)}, so a large part of any gene's agreement is shared tissue structure, and ${fmt(100 * RP.frac_above_null_q95, 1)} % of genes exceed the null's 95th percentile. `,
    `The cleaner statement is about the RNA panel: of the ${RP.n_testable} panel genes whose marker tissue is one of the 7 proteomics tissues and whose protein is quantified, ${RP.n_same_marker} keep the same marker tissue at the protein level (chance 1/7). The protein selector, however, picks other genes for the same tissues: the two k20 panels share ${M.b1.k20_intersection} gene(s), and ${M.b1.genes_in_both_candidate_lists} gene(s) appear in both layers' wider candidate lists.`,
  );
  await figure(document.getElementById("fig-markers"), {
    title: "RNA and protein agree across tissues for most genes, but so does the mismatched-pair null",
    subtitle: `Distribution of the cross-tissue Spearman correlation over ${RP.n_genes} genes (7 tissues, tissue means over the same animals); vertical marks: observed median, null median and null 95th percentile.`,
    build: () => {
      const t = tokens(); const p = palette();
      const centers = RP.hist_edges.slice(0, -1).map((e, i) => (e + RP.hist_edges[i + 1]) / 2);
      const traces = [{ type: "bar", x: centers, y: RP.hist_counts, width: 0.095, name: "genes", marker: { color: p[0], line: { color: t.surface, width: 1 } }, hovertemplate: "ρ ≈ %{x:.2f}: %{y} genes<extra></extra>" }];
      const vline = (x, text, color, dash) => ({ shape: { type: "line", x0: x, x1: x, y0: 0, y1: 1, xref: "x", yref: "paper", line: { color, width: 2, dash } }, ann: { x, y: 1, xref: "x", yref: "paper", text, showarrow: false, yanchor: "bottom", font: { size: 11, color: t.ink2 } } });
      const marks = [vline(RP.spearman_median, "median " + fmt(RP.spearman_median, 2), p[0], "solid"), vline(RP.null_median, "null median " + fmt(RP.null_median, 2), t.ink2, "dot"), vline(RP.null_q95, "null 95 % " + fmt(RP.null_q95, 2), t.ink2, "dash")];
      return { traces, layout: { xaxis: { title: { text: "Spearman correlation across the 7 tissues" }, range: [-1, 1] }, yaxis: { title: { text: "genes" } }, shapes: marks.map((m) => m.shape), annotations: marks.map((m) => m.ann), showlegend: false, margin: { t: 40 } },
               table: { columns: ["bin_low", "bin_high", "genes"], rows: RP.hist_counts.map((c, i) => ({ bin_low: RP.hist_edges[i], bin_high: RP.hist_edges[i + 1], genes: c })) } };
    },
    source: "results_multiomic/01_rii/rna_protein_correlation.csv (binned here), rna_protein_correlation_summary.csv",
    notShow: "gene-level significance: with 7 tissues a single correlation is coarse; the panel-gene table below is the per-gene statement.", height: "short",
  });
  document.getElementById("tbl-markers").replaceChildren(el("details", { class: "fig-notes", open: "" }, [el("summary", {}, `The ${RP.n_testable} testable RNA panel genes at the protein level`),
    tableFrom({ columns: ["gene_symbol", "lists", "rna_marker_tissue", "protein_marker_tissue", "agree"], rows: RP.panel_genes.map((r) => ({ ...r, agree: r.agree ? "yes" : "no" })) }),
    el("p", {}, [el("b", {}, "Source: "), "results_multiomic/01_rii/rna_protein_panel_genes.csv (testable_c rows); B1 overlap: results_multiomic/09_extensions/overlap_summary.csv"])]));

  // ---- 7. metabolites --------------------------------------------------------------------------------
  const stopped = M.metabolites_stopped;
  document.getElementById("p-metab").replaceChildren(
    `Metabolites are matched by RefMet name. The HILIC+ platform runs in all 19 MoTrPAC tissues but only 129 named metabolites are common to every tissue, and just ${stopped.matched} of them are in the mouse aging atlas, so that leg stopped under the pre-registered rule (< ${stopped.min_overlap}). The six-platform source on the 9 core tissues matches ${mw.matched} names and names the organ of ${fmt(mw.acc_k20)}${ci(mw.acc_k20_ci)} of ${mw.n_mapped} samples from ${mw.n_individuals} mice (chance ${fmt(mw.chance)}); the Sato 2022 legs are lower. `,
    `Exercise: Sato's mice ran a single acute bout (tissues collected right after; ${SA.n_sedentary_mice_union} sedentary and ${SA.n_exercised_mice_union} exercised mice, ${SA.n_sedentary_samples_total} and ${SA.n_exercised_samples_total} samples). A native metabolite panel fit and calibrated on the sedentary mice keeps accuracy ${fmt(SA.invariance_full_accuracy)} and coverage ${fmt(SA.invariance_full_coverage)} at ${perSet(SA.invariance_full_set_size, SA.invariance_n_classes)} on the exercised mice (all features; k20 ${fmt(SA.invariance_k20_accuracy)} / ${fmt(SA.invariance_k20_coverage)} at ${fmt(SA.invariance_k20_set_size, 2)}). The mice are exercised, not trained: this is not the same shift as the RNA controls → 8-week-trained rung.`,
  );
  document.getElementById("tbl-metab").replaceChildren(el("details", { class: "fig-notes", open: "" }, [el("summary", {}, "Metabolite transfer legs (k20 unless stated), with set sizes and label-space sizes"),
    tableFrom({ columns: ["leg", "source_tissues", "matched", "n_mapped", "n_individuals", "chance", "acc_k20", "acc_full", "coverage_k20", "coverage_k20_recal5", "set_size_k20_recal5", "n_classes"], rows: M.metabolites.map((m) => ({ ...m, acc_k20: m.acc_k20 })) }),
    el("p", {}, [el("b", {}, "Source: "), "results_multiomic/04_metab_transfer/legs_summary.csv; results_multiomic/08_verification/sato_invariance_design.csv (exercise design and counts)"]),
    el("p", {}, [el("b", {}, "What it does not show: "), "identical analyte measurements — MoTrPAC HILIC+ / six platforms vs a triple-quad RP-negative panel vs Metabolon HD4; serum vs plasma, whole brain vs three rat regions and quadriceps vs gastrocnemius are imperfect maps; the recalibrated Sato coverages hold 5–12 of 9–19 classes per set."])]));

  document.getElementById("p-prov").replaceChildren(
    `Every value on this page is read from site/data/multiomic.json, exported by scripts/multiomic/export_site_data.py from the CSVs under results_multiomic/ (branch multiomic-overnight) and results_frozen/ (the RNA comparison rows), and recorded in site/data/provenance.json as mo_* entries (${M._meta.sources.length} source files; git ${M._meta.git_hash}, ${M._meta.generated}). `,
    "The full report with its pre-registration, log and per-phase READMEs is docs/MULTIOMIC_REPORT.md on the branch.",
  );
}

main().catch((e) => { console.error(e); const l = document.getElementById("lede"); if (l) l.textContent = "The page could not load its data: " + e.message; });
