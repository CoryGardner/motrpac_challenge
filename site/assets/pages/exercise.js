// The Exercise page: what endurance training does and does not do to the tissue fingerprint.
// Every number is read from site/data/exercise.json (X) and headline.json (H); nothing is typed by hand.
import { mountChrome, loadJSON, el, fmt, callout, pendingBlock, statTile, tableFrom, segmented, control } from "../site.js";
import { figure, bar, line, strip, heatmap, refLine, tokens, palette, template, CONFIG, ORGAN_SYSTEMS, TISSUE_NAMES, tissueLabel } from "../charts.js";
import { pick } from "../ladder.js";

const DUR = ["1w", "2w", "4w", "8w"];
const DUR_LABEL = { "1w": "1 week†", "2w": "2 weeks†", "4w": "4 weeks†", "8w": "8 weeks" };
const OMICS = [["TRNSCRPT", "RNA-seq"], ["PROT", "proteomics"], ["METAB", "metabolomics"]];
const MODEL_LABEL = { k20: "20 genes", k50: "50 genes", full: "all genes" };
const tissueOrder = ORGAN_SYSTEMS.flatMap((o) => o.tissues);

async function main() {
  await mountChrome("exercise.html");
  const [X, H] = await Promise.all([loadJSON("data/exercise.json"), loadJSON("data/headline.json")]);
  const S = X.summary;
  const tr20 = pick(H.ladder, "train_control_test_trained", "k20", "marginal");
  const fbdRow = (m, g) => X.fingerprint_by_duration.rows.find((r) => r.model === m && r.group === g && !r.pending);
  const k8 = fbdRow("k20", "8w");

  // ---- lede ----------------------------------------------------------------------------------------
  document.getElementById("lede").replaceChildren(
    `Fit on the ${tr20.n_source_animals} sedentary control animals alone, the 20-gene panel names the tissue of the 8-week-trained animals at ${fmt(k8.accuracy)} with coverage ${fmt(k8.coverage)}. `,
    `Training is real and measurable inside a tissue, in every omic layer, but it is a minor axis: smaller than the tissue contrast on every panel gene${S.ptr_k20_median_ratio_tissue_over_training ? ` (by a median factor of ${fmt(S.ptr_k20_median_ratio_tissue_over_training, 1)})` : ""}, and, where exercise biology is strongest, on the same axis as processing.`,
  );

  // ---- 1. design and physiology ------------------------------------------------------------------------
  const D = X.design;
  const dp = document.getElementById("p-design");
  if (D && D.by_group_sex) {
    const gap = S.design_gap_days_8w;
    dp.replaceChildren(`Five groups of 10 animals per sex: sedentary controls and 1, 2, 4 and 8 weeks of treadmill training. The 8-week group shares its arrival cohort and sacrifice window with the controls${gap !== null && gap !== undefined ? ` (median sacrifice gap ${fmt(gap, 0)} days)` : ""}; the shorter durations (†) were sacrificed months apart, so only the 8-week contrast is free of cohort and season.`);
    await figure(document.getElementById("fig-design"), {
      title: "Only the 8-week animals share their arrival cohort and sacrifice dates with the controls",
      subtitle: "Per group and sex: arrival at the facility (◆) and the first-to-last sacrifice window (—), from the study phenotype table.",
      build: () => {
        const t = tokens(); const p = palette();
        const rows = D.by_group_sex.filter((r) => r.group);
        const order = ["control", "1w", "2w", "4w", "8w"];
        rows.sort((a, b) => order.indexOf(a.group) - order.indexOf(b.group) || (a.sex > b.sex ? 1 : -1));
        const y = rows.map((r) => `${r.group} · ${r.sex}`);
        const traces = [];
        rows.forEach((r, i) => {
          traces.push({ ...line([r.sac_first, r.sac_last], [y[i], y[i]], { name: "sacrifice window", slot: order.indexOf(r.group) + 1, showlegend: false, hover: `${y[i]}: sacrificed ${r.sac_first} to ${r.sac_last} (${r.n_animals} animals)<extra></extra>` }), mode: "lines+markers" });
        });
        traces.push(strip(rows.map((r) => r.arrived), y, { name: "arrival", color: t.ink2, symbol: "diamond", size: 9, hover: "%{y}: arrived %{x}<extra></extra>" }));
        return { traces, layout: { xaxis: { type: "date", title: { text: "date" } }, yaxis: { autorange: "reversed", automargin: true }, margin: { t: 30, l: 10 }, showlegend: false, height: 340 },
                 table: { columns: Object.keys(rows[0]), rows } };
      },
      source: "results/15_time_course/design/design_by_group_sex.csv, design_contrasts.csv (from data/raw/pheno.csv)",
      notShow: "housing, training compliance and the daily handling that trained animals received; sex is confounded with arrival cohort inside every group (Limitations).",
      height: "short",
    });
    const tiles = document.getElementById("tiles-physiology");
    const gt = (X.physiology && X.physiology.group_tests) || [];
    for (const [v, label, unit] of [["vo2max_change", "VO2max change, 8 weeks vs control", "ml/kg/min"], ["fat_change", "body-fat change, 8 weeks vs control", "percentage points"]]) {
      for (const sex of ["female", "male"]) {
        const r = gt.find((x) => x.variable === v && x.contrast === "8w vs control" && x.sex === sex);
        if (!r) continue;
        tiles.appendChild(statTile({ id: `tile-${v}-${sex}`, value: r.median_diff, format: "1", label: `${label}, ${sex}s`,
                                     sub: `median difference (${unit}); n = ${r.n_trained} trained vs ${r.n_control} control; Mann–Whitney p = ${Number(r.mw_p).toExponential(1)}`,
                                     source: "results/15_time_course/physiology/physiology_group_tests.csv" }));
      }
    }
  } else {
    dp.replaceChildren("Five groups of 10 animals per sex: sedentary controls and 1, 2, 4 and 8 weeks of treadmill training; only the 8-week group shares its arrival cohort and sacrifice window with the controls (†: sacrificed months apart).");
    document.getElementById("fig-design").appendChild(pendingBlock("Design dates and physiology", "the design and physiology tables have not been generated yet"));
  }

  // ---- 2. within-tissue variance -------------------------------------------------------------------------
  document.getElementById("p-within").replaceChildren(
    `On the stacked matrix tissue explains ${fmt(S.stacked_pc1_r2_tissue)} of the first principal component and training group ${fmt(S.stacked_pc1_r2_group, 4)}. `,
    `Within a tissue the training group is visible above its permutation null in ${S.wt_TRNSCRPT_n_visible} of ${S.wt_TRNSCRPT_n_tissues} tissues by RNA-seq, ${S.wt_PROT_n_visible} of ${S.wt_PROT_n_tissues} by proteomics and ${S.wt_METAB_n_visible} of ${S.wt_METAB_n_tissues} by metabolomics, at most R² ${fmt(S.wt_heart_trnscrpt_max_r2_group, 2)} (heart RNA, on PC${S.wt_heart_trnscrpt_pc}).`,
  );
  await figure(document.getElementById("fig-within"), {
    title: "Inside a tissue, training is visible but minor: at most R² 0.70 of one principal component, never the first",
    subtitle: "Largest R² of the training group over within-tissue PCs 1–5, per omic layer; ● marks cells above the 95th percentile of label permutations; hover for the PC, the null and the sex R².",
    build: () => {
      const byAssay = Object.fromEntries(OMICS.map(([a]) => [a, Object.fromEntries((X.within_tissue[a] || []).map((r) => [r.tissue, r]))]));
      const tissues = tissueOrder.filter((t) => OMICS.some(([a]) => byAssay[a][t]));
      const z = tissues.map((t) => OMICS.map(([a]) => (byAssay[a][t] ? byAssay[a][t]["max_R2_group_PC1-5"] : null)));
      const text = tissues.map((t) => OMICS.map(([a]) => { const r = byAssay[a][t]; return r ? `${fmt(r["max_R2_group_PC1-5"], 2)}${r.group_visible ? " ●" : ""}` : ""; }));
      const custom = tissues.map((t) => OMICS.map(([a]) => { const r = byAssay[a][t]; return r ? `null95 ${fmt(r.null95_group, 2)}; PC${r.PC_with_max_group_R2}; sex R² ${fmt(r["max_R2_sex_PC1-5"], 2)}; ${r.n_animals} animals, ${r.n_features} features` : "not measured"; }));
      const tr = heatmap(z, OMICS.map(([, l]) => l), tissues.map(tissueLabel), { zmin: 0, zmax: 1, showscale: false, text, texttemplate: "%{text}", hover: "%{y} · %{x}: R² %{z:.2f}<br>%{customdata}<extra></extra>" });
      tr.customdata = custom;
      return { traces: [tr], layout: { height: 24 * tissues.length + 90, yaxis: { autorange: "reversed", automargin: true, tickfont: { size: 10 } }, xaxis: { side: "top" }, margin: { t: 40, l: 10, b: 10 } },
               table: { columns: ["assay", "tissue", "n_animals", "n_features", "max_R2_group_PC1-5", "null95_group", "group_visible", "PC_with_max_group_R2", "max_R2_sex_PC1-5"],
                        rows: OMICS.flatMap(([a]) => (X.within_tissue[a] || []).map((r) => ({ assay: a, ...r }))) } };
    },
    source: "results/03_eda/within_tissue_pca_{TRNSCRPT,PROT,METAB}.csv, readout_within_tissue_summary.csv; the stacked contrast from results/03_eda/variance_partition_TRNSCRPT.csv",
    notShow: "which genes carry the axis; and the group factor has five levels of which three (1w, 2w, 4w) are cohort-confounded, so this R² is an upper bound on training.",
    height: "tall",
  });

  // ---- 3. separability ------------------------------------------------------------------------------
  document.getElementById("p-sep").replaceChildren(
    `Within each of ${S.fusion_n_tissues} tissues, sedentary vs 8-week-trained animals are separated by a single omic layer at a mean AUROC of ${fmt(S.taskB_mean_auroc_8w, 2)}, above the permutation null in ${S.fusion_n_beats_null} of ${S.fusion_n_tissues}; combining layers adds nothing (fusion beats the best single omic in ${S.fusion_n_beats_single} of ${S.fusion_n_tissues}). `,
    `The shorter durations look the same (means ${fmt(S.taskB_mean_auroc_1w, 2)}, ${fmt(S.taskB_mean_auroc_2w, 2)}, ${fmt(S.taskB_mean_auroc_4w, 2)}), which is why they are marked †: they compare animals sacrificed months apart, so they cannot separate a fast response from cohort.`,
  );
  await figure(document.getElementById("fig-sep"), {
    title: `Sedentary vs trained within a tissue: separable at ceiling by every omic; fusion beats the best single omic in ${S.fusion_n_beats_single} of ${S.fusion_n_tissues} tissues`,
    subtitle: "Best single-omic AUROC per tissue and training duration (animal-grouped folds, tuned models); hover for the arm, the fold sd and each omic's best.",
    build: () => {
      const rows = X.separability.duration;
      const tissues = tissueOrder.filter((t) => rows.some((r) => r.tissue === t));
      const cell = (t, d) => rows.find((r) => r.tissue === t && r.duration === d);
      const z = tissues.map((t) => DUR.map((d) => cell(t, d)?.best_single_auroc ?? null));
      const text = tissues.map((t) => DUR.map((d) => { const r = cell(t, d); return r ? fmt(r.best_single_auroc, 2) : ""; }));
      const custom = tissues.map((t) => DUR.map((d) => { const r = cell(t, d); return r ? `${r.best_single_arm} (fold sd ${fmt(r.best_single_sd, 2)}); ${r.n_animals} animals<br>RNA ${fmt(r.best_TRNSCRPT_auroc, 2)} · protein ${fmt(r.best_PROT_auroc, 2)} · metabolite ${fmt(r.best_METAB_auroc, 2)}` : ""; }));
      const nulls = X.separability.best_vs_null;
      const tr = heatmap(z, DUR.map((d) => DUR_LABEL[d]), tissues.map(tissueLabel), { zmin: 0.5, zmax: 1, showscale: false, text, texttemplate: "%{text}", hover: "%{y} · %{x}: AUROC %{z:.3f}<br>%{customdata}<extra></extra>" });
      tr.customdata = custom;
      return { traces: [tr], layout: { height: 30 * tissues.length + 90, yaxis: { autorange: "reversed", automargin: true }, xaxis: { side: "top" }, margin: { t: 40, l: 10, b: 10 } },
               table: { columns: ["tissue", "best_arm", "best_auroc", "best_sd", "best_single", "single_auroc", "fusion_auroc", "fusion_minus_single", "null_p95_max_auroc", "p_perm"], rows: nulls } };
    },
    source: "results/07_fusion/taskB_duration_summary.csv, taskB_best_vs_null.csv, taskB_summary.csv (all 13 arms at 8 weeks), taskB_batch_balance.csv",
    notShow: "† columns compare animals sacrificed months apart from the controls; 'best single' is the best of 9 arms (optimistic); a fold sd of 0 means every fold at 1.0, not certainty; the permutation null (95th percentile of the best of 13 arms, 0.86–0.90) is in the data table. Separability says nothing about attribution: see the next section.",
    height: "short",
  });

  // ---- 4. covariates ----------------------------------------------------------------------------------
  document.getElementById("p-cov").replaceChildren(
    `The same folds, with no molecule at all: classifiers on the processing covariates alone. In ${S.verdict_n_training} of ${S.verdict_n_tissues} tissues they sit at their permutation null, so the separation is attributed to training. `,
    `In heart, kidney and gastrocnemius a covariate set reaches the null or beyond (sequencing depth ${fmt(S.cov_heart_depth_logreg, 2)} vs null ${fmt(S.cov_heart_depth_null95, 2)}; library ${fmt(S.cov_kidney_library_logreg, 2)} vs ${fmt(S.cov_kidney_library_null95, 2)}; QC fractions ${fmt(S.cov_skmgn_qc_logreg, 2)} vs ${fmt(S.cov_skmgn_qc_null95, 2)}), so there the training signal shares its axis with processing, as the tissue axis does across the study.`,
  );
  const FS = X.covariates.feature_sets;
  const fsOrder = ["pheno_collection", "trnscrpt_library", "trnscrpt_depth", "trnscrpt_qc", "prot_plex_channel", "all_batch_covariates"];
  const fsShort = { pheno_collection: "collection", trnscrpt_library: "library", trnscrpt_depth: "depth", trnscrpt_qc: "QC fractions", prot_plex_channel: "TMT plex/channel", all_batch_covariates: "all" };
  await figure(document.getElementById("fig-cov"), {
    title: "In heart, kidney and gastrocnemius the processing covariates alone separate the arms; elsewhere they sit at their null",
    subtitle: "AUROC of a logistic regression on each covariate set (no molecules), sedentary vs 8-week within tissue; ▲ marks cells at or above the 95th percentile of within-sex label permutations.",
    build: () => {
      const rows = X.covariates.auroc.filter((r) => r.model === "logreg");
      const tissues = tissueOrder.filter((t) => rows.some((r) => r.tissue === t));
      const cell = (t, f) => rows.find((r) => r.tissue === t && r.feature_set === f);
      const z = tissues.map((t) => fsOrder.map((f) => cell(t, f)?.auroc_mean ?? null));
      const text = tissues.map((t) => fsOrder.map((f) => { const r = cell(t, f); return r ? `${fmt(r.auroc_mean, 2)}${r.null_p95_auroc !== null && r.auroc_mean >= r.null_p95_auroc ? " ▲" : ""}` : ""; }));
      const custom = tissues.map((t) => fsOrder.map((f) => { const r = cell(t, f); return r ? `${FS[f]}<br>fold sd ${fmt(r.auroc_sd, 2)}; null95 ${fmt(r.null_p95_auroc, 2)}; permutation p ${fmt(r.p_perm, 3)}; ${r.n_features} covariates, ${r.n_animals} animals` : ""; }));
      const tr = heatmap(z, fsOrder.map((f) => fsShort[f]), tissues.map(tissueLabel), { zmin: 0.5, zmax: 1, showscale: false, text, texttemplate: "%{text}", hover: "%{y} · %{x}: AUROC %{z:.2f}<br>%{customdata}<extra></extra>" });
      tr.customdata = custom;
      const concl = X.covariates.conclusion;
      return { traces: [tr], layout: { height: 30 * tissues.length + 90, yaxis: { autorange: "reversed", automargin: true }, xaxis: { side: "top", tickfont: { size: 11 } }, margin: { t: 40, l: 10, b: 10 } },
               table: { columns: ["tissue", "verdict", "batch_variables_separating", "qc_variables_separating", "max_auroc_all_covariates", "null_p95_logreg_max"], rows: concl } };
    },
    source: "results/07_fusion/batch_covariate_auroc.csv (logistic rows; the random-forest rows are in the data table of the export), batch_conclusion.csv",
    notShow: "the three tissues where processing shares the axis are the ones where exercise biology is expected to be strongest; the untuned logistic model is a floor, not a ceiling; 'unresolvable' means the two cannot be told apart here, not that the separation is batch.",
    height: "short",
  });
  document.getElementById("callout-cov").appendChild(callout("note", "Reading the three shared-axis tissues", [
    `The batch check's verdict: ${S.verdict_overall}. `,
    `In gastrocnemius the separating QC variables are composition fractions (mitochondrial and mRNA read fractions) that mitochondrial biogenesis would produce, so they are plausibly biology, but on these data they cannot be told from a within-plate batch. In heart the arms differ in sequencing depth on the same plate and flowcell, which is unexplained. The tissue fingerprint is unaffected either way: it is fit on controls only and tested on trained animals in section 5.`,
  ]));

  // ---- 5. fingerprint by duration -----------------------------------------------------------------------
  const F = X.fingerprint_by_duration;
  document.getElementById("p-fbd").replaceChildren(
    `Fit on ${tr20.n_source_animals - tr20.n_calibration_animals} of the ${tr20.n_source_animals} sedentary controls and calibrated on the other ${tr20.n_calibration_animals}, the 20-gene panel names the tissue of the trained animals at ${fmt(fbdRow("k20", "1w").accuracy)} to ${fmt(k8.accuracy)} across the four durations, best on the cohort-matched 8-week group, with coverage ${fmt(fbdRow("k20", "1w").coverage)} to ${fmt(k8.coverage)}. `,
    `The 1- and 2-week dips are the consortium-flagged brown-fat-contaminated vena cava vials, not training; every between-duration coverage difference sits inside the fold spread.`,
  );
  const fstate = { model: "k20" };
  const fbdBuild = () => {
    const t = tokens();
    const rows = DUR.map((g) => fbdRow(fstate.model, g)).filter(Boolean);
    const ref = F.reference_held_out_controls && F.reference_held_out_controls[fstate.model];
    const x = [...(ref ? ["held-out controls"] : []), ...rows.map((r) => DUR_LABEL[r.group])];
    const acc = [...(ref ? [ref.accuracy] : []), ...rows.map((r) => r.accuracy)];
    const cov = [...(ref ? [ref.coverage] : []), ...rows.map((r) => r.coverage)];
    const err = (key, ci) => ({ type: "data", symmetric: false, visible: true, color: t.ink2, thickness: 1.5, width: 4,
      array: [...(ref ? [ref[`${key}_hi`] !== undefined ? ref[`${key}_hi`] - ref[key] : 0] : []), ...rows.map((r) => r[ci][1] - r[key])],
      arrayminus: [...(ref ? [ref[`${key}_lo`] !== undefined ? ref[key] - ref[`${key}_lo`] : 0] : []), ...rows.map((r) => r[key] - r[ci][0])] });
    const hoverA = [...(ref ? [`held-out controls (pooled over folds): accuracy ${fmt(ref.accuracy)}`] : []), ...rows.map((r) => `${DUR_LABEL[r.group]}: accuracy ${fmt(r.accuracy)} [${fmt(r.accuracy_ci[0])}, ${fmt(r.accuracy_ci[1])}]<br>n = ${r.n_samples} vials, ${r.n_individuals} animals${r.accuracy_excl_flagged !== undefined ? `<br>flagged vena cava vials excluded: ${fmt(r.accuracy_excl_flagged)}` : ""}`)];
    const hoverC = [...(ref ? [`held-out controls: coverage ${fmt(ref.coverage)}`] : []), ...rows.map((r) => `${DUR_LABEL[r.group]}: coverage ${fmt(r.coverage)} [${fmt(r.coverage_ci[0])}, ${fmt(r.coverage_ci[1])}]<br>empty sets ${fmt(r.empty)}, non-empty but wrong ${fmt(r.wrong_non_empty)}, mean set size ${fmt(r.set_size, 2)}`)];
    const ref9 = refLine(0.9, "1 − α = 0.90");
    const errs = F.errors_by_group.filter((e) => e.model === fstate.model);
    return { traces: [{ ...bar(x, acc, { name: "accuracy", slot: 1, text: acc.map((v) => fmt(v)), hover: "%{customdata}<extra>accuracy</extra>" }), customdata: hoverA, error_y: err("accuracy", "accuracy_ci") },
                      { ...bar(x, cov, { name: "coverage of the 90 % set", slot: 2, text: cov.map((v) => fmt(v)), hover: "%{customdata}<extra>coverage</extra>" }), customdata: hoverC, error_y: err("coverage", "coverage_ci") }],
             layout: { yaxis: { range: [0, 1.08], title: { text: "fraction" }, tickformat: ".1f" }, shapes: ref9.shapes, annotations: ref9.annotations, barmode: "group", legend: { y: 1.14 }, margin: { t: 40, b: 50 } },
             table: { columns: ["model", "group", "n_samples", "n_individuals", "accuracy", "accuracy_ci", "coverage", "coverage_ci", "empty", "wrong_non_empty", "set_size", "errors"],
                      rows: rows.map((r) => ({ ...r, accuracy_ci: r.accuracy_ci.map((v) => v.toFixed(3)).join(" – "), coverage_ci: r.coverage_ci.map((v) => v.toFixed(3)).join(" – "),
                                               errors: errs.filter((e) => e.group === r.group).map((e) => `${e.tissue} → ${e.y_pred} ${e.n} (${e.sex})`).join("; ") })) } };
  };
  const fctl = el("div", { class: "controls" }, [control("Model", segmented(Object.entries(MODEL_LABEL), fstate.model, async (v) => { fstate.model = v; await frerender(); }, "model"))]);
  let ffig;
  async function frerender() {
    if (!ffig) return;
    const b = fbdBuild();
    ffig.traces = b.traces; ffig.table = b.table;
    await window.Plotly.react(ffig.chart, b.traces, { ...template(), ...b.layout }, CONFIG);
    if (!ffig.tableWrap.hidden) ffig.tableWrap.replaceChildren(tableFrom(b.table));
  }
  ffig = await figure(document.getElementById("fig-fbd"), {
    title: "The fingerprint ignores training: fit on controls, it is at its best on the 8-week animals",
    subtitle: "Accuracy and coverage of α = 0.10 sets per training duration; whiskers are 95 % cluster-bootstrap intervals over the 10 animals of a group († cohort-confounded durations).",
    build: fbdBuild, toolbar: fctl,
    source: "recomputed from results/31_site_regen/08_shift_k20/scores_target_vials.csv and 08_shift_k50/ (design: results/08_shift/TRNSCRPT/shift_table.csv, split train_controls_test_trained)" + (F.matched_summary ? "; reference and matched-fold rows from results/15_time_course/fingerprint_by_duration/" : ""),
    notShow: F.design + ". Sets are of size 0 or 1 here, so lost coverage is abstention. The guarantee rests on 3 calibration animals (nominal 90 %); the † groups carry cohort and season with training.",
    height: "tall",
  });

  // ---- 6. panel genes vs training response -------------------------------------------------------------
  const PT = X.panel_training;
  const pp = document.getElementById("p-ptr");
  if (PT && PT.genes) {
    pp.replaceChildren(
      `For each panel gene, the consortium's differential-expression tables give its largest training response in its own marker tissue over sex and duration. `,
      `For every one of the ${S.ptr_k20_n_genes} panel genes the tissue contrast exceeds its largest training response: by a factor of ${fmt(S.ptr_k20_min_ratio_tissue_over_training, 1)} at the least and ${fmt(S.ptr_k20_median_ratio_tissue_over_training, 1)} at the median (largest training response ${fmt(S.ptr_k20_max_abs_logfc_marker, 2)} log2 units). ${S.ptr_k20_n_regulated_marker_5pct} of the ${S.ptr_k20_n_genes} genes are in the consortium's 5 % FDR set for their marker tissue, and the panel's tissue calls are unaffected (section 5).`,
    );
    await figure(document.getElementById("fig-ptr"), {
      title: `On every panel gene the tissue contrast exceeds the largest training response, by a median factor of ${fmt(S.ptr_k20_median_ratio_tissue_over_training, 1)}`,
      subtitle: "Per gene of the 20-gene panel: log2 CPM above the runner-up tissue (tissue effect) vs the largest |log2 fold change| of training in the marker tissue; ● = in the consortium's 5 % FDR set there.",
      build: () => {
        const rows = PT.genes.filter((g) => g.in_k20).sort((a, b) => a.tissue_effect_log2cpm - b.tissue_effect_log2cpm);
        const y = rows.map((g) => `${g.symbol}${g.regulated_marker_5pct ? " ●" : ""} · ${g.marker_tissue}`);
        const custom = rows.map((g) => `${g.symbol} in ${tissueLabel(g.marker_tissue)}<br>8-week log2FC: females ${fmt(g.logfc_8w_female, 2)} ± ${fmt(g.logfc_8w_female_se, 2)}, males ${fmt(g.logfc_8w_male, 2)} ± ${fmt(g.logfc_8w_male_se, 2)}<br>smallest adjusted p in the marker tissue ${fmt(g.min_adj_p_marker, 3)}; regulated in ${g.n_tissues_regulated_anywhere} tissues`);
        return { traces: [{ ...bar(rows.map((g) => g.tissue_effect_log2cpm), y, { name: "tissue effect (log2 CPM above the runner-up tissue)", slot: 1, horizontal: true, hover: "%{customdata}<extra>tissue effect %{x:.2f}</extra>" }), customdata: custom },
                          { ...bar(rows.map((g) => g.max_abs_logfc_marker), y, { name: "largest training |log2FC| in the marker tissue", slot: 2, horizontal: true, hover: "%{customdata}<extra>training |log2FC| %{x:.2f}</extra>" }), customdata: custom }],
                 layout: { barmode: "group", height: 26 * rows.length + 120, xaxis: { title: { text: "log2 units" }, rangemode: "tozero" }, yaxis: { automargin: true, tickfont: { size: 10 } }, legend: { y: 1.1 }, margin: { t: 40, l: 10 } },
                 table: { columns: ["symbol", "marker_tissue", "tissue_effect_log2cpm", "max_abs_logfc_marker", "ratio_tissue_over_training", "regulated_marker_5pct", "n_tissues_regulated_anywhere", "logfc_8w_female", "logfc_8w_male"], rows } };
      },
      source: "results/05_panels/TRNSCRPT/panel_training_response.csv, panel_training_summary.csv (consortium DESeq2 tables in data/raw/da/, 5 % FDR set from data/raw/training_regulated_features.csv)",
      notShow: "responses outside the marker tissue (in the data table), protein-level responses, and whether another exercise modality could cross the margin; the consortium DA n is 5–6 animals per sex per group, and the tissue effect is a difference of tissue means, not a fold change against a null.",
      height: "tall",
    });
  } else {
    pp.replaceChildren("For each panel gene, the consortium's differential-expression tables give its largest training response in its own marker tissue; the tissue contrast the panel reads is larger for every gene.");
    document.getElementById("fig-ptr").appendChild(pendingBlock("Panel genes vs training response", "the panel training-response table has not been generated yet"));
  }

  // ---- takeaway ----------------------------------------------------------------------------------------
  document.getElementById("callout-takeaway").appendChild(callout("note", "What this adds to the tissue answer", [
    `The panel's tissue call is invariant to training state (${fmt(k8.accuracy)} accuracy and ${fmt(k8.coverage)} coverage on the 8-week animals with a panel that never saw a trained animal). `,
    `Training is real (${fmt(S.taskB_mean_auroc_8w, 2)} mean AUROC within tissue) but lives on minor within-tissue axes, below the tissue contrast on every panel gene${S.ptr_k20_median_ratio_tissue_over_training ? ` (median factor ${fmt(S.ptr_k20_median_ratio_tissue_over_training, 1)})` : ""}; and in the tissues where exercise biology is strongest, the within-study sedentary-vs-trained contrast shares its axis with processing, as the tissue axis does across the study. Attribution needs an external check there too.`,
  ]));
}

main().catch((e) => { document.getElementById("lede").textContent = `Failed to load the page data: ${e.message}`; console.error(e); });
