import { mountChrome, loadJSON, el, statTile, fmt, pct } from "../site.js";
import { mountOverview } from "../overview.js";
import { figure, bar, line, band, tokens, palette } from "../charts.js";
import { mountLadder, pick } from "../ladder.js";

async function main() {
  await mountChrome("index.html");
  const [H, PC, N, BM] = await Promise.all([loadJSON("data/headline.json"), loadJSON("data/panel_curve.json"), loadJSON("data/nesting.json"), loadJSON("data/bodymap.json")]);
  const ex = H.extras;
  const tile = Object.fromEntries(H.tiles.map((t) => [t.id, t]));
  const idk = pick(H.ladder, "in_distribution", "k20", "marginal", "pooled");
  const opa = pick(H.ladder, "in_distribution", "k20", "marginal", "one_per_animal");
  const tr = pick(H.ladder, "train_control_test_trained", "k20", "marginal");
  const mf = pick(H.ladder, "train_male_test_female", "k20", "marginal");
  const bm = pick(H.ladder, "different_lab", "k20", "marginal");
  const gt = pick(H.ladder, "different_species", "k20", "marginal");
  // BodyMap organs with a MoTrPAC counterpart: organ_map is organ → rat tissues, or null
  const nOrgans = Object.keys(BM.organ_map).length, nMapped = Object.values(BM.organ_map).filter(Boolean).length;
  const bmAcc = tile.tile_bodymap_k20.value;
  const perSet = (v) => `${fmt(v, 2)} tissue${fmt(v, 2) === "1.00" ? "" : "s"} per set`;

  // lede
  document.getElementById("lede").replaceChildren(
    el("b", {}, "Yes."),
    ` A ${H.design.k_panel}-gene panel identifies ${ex.n_tissues} rat tissues at ${fmt(tile.tile_acc_k20.value)} balanced accuracy and, fit on all animals, names ${bmAcc === 1 ? "every mapped adult organ" : fmt(bmAcc) + " of the mapped adult organs"} in another laboratory's rats (${fmt(bmAcc)}, ${nMapped} of ${nOrgans} organs). `,
    "Its 90 % guarantee ", el("a", { href: "#sec-guarantee" }, "holds within the study"), ", abstains rather than guesses under shift, and is restored beyond the study by recalibrating on three animals ",
    `(${fmt(bm?.recal_n3)}, within species).`,
  );

  // tiles
  // ---- overview: what we did, in one picture (numbers from headline.json only) --------------------------------
  {
    const d = H.design;
    const tr20 = pick(H.ladder, "train_control_test_trained", "k20", "marginal");
    const bm20 = pick(H.ladder, "different_lab", "k20", "marginal");
    const gt20 = pick(H.ladder, "different_species", "k20", "marginal");
    const bridge = typeof ex.bridge_sum_ratio_all_genes_pool99 === "number" ? ex.bridge_sum_ratio_all_genes_pool99 : null;
    document.getElementById("overview-summary").replaceChildren(
      `MoTrPAC's rat endurance-training study is the asset here: ${ex.n_tissues} tissues from the same ${ex.n_animals} animals, ${ex.n_vials} RNA-seq vials, a training time course, several omic layers, and reference-standard RNA pools on every extraction plate. `,
      `Inside animal-grouped folds we selected a ${d.k_panel}-gene panel, gave it split-conformal sets that promise the true tissue ${pct(1 - d.alpha)} of the time, carried both up a ladder of shifts (trained animals, the other sex, another laboratory's rats, human GTEx), and measured directly, on the consortium's bridging standards, how much processing contributes to the within-study signal.`,
    );
    const nodes = [
      { label: `MoTrPAC RNA-seq, ${ex.n_tissues} tissues`, num: String(ex.n_vials), unit: "vials" },
      { label: "Animal-grouped folds, whole animals held out", num: String(d.n_outer_folds), unit: "folds" },
      { label: "Panel selected inside each fold", num: fmt(H.accuracy.k20.mean), unit: `balanced accuracy, ${d.k_panel} genes` },
      { label: `${pct(1 - d.alpha)} sets calibrated on held-out animals`, num: String(d.n_cal_animals), unit: "calibration animals" },
      { label: "Shift ladder: state, sex, laboratory, species", num: fmt(bm20.accuracy), unit: "mapped organs named, other lab" },
      bridge !== null ? { label: "Batch measured on bridging standards", num: pct(bridge, 1), unit: "of tissue-separating variance" }
                      : { label: "Batch nested in tissue, by design", num: `${tile.tile_estimable.value} of ${tile.tile_estimable.total}`, unit: "pairs contrastable in a batch" },
      { label: "Verdict: sufficient, transferable, calibratable", num: null, unit: null, verdict: true },
    ];
    const alt = `The study in one picture: MoTrPAC RNA-seq of ${ex.n_tissues} tissues (${ex.n_vials} vials) → ${d.n_outer_folds} animal-grouped folds → a ${d.k_panel}-gene panel selected inside each fold (${fmt(H.accuracy.k20.mean)} balanced accuracy) → ${pct(1 - d.alpha)} conformal sets calibrated on ${d.n_cal_animals} held-out animals → a shift ladder (training state, sex, laboratory, species; ${fmt(bm20.accuracy)} of mapped organs named in another laboratory) → batch measured directly on bridging standards${bridge !== null ? ` (${pct(bridge, 1)} of tissue-separating variance)` : ""} → verdict: sufficient, transferable, calibratable.`;
    mountOverview(document.getElementById("overview-svg"), nodes, { alt, desc: nodes.map((n) => n.label + (n.num ? ` ${n.num} ${n.unit || ""}` : "")).join("; ") });
    document.getElementById("overview-caption").replaceChildren(
      `Coverage of the ${pct(1 - d.alpha)} sets: ${fmt(tr20.coverage)} on trained animals, ${fmt(bm20.coverage)} in another laboratory, ${fmt(gt20.coverage)} in human; three target animals restore ${fmt(bm20.recal_n3)} within species.`,
    );
    const li = (parts) => el("li", {}, parts);
    document.getElementById("overview-new").replaceChildren(
      li([`A class-aware round-robin selector makes ${d.k_panel} genes sufficient for ${ex.n_tissues} tissues (${fmt(H.accuracy.k20.mean)} ± ${fmt(H.accuracy.k20.sd)}); a univariate F-test at the same size reaches ${fmt(ex.acc_fclassif_k20)}.`]),
      li([`The ${pct(1 - d.alpha)} guarantee, tested under shift: it holds within the study, abstains rather than guesses beyond it, and three animals restore it within species (${fmt(bm20.recal_n3)}).`]),
      li([bridge !== null ? `Batch was measured directly on the consortium's bridging standards (${pct(bridge, 1)} of tissue-separating variance) and checked against an independent laboratory (${fmt(bm20.accuracy)} of mapped organs named).`
                          : `The processing design was audited layer by layer and the fingerprint checked against an independent laboratory (${fmt(bm20.accuracy)} of mapped organs named).`]),
    );
    document.getElementById("overview-why").replaceChildren(
      li([el("b", {}, "For the challenge question: "), "a compact, interpretable signature is sufficient, and “reliably” is measurable: a calibrated guarantee whose behaviour under shift is known and whose repair costs three animals."]),
      li([el("b", {}, "For MoTrPAC analysts: "), `each tissue was processed as a unit, as in every multi-tissue design (${tile.tile_estimable.value} of ${tile.tile_estimable.total} pairs contrastable within a batch), so pair within-study accuracy with an external check. Training state, the contrast the design supports directly, `, el("a", { href: "exercise.html" }, "has its own page"), `.`]),
      li([el("b", {}, "For anyone deploying a signature: "), `calibrate on the target. Three animals restore ${fmt(bm20.recal_n3)} coverage at ${fmt(bm20.recal_n3_size, 2)} tissue per set (within species); a zero-error certificate at α = δ = ${fmt(d.alpha, 2)} needs ${ex.n_zero_error_1010} animals.`]),
    );
  }

  const tiles = document.getElementById("tiles");
  for (const t of H.tiles) tiles.appendChild(statTile(t));

  // ladder with controls (shared with the Transfer page), and the pointer to the Exercise page
  await mountLadder(document.getElementById("fig-ladder"), H, { full: false });
  document.getElementById("p-exercise").replaceChildren(
    "Training state, the one contrast this design supports directly, has its own page → ",
    el("a", { href: "exercise.html" }, "Exercise"),
    ` (${fmt(tr?.accuracy)} accuracy, ${fmt(tr?.coverage)} coverage, panel fit on sedentary controls only).`,
  );

  // section 1: accuracy paragraph + panel curve mini
  document.getElementById("p-accuracy").replaceChildren(
    `Under five animal-grouped folds the round-robin selector reaches ${fmt(H.accuracy.k20.mean)} ± ${fmt(H.accuracy.k20.sd)} balanced accuracy at 20 genes, `,
    `${fmt(H.accuracy.k50.mean)} at 50 and ${fmt(H.accuracy.full.mean)} with all genes. `,
    `A univariate F-test at the same size reaches ${fmt(ex.acc_fclassif_k20)}: the selector, not the classifier, is the result.`,
  );
  await figure(document.getElementById("fig-curve"), {
    title: `The curve saturates by 15–20 genes with a class-aware selector; a univariate F-test at the same size reaches ${fmt(ex.acc_fclassif_k20)}`,
    subtitle: "Mean ± sd balanced accuracy over 5 animal-grouped folds vs panel size (log scale); logreg_l2 on the selected genes.",
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
    notShow: "which genes: the selection changes per fold (see the Fingerprint page for the stable core and the Explorer's panel builder).",
  });

  // section 2: guarantee paragraph + empty-set stack
  document.getElementById("p-guarantee").replaceChildren(
    `The 20-gene sets cover ${fmt(idk?.coverage)} of held-out vials (pooled calibration; ${fmt(opa?.coverage)} with one vial per animal) and ${fmt(tr?.coverage)} on trained animals when the panel is fit on sedentary controls only. `,
    `Beyond the study they abstain rather than guess — coverage ${fmt(bm?.coverage)} in another laboratory and ${fmt(gt?.coverage)} in human, wrong non-empty sets ${fmt(bm?.wrong_non_empty)} and ${fmt(gt?.wrong_non_empty)} — and three target animals restore ${fmt(bm?.recal_n3)} within species at ${perSet(bm?.recal_n3_size)}.`,
  );
  await figure(document.getElementById("fig-empty"), {
    title: "Under shift the sets abstain rather than guess: the shortfall is empty sets; wrong confident sets stay rare",
    subtitle: "Fraction of test samples whose α = 0.10 set (20-gene panel, marginal, source-calibrated) holds the true tissue, is wrong, or is empty.",
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

  // section 3: identifiability paragraph + nesting mini
  const est = N.estimable_pairs.find((r) => r.assay === "TRNSCRPT");
  const nest = Object.fromEntries(N.nesting.TRNSCRPT.map((r) => [r.variable, r]));
  const pairText = est.estimable_pairs.split(";").map((p) => p.split("|").map((t) => t.toLowerCase()).join(" and ")).join("; ");
  const facts = [`${fmt(bmAcc)} of mapped adult organs named in a laboratory with none of these batches`];
  if (typeof ex.bridge_sum_ratio_all_genes_pool99 === "number") {
    facts.push(`batch measured at ${(100 * ex.bridge_sum_ratio_all_genes_pool99).toFixed(1)} % of tissue-separating variance on bridging pools run on ${ex.bridge_n_plates_pool99} plates at both sites`);
  }
  document.getElementById("p-identifiability").replaceChildren(
    `As in any multi-tissue design, each tissue was processed as a unit: ${ex.n_plates} plates, ${ex.n_lib_batches} library batches and ${ex.n_flowcells} flowcells hold whole tissues, ${est.n_pairs_estimable} of ${est.n_pairs_total} pairs ${est.n_pairs_estimable === 1 ? "is" : "are"} contrastable inside a batch, and library QC numbers alone reach ${fmt(ex.qc_all)}. `,
    `${facts.length === 2 ? "Two external facts make" : "One external fact makes"} the fingerprint credible as biology: ${facts.join(", and ")}.`,
  );
  await figure(document.getElementById("fig-nesting"), {
    title: `As in any multi-tissue design each tissue was processed as a unit: ${est.n_pairs_estimable} of ${est.n_pairs_total} pairs contrastable inside a batch (${pairText})`,
    subtitle: `Of the ${est.n_pairs_total} RNA-seq tissue pairs, those sharing a level of each processing variable, and of all three.`,
    build: () => {
      const vars = [["RNA_extr_plate_ID", "RNA extraction plate"], ["Lib_batch_ID", "library batch"], ["Seq_flowcell_ID", "flowcell"]];
      const y = [...vars.map(([, l]) => l), `all three (estimable: ${pairText})`];
      const x = [...vars.map(([v]) => nest[v].n_pairs_sharing_level), est.n_pairs_estimable];
      return {
        traces: [bar(x, y, { horizontal: true, name: "tissue pairs sharing a level", slot: 1, text: x.map((v) => `${v} of ${est.n_pairs_total}`), textposition: "outside", hover: "%{y}: %{x} of " + est.n_pairs_total + " pairs<extra></extra>" })],
        layout: { xaxis: { range: [0, est.n_pairs_total], title: { text: `tissue pairs (of ${est.n_pairs_total})` } }, yaxis: { autorange: "reversed", automargin: true }, margin: { t: 20, l: 10 }, showlegend: false },
        table: { columns: ["variable", "n_levels", "max_tissues_per_level", "n_levels_shared", "tissues_in_one_level", "cramers_v", "n_pairs_sharing_level", "n_pairs_total"], rows: N.nesting.TRNSCRPT },
      };
    },
    source: "results/16_identifiability/nesting_TRNSCRPT.csv, results/16_identifiability/estimable_pairs.csv",
    notShow: "the other omic layers (Identifiability page): TMT proteomics is nested by construction, immunoassay plates are the one layer that mixes tissues.",
    height: "short",
  });
}

main().catch((e) => { console.error(e); document.getElementById("lede").textContent = "Failed to load site data: " + e.message; });
