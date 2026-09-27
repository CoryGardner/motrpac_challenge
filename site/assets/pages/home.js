import { mountChrome, loadJSON, el, statTile, fmt } from "../site.js";
import { figure, bar, line, band, tokens, palette } from "../charts.js";
import { mountLadder, pick } from "../ladder.js";

async function main() {
  await mountChrome("index.html");
  const [H, PC, N, Q] = await Promise.all([loadJSON("data/headline.json"), loadJSON("data/panel_curve.json"), loadJSON("data/nesting.json"), loadJSON("data/qc_baseline.json")]);
  const ex = H.extras;
  const tile = Object.fromEntries(H.tiles.map((t) => [t.id, t]));

  // lede
  document.getElementById("lede").replaceChildren(
    el("b", {}, "Yes"),
    ` — a 20-gene panel identifies 19 rat tissues at ${fmt(tile.tile_acc_k20.value)} balanced accuracy and, fit on all animals, names ${fmt(tile.tile_bodymap_k20.value, 3)} of the mapped adult organs in another laboratory's rats (9 of 11 organs, muscle and brain scored as super-classes). `,
    "But “reliably” has three parts, and only the first survives on its own: the ",
    el("a", { href: "#sec-guarantee" }, "coverage guarantee"), " does not travel, and within one study the ",
    el("a", { href: "#sec-identifiability" }, "tissue axis cannot be separated from processing"), ".",
  );

  // tiles
  const tiles = document.getElementById("tiles");
  for (const t of H.tiles) tiles.appendChild(statTile(t));

  // ladder with controls (shared with the Transfer page)
  await mountLadder(document.getElementById("fig-ladder"), H, { full: false });

  // section 1: accuracy paragraph + panel curve mini
  document.getElementById("p-accuracy").replaceChildren(
    `Under five animal-grouped folds the round-robin selector reaches ${fmt(H.accuracy.k20.mean)} ± ${fmt(H.accuracy.k20.sd)} balanced accuracy at 20 genes, `,
    `${fmt(H.accuracy.k50.mean)} at 50 and ${fmt(H.accuracy.full.mean)} with all genes. The selector, not the classifier, was the hard part: `,
    `a univariate F-test at the same k picks markers of the same easy tissues and reaches ${fmt(ex.acc_fclassif_k20)}.`,
  );
  await figure(document.getElementById("fig-curve"), {
    title: "The panel curve saturates by 15–20 genes with a class-aware selector; the F-test never gets there",
    subtitle: "Mean balanced accuracy ± sd over 5 animal-grouped folds vs panel size (log scale), round-robin vs F-test selection, logreg_l2 on the selected genes.",
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
  const idk = pick(H.ladder, "in_distribution", "k20", "marginal", "pooled");
  const bm = pick(H.ladder, "different_lab", "k20", "marginal");
  const gt = pick(H.ladder, "different_species", "k20", "marginal");
  const mf = pick(H.ladder, "train_male_test_female", "k20", "marginal");
  document.getElementById("p-guarantee").replaceChildren(
    `A split-conformal prediction set promises that the true tissue is in the set 90 % of the time. On held-out MoTrPAC animals the 20-gene sets cover ${fmt(idk?.coverage)} of vials. `,
    `With the same calibration they cover ${fmt(bm?.coverage)} of adult BodyMap organs and ${fmt(gt?.coverage)} of human GTEx samples, and the shortfall is almost entirely `,
    el("b", {}, "empty sets"), ` (${fmt(bm?.empty)} and ${fmt(gt?.empty)} of samples): the model names the tissue correctly but is less confident on another laboratory's libraries. `,
    `Three target animals repair it within species (${fmt(bm?.recal_n3)} coverage at ${fmt(bm?.recal_n3_size, 2)} tissues per set); across species the repaired sets hold ${fmt(gt?.recal_n3_size, 1)} of 19 tissues.`,
  );
  await figure(document.getElementById("fig-empty"), {
    title: "What replaces coverage is abstention, not confident error",
    subtitle: "Fraction of test samples whose α = 0.10 set (20-gene panel, marginal, source-calibrated) contains the true tissue, is non-empty but wrong, or is empty.",
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
  const qc = Object.fromEntries(Q.summary.map((r) => [r.features, r]));
  document.getElementById("p-identifiability").replaceChildren(
    `Each tissue sits entirely inside one RNA extraction plate, one library batch and one flowcell: ${ex.n_plates} plates and ${ex.n_lib_batches} library batches holding at most ${Math.max(nest.RNA_extr_plate_ID.max_tissues_per_level, nest.Lib_batch_ID.max_tissues_per_level)} tissues each, ${ex.n_flowcells} flowcells holding up to ${nest.Seq_flowcell_ID.max_tissues_per_level}. `,
    `Only ${est.n_pairs_estimable} of ${est.n_pairs_total} tissue pairs share a level of all three (${est.estimable_pairs.replace("|", " vs ").toLowerCase()}, which is also the sex contrast). `,
    `Library QC numbers alone, with no gene, classify the tissue at ${fmt(qc.all.acc_mean)} ± ${fmt(qc.all.acc_sd)} (technical numbers ${fmt(qc.technical.acc_mean)}, composition fractions ${fmt(qc.composition.acc_mean)}). `,
    "So within-study accuracy is not evidence that the signature is biology. The evidence is external: the panel transfers to a laboratory where none of these batches exist. ",
    ...(typeof ex.bridge_sum_ratio_all_genes_pool99 === "number"
      ? [`And where batch could be measured directly, on a reference RNA pool run on ${ex.bridge_n_plates_pool99} plates at both sites, it was ${(100 * ex.bridge_sum_ratio_all_genes_pool99).toFixed(1)} % of the variance that separates tissues.`]
      : []),
  );
  await figure(document.getElementById("fig-nesting"), {
    title: `${est.n_pairs_estimable} tissue pair${est.n_pairs_estimable === 1 ? "" : "s"} in ${est.n_pairs_total} can be contrasted inside a processing batch`,
    subtitle: `Number of the ${est.n_pairs_total} RNA-seq tissue pairs that share at least one level of each processing variable, and the pairs sharing a level of all three.`,
    build: () => {
      const vars = [["RNA_extr_plate_ID", "RNA extraction plate"], ["Lib_batch_ID", "library batch"], ["Seq_flowcell_ID", "flowcell"]];
      const y = [...vars.map(([, l]) => l), "all three (estimable)"];
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
