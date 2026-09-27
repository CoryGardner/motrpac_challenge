import { mountChrome, loadJSON, fmt } from "../site.js";
import { mountPipeline } from "../overview.js";

async function main() {
  await mountChrome("methods.html");
  const H = await loadJSON("data/headline.json");
  const d = H.design;
  mountPipeline(document.getElementById("fig-pipeline"), { nTrain: d.n_train_animals, nTest: d.n_test_animals, nTotal: H.extras.n_animals, prefilter: d.variance_prefilter,
                                                          cGrid: d.C_grid, nFit: d.n_fit_animals, nCal: d.n_cal_animals, nFolds: d.n_outer_folds });
  const row = (model, calib) => H.ladder.find((r) => r.rung_id === "in_distribution" && r.model === model && r.variant === "marginal" && r.calibration === calib);
  const opa = row("full", "one_per_animal"), pooled = row("full", "pooled"), opa20 = row("k20", "one_per_animal"), pooled20 = row("k20", "pooled");
  document.getElementById("p-coverage").replaceChildren(
    `On this design the α = 0.10 sets of the all-gene model cover ${fmt(opa?.coverage)} of held-out vials with one vial per animal and ${fmt(pooled?.coverage)} with pooled vials (the 20-gene panel, the model on the Home page: ${fmt(opa20?.coverage)} and ${fmt(pooled20?.coverage)} on the same design).`,
  );
}

async function multiomic() {
  const M = await loadJSON("data/multiomic.json");
  const H = M.how, S = M.scales;
  const mw = M.metabolites.find((m) => m.leg === "deep_mw");
  const jr = M.ladder_species.find((r) => r.target.startsWith("protein → Jiang 2020 (cleaned") && r.model === "k20");
  document.getElementById("p-mo-matrix").replaceChildren(
    `Proteins: MoTrPAC ran global proteomics on ${H.n_tissues_prot} of its ${H.n_tissues_rna} tissues. From the portal's peptide-level reporter-ion intensities (release ${H.release}), each TMT plex is rolled up to proteins (at least ${H.min_peptides} quantified peptides), its reference-pool channel dropped, each channel divided by its total and expressed as log2 ppm; proteins quantified in every tissue are kept (${S.n_proteins_inner} proteins, ${S.n_vials} vials, ${S.n_animals} animals). Metabolites use the distributed MoTrPAC tables.`);
  document.getElementById("p-mo-code").replaceChildren(
    "The same code as the RNA fingerprint: the round-robin selector, the tuned L2 logistic regression, animal-grouped folds and the split-conformal sets with recalibration on a few target individuals (src/tfp/models.py, conformal.py, transfer.py). The only new loader is src/tfp/rii.py, which reads the reporter-ion files.");
  document.getElementById("p-mo-transfer").replaceChildren(
    `Transfer: proteins map to the human atlas through the same 1:1 rat–human orthologs as GTEx (${jr.n_samples} Jiang 2020 samples from ${jr.n_individuals} donors map to MoTrPAC tissues); metabolites map by RefMet name (${mw.matched} names shared with the mouse aging atlas, ${mw.n_mapped} samples from ${mw.n_individuals} mice). Every feature is z-scored within its own dataset and the classifier's coefficients are applied unchanged.`);
}

main().then(multiomic).catch((e) => { console.error(e); });
