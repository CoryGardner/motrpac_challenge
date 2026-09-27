// Shared sentences built from site/data/product.json (pv_* provenance): the recalibration caveat (Methods, Limitations)
// and the composition sensitivity of within-set scaling (Limitations, Check samples). No number is typed here.
const f2 = (v) => Number(v).toFixed(2), f3 = (v) => Number(v).toFixed(3), pc = (v) => `${Math.round(100 * v)} %`;

export function recalNote(PR, design, nNeeded) {
  const alpha = design.alpha;
  const b = PR.recal_draws["bodymap.3"];
  return `The finite-sample conformal guarantee needs the calibration and test scores to be exchangeable. Recalibrating on a few target animals pools every tissue of those animals (${b.min_cal_scores}–${b.max_cal_scores} BodyMap scores from ${b.n_recal} animals per draw), so the recalibrated coverage (${f3(b.mean_coverage_all_draws)} on average, ${f2(b.min_coverage)}–${f2(b.max_coverage)} per draw, ${b["n_draws_below_0.90"]} of ${b.draws} draws below 0.90) is an observed coverage across draws, not a guarantee over new animals. A per-animal guarantee at α = ${alpha.toFixed(2)} needs at least ${nNeeded} calibration animals with one sample each (the ⌈(n + 1)(1 − α)⌉-th score must exist); with ${b.n_recal} it cannot exist. The BodyMap animal identities are themselves inferred from the replicate index of each sample name. The in-study 90 % sets keep the word guarantee: they are calibrated on ${design.n_cal_animals} held-out animals and reported both pooled and with one vial per animal.`;
}

export function compositionNote(PR) {
  const c = PR.composition;
  return `Measured on ${c.n_subsets} random subsets (${c.size_min}–${c.size_max} samples, at least three organs) of the rat BodyMap 21-week adults: with within-set scaling the coverage of the 90 % sets ranged ${pc(c["within.coverage.p05"])}–${pc(c["within.coverage.p95"])} and the accuracy ${pc(c["within.accuracy.p05"])}–${pc(c["within.accuracy.p95"])} (5th–95th percentile), depending only on which other samples were uploaded; reference scaling scores each sample alone, so its calls do not change with the upload (coverage ${pc(c["reference.coverage.p05"])}–${pc(c["reference.coverage.p95"])} across the same subsets, varying only with which samples each subset contains).`;
}

/** Which 20 genes: the panel card (all animals), the per-fold panels (the reported accuracy) and the browser model. */
export function modelNote(model, nAnimals) {
  const nCal = model.calibration.n_animals, shared = model.n_genes_shared_with_all_animal_panel, k = model.genes.length;
  const own = model.genes.filter((g) => !g.in_all_animal_panel).map((g) => g.symbol);
  return `One method, three fits of it. The balanced accuracy is for ${k}-gene panels re-selected inside each animal-grouped fold; the panel card on the Panel page is the ${k} genes selected on all ${nAnimals} MoTrPAC animals; Check samples runs the transfer model fit on ${nAnimals - nCal} animals and calibrated on the other ${nCal}, which shares ${shared} of its ${k} genes with the panel card (its own: ${own.join(", ")}). Upload the browser model's genes (listed in the template and in site/data/panel_model.json); a gene it needs but your table lacks is scored at the reference mean.`;
}
