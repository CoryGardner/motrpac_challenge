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

main().catch((e) => { console.error(e); });
