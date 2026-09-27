# Molecular Tissue Fingerprints — analysis phases, site export and checks.  `make help` lists the targets.
# Every analysis phase writes to results/ (and appends to results/REPORT.md unless MOTRPAC_NO_REPORT=1).
# The site reads results/ when a complete run is present, else the committed snapshot results_frozen/.
PY      ?= python
RAW     ?= data/raw
SYN     ?= data/synthetic
ASSAY   ?= TRNSCRPT
QUICK   ?= 0
Q       := $(if $(filter 1,$(QUICK)),--quick,)
MODELS  ?= centroid,logreg_l2,rf
CORE9   := SKM-GN,HEART,KIDNEY,LIVER,LUNG,BAT,WAT-SC,HIPPOC,PLASMA
METAB_FLAGS := --drop-incomplete-samples 0.2 --complete-features
GTEX_DIR   := data/external/gtex
GTEX_TPM   := $(GTEX_DIR)/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz
GTEX_READS := $(GTEX_DIR)/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz
GTEX_ATTRS := $(GTEX_DIR)/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt
REGEN      := results/31_site_regen
FROZEN     := results_frozen
RESULTS    ?=
MOTRPAC_PORTAL ?=
SITE_PORT  ?= 8765
export PYTHONPATH := src
export MOTRPAC_PORTAL
SMOKE_ENV := MOTRPAC_RAW=$(SYN) MOTRPAC_RESULTS=results_smoke MOTRPAC_NO_REPORT=1

.PHONY: help env synthetic smoke export check inventory eda eda-readout baselines prot-diagnostic panels annotate-panel panel-training conformal fusion batch-check shift time-course bodymap gtex transfer external test all clean-synthetic clean-results \
        identifiability portal-check regen-scores panel-model freeze-results verify-frozen site-data site-test site screenshots linkcheck summary-figure home-figure figures brand

help: ## list the targets
	@grep -E '^[a-zA-Z0-9_-]+:.*?## ' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "  %-18s %s\n", $$1, $$2}'
	@echo "variables: QUICK=1 MODELS=$(MODELS) RAW=$(RAW) RESULTS=<dir> MOTRPAC_PORTAL=<portal download root>"

env: ## install the Python dependencies (pinned) into the active environment
	$(PY) -m pip install -r requirements.txt

# ---- data ---------------------------------------------------------------------------------------------
synthetic: ## write a small synthetic data set (refuses to overwrite real data)
	@if [ -f $(RAW)/pheno.csv ] && [ ! -f $(RAW)/_SYNTHETIC ]; then echo "$(RAW) holds real data; use RAW=$(SYN)"; exit 1; fi
	$(PY) scripts/make_synthetic_data.py --out $(RAW)

export: ## export the MoTrPAC R package to data/raw (R env; ~15 min)
	Rscript R/install_deps.R
	Rscript R/export_motrpac.R $(RAW)

check: ## check the environment and the data export
	$(PY) scripts/00_check_env.py

# ---- analysis phases (need data/raw) ---------------------------------------------------------------------
inventory: ## phase 02: data inventory
	$(PY) scripts/02_inventory.py

eda: ## phase 03: variance partition, PCA, batch partition
	$(PY) scripts/03_eda.py $(Q)
	$(PY) scripts/03_eda.py --join inner --out results/03_eda_inner $(Q)
	$(PY) scripts/03_eda.py --assays METAB --join inner --tissues $(CORE9) --out results/03_eda_metab_core9 $(Q)
	$(PY) scripts/03_eda_readout.py

eda-readout:
	$(PY) scripts/03_eda_readout.py

baselines: ## phase 04: tuned simple baselines per omic
	$(PY) scripts/04_fingerprint_baselines.py --assay TRNSCRPT --models $(MODELS) $(Q)
	$(PY) scripts/04_fingerprint_baselines.py --assay METAB $(METAB_FLAGS) --models $(MODELS) $(Q)
	$(PY) scripts/04_fingerprint_baselines.py --assay METAB --tissues $(CORE9) $(METAB_FLAGS) --models $(MODELS) --out results/04_baselines/METAB_core9 $(Q)
	$(PY) scripts/04_fingerprint_baselines.py --assay PROT --models $(MODELS) $(Q)
	$(MAKE) prot-diagnostic PY=$(PY) QUICK=$(QUICK)

prot-diagnostic:
	$(PY) scripts/04_prot_diagnostic.py --assay PROT $(Q)
	$(PY) scripts/04_prot_diagnostic.py --assay METAB $(Q)
	$(PY) scripts/04_prot_diagnostic.py --assay METAB --tissues $(CORE9) --out results/04_baselines/METAB_core9 $(Q)

panels: ## phase 05: compact panels, selector comparison, stability, annotation
	$(PY) scripts/05_compact_panels.py --assay TRNSCRPT --compare-selectors --stability-k 20 $(Q)
	$(PY) scripts/05_compact_panels.py --assay METAB --tissues $(CORE9) $(METAB_FLAGS) --stability-k 10 --confusion-k 10,20 --out results/05_panels/METAB_core9 $(Q)
	$(MAKE) annotate-panel PY=$(PY)
	$(PY) scripts/05_readout.py

annotate-panel:
	$(PY) scripts/05_annotate_panel.py --assay TRNSCRPT
	$(PY) scripts/05_annotate_panel.py --assay TRNSCRPT --stability-file results/05_panels/TRNSCRPT/stability_k20.csv
	$(PY) scripts/05_annotate_panel.py --assay METAB --panel-file results/05_panels/METAB_core9/candidate_panel.csv --out results/05_panels/METAB_core9

panel-training: ## phase 05b: training response of the panel genes in their marker tissues (consortium DA tables)
	MOTRPAC_NO_REPORT=1 $(PY) scripts/05_panel_training_response.py

conformal: ## phase 06: split-conformal sets, certificate, sizing
	$(PY) scripts/06_conformal_certify.py --assay TRNSCRPT --alpha 0.1 --delta 0.1 --one-per-animal --cal-animals 22 --n-repeats 20 --conditional --grid 10,15,20,30,50,100 $(Q)
	$(PY) scripts/06_conformal_certify.py --assay TRNSCRPT --alpha 0.05 --delta 0.05 --n-repeats 20 --out results/06_conformal/TRNSCRPT_pooled_a05 $(Q)

fusion: ## phase 07: control vs trained within tissue, single omic vs fusion, batch covariates
	$(PY) scripts/07_fusion_vs_baselines.py --models $(MODELS) --durations 1w,2w,4w,8w $(Q)
	$(MAKE) batch-check PY=$(PY) QUICK=$(QUICK)

batch-check:
	$(PY) scripts/07_batch_check.py $(Q)

shift: ## phase 08: shift tests (training state, held-out sex, time point)
	$(PY) scripts/08_shift_tests.py --assay TRNSCRPT --k 20 --recalibrate-target 3,5 $(Q)
	$(PY) scripts/08_shift_tests.py --assay METAB --k 10 --tissues $(CORE9) $(METAB_FLAGS) --recalibrate-target 3,5 --out results/08_shift/METAB_core9 $(Q)
	$(PY) scripts/08_readout.py

time-course: ## phase 15: study-design dates, physiology, the fingerprint by training duration
	MOTRPAC_NO_REPORT=1 $(PY) scripts/15_time_course_design.py
	MOTRPAC_NO_REPORT=1 $(PY) scripts/15_fingerprint_by_duration.py $(Q)

bodymap: ## phase 12: rat BodyMap transfer (needs data/external)
	@test -f data/external/bodymap_counts.csv || Rscript R/export_bodymap.R data/external 100
	$(PY) scripts/12_bodymap_validate.py $(Q)

gtex: ## phases 11, 13, 14: GTEx preparation and transfer (needs the three GTEx files)
	$(PY) scripts/11_gtex_prepare.py --tpm $(GTEX_TPM) --reads $(GTEX_READS) --attrs $(GTEX_ATTRS)
	$(PY) scripts/13_gtex_transfer.py $(Q)
	$(MAKE) transfer PY=$(PY) QUICK=$(QUICK)

transfer:
	$(PY) scripts/14_transfer_representations.py $(Q)
	$(PY) scripts/14_transfer_representations.py --targets gtex --gtex-matrix cpm --skip-cv --out results/14_transfer_cpm $(Q)

external: bodymap gtex ## the two external transfers

identifiability: ## phase 16: batch nesting, estimable pairs, QC-only baseline, bridge (portal files via MOTRPAC_PORTAL)
	MOTRPAC_NO_REPORT=1 $(PY) scripts/16_identifiability.py --bridge $(if $(MOTRPAC_PORTAL),--portal $(MOTRPAC_PORTAL),)

portal-check: ## list the portal files phase 16 needs under MOTRPAC_PORTAL
	$(PY) scripts/16_identifiability.py --check-portal $(if $(MOTRPAC_PORTAL),--portal $(MOTRPAC_PORTAL),)

all: check inventory eda baselines panels conformal fusion shift ## phases 02–08 (about 45 min on 20 cores)

# ---- no-data paths --------------------------------------------------------------------------------------
smoke: ## synthetic data → phases 02, 04, 05, 06 in --quick mode → results_smoke/ (no real data needed)
	$(PY) scripts/make_synthetic_data.py --out $(SYN)
	$(SMOKE_ENV) $(PY) scripts/02_inventory.py
	$(SMOKE_ENV) $(PY) scripts/04_fingerprint_baselines.py --assay TRNSCRPT --models centroid,logreg_l2 --quick
	$(SMOKE_ENV) $(PY) scripts/05_compact_panels.py --assay TRNSCRPT --quick
	$(SMOKE_ENV) $(PY) scripts/06_conformal_certify.py --assay TRNSCRPT --quick

test: ## Python tests (results-dependent ones read results/ or results_frozen/), the JS conformal test, the snapshot check
	$(PY) -m pytest -q
	@if command -v node >/dev/null 2>&1 && test -f site/data/conformal_fixtures.json; then \
	  node tests/test_site_conformal.js && node tests/test_score_tool.js && node tests/test_check_core.js; \
	else echo "(node or site/data absent: JS tests skipped)"; fi
	$(MAKE) verify-frozen PY=$(PY)

# ---- the site ------------------------------------------------------------------------------------------
# the --save-scores reruns of phases 06, 08, 12 and 13 into results/31_site_regen (results/ itself is not touched;
# tests/test_regen_scores.py checks that they reproduce the published tables cell for cell)
regen-scores: ## per-sample score regenerations of phases 06, 08 (k = 20, 50), 12, 13 (needs data/)
	mkdir -p $(REGEN)/logs
	MOTRPAC_NO_REPORT=1 $(PY) scripts/06_conformal_certify.py --assay TRNSCRPT --alpha 0.1 --delta 0.1 --one-per-animal --cal-animals 22 --n-repeats 20 --conditional --grid 10,15,20,30,50,100 --out $(REGEN)/06_conformal/TRNSCRPT --save-scores
	MOTRPAC_NO_REPORT=1 $(PY) scripts/08_shift_tests.py --assay TRNSCRPT --k 20 --recalibrate-target 3,5 --out $(REGEN)/08_shift_k20 --save-scores
	MOTRPAC_NO_REPORT=1 $(PY) scripts/08_shift_tests.py --assay TRNSCRPT --k 50 --recalibrate-target 3,5 --out $(REGEN)/08_shift_k50 --save-scores
	MOTRPAC_NO_REPORT=1 $(PY) scripts/12_bodymap_validate.py --out $(REGEN)/12_bodymap --save-scores
	MOTRPAC_NO_REPORT=1 $(PY) scripts/13_gtex_transfer.py --out $(REGEN)/13_gtex --save-scores
	$(PY) -m pytest -q tests/test_regen_scores.py

panel-model: ## phase 34: export the 20-gene transfer model for the Check samples page (needs data/ and regen-scores)
	MOTRPAC_NO_REPORT=1 $(PY) scripts/34_panel_model.py

freeze-results: ## copy every result file the site reads into results_frozen/ (committed snapshot)
	$(PY) scripts/33_freeze_results.py

verify-frozen: ## check results_frozen/ against its manifest
	$(PY) scripts/33_freeze_results.py --verify

product-validation: ## phase 40: the Check samples page's validation numbers → results_product/40_product (needs data/external)
	$(PY) scripts/40_product_validation.py

product-data: ## site/data/product.json and the pv_* provenance entries from results_product/
	$(PY) scripts/41_export_product_data.py

site-data: ## export site/data with provenance from results/ (or RESULTS=<dir>, or the snapshot); check the anchors
	$(PY) scripts/30_export_site_data.py --check-anchors --reconciliation $(if $(RESULTS),--results $(RESULTS),)
	$(PY) scripts/multiomic/export_site_data.py   # re-appends the mo_* provenance entries the main export rewrites
	$(PY) scripts/41_export_product_data.py       # and the pv_* entries of the Check samples page

site-test: ## the site tests: provenance, regeneration, JS conformal port, links
	$(PY) -m pytest -q tests/test_site_data.py tests/test_regen_scores.py tests/test_frozen_results.py
	node tests/test_site_conformal.js
	node tests/test_score_tool.js
	node tests/test_check_core.js
	$(PY) tools/linkcheck.py

site: site-data site-test ## export the site data and run the site tests
	@echo "site ready: python -m http.server -d site 8000"

screenshots: ## render every page in both themes and widths (needs Chrome + playwright: cd tools && npm install)
	cd tools && npm install --no-audit --no-fund >/dev/null && cd .. && node tools/screenshot.js --base $${SITE_BASE:-http://localhost:8000}

linkcheck: ## static link, asset and size check of site/
	$(PY) tools/linkcheck.py

summary-figure: ## figures/summary_figure.png from site/data (no data needed)
	$(PY) scripts/32_summary_figure.py

home-figure: ## figures/home.png: a render of the home page (needs Chrome + playwright)
	cd tools && npm install --no-audit --no-fund >/dev/null
	@$(PY) -m http.server -d site $(SITE_PORT) >/dev/null 2>&1 & echo $$! > .server.pid; sleep 1; \
	  node tools/screenshot.js --base http://localhost:$(SITE_PORT) --pages=index; status=$$?; kill $$(cat .server.pid); rm -f .server.pid; exit $$status
	$(PY) -c "from PIL import Image; im = Image.open('site/_screenshots/index-light-desktop.png'); im = im.crop((0, 0, im.width, min(im.height, 1800))); im.thumbnail((1000, 100000)); im.save('figures/home.png', optimize=True); print('wrote figures/home.png', im.size)"

figures: summary-figure home-figure ## both README figures

brand: ## the badge, favicon, logo and social-preview images from branding/ratpac-logo-1600x1800.png
	$(PY) tools/make_logo_assets.py

# ---- cleanup -------------------------------------------------------------------------------------------
clean-synthetic: ## remove synthetic data (only if the _SYNTHETIC marker is present)
	@if [ -f $(RAW)/_SYNTHETIC ]; then find $(RAW) -mindepth 1 ! -name .gitkeep -delete && echo "removed synthetic data"; else echo "$(RAW) is not synthetic; not touching it"; fi
	rm -rf $(SYN) results_smoke

clean-results: ## delete results/ (never results_frozen/)
	find results -mindepth 1 ! -name .gitkeep ! -name ABSTRACT.md -delete
