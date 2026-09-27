# make <phase> [QUICK=1]   — every phase writes to results/ and appends to results/REPORT.md
# The targets encode the final configuration of every phase (the one results/SUMMARY.md is built from), so
# `make clean-results all bodymap gtex report` reproduces the whole pre-hackathon run from scratch.
# Site: `make site` regenerates the per-sample score exports, site/data and runs the site tests.
PY      ?= python
RAW     ?= data/raw
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
export PYTHONPATH := src

.PHONY: help env synthetic export check inventory eda eda-readout baselines prot-diagnostic panels annotate-panel conformal fusion batch-check shift discordance bodymap gtex transfer external report test all clean-synthetic clean-results \
        identifiability regen-scores site-data site-test site screenshots linkcheck summary-figure

help:
	@echo "targets: env synthetic export check inventory eda baselines panels conformal fusion shift discordance report test all"
	@echo "         bodymap gtex transfer external (need data/external/); clean-synthetic clean-results"
	@echo "         identifiability (phase 16) · regen-scores (phases 06/12/13 with --save-scores) · site-data · site-test · site · screenshots · linkcheck · summary-figure"
	@echo "variables: QUICK=1 MODELS=$(MODELS) RAW=data/raw"

env:
	$(PY) -m pip install -r requirements.txt

synthetic:
	$(PY) scripts/make_synthetic_data.py --out $(RAW)

export:
	Rscript R/install_deps.R
	Rscript R/export_motrpac.R $(RAW)

check:
	$(PY) scripts/00_check_env.py

inventory:
	$(PY) scripts/02_inventory.py

eda:
	$(PY) scripts/03_eda.py $(Q)
	$(PY) scripts/03_eda.py --join inner --out results/03_eda_inner $(Q)
	$(PY) scripts/03_eda.py --assays METAB --join inner --tissues $(CORE9) --out results/03_eda_metab_core9 $(Q)
	$(PY) scripts/03_eda_readout.py

eda-readout:
	$(PY) scripts/03_eda_readout.py

baselines:
	$(PY) scripts/04_fingerprint_baselines.py --assay TRNSCRPT --models $(MODELS) $(Q)
	$(PY) scripts/04_fingerprint_baselines.py --assay METAB $(METAB_FLAGS) --models $(MODELS) $(Q)
	$(PY) scripts/04_fingerprint_baselines.py --assay METAB --tissues $(CORE9) $(METAB_FLAGS) --models $(MODELS) --out results/04_baselines/METAB_core9 $(Q)
	$(PY) scripts/04_fingerprint_baselines.py --assay PROT --models $(MODELS) $(Q)
	$(MAKE) prot-diagnostic PY=$(PY) QUICK=$(QUICK)

prot-diagnostic:
	$(PY) scripts/04_prot_diagnostic.py --assay PROT $(Q)
	$(PY) scripts/04_prot_diagnostic.py --assay METAB $(Q)
	$(PY) scripts/04_prot_diagnostic.py --assay METAB --tissues $(CORE9) --out results/04_baselines/METAB_core9 $(Q)

panels:
	$(PY) scripts/05_compact_panels.py --assay TRNSCRPT --compare-selectors --stability-k 20 $(Q)
	$(PY) scripts/05_compact_panels.py --assay METAB --tissues $(CORE9) $(METAB_FLAGS) --stability-k 10 --confusion-k 10,20 --out results/05_panels/METAB_core9 $(Q)
	$(MAKE) annotate-panel PY=$(PY)
	$(PY) scripts/05_readout.py

annotate-panel:
	$(PY) scripts/05_annotate_panel.py --assay TRNSCRPT
	$(PY) scripts/05_annotate_panel.py --assay TRNSCRPT --stability-file results/05_panels/TRNSCRPT/stability_k20.csv
	$(PY) scripts/05_annotate_panel.py --assay METAB --panel-file results/05_panels/METAB_core9/candidate_panel.csv --out results/05_panels/METAB_core9

conformal:
	$(PY) scripts/06_conformal_certify.py --assay TRNSCRPT --alpha 0.1 --delta 0.1 --one-per-animal --cal-animals 22 --n-repeats 20 --conditional --grid 10,15,20,30,50,100 $(Q)
	$(PY) scripts/06_conformal_certify.py --assay TRNSCRPT --alpha 0.05 --delta 0.05 --n-repeats 20 --out results/06_conformal/TRNSCRPT_pooled_a05 $(Q)

fusion:
	$(PY) scripts/07_fusion_vs_baselines.py --models $(MODELS) --durations 1w,2w,4w,8w $(Q)
	$(MAKE) batch-check PY=$(PY) QUICK=$(QUICK)

batch-check:
	$(PY) scripts/07_batch_check.py $(Q)

shift:
	$(PY) scripts/08_shift_tests.py --assay TRNSCRPT --k 20 --recalibrate-target 3,5 $(Q)
	$(PY) scripts/08_shift_tests.py --assay METAB --k 10 --tissues $(CORE9) $(METAB_FLAGS) --recalibrate-target 3,5 --out results/08_shift/METAB_core9 $(Q)
	$(PY) scripts/08_readout.py

discordance:
	$(PY) scripts/09_discordance.py $(Q)

bodymap:
	@test -f data/external/bodymap_counts.csv || Rscript R/export_bodymap.R data/external 100
	$(PY) scripts/12_bodymap_validate.py $(Q)

gtex:
	$(PY) scripts/11_gtex_prepare.py --tpm $(GTEX_TPM) --reads $(GTEX_READS) --attrs $(GTEX_ATTRS)
	$(PY) scripts/13_gtex_transfer.py $(Q)
	$(MAKE) transfer PY=$(PY) QUICK=$(QUICK)

transfer:
	$(PY) scripts/14_transfer_representations.py $(Q)
	$(PY) scripts/14_transfer_representations.py --targets gtex --gtex-matrix cpm --skip-cv --out results/14_transfer_cpm $(Q)

external: bodymap gtex

report:
	$(PY) scripts/10_make_report.py

test:
	$(PY) -m pytest -q
	@command -v node >/dev/null 2>&1 && test -f site/data/conformal_fixtures.json && node tests/test_site_conformal.js || echo "(node or site/data absent: JS conformal test skipped)"

all: check inventory eda baselines panels conformal fusion shift discordance report

# ---- identifiability (phase 16) and the site ---------------------------------------------------------
identifiability:
	MOTRPAC_NO_REPORT=1 $(PY) scripts/16_identifiability.py

# the --save-scores reruns of phases 06, 12 and 13 into results/31_site_regen (results/ itself is not touched;
# tests/test_regen_scores.py checks that they reproduce the published tables cell for cell)
regen-scores:
	mkdir -p $(REGEN)/logs
	MOTRPAC_NO_REPORT=1 $(PY) scripts/06_conformal_certify.py --assay TRNSCRPT --alpha 0.1 --delta 0.1 --one-per-animal --cal-animals 22 --n-repeats 20 --conditional --grid 10,15,20,30,50,100 --out $(REGEN)/06_conformal/TRNSCRPT --save-scores
	MOTRPAC_NO_REPORT=1 $(PY) scripts/12_bodymap_validate.py --out $(REGEN)/12_bodymap --save-scores
	MOTRPAC_NO_REPORT=1 $(PY) scripts/13_gtex_transfer.py --out $(REGEN)/13_gtex --save-scores
	$(PY) -m pytest -q tests/test_regen_scores.py

site-data:
	$(PY) scripts/30_export_site_data.py --check-anchors --reconciliation

site-test:
	$(PY) -m pytest -q tests/test_site_data.py tests/test_regen_scores.py
	node tests/test_site_conformal.js
	$(PY) tools/linkcheck.py

# export data → run the site tests; screenshots are optional (need playwright: cd tools && npm install)
site: site-data site-test
	@echo "site ready: python -m http.server -d site 8000"

screenshots:
	cd tools && npm install --no-audit --no-fund >/dev/null && cd .. && node tools/screenshot.js --base $${SITE_BASE:-http://localhost:8000}

linkcheck:
	$(PY) tools/linkcheck.py

summary-figure:
	$(PY) scripts/32_summary_figure.py

clean-synthetic:
	@if [ -f $(RAW)/_SYNTHETIC ]; then find $(RAW) -mindepth 1 ! -name .gitkeep -delete && echo "removed synthetic data"; else echo "$(RAW) is not synthetic; not touching it"; fi

clean-results:
	find results -mindepth 1 ! -name .gitkeep ! -name ABSTRACT.md -delete
