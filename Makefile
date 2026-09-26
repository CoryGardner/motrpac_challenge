# make <phase> [QUICK=1]   — every phase writes to results/ and appends to results/REPORT.md
# The targets encode the final configuration of every phase (the one results/SUMMARY.md is built from), so
# `make clean-results all bodymap gtex report` reproduces the whole pre-hackathon run from scratch.
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

.PHONY: help env synthetic export check inventory eda eda-readout baselines prot-diagnostic panels annotate-panel conformal fusion batch-check shift discordance bodymap gtex transfer external report test all clean-synthetic clean-results

help:
	@echo "targets: env synthetic export check inventory eda baselines panels conformal fusion shift discordance report test all"
	@echo "         bodymap gtex transfer external (need data/external/); clean-synthetic clean-results"
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

all: check inventory eda baselines panels conformal fusion shift discordance report

clean-synthetic:
	@if [ -f $(RAW)/_SYNTHETIC ]; then find $(RAW) -mindepth 1 ! -name .gitkeep -delete && echo "removed synthetic data"; else echo "$(RAW) is not synthetic; not touching it"; fi

clean-results:
	find results -mindepth 1 ! -name .gitkeep ! -name ABSTRACT.md -delete
