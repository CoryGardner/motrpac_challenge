#!/usr/bin/env Rscript
# Export the rat BodyMap (bodymapRat, Bioconductor ExperimentHub) to CSV for scripts/12_bodymap_validate.py.
# Installs bodymapRat if needed (the same pinned path as R/install_deps.R: Bioconductor's version, checked against
# the 1.28.0 of the reference export), checks the hub file size (must be < max_mb, default 100) before downloading,
# sums technical runs per biological sample, writes:
#   data/external/bodymap_counts.csv       genes x samples (feature_ID = Ensembl gene id)
#   data/external/bodymap_meta.csv         sample, organ, sex, stage_weeks, replicate, animal_id, n_runs
#   data/external/bodymap_coldata.csv      the raw colData (one row per run)
#   data/external/bodymap_provenance.json  package version, ExperimentHub record(s), R / Bioconductor versions, date
# Usage: Rscript R/export_bodymap.R [out_dir] [max_mb]
args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1) args[1] else "data/external"
max_mb <- if (length(args) >= 2) as.numeric(args[2]) else 100
options(repos = c(CRAN = "https://cloud.r-project.org"), timeout = 1e5)
for (p in c("jsonlite", "BiocManager")) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
# ---- bodymapRat, as in R/install_deps.R ----------------------------------------------------------------------
BODYMAP_VERSION <- "1.28.0"  # documented for the reference export of 2026-09-17 (Bioconductor serves 1.26.0 under 3.22, 1.28.0 under 3.23)
if (!requireNamespace("bodymapRat", quietly = TRUE)) {
  message("Installing bodymapRat from Bioconductor ", as.character(BiocManager::version()))
  BiocManager::install("bodymapRat", update = FALSE, ask = FALSE)
}
if (!requireNamespace("bodymapRat", quietly = TRUE)) stop("bodymapRat is not installed (BiocManager::install(\"bodymapRat\") failed)")
bm_version <- as.character(packageVersion("bodymapRat"))
if (bm_version != BODYMAP_VERSION) {
  message("NOTE: bodymapRat ", bm_version, " is installed, not the ", BODYMAP_VERSION, " documented for the reference export; Bioconductor ",
          as.character(BiocManager::version()), " serves this version (3.22 serves 1.26.0, 3.23 serves 1.28.0). The ExperimentHub data are ",
          "the same GSE53960 counts; the version used is recorded in bodymap_provenance.json.")
}
suppressPackageStartupMessages({ library(ExperimentHub); library(SummarizedExperiment) })
eh <- ExperimentHub()
q <- query(eh, "bodymapRat")
cat("hub records:", paste(names(q), collapse = ","), "\n")
print(mcols(q)[, c("title", "rdataclass", "rdatapath")])
id <- names(q)[1]
url <- paste0(hubUrl(eh), "/fetch/", sub("^EH", "", id))
h <- tryCatch(system2("curl", c("-sI", "-L", "-m", "60", shQuote(url)), stdout = TRUE), error = function(e) character(0))
len <- suppressWarnings(as.numeric(sub(".*: *", "", grep("(?i)^content-length", h, value = TRUE, perl = TRUE))))
len <- if (length(len)) max(len, na.rm = TRUE) else NA
cat(sprintf("hub resource %s size: %s MB\n", id, if (is.na(len)) "unknown" else round(len / 1e6, 1)))
if (is.na(len) || len / 1e6 > max_mb) stop(sprintf("resource size %s MB is unknown or above %s MB; not downloading", round(len / 1e6, 1), max_mb))
se <- eh[[id]]
cat("SummarizedExperiment:", paste(dim(se), collapse = " x "), "\n")
cd <- as.data.frame(colData(se))
cat("colData columns:", paste(names(cd), collapse = ", "), "\n")
write.csv(cd, file.path(out_dir, "bodymap_coldata.csv"), row.names = TRUE)
counts <- as.matrix(assay(se, 1))
# biological sample = organ x sex x stage x replicate index; technical runs summed
pick <- function(cands) { hit <- cands[cands %in% names(cd)]; if (length(hit)) hit[1] else NA }
organ_col <- pick(c("organ", "tissue", "colOrgan")); sex_col <- pick(c("sex", "Sex")); stage_col <- pick(c("stage", "age", "developmentalStage"))
rep_col <- pick(c("bioRep", "replicate", "rep", "sample", "title"))
cat("using columns:", organ_col, sex_col, stage_col, rep_col, "\n")
key <- paste(cd[[organ_col]], cd[[sex_col]], cd[[stage_col]], cd[[rep_col]], sep = "|")
samples <- unique(key)
M <- sapply(samples, function(k) rowSums(counts[, key == k, drop = FALSE]))
colnames(M) <- make.names(samples)
meta <- data.frame(sample = colnames(M), key = samples, organ = cd[[organ_col]][match(samples, key)],
                   sex = cd[[sex_col]][match(samples, key)], stage_weeks = cd[[stage_col]][match(samples, key)],
                   replicate = cd[[rep_col]][match(samples, key)], n_runs = as.integer(table(key)[samples]), stringsAsFactors = FALSE)
write.csv(cbind(feature_ID = rownames(M), as.data.frame(M, check.names = FALSE)), file.path(out_dir, "bodymap_counts.csv"), row.names = FALSE)
write.csv(meta, file.path(out_dir, "bodymap_meta.csv"), row.names = FALSE)
cat("wrote", nrow(M), "genes x", ncol(M), "biological samples (", ncol(counts), "runs )\n")
# ---- provenance ----------------------------------------------------------------------------------------------
prov <- list(
  package = "bodymapRat",
  package_version = bm_version,
  reference_version = BODYMAP_VERSION,
  license = "CC BY 4.0 (bodymapRat); data: rat BodyMap, GEO GSE53960 (Yu et al., Nat Commun 2014)",
  experimenthub_records = as.list(names(q)),
  experimenthub_record_used = id,
  experimenthub_title = as.character(mcols(q)$title[1]),
  experimenthub_rdatapath = as.character(mcols(q)$rdatapath[1]),
  hub_url = url,
  hub_file_bytes = if (is.na(len)) NULL else len,
  r_version = R.version.string,
  bioconductor_version = as.character(BiocManager::version()),
  experimenthub_version = as.character(packageVersion("ExperimentHub")),
  summarizedexperiment_version = as.character(packageVersion("SummarizedExperiment")),
  dims = list(genes = nrow(se), runs = ncol(se), biological_samples = ncol(M)),
  columns_used = list(organ = organ_col, sex = sex_col, stage = stage_col, replicate = rep_col),
  outputs = c("bodymap_counts.csv", "bodymap_meta.csv", "bodymap_coldata.csv"),
  exported_on = format(Sys.time(), "%Y-%m-%d %H:%M:%S UTC", tz = "UTC")
)
jsonlite::write_json(prov, file.path(out_dir, "bodymap_provenance.json"), auto_unbox = TRUE, pretty = TRUE, null = "null")
cat("wrote", file.path(out_dir, "bodymap_provenance.json"), "\n")
