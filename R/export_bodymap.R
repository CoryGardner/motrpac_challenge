#!/usr/bin/env Rscript
# Export the rat BodyMap (bodymapRat, Bioconductor ExperimentHub) to CSV for scripts/12_bodymap_validate.py.
# Installs bodymapRat if needed, checks the hub file size (must be < max_mb, default 100) before
# downloading, sums technical runs per biological sample, writes:
#   data/external/bodymap_counts.csv   genes x samples (feature_ID = Ensembl gene id)
#   data/external/bodymap_meta.csv     sample, organ, sex, stage_weeks, replicate, animal_id, n_runs
#   data/external/bodymap_coldata.csv  the raw colData (one row per run)
args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1) args[1] else "data/external"
max_mb <- if (length(args) >= 2) as.numeric(args[2]) else 100
options(repos = c(CRAN = "https://cloud.r-project.org"), timeout = 1e5)
if (!requireNamespace("bodymapRat", quietly = TRUE)) {
  message("installing bodymapRat from Bioconductor")
  BiocManager::install("bodymapRat", update = FALSE, ask = FALSE)
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
