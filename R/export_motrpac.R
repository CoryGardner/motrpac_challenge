#!/usr/bin/env Rscript
# Export the MotrpacRatTraining6moData objects to flat CSVs that the Python side reads.
#
# Usage:
#   Rscript R/export_motrpac.R [out_dir] [assays]
#     out_dir  default "data/raw"
#     assays   optional comma-separated subset, e.g. "TRNSCRPT,PROT,METAB" (default: all)
#
# Output layout (see docs/DATA_GUIDE.md):
#   out_dir/pheno.csv
#   out_dir/norm/<ASSAY>__<TISSUE>.csv        sample-level normalized data (features x samples)
#   out_dir/counts/<ASSAY>__<TISSUE>.csv      raw counts (RNA-seq)
#   out_dir/da/<ASSAY>__<TISSUE>.csv          differential analysis tables
#   out_dir/meta/<ASSAY>.csv                  assay-level sample metadata
#   out_dir/meta/<ASSAY>_VIALS.csv            METAB/IMMUNO: vial -> platform -> bid -> exported column
#   out_dir/meta/<ASSAY>_FEATURES.csv         METAB/IMMUNO: exported feature rows -> platform
#   out_dir/codes/{tissue_abbrev,assay_abbrev,group_colors}.csv
#   out_dir/{feature_to_gene,rat_to_human_gene,metab_feature_id_map,
#            training_regulated_features,outliers}.csv
#   out_dir/manifest.json
#
# Objects are discovered by regex so minor naming differences between package versions
# (e.g. HEART_PROT_DA vs PROT_HEART_DA) are tolerated; assay and tissue tokens are validated
# against ASSAY_ABBREV / TISSUE_ABBREV so nothing that merely looks like a pattern is exported.
# Metabolomics and immunoassays are stored as nested lists ([[platform]][[tissue]]) in the
# package, not as per-tissue objects; they are flattened to one table per tissue with one
# column per biospecimen (bid), named by a representative viallabel so that
# bid == viallabel[:5] keeps holding downstream. Anything unrecognized is listed in the
# manifest under "skipped" with a reason, so nothing disappears silently.

args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1) args[1] else "data/raw"
assay_filter <- if (length(args) >= 2) strsplit(args[2], ",")[[1]] else NULL

suppressPackageStartupMessages({
  library(MotrpacRatTraining6moData)
  library(data.table)
  library(jsonlite)
})

pkg <- "MotrpacRatTraining6moData"
pkg_version <- as.character(packageVersion(pkg))
t_start <- Sys.time()
message("Exporting ", pkg, " ", pkg_version, " -> ", out_dir)

for (d in c("", "norm", "counts", "da", "meta", "codes")) {
  dir.create(file.path(out_dir, d), recursive = TRUE, showWarnings = FALSE)
}

# ---- helpers -------------------------------------------------------------------------------
env <- new.env()
load_obj <- function(name) {
  # Lazy-loaded datasets: data(list=...) materializes them into env
  suppressWarnings(data(list = name, package = pkg, envir = env))
  if (!exists(name, envir = env, inherits = FALSE)) return(NULL)
  get(name, envir = env, inherits = FALSE)
}
drop_obj <- function(name) if (exists(name, envir = env, inherits = FALSE)) rm(list = name, envir = env)

to_df <- function(obj) {
  if (is.data.frame(obj)) return(obj)
  if (is.matrix(obj)) return(as.data.frame(obj, stringsAsFactors = FALSE))
  if (is.atomic(obj) && !is.null(names(obj))) {
    return(data.frame(name = names(obj), value = unname(as.character(obj)), stringsAsFactors = FALSE))
  }
  if (is.atomic(obj)) return(data.frame(value = as.character(obj), stringsAsFactors = FALSE))
  NULL
}

# Identifier columns are written as plain digit strings: fwrite would otherwise print doubles
# such as 90217015305 as 9.0217e+10 (default scipen), which breaks every viallabel join.
ID_COLS <- c("viallabel", "pid", "bid", "labelid", "PID", "BID")
format_id <- function(x) {
  out <- if (is.double(x)) sprintf("%.0f", x) else as.character(x)
  out[is.na(x)] <- NA_character_
  out
}
write_csv <- function(df, path) {
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  for (j in seq_along(df)) {
    if (is.factor(df[[j]])) {
      df[[j]] <- as.character(df[[j]])          # factors -> character so nothing turns into integer codes
    } else if (is.list(df[[j]])) {
      df[[j]] <- vapply(df[[j]], function(v) paste(as.character(v), collapse = "|"), "")
      message("    (list column '", names(df)[j], "' collapsed with '|')")
    }
  }
  for (cn in intersect(ID_COLS, names(df))) if (!is.character(df[[cn]])) df[[cn]] <- format_id(df[[cn]])
  fwrite(df, path, na = "NA")
  invisible(dim(df))
}

# Abbreviation tables (validate object-name tokens; hyphenated labels stay index-aligned with tokens)
tissue_abbrev <- load_obj("TISSUE_ABBREV")
assay_abbrev  <- load_obj("ASSAY_ABBREV")
FALLBACK_TISSUES <- c("ADRNL","BAT","BLOOD","COLON","CORTEX","HEART","HIPPOC","HYPOTH","KIDNEY","LIVER",
                      "LUNG","OVARY","PLASMA","SKM-GN","SKM-VL","SMLINT","SPLEEN","TESTES","VENACV","WAT-SC")
FALLBACK_ASSAYS <- c("TRNSCRPT","PROT","PHOSPHO","ACETYL","UBIQ","METAB","IMMUNO","ATAC","METHYL")
tissue_abbrev_vec <- unique(c(as.character(unlist(tissue_abbrev)), FALLBACK_TISSUES))
tissue_tokens <- gsub("-", "", toupper(tissue_abbrev_vec))
assay_tokens  <- unique(c(as.character(unlist(assay_abbrev)), FALLBACK_ASSAYS))
tissue_label_for <- function(tok) tissue_abbrev_vec[match(tok, tissue_tokens)]

keep_assay <- function(a) is.null(assay_filter) || a %in% assay_filter

# PHENO identifier sets, used to report (not fix) samples the Python loader will not match
pheno <- load_obj("PHENO")
pheno_vl  <- as.character(pheno$viallabel)
pheno_bid <- format_id(pheno$bid)
n_bad_bid <- sum(substr(pheno_vl, 1, 5) != pheno_bid, na.rm = TRUE)
message("PHENO: ", nrow(pheno), " vials, ", length(unique(pheno$pid)), " animals; rows where bid != substr(viallabel,1,5): ", n_bad_bid)
if (n_bad_bid > 0) warning("PHENO bid/viallabel mismatch: the Python bid fallback join will be wrong for these rows")
drop_obj("PHENO")

# ---- enumerate objects ----------------------------------------------------------------------
items <- data(package = pkg)$results[, "Item"]
items <- unique(sub(" .*$", "", items))  # "X (Y)" -> "X"
message(length(items), " data objects found")

manifest <- list(package = pkg, version = pkg_version, exported_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
                 files = list(), skipped = character(0), skip_reasons = list())
add_file <- function(kind, name, path, dims, assay = NA, tissue = NA) {
  manifest$files[[length(manifest$files) + 1]] <<- list(
    kind = kind, object = name, path = path, n_rows = dims[1], n_cols = dims[2],
    assay = assay, tissue = tissue)
}
skip <- function(name, reason) {
  manifest$skipped <<- unique(c(manifest$skipped, name))
  manifest$skip_reasons[[name]] <<- reason
}
unskip <- function(name) {
  manifest$skipped <<- setdiff(manifest$skipped, name)
  manifest$skip_reasons[[name]] <<- NULL
}

feature_ids <- character(0)  # collected from norm/count files to subset FEATURE_TO_GENE

KIND_OF <- c(NORM_DATA = "norm", RAW_COUNTS = "counts", DA = "da")
for (name in sort(items)) {
  m1 <- regmatches(name, regexec("^([A-Z]+)_([A-Z]+)_(NORM_DATA|RAW_COUNTS|DA)$", name))[[1]]
  m2 <- regmatches(name, regexec("^([A-Z]+)_META$", name))[[1]]
  if (length(m1) == 4) {
    kind <- KIND_OF[[m1[4]]]
    a <- m1[2]; b <- m1[3]
    # <ASSAY>_<TISSUE>_...; DA tables also tolerate <TISSUE>_<ASSAY>_DA
    assay <- if (a %in% assay_tokens) a else if (kind == "da" && b %in% assay_tokens) b else NA_character_
    tissue <- if (identical(assay, a)) b else a
    if (is.na(assay) || !(tissue %in% tissue_tokens)) {
      skip(name, sprintf("unrecognized assay/tissue token (%s, %s)", a, b)); next
    }
    if (!keep_assay(assay)) { skip(name, "excluded by assay filter"); next }
    obj <- load_obj(name); df <- to_df(obj); drop_obj(name)
    if (is.null(df)) { skip(name, "not tabular"); next }
    path <- file.path(out_dir, kind, sprintf("%s__%s.csv", assay, tissue))
    dims <- write_csv(df, path)
    if (kind != "da" && "feature_ID" %in% names(df)) feature_ids <- c(feature_ids, as.character(df$feature_ID))
    add_file(kind, name, path, dims, assay, tissue)
    message(sprintf("  %-6s %s  %d x %d", kind, name, dims[1], dims[2]))
  } else if (length(m2) == 2) {
    assay <- m2[2]
    if (!(assay %in% assay_tokens)) { skip(name, "unrecognized assay token"); next }
    if (!keep_assay(assay)) { skip(name, "excluded by assay filter"); next }
    obj <- load_obj(name); df <- to_df(obj); drop_obj(name)
    if (is.null(df)) { skip(name, "not tabular"); next }
    path <- file.path(out_dir, "meta", sprintf("%s.csv", assay))
    dims <- write_csv(df, path)
    add_file("meta", name, path, dims, assay)
    message(sprintf("  %-6s %s  %d x %d", "meta", name, dims[1], dims[2]))
  } else {
    skip(name, "no pattern")  # singletons are handled below (and un-skipped there)
  }
}

# ---- nested sample-level objects (METAB, IMMUNO) ------------------------------------------
# Leaf -> numeric matrix features x viallabels. METAB leaves are features x viallabels with the
# feature ID in rownames; IMMUNO leaves are samples x analytes with a leading `viallabel` column.
leaf_to_matrix <- function(leaf, key) {
  if (is.null(leaf)) return(NULL)
  if (is.data.frame(leaf) && "viallabel" %in% names(leaf)) {
    vl <- format_id(leaf$viallabel)
    num <- leaf[, vapply(leaf, is.numeric, TRUE), drop = FALSE]
    m <- t(as.matrix(num)); colnames(m) <- vl
  } else if (is.data.frame(leaf)) {
    fid <- if ("feature_ID" %in% names(leaf)) as.character(leaf$feature_ID) else rownames(leaf)
    num <- leaf[, vapply(leaf, is.numeric, TRUE), drop = FALSE]
    m <- as.matrix(num); rownames(m) <- fid
  } else if (is.matrix(leaf)) {
    m <- leaf
  } else {
    message("    (", key, ": leaf is ", class(leaf)[1], ", skipped)"); return(NULL)
  }
  if (is.null(rownames(m)) || identical(rownames(m), as.character(seq_len(nrow(m))))) {
    message("    (", key, ": no feature ids, skipped)"); return(NULL)
  }
  cn <- sub("^X(?=[0-9]{11}$)", "", colnames(m), perl = TRUE)
  ok <- grepl("^[0-9]{11}$", cn)
  if (any(!ok)) message("    (", key, ": dropping ", sum(!ok), " non-viallabel columns: ", paste(head(cn[!ok], 5), collapse = ","), ")")
  m <- m[, ok, drop = FALSE]; colnames(m) <- cn[ok]
  storage.mode(m) <- "double"
  m
}

# Within one platform, average any vials that share a biospecimen (bid = first 5 digits).
collapse_to_bid <- function(m) {
  bid <- substr(colnames(m), 1, 5)
  if (!anyDuplicated(bid)) { colnames(m) <- bid; return(list(m = m, n_collapsed = 0L)) }
  groups <- split(seq_along(bid), bid)
  out <- vapply(groups, function(j) rowMeans(m[, j, drop = FALSE], na.rm = TRUE), numeric(nrow(m)))
  out <- matrix(out, nrow = nrow(m), dimnames = list(rownames(m), names(groups)))
  out[is.nan(out)] <- NA_real_
  list(m = out, n_collapsed = length(bid) - length(groups))
}

export_nested <- function(obj_name, assay, feature_lookup) {
  if (!(obj_name %in% items)) { message("  (absent) ", obj_name); return(invisible(NULL)) }
  if (!keep_assay(assay)) { skip(obj_name, "excluded by assay filter"); return(invisible(NULL)) }
  obj <- load_obj(obj_name)
  if (!is.list(obj) || is.data.frame(obj)) { skip(obj_name, "not a nested list"); return(invisible(NULL)) }
  platforms <- sort(names(obj))   # fixed order -> deterministic row order for repeated features
  tissue_keys <- sort(unique(unlist(lapply(obj, names))))
  message("  nested ", obj_name, ": ", length(platforms), " platforms, ", length(tissue_keys), " tissue keys")
  vial_map <- list(); feat_map <- list(); n_files <- 0L
  for (key in tissue_keys) {
    tok <- gsub("-", "", toupper(key))
    if (!(tok %in% tissue_tokens)) { message("    (unknown tissue key '", key, "', skipped)"); next }
    tissue_label <- tissue_label_for(tok)
    parts <- list(); n_collapsed <- 0L; vm_t <- list()
    for (p in platforms) {
      m <- leaf_to_matrix(obj[[p]][[key]], paste(p, key))
      if (is.null(m) || nrow(m) == 0 || ncol(m) == 0) next
      vm_t[[p]] <- data.table(tissue = tissue_label, dataset = p, viallabel = colnames(m),
                              bid = substr(colnames(m), 1, 5), n_features = nrow(m))
      cm <- collapse_to_bid(m); n_collapsed <- n_collapsed + cm$n_collapsed
      parts[[p]] <- data.table(feature_ID = rownames(m), dataset = p, cm$m)
    }
    if (!length(parts)) next
    vm <- rbindlist(vm_t)
    vm[, in_pheno := viallabel %in% pheno_vl]
    # one exported column per bid, named by a representative real viallabel (present in PHENO if any)
    rep <- vm[order(bid, -in_pheno, viallabel), .(representative_viallabel = viallabel[1]), by = bid]
    rep_vec <- setNames(rep$representative_viallabel, rep$bid)
    vm <- merge(vm, rep, by = "bid", all.x = TRUE)
    dt <- rbindlist(parts, use.names = TRUE, fill = TRUE)
    bid_cols <- setdiff(names(dt), c("feature_ID", "dataset"))
    stopifnot(all(bid_cols %in% names(rep_vec)), all(substr(rep_vec[bid_cols], 1, 5) == bid_cols))
    setnames(dt, bid_cols, unname(rep_vec[bid_cols]))
    feat <- feature_lookup(tissue_label, dt$dataset, dt$feature_ID)
    repeated <- duplicated(dt$feature_ID) | duplicated(dt$feature_ID, fromLast = TRUE)
    feat_map[[key]] <- data.table(tissue = tissue_label, dataset = dt$dataset, feature_ID = dt$feature_ID,
                                  feature = feat, repeated = repeated)
    dt[, dataset := NULL]
    dt[, `:=`(feature = feat, tissue = tissue_label, assay = assay)]
    setcolorder(dt, c("feature", "feature_ID", "tissue", "assay", sort(unname(rep_vec[bid_cols]))))
    path <- file.path(out_dir, "norm", sprintf("%s__%s.csv", assay, tok))
    dims <- write_csv(dt, path)
    add_file("norm", obj_name, path, dims, assay, tok)
    n_files <- n_files + 1L
    message(sprintf("  norm   %s[%s]  %d x %d  (%d platforms, %d bids, %d vials averaged within platform, %d repeated feature rows, %d regulated-flagged, %d bids not in PHENO)",
                    assay, key, dims[1], dims[2], length(parts), length(bid_cols), n_collapsed,
                    sum(repeated), sum(!is.na(feat)), sum(!(bid_cols %in% pheno_bid))))
    vial_map[[key]] <- vm
  }
  if (length(vial_map)) {
    vm_all <- rbindlist(vial_map)
    setcolorder(vm_all, c("tissue", "dataset", "viallabel", "bid", "representative_viallabel", "n_features", "in_pheno"))
    path <- file.path(out_dir, "meta", sprintf("%s_VIALS.csv", assay))
    dims <- write_csv(vm_all, path); add_file("meta", obj_name, path, dims, assay)
    fm <- rbindlist(feat_map)
    path2 <- file.path(out_dir, "meta", sprintf("%s_FEATURES.csv", assay))
    dims2 <- write_csv(fm, path2); add_file("meta", obj_name, path2, dims2, assay)
    message(sprintf("  meta   %s: %d tissue files; %d vials on %d platforms (%d vials not in PHENO, %d bids not in PHENO); %d feature rows (%d repeated)",
                    assay, n_files, nrow(vm_all), length(unique(vm_all$dataset)), sum(!vm_all$in_pheno),
                    length(setdiff(unique(vm_all$bid), pheno_bid)), nrow(fm), sum(fm$repeated)))
  }
  unskip(obj_name)
  drop_obj(obj_name)
  invisible(NULL)
}

reg <- load_obj("TRAINING_REGULATED_FEATURES")
reg_features <- if (!is.null(reg) && "feature" %in% names(reg)) unique(as.character(reg$feature)) else character(0)
drop_obj("reg")
message(length(reg_features), " training-regulated feature strings")

metab_lookup <- local({
  map <- load_obj("METAB_FEATURE_ID_MAP")
  if (is.null(map)) return(function(tissue_label, dataset, feature_ID) {
    f <- sprintf("METAB;%s;%s", tissue_label, feature_ID); ifelse(f %in% reg_features, f, NA_character_) })
  norm_t <- function(x) gsub("-", "", toupper(as.character(x)))
  f1 <- setNames(as.character(map$feature), paste(norm_t(map$tissue), map$dataset, map$feature_ID_sample_data, sep = "\t"))
  f2 <- setNames(as.character(map$feature), paste(norm_t(map$tissue), map$feature_ID_sample_data, sep = "\t"))
  function(tissue_label, dataset, feature_ID) {
    tok <- norm_t(tissue_label)
    f <- unname(f1[paste(tok, dataset, feature_ID, sep = "\t")])
    miss1 <- is.na(f)
    f[miss1] <- unname(f2[paste(tok, feature_ID[miss1], sep = "\t")])
    miss2 <- is.na(f)
    f[miss2] <- sprintf("METAB;%s;%s", tissue_label, feature_ID[miss2])
    message(sprintf("    feature strings: %d via (tissue,platform,id), %d via (tissue,id), %d string-built",
                    sum(!miss1), sum(miss1 & !miss2), sum(miss2)))
    ifelse(f %in% reg_features, f, NA_character_)
  }
})
immuno_lookup <- function(tissue_label, dataset, feature_ID) {
  f <- sprintf("IMMUNO;%s;%s", tissue_label, feature_ID)
  ifelse(f %in% reg_features, f, NA_character_)
}
export_nested("METAB_NORM_DATA_NESTED", "METAB", metab_lookup)
export_nested("IMMUNO_NORM_DATA_NESTED", "IMMUNO", immuno_lookup)

# ---- singletons -----------------------------------------------------------------------------
export_single <- function(name, path, transform = identity) {
  obj <- load_obj(name)
  if (is.null(obj)) { message("  (absent) ", name); return(invisible(NULL)) }
  df <- to_df(transform(obj))
  if (is.null(df)) { skip(name, "not tabular"); message("  (skip, not tabular) ", name); return(invisible(NULL)) }
  dims <- write_csv(df, path)
  add_file("table", name, path, dims)
  unskip(name)
  message("  table  ", name, "  ", dims[1], " x ", dims[2])
  drop_obj(name)
}

# PHENO: make sure viallabel is a column (rownames are viallabels in some package versions)
export_single("PHENO", file.path(out_dir, "pheno.csv"), function(x) {
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  if (!"viallabel" %in% names(x)) x <- cbind(viallabel = rownames(x), x)
  x$viallabel <- as.character(x$viallabel)
  x
})

# FEATURE_TO_GENE is ~4M rows; keep only features we exported
feature_ids <- unique(feature_ids)
message(length(feature_ids), " unique feature IDs collected from norm/counts tables")
export_single("FEATURE_TO_GENE", file.path(out_dir, "feature_to_gene.csv"), function(x) {
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  if (length(feature_ids) > 0 && "feature_ID" %in% names(x)) {
    n0 <- nrow(x); x <- x[x$feature_ID %in% feature_ids, ]
    message("  FEATURE_TO_GENE: ", n0, " -> ", nrow(x), " rows after subsetting to exported features")
  }
  x
})

export_single("RAT_TO_HUMAN_GENE", file.path(out_dir, "rat_to_human_gene.csv"))
export_single("METAB_FEATURE_ID_MAP", file.path(out_dir, "metab_feature_id_map.csv"))
export_single("TRAINING_REGULATED_FEATURES", file.path(out_dir, "training_regulated_features.csv"))
export_single("OUTLIERS", file.path(out_dir, "outliers.csv"))
export_single("REPEATED_FEATURES", file.path(out_dir, "repeated_features.csv"))
export_single("TISSUE_ABBREV", file.path(out_dir, "codes", "tissue_abbrev.csv"))
export_single("ASSAY_ABBREV", file.path(out_dir, "codes", "assay_abbrev.csv"))
export_single("TISSUE_CODE_TO_ABBREV", file.path(out_dir, "codes", "tissue_code_to_abbrev.csv"))
export_single("ASSAY_CODE_TO_ABBREV", file.path(out_dir, "codes", "assay_code_to_abbrev.csv"))
export_single("GROUP_COLORS", file.path(out_dir, "codes", "group_colors.csv"))
export_single("TISSUE_COLORS", file.path(out_dir, "codes", "tissue_colors.csv"))
export_single("ASSAY_COLORS", file.path(out_dir, "codes", "assay_colors.csv"))
for (nm in c("TRNSCRPT_FEATURE_ANNOT", "PROT_FEATURE_ANNOT", "PHOSPHO_FEATURE_ANNOT",
             "ACETYL_FEATURE_ANNOT", "UBIQ_FEATURE_ANNOT")) {
  if (nm %in% items) export_single(nm, file.path(out_dir, "codes", paste0(tolower(nm), ".csv")))
}

manifest$skipped <- I(manifest$skipped)  # always a JSON array, even with one element
write_json(manifest, file.path(out_dir, "manifest.json"), auto_unbox = TRUE, pretty = TRUE, null = "null")
kinds <- table(vapply(manifest$files, function(f) f$kind, ""))
message("Done in ", round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 1), " min. ",
        length(manifest$files), " files written (", paste(names(kinds), kinds, sep = "=", collapse = ", "), "); ",
        length(manifest$skipped), " objects skipped (listed with reasons in manifest.json).")
