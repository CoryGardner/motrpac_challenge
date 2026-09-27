#!/usr/bin/env Rscript
# One-time: install the R dependencies of the export scripts, the MoTrPAC data package (pinned to one commit)
# and the rat BodyMap data package (Bioconductor ExperimentHub).
# Usage: Rscript R/install_deps.R
# Reference environment: R 4.5.3, Bioconductor 3.22 (environment-r.yml). The MoTrPAC package install is slow
# (~1.3 GB download, then ~5 min "moving datasets to lazyload DB"); that is normal.

options(repos = c(CRAN = "https://cloud.r-project.org"))
options(timeout = 1e5)  # the package tarball is large; extend the download timeout

# ---- CRAN -----------------------------------------------------------------------------------------------------
for (p in c("data.table", "jsonlite", "remotes", "BiocManager")) {
  if (!requireNamespace(p, quietly = TRUE)) {
    message("Installing ", p)
    install.packages(p)
  }
}

# ---- MotrpacRatTraining6moData, pinned ------------------------------------------------------------------------
# Pinned to commit f831a4fe421ec11687640452484a8247137aa74a (GitHub tag v2.1.0, 2025-08-13). The package's
# DESCRIPTION reads Version 2.0.0 at both tags v2.0.0 and v2.1.0, so packageVersion() cannot tell the two apart;
# the commit can, and it is what the reference export of 2026-09-17 used.
MOTRPAC_REPO    <- "MoTrPAC/MotrpacRatTraining6moData"
MOTRPAC_SHA     <- "f831a4fe421ec11687640452484a8247137aa74a"
MOTRPAC_VERSION <- "2.0.0"
MOTRPAC_TARBALL <- sprintf("https://api.github.com/repos/%s/tarball/%s", MOTRPAC_REPO, MOTRPAC_SHA)

motrpac_installed <- function() requireNamespace("MotrpacRatTraining6moData", quietly = TRUE)
motrpac_sha <- function() {
  # RemoteSha is written by remotes::install_github; a tarball install has none (NA).
  d <- tryCatch(packageDescription("MotrpacRatTraining6moData"), error = function(e) NULL, warning = function(w) NULL)
  if (is.null(d) || is.null(d$RemoteSha)) NA_character_ else d$RemoteSha
}

need_install <- !motrpac_installed() || (!is.na(motrpac_sha()) && !identical(motrpac_sha(), MOTRPAC_SHA))
if (need_install) {
  if (motrpac_installed()) message("MotrpacRatTraining6moData is installed from commit ", motrpac_sha(), "; reinstalling the pinned one")
  message("Installing MotrpacRatTraining6moData from GitHub commit ", substr(MOTRPAC_SHA, 1, 7), " (several minutes)")
  ok <- tryCatch({
    remotes::install_github(MOTRPAC_REPO, ref = MOTRPAC_SHA, upgrade = "never")
    TRUE
  }, error = function(e) {
    message("install_github failed: ", conditionMessage(e))
    FALSE
  })
  if (!ok) {
    # Fallback when install_github fails (seen intermittently on some networks): the tarball of the same commit.
    # Equivalent by hand:
    #   1. download MOTRPAC_TARBALL (redirects to codeload.github.com/.../legacy.tar.gz/<sha>)
    #   2. install.packages("<file>.tar.gz", repos = NULL, type = "source")
    message("Falling back to the tarball of commit ", substr(MOTRPAC_SHA, 1, 7), ": ", MOTRPAC_TARBALL)
    tgz <- file.path(tempdir(), paste0("MotrpacRatTraining6moData-", substr(MOTRPAC_SHA, 1, 7), ".tar.gz"))
    download.file(MOTRPAC_TARBALL, tgz, mode = "wb")
    install.packages(tgz, repos = NULL, type = "source")
  }
}

if (!motrpac_installed()) stop("MotrpacRatTraining6moData is not installed; see the fallback note in this file")
stopifnot(packageVersion("MotrpacRatTraining6moData") == MOTRPAC_VERSION)
sha <- motrpac_sha()
if (is.na(sha)) {
  message("OK: MotrpacRatTraining6moData ", MOTRPAC_VERSION, " (tarball install: the DESCRIPTION records no commit; pinned ",
          substr(MOTRPAC_SHA, 1, 7), ")")
} else if (!identical(sha, MOTRPAC_SHA)) {
  stop("MotrpacRatTraining6moData is installed from commit ", sha, ", not the pinned ", MOTRPAC_SHA)
} else {
  message("OK: MotrpacRatTraining6moData ", MOTRPAC_VERSION, " from commit ", substr(sha, 1, 7))
}

# ---- bodymapRat (rat BodyMap, GEO GSE53960; Bioconductor ExperimentHub) ---------------------------------------
# The version documented for the reference export (2026-09-17) is bodymapRat 1.28.0 (CC BY 4.0). Bioconductor serves
# one version per release (3.22 installs 1.26.0, 3.23 installs 1.28.0), so the version depends on the Bioconductor in
# use; the ExperimentHub resource behind them is the same GSE53960 STAR gene-count matrix. R/export_bodymap.R records
# the version actually used in data/external/bodymap_provenance.json.
BODYMAP_VERSION <- "1.28.0"
if (!requireNamespace("bodymapRat", quietly = TRUE)) {
  message("Installing bodymapRat from Bioconductor ", as.character(BiocManager::version()))
  BiocManager::install("bodymapRat", update = FALSE, ask = FALSE)
}
if (!requireNamespace("bodymapRat", quietly = TRUE)) stop("bodymapRat is not installed (BiocManager::install(\"bodymapRat\") failed)")
bm <- as.character(packageVersion("bodymapRat"))
if (bm == BODYMAP_VERSION) {
  message("OK: bodymapRat ", bm)
} else {
  message("NOTE: bodymapRat ", bm, " is installed, not the ", BODYMAP_VERSION, " documented for the reference export. ",
          "Bioconductor ", as.character(BiocManager::version()), " serves this version (3.22 serves 1.26.0, 3.23 serves 1.28.0); ",
          "the ExperimentHub data are the same GSE53960 counts. The version used is recorded in data/external/bodymap_provenance.json.")
}
