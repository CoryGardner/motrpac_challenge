#!/usr/bin/env Rscript
# One-time: install CRAN deps and the MoTrPAC data package from GitHub.
# Usage: Rscript R/install_deps.R
# The MoTrPAC package install is slow (~1.3 GB download, then ~5 min "moving datasets to
# lazyload DB"); that's normal.

options(repos = c(CRAN = "https://cloud.r-project.org"))
options(timeout = 1e5)  # the package tarball is large; extend download timeout

# `remotes` is all install_github needs; `devtools` (much heavier) is only a fallback.
cran_pkgs <- c("data.table", "jsonlite")
if (!requireNamespace("remotes", quietly = TRUE) && !requireNamespace("devtools", quietly = TRUE)) {
  cran_pkgs <- c(cran_pkgs, "remotes")
}
for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    message("Installing ", p)
    install.packages(p)
  }
}

if (!requireNamespace("MotrpacRatTraining6moData", quietly = TRUE)) {
  message("Installing MotrpacRatTraining6moData from GitHub (this takes several minutes)")
  install_github <- if (requireNamespace("remotes", quietly = TRUE)) remotes::install_github else devtools::install_github
  install_github("MoTrPAC/MotrpacRatTraining6moData", upgrade = "never")
}

# Last resort if install_github fails (seen intermittently on macOS):
#   1. download https://api.github.com/repos/MoTrPAC/MotrpacRatTraining6moData/tarball/HEAD
#   2. install.packages("path/to/MoTrPAC-MotrpacRatTraining6moData-<hash>.tar.gz", repos = NULL, type = "source")

if (requireNamespace("MotrpacRatTraining6moData", quietly = TRUE)) {
  message("OK: MotrpacRatTraining6moData ", as.character(packageVersion("MotrpacRatTraining6moData")))
} else {
  stop("MotrpacRatTraining6moData is not installed; see the last-resort note in this file")
}
