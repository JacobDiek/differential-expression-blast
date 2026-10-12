# Installs the R packages the pipeline needs (only the ones that are missing).
#   Rscript R/install_packages.R              # pipeline packages
#   Rscript R/install_packages.R --analysis   # also the brook stickleback analysis script's packages

args <- commandArgs(trailingOnly = TRUE)
repos <- "https://cloud.r-project.org"

cran_packages <- c("ggplot2", "jsonlite")
bioc_packages <- c("DESeq2")

if ("--analysis" %in% args) {
  cran_packages <- c(cran_packages, "ashr", "ggrepel")
  bioc_packages <- c(bioc_packages, "apeglm")
}

missing <- function(packages) {
  packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
}

to_install <- missing(cran_packages)
if (length(to_install) > 0) {
  cat("Installing from CRAN:", to_install, "\n")
  install.packages(to_install, repos = repos)
}

to_install <- missing(bioc_packages)
if (length(to_install) > 0) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = repos)
  }
  cat("Installing from Bioconductor:", to_install, "\n")
  BiocManager::install(to_install, update = FALSE, ask = FALSE)
}

still_missing <- missing(c(cran_packages, bioc_packages))
if (length(still_missing) > 0) {
  stop("Could not install: ", paste(still_missing, collapse = ", "))
}
cat("R packages are installed.\n")
