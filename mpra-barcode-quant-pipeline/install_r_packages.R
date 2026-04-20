#!/usr/bin/env Rscript
# Install all R package dependencies for the MPRA barcode pipeline
packages <- c(
  "dplyr", "tidyverse", "ggplot2", "cowplot", "reshape2",
  "gtools", "patchwork", "stringr", "ggrastr", "DescTools",
  "castor", "ape", "phytools"
)

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing:", pkg, "\n")
    install.packages(pkg, repos = "https://cloud.r-project.org")
  } else {
    cat("Already installed:", pkg, "\n")
  }
}
