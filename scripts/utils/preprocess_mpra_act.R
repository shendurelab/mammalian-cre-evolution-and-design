#!/usr/bin/env Rscript
# Preprocess the raw per-BC MPRA count table into a slim, pre-winsorized
# per-CRE activity table for the figure scripts.
#
# The raw data/oCRE_and_CRE_optimization_MPRA_count_table.txt.gz (~62 MB) is
# too large to ship on GitHub, but virtually every figure script only needs
# the derived df_mpra_act (per CRE, per biological replicate, winsorized
# MPRA activity + per-CRE summary stats). This script computes that table
# once and writes it to data/, so the main figure scripts can load it
# directly without touching the raw BC-level counts.
#
# Run once from the repo root (after restoring the raw count table):
#   Rscript scripts/utils/preprocess_mpra_act.R
#
# Output (written to data/):
#   oCRE_and_CRE_optimization_MPRA_act.txt.gz

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(DescTools)
})

source("scripts/utils/mpra_tablemaker.R")

IN  <- "data/oCRE_and_CRE_optimization_MPRA_count_table.txt.gz"
OUT <- "data/oCRE_and_CRE_optimization_MPRA_act.txt.gz"

if (!file.exists(IN)) {
  stop("Missing input: ", IN, "\n",
       "Restore the raw count table from external storage before running this preprocessor.")
}

cat("Loading", IN, "...\n")
df_counts_w_CREs2 <- read.table(IN, header = TRUE)

cat("calculate_MPRA_coverage ...\n")
coverage_lib <- calculate_MPRA_coverage(df_counts_w_CREs2)

cat("generate_MPRA_wide_table ...\n")
df_wide <- generate_MPRA_wide_table(df_counts_w_CREs2, coverage_lib)

cat("generate_winsorize_MPRA_table ...\n")
df_mpra_act <- generate_winsorize_MPRA_table(df_wide)

cat("Writing", OUT, "(", nrow(df_mpra_act), "rows)\n")
write.table(df_mpra_act, gzfile(OUT), sep = "\t", row.names = FALSE, quote = FALSE)
cat("Done.\n")
