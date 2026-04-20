#!/usr/bin/env Rscript
# Preprocess the CRE optimization ProBound TFBS prediction table into a slim
# version for committing to the repo.
#
# Mirrors scripts/utils/preprocess_oCRE_TFBS.R (see that file for the full
# rationale). Keeps only the 6 columns and 7 TFs that load_figure5_data.R
# actually uses and rounds norm_affinity to 5 significant figures.
#
# Run once from the repo root (after restoring the raw file):
#   Rscript scripts/utils/preprocess_CRE_optimization_TFBS.R
#
# Output:
#   data/CRE_optimization_ProBound_ParEndo_TFBS_predictions_slim.txt.gz

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

IN  <- "data/CRE_optimization_ProBound_ParEndo_TFBS_predictions.txt.gz"
OUT <- "data/CRE_optimization_ProBound_ParEndo_TFBS_predictions_slim.txt.gz"
KEEP_TFS <- c("Jun_Atf3", "Foxa2", "Gata4", "Gata6", "Sox17", "Klf4", "Hnf1b")

if (!file.exists(IN)) {
  stop("Missing input: ", IN, "\n",
       "Restore the raw TFBS prediction table from external storage first.")
}

cat("Loading", IN, "...\n")
df <- read_tsv(IN, show_col_types = FALSE)

cat("  rows:", nrow(df), "TFs:", length(unique(df$TF_name)), "\n")

slim <- df %>%
  filter(TF_name %in% KEEP_TFS) %>%
  select(start_pos, end_pos, TF_name, CRE, TFBS_orientation, norm_affinity) %>%
  mutate(norm_affinity = signif(norm_affinity, 5))

cat("  slim rows:", nrow(slim), "\n")
cat("Writing", OUT, "...\n")
write_tsv(slim, gzfile(OUT))
cat("Done:", file.info(OUT)$size / 1024^2, "MB\n")
