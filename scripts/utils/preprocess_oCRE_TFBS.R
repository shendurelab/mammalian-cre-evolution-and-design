#!/usr/bin/env Rscript
# Preprocess the full oCRE ProBound TFBS prediction table into a slim version
# for committing to the repo.
#
# The raw data/oCRE_ProBound_ParEndo_TFBS_predictions.txt.gz (~264 MB) has
# 12 columns across ~6.1 M rows for 12 TFs. The figure scripts only use
# 6 of those columns, and filter to 5 TF families (merging Gata4 + Gata6
# into Gata4/6). This script:
#   - drops unused columns (kmer_seq, affinity_FOR, affinity_REV,
#     substring_l, norm_affinity_FOR, norm_affinity_REV)
#   - filters to TFs in the final set (Jun_Atf3, Foxa2, Gata4, Gata6,
#     Sox17, Klf4, Hnf1b)
#   - rounds norm_affinity to 5 significant figures (more than enough
#     for downstream filters, Gata4/6 slice_max, and binding_div)
#
# Result: ~30 MB file (vs 264 MB raw), which fits in the repo.
#
# Run once from the repo root (after restoring the raw file):
#   Rscript scripts/utils/preprocess_oCRE_TFBS.R
#
# Output:
#   data/oCRE_ProBound_ParEndo_TFBS_predictions_slim.txt.gz

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

IN  <- "data/oCRE_ProBound_ParEndo_TFBS_predictions.txt.gz"
OUT <- "data/oCRE_ProBound_ParEndo_TFBS_predictions_slim.txt.gz"
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
