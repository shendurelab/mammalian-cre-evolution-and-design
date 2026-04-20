#!/usr/bin/env Rscript
# Step 3: Filter BC-CRE pileup into a clean dictionary
#
# Reads pileup file(s), removes polyG barcodes, filters by minimum read count,
# plots read-count distribution, and writes final BC list.
#
# Usage: Rscript 03_filter_dictionary.R <pileup_file.txt.gz> <min_reads> <output_BC_list.txt>

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript 03_filter_dictionary.R <pileup.txt.gz> <min_reads> <output.txt>")
}

pileup_file <- args[1]
min_reads   <- as.numeric(args[2])
output_file <- args[3]

cat("Reading pileup:", pileup_file, "\n")
df <- read.delim(pileup_file, sep = '\t') %>%
  group_by(CRE_id, BC1) %>%
  summarise(n_count = sum(n_count), .groups = "drop") %>%
  group_by(BC1) %>%
  mutate(n_prct = n_count / sum(n_count))

# Remove polyG barcodes
df <- df %>% filter(!grepl("^G+$", BC1))

cat("Total unique BCs:", length(unique(df$BC1)), "\n")

# Plot read count distribution
pdf_out <- sub("\\.[^.]+$", "_read_dist.pdf", output_file)
bc_dist <- ggplot(df) +
  stat_bin(aes(x = n_count), geom = "step", bins = 50) +
  scale_x_log10() + scale_y_log10() +
  annotation_logticks() + theme_classic() +
  labs(x = "Reads per BC-CRE pair", y = "Count",
       title = paste0("Min reads filter: ", min_reads))
ggsave(pdf_out, plot = bc_dist, height = 5, width = 5)
cat("Saved distribution plot:", pdf_out, "\n")

# Filter
filtered <- df %>% filter(n_count >= min_reads)
cat("After filtering (>=", min_reads, "reads):\n")
cat("  Unique CREs:", length(unique(filtered$CRE_id)), "\n")
cat("  Unique BCs: ", length(unique(filtered$BC1)), "\n")

# Write BC list
write.table(filtered[, c("BC1", "CRE_id")], output_file,
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("Wrote:", output_file, "\n")
