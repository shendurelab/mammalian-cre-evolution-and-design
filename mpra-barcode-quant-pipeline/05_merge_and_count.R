#!/usr/bin/env Rscript
# Step 5: Merge BC dictionaries and join with per-sample quantification
#
# Customized for PYS2 oCREv2 / CRE trajectory MPRA.
# Merges 6 BC dictionary sources (oCRE, CRE_traj, controls, full_CRE)
# and joins with per-sample BC quantification from Step 4.
#
# PYS2 sample naming: {lib}_{material}_{rep}_S{n}
#   e.g., oCREv2_DNA_1_S4 -> material=DNA, rep=1
#
# Usage: Rscript 05_merge_and_count.R <bc_quant_dir> <output_prefix>
#
# Outputs:
#   {output_prefix}_merged_dictionary.txt.gz     - deduplicated BC dictionary
#   {output_prefix}_count_table.txt.gz           - BC x sample count table
#   {output_prefix}_count_table_wo_zeroes.txt.gz - same, excluding undetected BCs
#   {output_prefix}_recovery_per_lib.txt.gz      - BC recovery stats per library

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript 05_merge_and_count.R <bc_quant_dir> <output_prefix>")
}

bc_quant_dir <- args[1]
out_prefix   <- args[2]

# --- 1. Load and merge dictionaries ---
cat("Merging BC dictionaries\n")

# Dictionary files are stored locally in data/dictionaries/
dict_dir <- file.path(dirname(sub("--file=", "", commandArgs()[grep("--file=", commandArgs())])), "data", "dictionaries")
if (!dir.exists(dict_dir)) dict_dir <- "data/dictionaries"

df_traj <- read.table(file.path(dict_dir, "CRE_traj_final_BC_list.txt.gz"),
                       header = TRUE)
df_oCRE <- read.table(file.path(dict_dir, "oCRE_v2_final_BC_list.txt.gz"),
                       header = TRUE)
df_engreitz <- read.table(file.path(dict_dir, "Engreitz_controls_final_BC_list.txt"),
                           header = TRUE)
df_min <- read.table(file.path(dict_dir, "minP_final_BC_list.txt.gz"),
                      header = TRUE)
df_EEF1A <- read.table(file.path(dict_dir, "EEF1aP_final_BC_list.txt.gz"),
                        header = TRUE)
df_full <- read.table(file.path(dict_dir, "fullCRE_BC_asso_10percent_complexity_lib.txt.gz"),
                       header = TRUE)

# Normalize columns and assign class labels
df_min2 <- df_min %>% transform(CRE_id = "minP", class = "neg_control") %>% select(BC, CRE_id, class)
df_EEF1A2 <- df_EEF1A %>% transform(CRE_id = "EEF1A1p", class = "pos_control") %>% select(BC, CRE_id, class)
df_full2 <- df_full %>% select(BC, CRE_id) %>% transform(class = "full_CRE")
df_engreitz2 <- df_engreitz %>% select(BC = BC1, CRE_id = CRE_id) %>% transform(class = "Engreitz_control")
df_traj2 <- df_traj %>% select(BC = BC1, CRE_id = CRE_id) %>% transform(class = "CRE_traj")
df_oCRE2 <- df_oCRE %>% select(BC = BC1, CRE_id) %>% transform(class = "oCRE")

df_all_CREs <- rbind(df_min2, df_EEF1A2, df_full2, df_engreitz2, df_oCRE2, df_traj2)

cat("Total BC entries:", nrow(df_all_CREs), "\n")
cat("Unique BCs:", length(unique(df_all_CREs$BC)), "\n")

# Remove BCs that map to multiple CREs (duplicates)
df_dedup <- df_all_CREs %>%
  transform(is_duplicate = (duplicated(BC) | duplicated(BC, fromLast = TRUE))) %>%
  filter(!is_duplicate)

cat("After deduplication:", nrow(df_dedup), "BCs\n")

cat("Writing out BC dictionaries\n")
dict_out <- paste0(out_prefix, "_merged_dictionary.txt.gz")
write.table(df_dedup, gzfile(dict_out), sep = "\t", row.names = FALSE, quote = FALSE)
cat("Wrote:", dict_out, "\n")

# --- 2. Load BC quantification and join ---
# PYS2 sample naming: oCREv2_DNA_1_S4 -> s[2]=material, s[3]=rep
quant_files <- Sys.glob(file.path(bc_quant_dir, "*_BC_quant_*.txt.gz"))
cat("\nFound", length(quant_files), "quantification file(s)\n")

df_all_CREs_w_counts <- data.frame()
for (name_oi in quant_files) {
  cat("  ", basename(name_oi), "\n")
  s <- strsplit(basename(name_oi), "_")[[1]]
  df_BC_counts <- read.table(name_oi, header = TRUE)
  nuc_type <- s[2]
  rep_type <- s[3]
  df_BC_counts2 <- df_BC_counts %>%
    select(BC = mBC, UMI = n_UMI_per_mBC, reads = n_reads_per_mBC) %>%
    transform(lib_name = paste0(nuc_type, "_", rep_type),
              material = nuc_type,
              biol_rep = rep_type)
  df_all_CREs_w_counts <- rbind(df_all_CREs_w_counts,
                                 df_dedup %>% left_join(df_BC_counts2, by = "BC"))
}

df_all_CREs_w_counts <- df_all_CREs_w_counts[complete.cases(df_all_CREs_w_counts), ]
cat("\nTotal rows in count table:", nrow(df_all_CREs_w_counts), "\n")

cat("Writing out BC count summaries...\n")

count_out <- paste0(out_prefix, "_count_table.txt.gz")
write.table(df_all_CREs_w_counts, gzfile(count_out),
            sep = "\t", row.names = FALSE, quote = FALSE)

count_nz <- paste0(out_prefix, "_count_table_wo_zeroes.txt.gz")
write.table(df_all_CREs_w_counts %>% filter(!is.na(UMI)), gzfile(count_nz),
            sep = "\t", row.names = FALSE, quote = FALSE)

# --- 3. Recovery stats ---
recovery <- df_all_CREs_w_counts %>%
  group_by(class, lib_name) %>%
  summarize(
    frac_BC_recovered = mean(!is.na(UMI)),
    frac_BC_recovered_2plus = sum(UMI >= 2, na.rm = TRUE) / length(UMI),
    .groups = "drop"
  ) %>% arrange(lib_name)

recovery_out <- paste0(out_prefix, "_recovery_per_lib.txt.gz")
write.table(recovery, gzfile(recovery_out), sep = "\t", row.names = FALSE, quote = FALSE)

cat("\nWrote:\n")
cat("  ", count_out, "\n")
cat("  ", count_nz, "\n")
cat("  ", recovery_out, "\n")
