#!/usr/bin/env Rscript
# Step 6: First-pass MPRA analysis for PYS2 oCREv2
#
# QC, controls, replicate correlation, oCRE activity, phylogenetic analysis.
#
# Usage: Rscript 06_analyze_mpra.R <count_table_wo_zeroes_w_dups.txt.gz> <output_dir>
#
# Outputs:
#   {output_dir}/lib_saturation_calculations.txt  - sequencing saturation per library
#   {output_dir}/dna_barcode_recovery.pdf         - DNA BC recovery histogram
#   {output_dir}/all_class_expression_hex_v2.pdf  - DNA vs RNA hex plots by class
#   {output_dir}/control_expression.pdf           - control scatter plots
#   {output_dir}/control_expression_hist_w_on_target_norm.pdf - normalized control histograms
#   {output_dir}/cactus_biol_rep_corr_mpra_act.pdf - biological replicate correlation
#   {output_dir}/oCRE_act_hist_replicates.pdf     - oCRE activity distributions
#   {output_dir}/oCRE_cre_hist.pdf                - activity by CRE
#   {output_dir}/mouse_expression_bxplt.pdf       - mouse element activity boxplot
#   {output_dir}/PYS2_oCREv2_cactus_tiles_MPRA_mean_act.txt - mean oCRE activity table
#   {output_dir}/HepG2_MPRA_oCRE_act_table_wo_zeros_w_dups.txt - per-tile activity table

### module load pcre2/10.39; module load R/4.3.1

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(cowplot)
  library(ggrastr)
  library(castor)
  library(ape)
  library(phytools)
  library(reshape2)
  library(DescTools)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript 06_analyze_mpra.R <counts.txt.gz> <output_dir>")
}

count_file <- args[1]
out_dir    <- args[2]

# Source helper functions
script_dir <- file.path(dirname(sub("--file=", "", commandArgs()[grep("--file=", commandArgs())])), "scripts")
if (file.exists(file.path(script_dir, "mpra_tablemaker.R"))) {
  source(file.path(script_dir, "mpra_tablemaker.R"))
} else {
  source("scripts/mpra_tablemaker.R")
}

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
p_count <- 0.05

# --- Load data ---
cat("Loading:", count_file, "\n")
df_counts_w_CREs2 <- read.table(count_file, header = TRUE)

# ========================================
# Saturation calculations
# ========================================
coverage_lib <- calculate_MPRA_coverage(df_counts_w_CREs2)
write.table(coverage_lib, file.path(out_dir, "lib_saturation_calculations.txt"),
            sep = '\t', row.names = FALSE, quote = FALSE)
cat("Saturation:\n")
print(coverage_lib %>% select(lib_name, saturation, frac_BC_recovered))

# ========================================
# Wide table generation
# ========================================
df_wide <- generate_MPRA_wide_table(df_counts_w_CREs2, coverage_lib)

# ========================================
# BC recovery analysis
# ========================================
DNA_umi_thresh <- 2
df_BC_recovered <- df_counts_w_CREs2 %>%
  group_by(class, CRE_id, material, biol_rep) %>%
  summarize(n_BC_recovered = length(BC[UMI >= DNA_umi_thresh]), .groups = "drop")

bc_hist <- ggplot(df_BC_recovered %>% filter(material == "DNA")) +
  stat_bin(aes(x = n_BC_recovered, color = class),
           geom = "step", position = "identity", bins = 50) +
  facet_wrap(~biol_rep) +
  scale_x_log10() + scale_y_log10() + theme_cowplot()
ggsave(file.path(out_dir, "dna_barcode_recovery.pdf"),
       plot = bc_hist, height = 6, width = 10, dpi = 300, device = "pdf")

# ========================================
# DNA vs RNA expression hex plots
# ========================================
all_hex <- ggplot(df_wide) +
  stat_bin2d(aes(x = UMI_DNA + p_count, y = UMI_RNA + p_count)) +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  geom_abline(intercept = log10(30)) +
  geom_abline(intercept = log10(1/20), linetype = "dashed") +
  scale_fill_gradient() +
  facet_grid(biol_rep ~ class) + theme_cowplot() +
  annotation_logticks()
ggsave(file.path(out_dir, "all_class_expression_hex_v2.pdf"),
       plot = all_hex, height = 8, width = 16, dpi = 300, device = "pdf")

# ========================================
# Controls
# ========================================
control_oi_wide <- df_wide %>% filter(class %in% c("pos_control", "neg_control"))

controls_plot <- ggplot(control_oi_wide) +
  geom_point_rast(aes(x = UMI_DNA + p_count, y = UMI_RNA + p_count, color = class), alpha = 0.7) +
  scale_color_manual(values = c("#006ab5", "#e4181b")) +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  geom_abline(intercept = log10(30)) +
  geom_abline(intercept = log10(1/20), linetype = "dashed") +
  facet_wrap(~biol_rep) + theme_cowplot() +
  annotation_logticks()
ggsave(file.path(out_dir, "control_expression.pdf"),
       plot = controls_plot, height = 4, width = 10, dpi = 300, device = "pdf")

options(scipen = 999)
controls_hist <- ggplot(control_oi_wide) +
  stat_bin(aes(x = ((UMI_RNA / on_target_UMI_RNA) + p_count) / ((UMI_DNA / on_target_UMI_DNA) + p_count),
               color = class),
           geom = "step", position = "identity", bins = 60) +
  scale_color_manual(values = c("#006ab5", "#e4181b")) +
  scale_x_log10(limits = c(0.005, 200), expand = expansion(0, 0)) +
  annotation_logticks(sides = "b") +
  scale_y_continuous(expand = expansion(0, 0.1)) +
  facet_wrap(~biol_rep) + theme_cowplot()
ggsave(file.path(out_dir, "control_expression_hist_w_on_target_norm.pdf"),
       plot = controls_hist, height = 6, width = 12, dpi = 300, device = "pdf")

# ========================================
# Winsorized MPRA activity
# ========================================
df_mpra_act <- generate_winsorize_MPRA_table(df_wide)
oCRE_winsor_act <- df_mpra_act %>% filter(class == "oCRE")

# ========================================
# Biological replicate correlation
# ========================================
oCRE_winsor_cast <- dcast(oCRE_winsor_act, CRE_id ~ biol_rep, value.var = "MPRA_act")
oCRE_winsor_cast <- oCRE_winsor_cast[, -1]
colnames(oCRE_winsor_cast) <- paste0("rep", colnames(oCRE_winsor_cast))
rep_corrs <- data.frame(
  replicates = c("rep1-rep2", "rep2-rep3", "rep3-rep1"),
  corr = c(cor(oCRE_winsor_cast$rep1, oCRE_winsor_cast$rep2, method = "spearman"),
           cor(oCRE_winsor_cast$rep2, oCRE_winsor_cast$rep3, method = "spearman"),
           cor(oCRE_winsor_cast$rep1, oCRE_winsor_cast$rep3, method = "spearman"))
)

plot_pair <- function(replicate1, replicate2, corr) {
  rep_comp <- paste0(replicate1, "-", replicate2)
  corr <- corr %>% filter(replicates == rep_comp)
  ggplot(oCRE_winsor_cast, aes_string(x = replicate1, y = replicate2)) +
    geom_point_rast(alpha = 0.5) +
    geom_smooth(method = "lm", col = "blue", se = FALSE) +
    labs(title = paste(replicate1, "vs", replicate2),
         x = paste(replicate1, "\nMPRA act."),
         y = paste(replicate2, "\nMPRA act."),
         subtitle = paste("Correlation:", round(corr$corr, 3))) +
    scale_x_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) + annotation_logticks() +
    theme_classic() + theme(text = element_text(size = 16))
}

p1 <- plot_pair("rep1", "rep2", rep_corrs)
p2 <- plot_pair("rep2", "rep3", rep_corrs)
p3 <- plot_pair("rep3", "rep1", rep_corrs)

corr_plots <- plot_grid(p1, p2, p3, ncol = 3, align = "vh")
ggsave(file.path(out_dir, "cactus_biol_rep_corr_mpra_act.pdf"),
       height = 4, width = 12, plot = corr_plots, device = "pdf", dpi = 300)

# Mean activity per oCRE tile
oCRE_mean_act <- oCRE_winsor_act %>%
  group_by(CRE_id) %>%
  summarise(mean_MPRA_act = mean(MPRA_act),
            sd_MPRA_act = sd(MPRA_act))
write.table(oCRE_mean_act, file.path(out_dir, "PYS2_oCREv2_cactus_tiles_MPRA_mean_act.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# ========================================
# oCRE tile-level analysis with phylogenetic distances
# ========================================
cactus_tree <- read.tree("/net/shendure/vol8/projects/tli/ucsc_cactus/241-mammalian-2020v2.phast-242.nh")
cactus_dist <- cophenetic(cactus_tree)
cactus_dist <- as.data.frame(cactus_dist[, c("Mus_musculus")])
colnames(cactus_dist) <- c("distance")

winsor_cut <- 0.01
DNA_UMI_thresh <- 2

# Filter to oCRE class
df_oi <- df_counts_w_CREs2 %>% filter(class == "oCRE") %>% select(-lib_name)

df_oi_wide <- df_oi %>%
  pivot_wider(id_cols = c(BC, CRE_id, class, biol_rep),
              values_from = c(UMI, reads),
              names_from = material)
df_oi_wide2 <- df_oi_wide %>% filter(!is.na(CRE_id))
df_oi_wide2[is.na(df_oi_wide2)] <- 0

# Parse species from CRE_id (format: {full_CRE_id}__{species}__{coords})
df_oi_wide2 <- df_oi_wide2 %>%
  mutate(
    full_CRE_id = sub("__.*", "", CRE_id),
    species = sub("_20240416\\.fa$", "", sub("^[^_]+__", "", sub("__[^_]+$", "", CRE_id)))
  )

# Tile-level winsorized activity
df_tiled_mean_act <- df_oi_wide2 %>%
  filter(UMI_DNA >= DNA_UMI_thresh) %>%
  group_by(CRE_id, full_CRE_id, species, biol_rep) %>%
  summarize(
    n_BC = length(BC),
    sum_winsor_act = sum(Winsorize(UMI_RNA, probs = c(0, 1 - winsor_cut))) /
                     sum(Winsorize(UMI_DNA, probs = c(0, 1 - winsor_cut))),
    .groups = "drop"
  )

# Add phylogenetic distances
mm <- match(df_tiled_mean_act$species, rownames(cactus_dist))
df_tiled_mean_act$distance <- cactus_dist$distance[mm]

write.table(df_tiled_mean_act,
            file.path(out_dir, "HepG2_MPRA_oCRE_act_table_wo_zeros_w_dups.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# ========================================
# oCRE activity histograms by replicate
# ========================================
oCRE_act <- ggplot(df_tiled_mean_act) +
  stat_bin(aes(x = sum_winsor_act + p_count, color = biol_rep),
           geom = "step", position = "identity", bins = 50) +
  facet_wrap(~full_CRE_id, ncol = 4) +
  scale_x_log10() + scale_y_log10() + theme_cowplot() + annotation_logticks()
ggsave(file.path(out_dir, "oCRE_cre_hist.pdf"),
       plot = oCRE_act, height = 6, width = 16, dpi = 300, device = "pdf")

# ========================================
# Mouse element activity boxplot
# ========================================
mouse_mean_act <- df_tiled_mean_act %>% filter(species == "Mus_musculus")

mouse_bxplot <- ggplot(mouse_mean_act,
                       aes(x = full_CRE_id, y = sum_winsor_act + p_count, fill = full_CRE_id)) +
  geom_boxplot() +
  labs(x = "CRE", y = "RNA/DNA") +
  scale_y_log10() + annotation_logticks(sides = "l") +
  geom_hline(yintercept = 10, linetype = "dashed", color = "#e4181b") +
  geom_hline(yintercept = 0.02, linetype = "dashed", color = "#006ab5") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "none")
ggsave(file.path(out_dir, "mouse_expression_bxplt.pdf"),
       plot = mouse_bxplot, height = 5, width = 5, dpi = 300, device = "pdf")

# ========================================
# Activity vs phylogenetic distance
# ========================================
oCRE_dist <- ggplot(df_tiled_mean_act) +
  geom_point_rast(aes(x = distance, y = sum_winsor_act + p_count, color = biol_rep)) +
  labs(x = "Phylogenetic distance (from mouse)") +
  facet_grid(rows = vars(biol_rep), cols = vars(full_CRE_id)) +
  scale_y_log10() + theme_cowplot() + annotation_logticks(sides = "l")
ggsave(file.path(out_dir, "oCRE_cre_dist.pdf"),
       plot = oCRE_dist, height = 8, width = 16, dpi = 300, device = "pdf")

cat("\nDone. Outputs in:", out_dir, "\n")
