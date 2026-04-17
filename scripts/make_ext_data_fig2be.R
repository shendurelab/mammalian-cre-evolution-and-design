### make_ext_data_fig2be.R
###
### Reproducible generation of extended-data figure 2b/e for paper 1.
###
### Usage (from the repo root):
###   Rscript scripts/make_ext_data_fig2be.R
###
### Inputs are resolved relative to REPO_ROOT. Two large data files are
### NOT included in the repo (see .gitignore) and must be restored to
### data/ before running:
###     data/oCRE_and_CRE_optimization_MPRA_count_table.txt.gz  (~62 MB)
###     data/oCRE_ProBound_ParEndo_TFBS_predictions.txt.gz      (~264 MB)
###
### Outputs are written to REPO_ROOT/figures/.

library(stringr)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(ggrastr)
library(castor)
library(ape)
library(phytools)
library(reshape2)
library(ggtree)
library(readxl)
library(tidytree)
library(ggtreeExtra)
library(IRanges)
library(seqinr)
library(scales)
library(ggnewscale)
library(patchwork)
library(ggplotify)

# ---- Resolve repo paths robustly (runnable from REPO_ROOT or scripts/) ----
find_repo_root <- function() {
  for (d in c(".", "..", dirname(sys.frame(1)$ofile %||% ""))) {
    if (dir.exists(file.path(d, "data")) && dir.exists(file.path(d, "scripts"))) {
      return(normalizePath(d))
    }
  }
  stop("Could not locate repo root (expected data/ and scripts/ directories).")
}
`%||%` <- function(x, y) if (is.null(x)) y else x
REPO_ROOT <- Sys.getenv("REPO_ROOT", unset = find_repo_root())
# The util scripts below use relative paths like "../data/...", so run from scripts/.
setwd(file.path(REPO_ROOT, "scripts"))
DATA_DIR  <- file.path(REPO_ROOT, "data")
UTIL_DIR  <- file.path(REPO_ROOT, "scripts", "utils")
OUT_DIR   <- file.path(REPO_ROOT, "figures")
dir.create(OUT_DIR, showWarnings = FALSE)

source(file.path(UTIL_DIR, "mpra_tablemaker.R"))
source(file.path(UTIL_DIR, "make_tree_functions.R"))
source(file.path(UTIL_DIR, "trajectory_utils.R"))
# figure1_func.R and figure2_func.R are identical, so we reuse figure1_func.R
source(file.path(UTIL_DIR, "figure1_func.R"))
source(file.path(UTIL_DIR, "load_figure1_2_data.R"))

# https://joeystanley.com/blog/custom-themes-in-ggplot2
pub_theme <- function() {
  theme_cowplot(font_size = 6) %+replace%
    theme(
      strip.background = element_rect(color='black', 
                                      fill='grey90',
                                      size=0.5, 
                                      linetype="solid",
                                      linewidth = 0.25),
      text = element_text(size=6, family='Helvetica'),
      strip.text = element_text(family = "Helvetica", size=6, face='bold', margin=margin(b=2.5, t=2.5)),
      axis.text.x = element_text(size=6, color = "black"),
      axis.text.y = element_text(size=6, color = "black"),
      legend.text = element_text(size=6),
      legend.title = element_text(size=6),
      axis.line = element_line(linewidth = 0.25),
      axis.title = element_text(size=7),
      legend.position = 'bottom',
      legend.title.align = 0.5,
      plot.tag = element_text(size = 10, face='bold'),
      plot.title = element_text(hjust = 0.5, face='plain', size=9,
                                margin=margin(b=10))
    )
}
theme_set(pub_theme())


cor.pred <- oCRE_mean_act %>% group_by(full_CRE_id) %>% summarise(rho = cor(mean_MPRA_act, log2_norm_score, method = 'spearman'))
cor.pred$full_CRE_id <- factor(cor.pred$full_CRE_id, levels = c('Gata4_chr14_5729','Epas1_chr17_10063','Lama1_chr17_7784','Sparc_chr11_7211','Bend5_chr4_8201'))
cor.pred
# full_CRE_id       rho
#   <chr>                <dbl>
# 1 Bend5_chr4_8201      0.353
# 2 Epas1_chr17_10063    0.643
# 3 Gata4_chr14_5729     0.364
# 4 Lama1_chr17_7784     0.615
# 5 Sparc_chr11_7211     0.368
big.tree <- as_tibble(tree)
big.tree <- groupClade(big.tree, c(293, 247,353,383,429,346,459))
group.order <- data.frame(group = seq(1:7),
                         Order = c('RODENTIA','PRIMATES',
                                  'CHIROPTERA','CETARTIODACTYLA','CARNIVORA',
                                  'EULIPOTYPHLA','PERISSODACTYLA'))
group.order <- bind_rows(group.order, data.frame(group = 0, Order = 'OTHERS'))
group.order$group <- factor(group.order$group)
big.tree <- left_join(big.tree, group.order, by = c('group'))

clade_oCRE_mean_act <- left_join(oCRE_mean_act, as.data.frame(big.tree) %>% dplyr::select(species = label, Order), by = 'species')               
clade_oCRE_mean_act <- clade_oCRE_mean_act %>% mutate(genome = case_when(
                        startsWith(species, 'full') ~ 'ancestor',
                        TRUE ~ 'extant'))
clade_oCRE_mean_act$full_CRE_id <- factor(clade_oCRE_mean_act$full_CRE_id, levels = c('Gata4_chr14_5729','Epas1_chr17_10063','Lama1_chr17_7784','Sparc_chr11_7211','Bend5_chr4_8201'))
df_mouse_tiles_mpra$full_CRE_id <- factor(df_mouse_tiles_mpra$full_CRE_id, levels = c('Gata4_chr14_5729','Epas1_chr17_10063','Lama1_chr17_7784','Sparc_chr11_7211','Bend5_chr4_8201'))

clade_oCRE_mean_act <- clade_oCRE_mean_act %>% 
    left_join(cor.pred, by = 'full_CRE_id') %>% 
    mutate(pred_label = paste0(full_CRE_id, "\nrho = ", round(rho, 3)))  # use for labeling only

corr.plot <- ggplot() + 
    geom_point_rast(data = clade_oCRE_mean_act, aes(x = log2_norm_score, y = mean_MPRA_act, color = Order, shape = genome, alpha = genome), 
                                        size = 1, stroke = 1.1) +
    scale_alpha_manual(values = c(0.3, 0.8)) + 
    geom_rect(data = control_mean_act %>% filter(CRE_id == 'minP'), aes(ymin = mean_MPRA_act-sd_MPRA_act, ymax = mean_MPRA_act+sd_MPRA_act), xmin = -Inf, xmax = Inf,
           alpha = .4,fill = "blue")  +
    geom_hline(yintercept = control_mean_act %>% filter(CRE_id == 'minP') %>% pull(mean_MPRA_act) ,color = "blue", linetype = 'dashed') + 
    geom_rect(data = df_mouse_tiles_mpra, aes(ymin = mean_MPRA_act-sd_MPRA_act, ymax = mean_MPRA_act+sd_MPRA_act), xmin = -Inf, xmax = Inf,
           alpha = .4,fill = "#008000") + 
    geom_hline(data = df_mouse_tiles_mpra, aes(yintercept = mean_MPRA_act), color = "#008000", linetype = 'dashed') + 
#     scale_x_continuous(limits = c(40, 100), breaks = c(40,60, 80, 100)) +
    scale_y_log10(
       breaks = c(10^-2,10^-1, 10^0, 10^1), limits = c(0.01, 20),
       labels = scales::trans_format("log10", scales::math_format(10^.x))
     ) + annotation_logticks(sides="l") + scale_color_manual(values = order.colors) + theme_classic() + 
#     geom_text(data = cor.pred, 
#                           aes(label = paste0("rho = ", round(rho, 3))),
#                           x= Inf, y = Inf,
#                           hjust = 1,  # Align to the right
#                         vjust = 1.5, size = 5) +   # Align to the bottom  
    facet_wrap(~full_CRE_id, ncol = 5, labeller = as_labeller(setNames(clade_oCRE_mean_act$pred_label, clade_oCRE_mean_act$full_CRE_id))) + 
    scale_shape_manual(values = c(23, 21), guide = guide_legend(override.aes = list(size = 2, color = 'black'))) +
    labs(y = 'MPRA activity', x = 'predicted log2(accessibility)') + 
    theme(legend.position = 'none')
ggsave(file.path(OUT_DIR, 'ext_data_fig2b.pdf'), plot = corr.plot, width = 210, height = 50, useDingbats = F, units = 'mm', device = 'pdf')

### load DMS data
mouse_DMS_MPRA$mut <- toupper(mouse_DMS_MPRA$mut)

CRE_oi <- "Epas1_chr17_10063"
fasta_file <- file.path(DATA_DIR, "oCRE_fasta", "Epas1_chr17_10063_oCRE.fasta")

keep.tips <- c('Mus_musculus','Lycaon_pictus')
epas1_msa_df <- make_alignment_dataframe(fasta_file, keep.tips, 'Mus_musculus',CRE_oi)
mouse_msa_df <- epas1_msa_df %>% filter(species == 'Mus_musculus')

epas1_msa_list <- make_msa_tfbs_traj_df(fasta_file, rev(keep.tips), 'Mus_musculus', CRE_oi, cactus.tfbs, mouse_func_tfbs)
labels <- collect_labels(epas1_msa_list$msa_df_TFBS, side = T)

epas1_func_tfbs <- mouse_func_tfbs %>% filter(CRE == CRE_oi)
mm.start <- match(epas1_func_tfbs$TFBS_start, mouse_msa_df$ref_pos)
mm.end <- match(epas1_func_tfbs$TFBS_end, mouse_msa_df$ref_pos)
epas1_func_tfbs$msa_start <- mouse_msa_df$curpos[mm.start]
epas1_func_tfbs$msa_end <- mouse_msa_df$curpos[mm.end]
epas1_msa_df <- left_join(epas1_msa_df, mouse_DMS_MPRA %>% filter(CRE == CRE_oi) %>%
                          dplyr::select(ref_pos = pos, mut, log2FC_MPRA, sd_log2FC_MPRA, MPRA_act, sd_MPRA_act, min_act,max_act, wilcox_BH),
                         by = c("ref_pos","mut"))
triplet_cluster_msa <- epas1_msa_df %>% filter(ref_pos >= 214) %>% filter(ref_pos <= 254)
triplet_cluster_msa$species <- factor(triplet_cluster_msa$species, levels = rev(keep.tips))

labels <- c()

for(species_oi in levels(triplet_cluster_msa$species)){
    if(species_oi%in%select.images$label){
        img_path <- select.images %>% filter(label == species_oi) %>% pull(uid)
        name_oi <- select.images %>% filter(label == species_oi) %>% pull(name)
        labels <- c(labels, paste0(
                      "<span style='vertical-align: middle;'>*", name_oi, "*</span> ",
                      "<img src='", img_path, "' width='25' style='vertical-align: middle;'/>"
                    ))
    }
}

triplet_cluster_msa <- triplet_cluster_msa %>% mutate(mut2 = case_when(
                        mut == 'DEL' ~ '-',
                        TRUE ~ mut))

triplet_cluster_map <- epas1_func_tfbs %>% filter(TFBS_start >= 214) %>% filter(TFBS_end <= 254) %>% dplyr::select(msa_start, msa_end, TF_name)
triplet_cluster_map$log2FC_MPRA <- 0
triplet_cluster_map$mut2 <- NA
triplet_cluster_map$species <- NA
triplet_cluster_map$curpos <- NA

library(ggtext)
species_levels <- rev(keep.tips)  # one source of truth


triplet_cluster_plot <- ggplot(triplet_cluster_msa, aes(x = curpos, y = species, fill = log2FC_MPRA, label = mut2)) +
  geom_tile() +
  geom_text(aes(fontface = ifelse(mut_type == "Mismatch", "bold",
                                 ifelse(mut_type == "DEL", "bold", "plain")), alpha = ifelse(mut_type == "Mismatch", 1,
                             ifelse(mut_type == "DEL", 1, 0.9))),color = 'black',
           show.legend = FALSE) +
  scale_fill_gradient2(low = "dodgerblue2",
                             mid = "white",
                             high = "firebrick1",
                             midpoint = 0,
                             space = "Lab",
                             na.value = "gainsboro",
                             limits=c(-1.5,1.5),
                             breaks = c(-1.5, 0, 1.5),
                             oob=squish, name = expression(log[2] ~ "FC (MPRA)")) +
  theme_classic() + scale_y_discrete(limits = species_levels, labels = labels, expand = c(0,0)) +
  theme(
    panel.grid = element_blank(),
#     axis.text.x = element_blank(),
    axis.title = element_blank(),
    axis.text.y = element_markdown(color = "black", size = 9),
    axis.ticks.y = element_blank(), # Optionally remove ticks if necessary
    legend.position = 'bottom'
  ) +
    new_scale_color() +
  # Adding lines spanning the top of the plot based on msa_start and msa_end
  geom_segment(data = triplet_cluster_map, aes(x = msa_start, xend = msa_end, y = 3.7, yend = 3.7, color = TF_name),
               size = 1) +
  # Optionally add labels for TF names
  geom_text(data = triplet_cluster_map, aes(x = (msa_start + msa_end) / 2, y = Inf, label = TF_name, color = TF_name),
            vjust = -1, size = 4) +
  scale_color_manual(values = col_TFs, guide = 'none')


### load entire affinity normalized across species tested
source(file.path(UTIL_DIR, "analyze_figure3_data.R"))


triplet.binding <- aff.mat %>% filter(species%in%keep.tips) %>%
    filter(TF_name2%in%c("Foxa2_3","Sox17_1","Gata4/6_4")) #%>% filter(msa_start >= 214) %>% filter(msa_end <= 254)
triplet.binding$species <- factor(triplet.binding$species, levels = rev(keep.tips))
triplet.binding$TF_name <- factor(triplet.binding$TF_name, levels = c("Foxa2",'Sox17',"Gata4/6"))

triplet_cluster_aff <- ggplot(triplet.binding, aes(x = TF_name, y = species, fill = hue_color)) +
    geom_point(shape = 21, color = "black", size = 3, stroke = 0.5) +  # black border
#   geom_text(
#     aes(label = round(norm_affinity, 2)),
#     color = "black"
#   ) +
  scale_fill_identity() +  # Use precomputed fill colors
#   scale_y_discrete(limits = species_levels, expand = c(0,0)) +  # same y scale, no padding
  theme_classic() + labs(x = "predicted\naffinity", y = "") + 
  # add vertical borders between columns
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(), # Optionally remove ticks if necessary
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(), # Optionally remove ticks if necessary
      axis.line = element_blank(),
      plot.margin = margin(5.5, 5.5, 5.5, 2),   # left margin tighter
      legend.position = 'none'
  )

layout <- '
    AAAAAAAB
    AAAAAAAB
    '

p <- triplet_cluster_plot + triplet_cluster_aff + plot_layout(design = layout)
ggsave(file.path(OUT_DIR, "ext_data_fig2e.pdf"), plot = p, width = 170, height = 60, useDingbats = F, units = 'mm', device = 'pdf')

