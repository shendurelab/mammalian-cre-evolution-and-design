### make_ext_data_fig1.R
###
### Reproducible generation of extended-data figure panels for paper 1.
###
### Usage:
###   # from the repo root:
###   Rscript scripts/make_ext_data_fig1.R
###   # or from inside scripts/:
###   Rscript make_ext_data_fig1.R
###
### Required R packages: see library() calls below.
### Required data:
###   - Most inputs are resolved relative to REPO_ROOT/data/.
###   - Two large data files are NOT included in the repo (see .gitignore);
###     they must be restored to data/ before running:
###         data/oCRE_and_CRE_optimization_MPRA_count_table.txt.gz  (~62 MB)
###         data/oCRE_ProBound_ParEndo_TFBS_predictions.txt.gz      (~264 MB)
###   - One panel set (Gviz browser tracks -> ext_data_fig3_pII.pdf) depends
###     on large per-cell bigwigs that are also not in the repo. That panel
###     is skipped automatically if the bigwigs are missing; to regenerate
###     it, set the environment variable BIGWIG_DIR to a directory
###     containing, for cell in {endo, meso, ecto, pluri}:
###         {cell}.bw                     (observed scATAC signal)
###         {cell}_chrombpnet_nobias.bw   (chromBPNet predicted signal)

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
library(rhdf5)

# ---- Resolve repo paths robustly ----
# Allow running from either REPO_ROOT or REPO_ROOT/scripts.
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
EXTDATA_DIR <- file.path(REPO_ROOT, "extdata")
FIG_DATA  <- file.path(DATA_DIR, "ext_data_fig1")
UTIL_DIR  <- file.path(REPO_ROOT, "scripts", "utils")
OUT_DIR   <- file.path(REPO_ROOT, "figures")
dir.create(OUT_DIR, showWarnings = FALSE)

source(file.path(UTIL_DIR, "mpra_tablemaker.R"))
source(file.path(UTIL_DIR, "make_tree_functions.R"))
source(file.path(UTIL_DIR, "trajectory_utils.R"))
source(file.path(UTIL_DIR, "figure1_func.R"))
source(file.path(UTIL_DIR, "load_figure1_2_data.R"))
load(file.path(EXTDATA_DIR, "oCRE_func_tfbs_info.RData"))

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

### supp fig 8
cre_levels <- c('Gata4:chr14_5729','Epas1:chr17_10063','Lama1:chr17_7784','Sparc:chr11_7211','Bend5:chr4_8201')
filter_oCRE_func_tfbs.count <- oCRE_func_tfbs.count %>%
  filter(Order != 'OTHERS') %>% 
  mutate(
    n_func_aff_tfbs_scaled = pmin(n_func_aff_tfbs, 10),  # Cap at 95th percentile 
    genome = case_when(
      startsWith(species, 'full') ~ 'ancestor',
      TRUE ~ 'extant'
    ),
    CRE = sub("_", ":", CRE)  # Replace underscores with colons
  )

filter_oCRE_func_tfbs.count$CRE <- factor(filter_oCRE_func_tfbs.count$CRE, 
  levels = cre_levels)

filter_oCRE_func_tfbs.count$Order <- factor(filter_oCRE_func_tfbs.count$Order, 
  levels = c('RODENTIA','PRIMATES','CARNIVORA','CHIROPTERA',
             'CETARTIODACTYLA','PERISSODACTYLA','EULIPOTYPHLA'))

wt_mouse_tiles_mpra <- df_mouse_tiles_mpra %>%
  dplyr::select(CRE = full_CRE_id, species, mean_MPRA_act, sd_MPRA_act, 
                scores, norm_score, common_name, log2_norm_score) %>%
  mutate(CRE = sub("_", ":", CRE))  # Apply the same replacement here

wt_mouse_tiles_mpra$CRE <- factor(wt_mouse_tiles_mpra$CRE, 
  levels = cre_levels)

mouse.points <- oCRE_mean_act %>% 
  filter(species == 'Mus_musculus') %>% 
  dplyr::select(CRE = full_CRE_id, species, mean_MPRA_act, sd_MPRA_act, 
                scores, norm_score, common_name, log2_norm_score) %>%
  mutate(CRE = sub("_", ":", CRE))

mouse.points <- left_join(
  mouse.points, 
  oCRE_func_tfbs.count %>%
    dplyr::select(CRE, species, n_map_tfbs, n_map_aff_tfbs, n_func_tfbs, 
                  n_func_aff_tfbs, Order) %>%
    mutate(CRE = sub("_", ":", CRE)),
  by = c('CRE','species')
)

mouse.points$point_type <- 'WT-Mouse'
mouse.points$CRE <- factor(mouse.points$CRE, 
  levels = cre_levels)

count.corr <- filter_oCRE_func_tfbs.count %>%
  group_by(CRE) %>%
  summarise(rho = cor(mean_MPRA_act, n_func_aff_tfbs, method = 'spearman'), .groups = "drop")

count.corr$CRE <- factor(count.corr$CRE, levels = cre_levels)

filter_oCRE_func_tfbs.count <- filter_oCRE_func_tfbs.count %>%
  left_join(count.corr, by = "CRE") %>%
  mutate(count_label = paste0(CRE, "\nrho = ", round(rho, 3)))


func_tfbs_points <- ggplot() + 
        geom_quasirandom_rast(data = filter_oCRE_func_tfbs.count, aes(x = n_func_aff_tfbs, y = mean_MPRA_act, color = Order, shape = genome, alpha = genome), 
            dodge.width=1, size = 1, stroke = 1.1) +
        scale_alpha_manual(values = c(0.3, 0.8)) +
        scale_color_manual(values = order.colors) +
#         scale_fill_manual(values = c("ancestor" = "white", "extant" = "black")) +
        scale_shape_manual(values = c(23, 21), guide = guide_legend(override.aes = list(size = 2, color = 'black'))) +
#         scale_color_viridis_c(option = "cividis", breaks = c(0, 5, 10),
#                              labels = c("0","5","10+"), 
#                               name = expression("TFBS affinity ≥ WT TFBS affinity")) + # Blue to yellow (designed for colorblind users)
        geom_rect(data = wt_mouse_tiles_mpra, aes(ymin = mean_MPRA_act-sd_MPRA_act, ymax = mean_MPRA_act+sd_MPRA_act, color = NULL), xmin = -Inf, xmax = Inf,
                   alpha = .2, fill = "#008000", show.legend = FALSE) + 
        geom_hline(data = wt_mouse_tiles_mpra, aes(yintercept = mean_MPRA_act), color = '#008000', linetype = 'dashed', show.legend = FALSE) +
        geom_rect(data = control_mean_act %>% filter(CRE_id == 'minP'), aes(ymin = mean_MPRA_act-sd_MPRA_act, ymax = mean_MPRA_act+sd_MPRA_act, color = NULL), 
                  xmin = -Inf, xmax = Inf, show.legend = FALSE,alpha = .2,fill = "blue") +
        geom_hline(data = control_mean_act %>% filter(CRE_id == 'minP'), aes(yintercept = mean_MPRA_act), color = "blue", linetype = 'dashed', show.legend = FALSE) +
#         geom_label_repel(data = label.points, aes(label = common_name, x = n_func_aff_tfbs, mean_MPRA_act), 
#                         color = 'black', size = 3, box.padding = 1) +
        new_scale_fill() +
        geom_star(data = mouse.points, aes(x = n_func_aff_tfbs, y = mean_MPRA_act, fill = point_type, starshape = point_type), 
                  size = 3, show.legend = FALSE) +  
        scale_fill_manual(values = c('WT-Mouse' = '#cc0000')) + 
        scale_starshape_manual(values = c('WT-Mouse' = 1)) + 
        theme_classic() + 
        facet_wrap(~CRE, ncol = 5, scales = "free_x", labeller = as_labeller(setNames(filter_oCRE_func_tfbs.count$count_label, filter_oCRE_func_tfbs.count$CRE))) + 
        scale_y_log10(breaks = c(10^-1, 10^0, 10^1),
                   labels = scales::trans_format("log10", scales::math_format(10^.x)),
                          limits = c(0.05,15)) + annotation_logticks(sides = 'l') + 
        scale_x_continuous(breaks = scales::breaks_pretty(5)) + 
#         scale_x_continuous(labels = scales::number_format(accuracy = 1), limits = c(0,12)) +
        labs(x = "# of equivalent functional TFBS", y = "MPRA activity") + 
        theme(panel.spacing = unit(1, "lines"),
             legend.position = 'none')

library(ggbeeswarm)

order.species <- oCRE_sum_aff_tfbs.count %>%
  dplyr::select(species, common_name, Order) %>% distinct()

## build wide matrix; make CRE colon-style
aff_mat <- oCRE_sum_aff_tfbs.count %>%
  mutate(CRE = sub("_", ":", CRE)) %>%
  dcast(species + CRE ~ TF_name, value.var = 'sum_total_aff')

aff_mat[is.na(aff_mat)] <- 0

aff_mat <- left_join(
  aff_mat,
  oCRE_mean_act %>%
    dplyr::select(species, CRE = full_CRE_id, mean_MPRA_act, sd_MPRA_act) %>%
    mutate(CRE = sub("_", ":", CRE)),
  by = c("species","CRE")
)

aff_mat <- aff_mat[complete.cases(aff_mat),]

final_aff.mat <- data.frame()
glm.results   <- data.frame()

library(broom)
for (CRE_oi in unique(aff_mat$CRE)) {
  tmp_mat <- aff_mat %>% filter(CRE == CRE_oi)
  species.order <- tmp_mat %>% pull(species)
  cre.order     <- tmp_mat %>% pull(CRE)
  sd.order      <- tmp_mat %>% pull(sd_MPRA_act)

  tmp_mat <- tmp_mat %>% dplyr::select(-species, -CRE, -sd_MPRA_act, -Hnf1b)

  glm_model   <- glm(mean_MPRA_act ~ ., data = tmp_mat, family = Gamma(link = "log"))
  glm_summary <- tidy(glm_model)
  glm_summary$CRE <- CRE_oi
  glm.results <- bind_rows(glm.results, glm_summary)

  tmp_mat <- tmp_mat %>%
    mutate(
      pred_mod   = predict(glm_model),
      fitted_mod = fitted(glm_model)
    )

  tmp_mat$species     <- species.order
  tmp_mat$CRE         <- cre.order
  tmp_mat$sd_MPRA_act <- sd.order
  tmp_mat <- left_join(tmp_mat, order.species, by = c("species"))
  final_aff.mat <- bind_rows(final_aff.mat, tmp_mat)
}

aff.corr <- final_aff.mat %>%
  group_by(CRE) %>%
  summarise(rho = cor(log2(mean_MPRA_act), pred_mod, method = 'spearman'), .groups = "drop")

final_aff.mat <- final_aff.mat %>%
  mutate(genome = case_when(startsWith(species, 'full') ~ 'ancestor', TRUE ~ 'extant'))

final_aff.mat$CRE <- factor(final_aff.mat$CRE, levels = cre_levels)
aff.corr$CRE      <- factor(aff.corr$CRE,      levels = cre_levels)

oCRE_sum_aff_tfbs.count[is.na(oCRE_sum_aff_tfbs.count)] <- 0

subset_sum_aff_tfbs.count <- oCRE_sum_aff_tfbs.count %>%
  filter(Order != 'OTHERS') %>%
  mutate(CRE = sub("_", ":", CRE))

subset_sum_aff_tfbs.count$CRE <- factor(subset_sum_aff_tfbs.count$CRE, levels = cre_levels)

mouse_sum_aff_tfbs.count <- subset_sum_aff_tfbs.count %>%
  filter(species == 'Mus_musculus')

subset_total_sum_aff_tfbs.count <- subset_sum_aff_tfbs.count %>%
  group_by(species, CRE, Order, common_name) %>%
  summarise(
    sum_mapped_aff      = sum(sum_mapped_aff),
    sum_mapped_func_aff = sum(sum_mapped_func_aff),
    sum_total_aff       = sum(sum_total_aff),
    mean_MPRA_act       = max(mean_MPRA_act),
    sd_MPRA_act         = max(sd_MPRA_act),
    .groups = "drop"
  ) %>%
  mutate(genome = case_when(startsWith(species, 'full') ~ 'ancestor', TRUE ~ 'extant'))

mouse_total_sum_aff_tfbs.count <- subset_total_sum_aff_tfbs.count %>%
  filter(species == 'Mus_musculus')

final_aff.mat <- final_aff.mat %>%
  left_join(aff.corr, by = "CRE") %>%
  mutate(aff_label = paste0(CRE, "\nrho = ", round(rho, 3)))

mouse.points <- final_aff.mat %>% filter(species == 'Mus_musculus')
mouse.points$point_type <- 'WT-Mouse'
mouse.points$CRE <- factor(mouse.points$CRE, levels = cre_levels)

glm.results$CRE <- factor(glm.results$CRE, levels = cre_levels)
glm.results <- glm.results %>% filter(term != "(Intercept)")
glm.results <- glm.results %>%
  mutate(term = gsub("`", "", term))    # remove backticks
glm.results$term <- factor(glm.results$term, levels = final_TFs)
glm.results <- glm.results %>%
  mutate(est_sig = if_else(p.value < 0.05, "signf.", "not signf."))
glm.results <- glm.results %>% mutate(TF_name = case_when(
                term == 'Jun_Atf3' ~ 'AP-1',
                TRUE ~ term))
glm.results$TF_name <- factor(glm.results$TF_name, levels = c("AP-1","Foxa2","Gata4/6","Sox17","Klf4"))

glm.est <- ggplot(data=glm.results, aes(x=estimate, y=TF_name, fill = est_sig)) +
  geom_bar(stat="identity") +
  geom_vline(xintercept=0, color = 'black', linetype = 'dashed') + 
    scale_fill_manual(
    values = c("signf." = 'black', "not signf." = "grey80"),  # merge your TF colors with a gray fallback
      name = "effect significance") +
    facet_wrap(~CRE, ncol = 5) +
    labs(x = "Effect size (GLM coefficient, log scale)", y = "") + 
    theme_classic() + 
        theme(panel.spacing = unit(1, "lines"), legend.position = 'bottom')

glm.aff <- ggplot() + 
        geom_quasirandom_rast(data = final_aff.mat, aes(x = pred_mod, y = mean_MPRA_act, color = Order, shape = genome, alpha = genome),
                             size = 1, stroke = 1.1) + 
        scale_alpha_manual(values = c(0.3, 0.8)) + 
        geom_rect(data = mouse_total_sum_aff_tfbs.count, 
                  aes(ymin = mean_MPRA_act-sd_MPRA_act, ymax = mean_MPRA_act+sd_MPRA_act), xmin = -Inf, xmax = Inf,
               alpha = .2,fill = "#008000") +
        geom_hline(data = mouse_total_sum_aff_tfbs.count, aes(yintercept = mean_MPRA_act), color = "#008000", linetype = 'dashed') +

        geom_rect(data = control_mean_act %>% filter(CRE_id == 'minP'), aes(ymin = mean_MPRA_act-sd_MPRA_act, ymax = mean_MPRA_act+sd_MPRA_act, color = NULL), 
                  xmin = -Inf, xmax = Inf, show.legend = FALSE,alpha = .2,fill = "blue") +
        geom_hline(data = control_mean_act %>% filter(CRE_id == 'minP'), aes(yintercept = mean_MPRA_act), color = "blue", linetype = 'dashed', show.legend = FALSE) +
#         geom_text(data = aff.corr, 
#                           aes(label = paste0("rho = ", round(rho, 3))),
#                           x= -Inf, y = Inf,
#                           hjust = -0.3,  # Align to the right
#                         vjust = 1.5, size = 4) +   # Align to the bottom  
        geom_star(data = mouse.points, aes(x = pred_mod, y = mean_MPRA_act, fill = point_type, starshape = point_type), 
                  size = 3, show.legend = FALSE) +  
        scale_fill_manual(values = c('WT-Mouse' = '#cc0000')) + 
        scale_starshape_manual(values = c('WT-Mouse' = 1)) + 
        facet_wrap(~CRE, ncol = 5, labeller = as_labeller(setNames(final_aff.mat$aff_label, final_aff.mat$CRE))) + 
        labs(y = 'MPRA activity', x = expression("GLM predictions")) +
        scale_color_manual(values = order.colors, guide = guide_legend(override.aes = list(size = 3))) + 
        scale_shape_manual(values = c(23, 21), guide = guide_legend(override.aes = list(size = 2, color = 'black'))) +
        scale_y_log10(breaks = c(10^-1, 10^0, 10^1),
                   labels = scales::trans_format("log10", scales::math_format(10^.x)),
                          limits = c(0.05,15))+ annotation_logticks(sides = 'l')+ 
        theme_classic() + 
        theme(panel.spacing = unit(1, "lines"), legend.position = 'none')

#### gkmSVM predictions
load(file.path(EXTDATA_DIR, "oCRE_func_tfbs_info.RData"))
order.species <- oCRE_sum_aff_tfbs.count %>% dplyr::select(species, common_name, Order) %>% distinct()
oCRE_gkm_pred <- read.delim(file.path(DATA_DIR, "oCRE_gkmSVM_endo_prediction.txt"),
                           sep = '\t', header = F)
colnames(oCRE_gkm_pred) <- c("id", 'gkm_score')
oCRE_gkm_pred <- oCRE_gkm_pred %>% separate(id, into = c("full_CRE_id", 'species'), sep = '__')
oCRE_gkm_act <- oCRE_mean_act %>% left_join(oCRE_gkm_pred, by = c("full_CRE_id", 'species'))
mm <- match(oCRE_gkm_act$species, order.species$species)
oCRE_gkm_act$Order <- order.species$Order[mm]
oCRE_gkm_act.filter <- oCRE_gkm_act %>% filter(Order != "OTHERS")
oCRE_gkm_act.filter <- oCRE_gkm_act.filter %>% mutate(genome = case_when(
                        startsWith(species, 'full') ~ 'ancestor',
                        TRUE ~ 'extant'))

# After building oCRE_gkm_act but before gkm.corr:
oCRE_gkm_act.filter <- oCRE_gkm_act.filter %>%
  mutate(full_CRE_id = sub("_", ":", full_CRE_id))

# Now compute correlations on the colonized names
gkm.corr <- oCRE_gkm_act.filter %>%
  group_by(full_CRE_id) %>%
  summarise(rho = cor(mean_MPRA_act, gkm_score, method = 'spearman'), .groups = "drop")

# full_CRE_id           rho
#   <chr>               <dbl>
# 1 Bend5:chr4_8201   0.410  
# 2 Epas1:chr17_10063 0.153  
# 3 Gata4:chr14_5729  0.168  
# 4 Lama1:chr17_7784  0.574  
# 5 Sparc:chr11_7211  0.00771

# Label text
oCRE_gkm_act.filter <- oCRE_gkm_act.filter %>% 
  left_join(gkm.corr, by = "full_CRE_id") %>%
  mutate(gkm_label = paste0(full_CRE_id, "\nrho = ", round(rho, 3)))

mouse.points <- oCRE_gkm_act.filter %>% filter(species == 'Mus_musculus')
mouse.points$point_type <- 'WT-Mouse'
mouse.points$full_CRE_id <- factor(mouse.points$full_CRE_id, levels = cre_levels)

# Use same levels in the plot data too (optional but keeps facet order stable)
oCRE_gkm_act.filter$full_CRE_id <- factor(oCRE_gkm_act.filter$full_CRE_id, levels = cre_levels)

gkm.plot <- ggplot() + 
        geom_quasirandom_rast(data = oCRE_gkm_act.filter, aes(x = gkm_score, y = mean_MPRA_act, color = Order, shape = genome, alpha = genome),
                             size = 1, stroke = 1.1) + 
        scale_alpha_manual(values = c(0.3, 0.8)) + 
        geom_rect(data = mouse.points, 
                  aes(ymin = mean_MPRA_act-sd_MPRA_act, ymax = mean_MPRA_act+sd_MPRA_act), xmin = -Inf, xmax = Inf,
               alpha = .2,fill = "#008000") +
        geom_hline(data = mouse.points, aes(yintercept = mean_MPRA_act), color = "#008000", linetype = 'dashed') +

        geom_rect(data = control_mean_act %>% filter(CRE_id == 'minP'), aes(ymin = mean_MPRA_act-sd_MPRA_act, ymax = mean_MPRA_act+sd_MPRA_act, color = NULL), 
                  xmin = -Inf, xmax = Inf, show.legend = FALSE,alpha = .2,fill = "blue") +
        geom_hline(data = control_mean_act %>% filter(CRE_id == 'minP'), aes(yintercept = mean_MPRA_act), color = "blue", linetype = 'dashed', show.legend = FALSE) +
#         geom_text(data = aff.corr, 
#                           aes(label = paste0("rho = ", round(rho, 3))),
#                           x= -Inf, y = Inf,
#                           hjust = -0.3,  # Align to the right
#                         vjust = 1.5, size = 4) +   # Align to the bottom  
        geom_star(data = mouse.points, aes(x = gkm_score, y = mean_MPRA_act, fill = point_type, starshape = point_type), 
                  size = 3, show.legend = FALSE) +  
        scale_fill_manual(values = c('WT-Mouse' = '#cc0000')) + 
        scale_starshape_manual(values = c('WT-Mouse' = 1)) + 
        facet_wrap(~full_CRE_id, ncol = 5, labeller = as_labeller(setNames(oCRE_gkm_act.filter$gkm_label, oCRE_gkm_act.filter$full_CRE_id))) + 
        labs(y = 'MPRA activity', x = expression("gkm-SVM predictions")) +
        scale_color_manual(values = order.colors, guide = guide_legend(override.aes = list(size = 3))) + 
        scale_shape_manual(values = c(23, 21), guide = guide_legend(override.aes = list(size = 2, color = 'black'))) +
        scale_y_log10(breaks = c(10^-1, 10^0, 10^1),
                   labels = scales::trans_format("log10", scales::math_format(10^.x)),
                          limits = c(0.05,15))+ annotation_logticks(sides = 'l')+ theme_classic() + 
        theme(panel.spacing = unit(1, "lines"), legend.position = 'none')


p <- (func_tfbs_points + theme(legend.position = 'none')) / glm.aff / glm.est / gkm.plot + plot_annotation(tag_levels = 'a') &
  theme(plot.tag = element_text(size = 12, face='bold'))
ggsave(file.path(OUT_DIR, 'ext_data_fig1.pdf'), plot = p, width = 210, height = 200, useDingbats = F, units = 'mm', device = 'pdf')


### supp fig 9
# Load pre-extracted per-cell consensus & DA peak predictions.
# (Prepared from the large chromBPNet H5 files by
#  scripts/utils/preprocess_ext_data_fig1_h5.R.)

CELLS <- c("endo", "pluri", "meso", "ecto")

cell.corrs <- data.frame()
all.peaks  <- data.frame()
for (cell in CELLS) {
  f <- file.path(FIG_DATA, "consensus_peak_preds",
                 paste0(cell, "_cons_peaks_obs_pred.txt.gz"))
  obs.df <- read.delim(f, sep = '\t', header = TRUE)
  obs.df <- obs.df %>% mutate(cell_type = cell)
  all.peaks  <- bind_rows(all.peaks, obs.df)
  cell.corrs <- bind_rows(cell.corrs,
                          data.frame(cell_type = cell,
                                     rho = cor(obs.df$logcts, obs.df$pred_logcts, method = 'spearman')))
}

da.corrs     <- data.frame()
all.da_peaks <- data.frame()
for (cell in CELLS) {
  f <- file.path(FIG_DATA, "da_peak_preds",
                 paste0(cell, "_da_peaks_obs_pred.txt.gz"))
  obs.df <- read.delim(f, sep = '\t', header = TRUE)
  obs.df <- obs.df %>% mutate(cell_type = cell)
  all.da_peaks <- bind_rows(all.da_peaks, obs.df)
  da.corrs <- bind_rows(da.corrs,
                        data.frame(cell_type = cell,
                                   rho = cor(obs.df$logcts, obs.df$pred_logcts, method = 'spearman')))
}

shared_fill_scale <- scale_fill_viridis_c(
  trans = "log10", 
  limits = c(1, 1e5),  # ensure both plots use same range
  breaks = c(10^0,10^3, 10^5), 
  labels = scales::trans_format("log10", scales::math_format(10^.x)),
  name = "peak count"
)

full_cell_names <- data.frame(cell_type = c("endo","pluri","meso",'ecto'),
                    cell_name = c("Endoderm","Pluripotent","Mesoderm",'Ectoderm'))

all.peaks$cell_type <- factor(all.peaks$cell_type, levels = c("endo","pluri","meso",'ecto'))
mm <- match(all.peaks$cell_type, full_cell_names$cell_type)
all.peaks$cell_name <- full_cell_names$cell_name[mm]
all.peaks$cell_name <- factor(all.peaks$cell_name, levels = c("Endoderm","Pluripotent","Mesoderm",'Ectoderm'))
cell.corrs$cell_type <- factor(cell.corrs$cell_type, levels = c("endo","pluri","meso",'ecto'))
mm <- match(cell.corrs$cell_type, full_cell_names$cell_type)
cell.corrs$cell_name <- full_cell_names$cell_name[mm]
cell.corrs$cell_name <- factor(cell.corrs$cell_name, levels = c("Endoderm","Pluripotent","Mesoderm",'Ectoderm'))


all.da_peaks$cell_type <- factor(all.da_peaks$cell_type, levels = c("endo","pluri","meso",'ecto'))
mm <- match(all.da_peaks$cell_type, full_cell_names$cell_type)
all.da_peaks$cell_name <- full_cell_names$cell_name[mm]
all.da_peaks$cell_name <- factor(all.da_peaks$cell_name, levels = c("Endoderm","Pluripotent","Mesoderm",'Ectoderm'))
da.corrs$cell_type <- factor(da.corrs$cell_type, levels = c("endo","pluri","meso",'ecto'))
mm <- match(da.corrs$cell_type, full_cell_names$cell_type)
da.corrs$cell_name <- full_cell_names$cell_name[mm]
da.corrs$cell_name <- factor(da.corrs$cell_name, levels = c("Endoderm","Pluripotent","Mesoderm",'Ectoderm'))

all.peaks <- all.peaks %>% 
    left_join(cell.corrs, by = "cell_name") %>%
  mutate(peak_label = paste0(cell_name, "\nrho = ", round(rho, 3)))  # use for labeling only

all.da_peaks <- all.da_peaks %>% 
    left_join(da.corrs, by = "cell_name") %>%
  mutate(peak_label = paste0(cell_name, "\nrho = ", round(rho, 3)))  # use for labeling only


all_peaks.scatter <- ggplot() +
  stat_bin2d(data = all.peaks, aes(x = logcts, y = pred_logcts, fill = after_stat(count)), bins = 30) +
  shared_fill_scale + 
  theme_classic() +
  facet_wrap(~cell_name, ncol = 4, labeller = as_labeller(setNames(all.peaks$peak_label, all.peaks$cell_name))) +
  labs(x = "log2 (observed counts)", y = "log2 (predicted counts)") +
  theme(legend.position = 'bottom',
        panel.spacing = unit(1, "lines"))


all_da_peaks.scatter <- ggplot() +
  stat_bin2d(data = all.da_peaks, aes(x = logcts, y = pred_logcts, fill = after_stat(count)), bins = 30) +
  shared_fill_scale + 
  theme_classic() +
  facet_wrap(~cell_name, ncol = 4, labeller = as_labeller(setNames(all.da_peaks$peak_label, all.da_peaks$cell_name))) +
  labs(x = "log2 (observed counts)", y = "log2 (predicted counts)") +
  theme(legend.position = 'bottom',
        panel.spacing = unit(1, "lines"))


layout <- '
        AAAAAAAA
        BBBBBBBB
        '
p <- (all_peaks.scatter + ggtitle("Consensus peak predictions") + theme(plot.title = element_text(size = 10, face = "bold", hjust= 0.5))) + (all_da_peaks.scatter + ggtitle("Differential accessible (DA) peak predictions") + theme(plot.title = element_text(size = 10, face = "bold", hjust= 0.5))) + plot_layout(design = layout) + plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(size = 12, face='bold'))
ggsave(file.path(OUT_DIR, 'ext_data_fig3_partI.pdf'), plot  = p, width = 210, height = 160, useDingbats = F, units = 'mm', device = 'pdf')

# ============================================================================
# ext_data_fig3_pII: Gviz genome-browser tracks for 5 oCREs
# ============================================================================
# This section requires large per-cell bigwigs that are NOT in the repo.
# To regenerate it, set the BIGWIG_DIR environment variable to a directory
# containing, for cell in {endo, meso, ecto, pluri}:
#   {cell}.bw                     (observed merged scATAC signal)
#   {cell}_chrombpnet_nobias.bw   (chromBPNet predicted signal)
# Otherwise, this block is skipped with a warning.

BIGWIG_DIR <- Sys.getenv("BIGWIG_DIR", unset = "")
bw_cells <- c("endo","meso","ecto","pluri")
required_bw <- if (nzchar(BIGWIG_DIR)) {
  c(file.path(BIGWIG_DIR, paste0(bw_cells, ".bw")),
    file.path(BIGWIG_DIR, paste0(bw_cells, "_chrombpnet_nobias.bw")))
} else character(0)

if (nzchar(BIGWIG_DIR) && all(file.exists(required_bw))) {

library(GenomicRanges)
library(Gviz)
library(ggplotify)
library(rtracklayer)
library(BSgenome.Mmusculus.UCSC.mm10)
mm10 = BSgenome.Mmusculus.UCSC.mm10


cluster_config = read.csv(file.path(FIG_DATA, "clusters.csv"))
cluster_config$bigwig      <- file.path(BIGWIG_DIR, paste0(cluster_config$cluster, ".bw"))
cluster_config$pred_bigwig <- file.path(BIGWIG_DIR, paste0(cluster_config$cluster, "_chrombpnet_nobias.bw"))
cluster_config$colour = as.character(cluster_config$colour)
cluster_config$pred_colour <- c("#f6b26b","#93c47d","#e06666","#6fa8dc")
cluster_config$description = as.character(cluster_config$description)


cluster_config = cluster_config[!is.na(cluster_config$cluster), ]
rownames(cluster_config) = as.factor(cluster_config$cluster)
CLUSTERS <- rownames(cluster_config)

gene.bed <- read.delim(file.path(FIG_DATA, 'mm10_biomart_GeneExport.txt'), sep = '\t')
gene.bed <- gene.bed[gene.bed$Gene.type == 'protein_coding',]

refGene.mm10 <- import.gff(gzfile(file.path(DATA_DIR, 'references', 'mm10.refGene.gtf.gz')))
refGene.mm10 = keepStandardChromosomes(refGene.mm10, pruning.mode = "coarse")
gene_anno <- data.frame(refGene.mm10)
gene_anno <- gene_anno[!is.na(gene_anno$exon_number),]
gene_anno <- gene_anno[gene_anno$gene_name%in%gene.bed$Gene.name, ]

# rename some columns to match requirements
gene_anno <- gene_anno[,c('seqnames','start','end','strand','gene_id','gene_name','transcript_id','type','exon_id')]
colnames(gene_anno) <- c("chromosome", "start", "end", "strand", "gene", "symbol", "transcript", "feature", "exon")

gene_anno_subset = gene_anno[gene_anno$symbol=="Gata4",]
gene_anno_subset = gene_anno_subset[gene_anno_subset$feature=="start_codon",]
gene_anno_subset <- head(gene_anno_subset, n = 1)

all_peaks = read.table(file.path(FIG_DATA, "overlap.peaks.resolved.500.filtered_cmb_peaks.bed"))[,1:3]
colnames(all_peaks) = c("chr", "start", "end")
all_peaks = GRanges(all_peaks)

# one bigwig at a time, any number of peaks
get_matrix_from_bigwig <- function(bigwig_path, peak_set) {
    # ensure peak set has fixed width
    stopifnot(length(unique(width(peak_set)))==1)
    
    as.matrix(import(bigwig_path, 
      which=peak_set, as="NumericList"))
}

# calculate top and bottom percentiles for obs/pred for each state

obs.upper_lims = c()
obs.lower_lims = c()

pred.upper_lims = c()
pred.lower_lims = c()

NUM_SAMP = 1000

sampled_peaks = sample(all_peaks, NUM_SAMP)
OBS.BIGWIGS <- as.vector(cluster_config$bigwig)
PRED.BIGWIGS <- as.vector(cluster_config$pred_bigwig)

for (i in seq(1:length(OBS.BIGWIGS)))  {
    
    vals = as.vector(get_matrix_from_bigwig(OBS.BIGWIGS[[i]], resize(sampled_peaks, width=500, fix='center')))
    obs.upper_lims = c(obs.upper_lims, quantile(vals, 0.9999))
    obs.lower_lims = c(obs.lower_lims, quantile(vals, 0.001))
    
    vals = as.vector(get_matrix_from_bigwig(PRED.BIGWIGS[[i]], resize(sampled_peaks, width=500, fix='center')))
    
    pred.upper_lims = c(pred.upper_lims, quantile(vals, 0.997))
    pred.lower_lims = c(pred.lower_lims, quantile(vals, 0.001))
}

names(obs.upper_lims) <- CLUSTERS
names(obs.lower_lims) <- CLUSTERS
names(pred.upper_lims) <- CLUSTERS
names(pred.lower_lims) <- CLUSTERS

get_region_tracks <- function(chr, cluster_config, show_axis=T) {
    # NOTE: this does not perform any max aggregation
    
    bw_tracks = c()
    obs_bw_path = cluster_config["endo", "bigwig"]
    obs_bw_path <- as.character(obs_bw_path)
    
    obs_track = DataTrack(obs_bw_path, 
                             genome="mm10", 
                             chromosome = chr, 
                             name=sprintf("%s", 'parietal\nendo.'),
                             ylim=c(0,obs.upper_lims[['endo']]),
                             ##window=-1, 
                             type="hist",fontsize = 6,
                             col.histogram=cluster_config['endo', "colour"],
                             fill.histogram=cluster_config['endo', "colour"],
                             background.title = cluster_config['endo', "colour"], littleTicks=F)
    
    displayPars(obs_track)$showAxis = FALSE
    if (!show_axis) {
        displayPars(obs_track)$showTitle = FALSE
    } 
    bw_tracks = c(bw_tracks, obs_track)
    
    for (i in CLUSTERS) {
        pred_bw_path <- cluster_config[i, "pred_bigwig"]
        pred_bw_path <- as.character(pred_bw_path)
        
        pred_track = DataTrack(pred_bw_path, 
                             genome="mm10", 
                             chromosome = chr, 
                             name=sprintf("%s", paste0(cluster_config[i, "cluster"],"\npred.")),
                             ylim=c(0,pred.upper_lims[[i]]),
                             ##window=-1, 
                             type="hist",fontsize = 6,
                             col.histogram=cluster_config[i, "pred_colour"],
                             fill.histogram=cluster_config[i, "pred_colour"],
#                              background.title = "white",
                             background.title = cluster_config[i, "pred_colour"],
                               littleTicks=F)
        # don't show axis ticks for each plot
        displayPars(pred_track)$showAxis = FALSE
        if (!show_axis) {
            displayPars(pred_track)$showAxis = FALSE
        } 
        bw_tracks = c(bw_tracks, pred_track) 
    } 
#     bw_tracks = bw_tracks[order(cluster_config$cluster)] # reorder
    
    if (!show_axis) {
            displayPars(gene_track)$showAxis = F
    }
    
    return(bw_tracks)
}

get_gene_track <- function(chr, gene_anno, show_axis=T) {
    gene_track <- GeneRegionTrack(gene_anno, genome = "mm10", 
                                 chromosome = chr, 
                                 name = "", 
                                 collapseTranscripts="longest",
                                 transcriptAnnotation="symbol",
                                background.title = 'transparent',
                                  background.panel = "#FFFEDB",
                                fill='#000000',
                                stackHeight=0.5)
    if (!show_axis) {
            displayPars(gene_track)$showTitle = F
        }
    
    return(gene_track)
}
### coordinates for max mouse tiles
max.tiles <- read.delim(file.path(FIG_DATA, "DMS_300bp_tiles.bed"),
                       sep = '\t', header = F)
colnames(max.tiles) <- c("chr",'start','end','CRE_ID')

#### for Gata4 CRE
CHR = "chr14"
FROM = 63230939 - 500
TO = 63231783 + 500
HIGHLIGHT_COL = "black"
MAX_COL = "red"
MAX_COORDS = max.tiles %>% filter(CRE_ID == 'Gata4_chr14_5729')


bw_tracks = get_region_tracks(CHR, cluster_config, show_axis=T)
gene_track = get_gene_track(CHR, gene_anno, show_axis=T)
#### Gata4 CRE max tile
bw_highlighted = HighlightTrack(trackList = c(bw_tracks),
                               start=c(MAX_COORDS$start, 63230939), 
                               end=c(MAX_COORDS$end, 63231783), 
                               chromosome=CHR,
                               col=c(MAX_COL, HIGHLIGHT_COL), fill=c(MAX_COL, HIGHLIGHT_COL), alpha = c(0.2, 0.1))

GATA4 = as.ggplot(~plotTracks(c(bw_highlighted), from=FROM, to=TO,
                             sizes=c(rep(0.2,4),0.1)))+
       ggtitle("Gata4:chr14_5729") + theme(plot.title = element_text(size = 12, face = "bold"))

# ggsave('test.pdf', height = 4, width = 8, plot = GATA4, useDingbats = F)
# Rasterize only the plot elements from plotTracks()
# GATA4_rasterized <- rasterize(GATA4, dpi = 500) # Adjust DPI as needed

#### for Epas1 CRE
CHR = "chr17"
FROM = 86740868 - 500
TO = 86741862 + 500
HIGHLIGHT_COL = "black"
MAX_COL = "red"
MAX_COORDS = max.tiles %>% filter(CRE_ID == 'Epas1_chr17_10063')

bw_tracks = get_region_tracks(CHR, cluster_config, show_axis=T)
gene_track = get_gene_track(CHR, gene_anno, show_axis=T)
#### Gata4 CRE max tile
bw_highlighted = HighlightTrack(trackList = c(bw_tracks),
                               start=c(MAX_COORDS$start,86740868), 
                               end=c(MAX_COORDS$end, 86741862), 
                               chromosome=CHR,
                               col=c(MAX_COL, HIGHLIGHT_COL), fill=c(MAX_COL, HIGHLIGHT_COL), alpha = c(0.2, 0.1))

EPAS1 = as.ggplot(~plotTracks(c(bw_highlighted), from=FROM, to=TO,
                             sizes=c(rep(0.2,4),0.1)))+
       ggtitle("Epas1:chr17_10063") + theme(plot.title = element_text(size = 12, face = "bold"))

# Rasterize only the plot elements from plotTracks()
# EPAS1_rasterized <- rasterize(EPAS1, dpi = 500) # Adjust DPI as needed
#### for Lama1 CRE
CHR = "chr17"
FROM = 67653831 - 500
TO = 67654710 + 500
HIGHLIGHT_COL = "black"
MAX_COL = "red"
MAX_COORDS = max.tiles %>% filter(CRE_ID == 'Lama1_chr17_7784')

bw_tracks = get_region_tracks(CHR, cluster_config, show_axis=T)
gene_track = get_gene_track(CHR, gene_anno, show_axis=T)
#### Lama1 CRE max tile
bw_highlighted = HighlightTrack(trackList = c(bw_tracks),
                               start=c(MAX_COORDS$start, 67653831), 
                               end=c(MAX_COORDS$end, 67654710), 
                               chromosome=CHR,
                               col=c(MAX_COL, HIGHLIGHT_COL), fill=c(MAX_COL, HIGHLIGHT_COL), alpha = c(0.2, 0.1))

LAMA1 = as.ggplot(~plotTracks(c(bw_highlighted), from=FROM, to=TO,
                              sizes=c(rep(0.2,4),0.1)))+
       ggtitle("Lama1:chr17_7784") + theme(plot.title = element_text(size = 12, face = "bold"))

#### for Sparc CRE
CHR = "chr11"
FROM = 55428316 - 500
TO = 55429205 + 500
HIGHLIGHT_COL = "black"
MAX_COL = "red"
MAX_COORDS = max.tiles %>% filter(CRE_ID == 'Sparc_chr11_7211')

bw_tracks = get_region_tracks(CHR, cluster_config, show_axis=T)
gene_track = get_gene_track(CHR, gene_anno, show_axis=T)
#### Sparc CRE max tile
bw_highlighted = HighlightTrack(trackList = c(bw_tracks),
                               start=c(MAX_COORDS$start, 55428316), 
                               end=c(MAX_COORDS$end, 55429205), 
                               chromosome=CHR,
                               col=c(MAX_COL, HIGHLIGHT_COL), fill=c(MAX_COL, HIGHLIGHT_COL), alpha = c(0.2, 0.1))

SPARC = as.ggplot(~plotTracks(c(bw_highlighted), from=FROM, to=TO,
                             sizes=c(rep(0.2,4),0.1)))+
       ggtitle("Sparc:chr11_7211") + theme(plot.title = element_text(size = 12, face = "bold"))

#### for Bend5 CRE
CHR = "chr4"
FROM = 111483057 - 500
TO = 111483881 + 500
HIGHLIGHT_COL = "black"
MAX_COL = "red"
MAX_COORDS = max.tiles %>% filter(CRE_ID == 'Bend5_chr4_8201')


bw_tracks = get_region_tracks(CHR, cluster_config, show_axis=T)
gene_track = get_gene_track(CHR, gene_anno, show_axis=T)
#### Bend5 CRE max tile
bw_highlighted = HighlightTrack(trackList = c(bw_tracks),
                               start=c(MAX_COORDS$start,111483057), 
                               end=c(MAX_COORDS$end, 111483881), 
                               chromosome=CHR,
                               col=c(MAX_COL, HIGHLIGHT_COL), fill=c(MAX_COL, HIGHLIGHT_COL), alpha = c(0.2, 0.1))

BEND = as.ggplot(~plotTracks(c(bw_highlighted), from=FROM, to=TO,
                             sizes=c(rep(0.2,4),0.1)))+
       ggtitle("Bend5:chr4_8201") + theme(plot.title = element_text(size = 12, face = "bold"))

layout <- 'ABCDE'

p <- GATA4 + EPAS1 + LAMA1 + SPARC + BEND + plot_layout(design = layout)
ggsave(file.path(OUT_DIR, 'ext_data_fig3_pII.pdf'), plot = p, width = 210, height = 60, useDingbats = F, units = 'mm', device = 'pdf')

} else {
  # BIGWIG_DIR unset or required bigwigs missing
  message("Skipping ext_data_fig3_pII (Gviz browser tracks): ",
          "set BIGWIG_DIR to a directory containing the required bigwigs. ",
          "See the header of this script for details.")
}

#### holdout test set correlations
shared_fill_holdout_scale <- scale_fill_viridis_c(
  trans = "log10", 
  limits = c(1, 1e3),  # ensure both plots use same range
  breaks = c(10^0,10^1,10^2, 10^3), 
  labels = scales::trans_format("log10", scales::math_format(10^.x)),
  name = "peak count"
)

GKM <- file.path(FIG_DATA, "gkm_svm")

gkm.pos <- read.delim(file.path(GKM, 'gkm_length_pos_pred.txt'), sep = '\t', header = F)
colnames(gkm.pos) <- c("peak_id","pred")
gkm.pos$peak_id <- sub("::.*", "", gkm.pos$peak_id)
gkm.pos_acc <- read.delim(file.path(GKM, 'gkm_length_pos_test_acc.txt') , sep = '\t', header = F)
colnames(gkm.pos_acc) <- c("peak_id","size","covered","sum","mean0","mean")
gkm.pos <- left_join(gkm.pos, gkm.pos_acc, by = "peak_id")
gkm.DA <- read.delim(file.path(GKM, 'gkm_length_pos_test.DA_overlap.bed'), sep = '\t', header = F)
gkm.pos$class <- 'Positive'

gkm.neg <- read.delim(file.path(GKM, 'gkm_length_neg_pred.txt'), sep = '\t', header = F)
colnames(gkm.neg) <- c("peak_id","pred")
gkm.neg$peak_id <- sub("::.*", "", gkm.pos$peak_id)
gkm.neg_acc <- read.delim(file.path(GKM, 'gkm_length_neg_test_acc.txt') , sep = '\t', header = F)
colnames(gkm.neg_acc) <- c("peak_id","size","covered","sum","mean0","mean")
gkm.neg$peak_id <- gkm.neg_acc$peak_id
gkm.neg <- left_join(gkm.neg, gkm.neg_acc, by = "peak_id")
gkm.neg$class <- 'Negative'

gkm.total_counts <- bind_rows(gkm.pos, gkm.neg)
gkm.total_counts$log_sum <- log(gkm.total_counts$sum)
gkm.total_counts <- gkm.total_counts %>% mutate(DA_call = case_when(
                    peak_id%in%gkm.DA$V4 ~ TRUE,
                    TRUE ~ FALSE))
gkm.cor <- cor(gkm.total_counts$pred, gkm.total_counts$log_sum, method = 'spearman')
gkm_DA.cor <- cor(gkm.total_counts %>% filter(DA_call) %>% pull(pred), 
                  gkm.total_counts %>% filter(DA_call) %>% pull(log_sum), method = 'spearman')


pos.logcounts <- read.delim(file.path(GKM, 'endo_pos_chrombpnet_logcounts.txt.gz'), sep = '\t', header = TRUE)$pred
cmb.pos_acc <- read.delim(file.path(GKM, 'chrombpnet_length_pos_test_acc.txt'), sep = '\t', header = F)
colnames(cmb.pos_acc) <- c("peak_id","size","covered","sum","mean0","mean")
cmb.pos_acc$pred <- pos.logcounts
cmb.DA <- read.delim(file.path(GKM, 'chrombpnet_length_pos_test.DA_overlap.bed'), sep = '\t', header = F)
cmb.pos_acc$class <- 'Positive'

neg.logcounts <- read.delim(file.path(GKM, 'endo_neg_chrombpnet_logcounts.txt.gz'), sep = '\t', header = TRUE)$pred
cmb.neg_acc <- read.delim(file.path(GKM, 'chrombpnet_length_neg_test_acc.txt'), sep = '\t', header = F)
colnames(cmb.neg_acc) <- c("peak_id","size","covered","sum","mean0","mean")
cmb.neg_acc$pred <- neg.logcounts
cmb.neg_acc$class <- 'Negative'

cmb.total_counts <- bind_rows(cmb.pos_acc, cmb.neg_acc)
cmb.total_counts$log_sum <- log(cmb.total_counts$sum)
cmb.total_counts <- cmb.total_counts %>% mutate(DA_call = case_when(
                    peak_id%in%cmb.DA$V4 ~ TRUE,
                    TRUE ~ FALSE))
cmb.cor <- cor(cmb.total_counts$pred, cmb.total_counts$log_sum, method = 'spearman')
cmb_DA.cor <- cor(cmb.total_counts %>% filter(DA_call) %>% pull(pred), 
                  cmb.total_counts %>% filter(DA_call) %>% pull(log_sum), method = 'spearman')

cmb_counts.scatter <- ggplot() +
  stat_bin2d(data = cmb.total_counts, aes(x = pred, y = log_sum, fill = after_stat(count)), bins = 30) +
  annotate("text",
           x = Inf, y = -Inf,               # bottom-right corner
           label = paste0("rho = ", round(cmb.cor, 3)),
           hjust = 1.1, vjust = -0.5,      # fine-tune position inside margins
           color = "black") +
  shared_fill_holdout_scale + 
  theme_classic() +
  labs(x = "log (predicted counts)", y = "log (observed counts)") +
  theme(legend.position = 'none')

gkm_counts.scatter <- ggplot() +
  stat_bin2d(data = gkm.total_counts, aes(x = pred, y = log_sum, fill = after_stat(count)), bins = 30) +
  annotate("text",
           x = Inf, y = -Inf,               # bottom-right corner
           label = paste0("rho = ", round(gkm.cor, 3)),
           hjust = 1.1, vjust = -0.5,      # fine-tune position inside margins
           color = "black") +
  shared_fill_holdout_scale + 
  theme_classic() +
  labs(x = "gkm-SVM predictions", y = "log (observed counts)") +
  theme(legend.position = 'none')

cmb_DA_counts.scatter <- ggplot() +
  stat_bin2d(data = cmb.total_counts %>% filter(DA_call), aes(x = pred, y = log_sum, fill = after_stat(count)), bins = 30) +
  annotate("text",
           x = Inf, y = -Inf,               # bottom-right corner
           label = paste0("rho = ", round(cmb_DA.cor, 3)),
           hjust = 1.1, vjust = -0.5,      # fine-tune position inside margins
           color = "black") +
  shared_fill_holdout_scale + 
  theme_classic() +
  labs(x = "log (predicted counts)", y = "log (observed counts)") +
  theme(legend.position = 'none')

gkm_DA_counts.scatter <- ggplot() +
  stat_bin2d(data = gkm.total_counts %>% filter(DA_call), aes(x = pred, y = log_sum, fill = after_stat(count)), bins = 30) +
  annotate("text",
           x = Inf, y = -Inf,               # bottom-right corner
           label = paste0("rho = ", round(gkm_DA.cor, 3)),
           hjust = 1.1, vjust = -0.5,      # fine-tune position inside margins
           color = "black") +
  shared_fill_holdout_scale + 
  theme_classic() +
  labs(x = "gkm-SVM predictions", y = "log (observed counts)") +
  theme(legend.position = 'none')


layout <- '
        AAAABBBB
        AAAABBBB
        CCCCDDDD
        CCCCDDDD
        '
p <- (cmb_counts.scatter + ggtitle("ChromBPNet Holdout prediction") + theme(plot.title = element_text(size = 10, face = "bold", hjust= 0.5))) + (gkm_counts.scatter + ggtitle("gkm-SVM Holdout prediction") + theme(plot.title = element_text(size = 10, face = "bold", hjust= 0.5))) + cmb_DA_counts.scatter + gkm_DA_counts.scatter + plot_layout(design = layout) + plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(size = 12, face='bold'))
ggsave(file.path(OUT_DIR, 'ext_data_fig3_partIII.pdf'), plot  = p, width = 170, height = 100, useDingbats = F, units = 'mm', device = 'pdf')



# cmb_counts.scatter <- ggplot() +
#   geom_point_rast(data = cmb.total_counts, aes(x = pred, y = log_sum, color = class), alpha = 0.5) +
#   scale_color_manual(values = c("Positive" = 'red', "Negative" = "blue")) +
#   theme_classic() +
#   labs(x = "log (predicted counts)", y = "log (observed counts)") +
#   theme(legend.position = 'bottom')

# gkm_counts.scatter <- ggplot() +
#   geom_point_rast(data = gkm.total_counts, aes(x = pred, y = log_sum, color = class), alpha = 0.5) +
#   scale_color_manual(values = c("Positive" = 'red', "Negative" = "blue")) +
#   theme_classic() +
#   labs(x = "gkm-SVM predictions", y = "log (observed counts)") +
#   theme(legend.position = 'bottom')

# p <- (cmb_counts.scatter + ggtitle("ChromBPNet Holdout prediction") + theme(plot.title = element_text(size = 10, face = "bold", hjust= 0.5))) | (gkm_counts.scatter + ggtitle("gkm-SVM Holdout prediction") + theme(plot.title = element_text(size = 10, face = "bold", hjust= 0.5)))
# ggsave('ext_data_fig3_pII.pdf', plot  = p, width = 210, height = 80, useDingbats = F, units = 'mm', device = 'pdf')


##### extended figure 4

INSILICO_DIR <- file.path(FIG_DATA, "insilico_tf_profiles")
insilico_path <- function(cell) file.path(INSILICO_DIR, paste0('insilico_inj_', cell, 'TF_profiles_across_cell_models_random_seq.tsv.gz'))

all.profiles <- data.frame()
cell.profiles <- read.delim(gzfile(insilico_path('endo')), sep = '\t')
cell.profiles <- cell.profiles %>% filter(motifA%in%c('GATA','Foxa','Sox','Klf',"JUND"))
all.profiles <- bind_rows(cell.profiles, all.profiles)

cell.profiles <- read.delim(gzfile(insilico_path('pluri')), sep = '\t')
cell.profiles <- cell.profiles %>% filter(motifA%in%c('Pou5f1+Sox2','OCT4+SOX2'))
all.profiles <- bind_rows(cell.profiles, all.profiles)

cell.profiles <- read.delim(gzfile(insilico_path('meso')), sep = '\t')
cell.profiles <- cell.profiles %>% filter(motifA%in%c('TWIST1'))
all.profiles <- bind_rows(cell.profiles, all.profiles)

cell.profiles <- read.delim(gzfile(insilico_path('ecto')), sep = '\t')
cell.profiles <- cell.profiles %>% filter(motifA%in%c('OTX2'))
all.profiles <- bind_rows(cell.profiles, all.profiles)

all.profiles$cell_type <- factor(all.profiles$cell_type, levels = c('endo','pluri','ecto','meso'))
filter.profiles <- all.profiles %>% filter(motifA%in%c('GATA','Foxa','Sox','Klf',"JUND"))
filter.profiles$motifA <- factor(filter.profiles$motifA, levels = c('GATA','Foxa','Sox','Klf',"JUND"))

cell.plot <- ggplot(filter.profiles, mapping = aes(x = position, y = signal, color = type))+
        geom_line()+
        scale_x_continuous(name = 'Distance from injection (bp)', breaks = seq(-250,250,250),limits = c(-250, 250))+ 
        scale_y_log10(name = 'Average predicted signal', breaks = c(10^-3),labels = scales::trans_format("log10", scales::math_format(10^.x))) + annotation_logticks(sides = 'l')+
#          scale_y_continuous(name = 'Average predicted signal', breaks = seq(0,0.003,0.001),limits = c(0.0008,0.003) )+ 
        scale_color_manual(values = c('red3', '#000000'), labels = c('with motif', 'no motif'), name = "")+
        facet_grid(motifA~cell_type) + 
        theme_classic()+
        theme(legend.position = 'bottom',
             axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, color = 'black'))
#### plot instances
library(XML)
library(RCurl)
library(rlist)
library(universalmotif)
library(htmltools)

tables <- readHTMLTable(file.path(FIG_DATA, 'endo_modisco_motifs.html')) %>% as.data.frame()
colnames(tables) <- sub("^NULL\\.", "", colnames(tables))

motifs <- read_meme(file.path(DATA_DIR, 'motif_databases', 'CIS-BP_Mus_musculus.meme'))
motif_table <- data.frame(
  name = sapply(motifs, function(m) m@name),
  altname = sapply(motifs, function(m) m@altname),
  stringsAsFactors = FALSE
)
colnames(motif_table) <- c("match0","alt0")
tables <- left_join(tables, motif_table, by = c("match0"))
colnames(motif_table) <- c("match1","alt1")
tables <- left_join(tables, motif_table, by = c("match1"))
colnames(motif_table) <- c("match2","alt2")
tables <- left_join(tables, motif_table, by = c("match2"))
                
table2 <- tables %>% dplyr::select(pattern, num_seqlets, modisco_cwm_fwd, modisco_cwm_rev,
                           match0 = alt0, qval0, match0_logo, match1 = alt1,
                           qval1, match1_logo, match2 = alt2, qval2, match2_logo)
table2.filter <- table2 %>% filter(match0%in%c("(Gata2)_(Homo_sapiens)_(DBD_1.00)",
                                              "Klf4","Sox17","(Foxi1)_(Homo_sapiens)_(DBD_0.99)","(Nfe2l1)_(Homo_sapiens)_(DBD_1.00)")) %>% head(n = 5) ### only positive patterns
table2.filter$TF_name <- c("Gata4/6","Klf4","Sox17","Foxa2","AP-1")
table2.filter$num_seqlets <- as.numeric(table2.filter$num_seqlets)
                
seq.counts <- ggplot(table2.filter, aes(x = reorder(TF_name, -num_seqlets), y = 1, fill = num_seqlets)) +
  geom_tile(width = 0.95, height = 0.9, color = "white") +
  geom_text(aes(label = comma(num_seqlets)), fontface = "bold", color = "white") +
  scale_fill_viridis_c(name = "Total number of instances", 
                       breaks = c(10^2, 10^3,10^4),
                       limits = c(100, 11000),
                       labels = scales::trans_format("log10", scales::math_format(10^.x)),
                      trans = "log10") +
  labs(x = NULL, y = NULL) +
  theme_classic() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "bottom"
  )

layout <- 'AAAA######
           AAAA######
           AAAA######
           AAAA######
           AAAA######
           AAAA######
           BBBBBBBBBB'
p <- cell.plot + seq.counts + 
plot_layout(design = layout) + plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(size = 12, face='bold'))
ggsave(file.path(OUT_DIR, 'ext_data_fig4.pdf'), plot  = p, width = 210, height = 180, useDingbats = F, units = 'mm', device = 'pdf')