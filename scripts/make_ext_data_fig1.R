### module load pcre2/10.39; module load R/4.3.1
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
#### chrombpnet predictions ####
library(rhdf5)
source("utils/mpra_tablemaker.R")
source("utils/make_tree_functions.R")
source("utils/trajectory_utils.R")
source("utils/figure1_func.R")
source("utils/load_figure1_2_data.R")
load("../extdata/oCRE_func_tfbs_info.RData")

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
load("../extdata/oCRE_func_tfbs_info.RData")
order.species <- oCRE_sum_aff_tfbs.count %>% dplyr::select(species, common_name, Order) %>% distinct()
oCRE_gkm_pred <- read.delim("../data/oCRE_gkmSVM_endo_prediction.txt",
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
ggsave('ext_data_fig1.pdf', plot = p, width = 210, height = 200, useDingbats = F, units = 'mm', device = 'pdf')