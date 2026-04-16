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
library(ggplotify)
library(gtools)
library(ggstar)
library(purrr)
library(IRanges)
source("utils/mpra_tablemaker.R")
source("utils/make_tree_functions.R")
source("utils/trajectory_utils.R")
source("utils/figure1_func.R")
source("utils/analyze_sup_figure6_data.R")

# https://joeystanley.com/blog/custom-themes-in-ggplot2
pub_theme <- function() {
  theme_cowplot(font_size = 6) %+replace%
    theme(
      strip.background = element_rect(color='black', 
                                      fill='grey90',
                                      size=0.5, 
                                      linetype="solid",
                                      linewidth = 0.25),
      text = element_text(size=6, family='Helvetica', color = "black"),  # <--- add color
      strip.text = element_text(family = "Helvetica", size=6, face='bold', 
                                margin=margin(b=2.5, t=2.5)),
      axis.text.x = element_text(size=6, color = "black"),
      axis.text.y = element_text(size=6, color = "black"),
      legend.text = element_text(size=6, color = "black"),
      legend.title = element_text(size=6, color = "black"),
      axis.line = element_line(linewidth = 0.25),
      axis.title = element_text(size=7, color = "black"),
      legend.position = 'bottom',
      legend.title.align = 0.5,
      plot.tag = element_text(size = 10, face='bold', color = "black"),
      plot.title = element_text(hjust = 0.5, face='plain', size=9,
                                margin=margin(b=10), color = "black")
    )
}
theme_set(pub_theme())

### sup fig 6a
species.pair <- c("Mus_musculus",'fullTreeAnc239')
mouse.nodes <- make.nodepath(tree, species.pair)

epas1_mouse_msa_list$msa_df_TFBS <- epas1_mouse_msa_list$msa_df_TFBS %>% 
                    mutate(species = factor(species, levels = keep.tips))

labels <- collect_labels(epas1_mouse_msa_list$msa_df_TFBS, side = T)
epas1_gg_tfbs <- plot_tfbs_trajectory(epas1_mouse_msa_list$msa_df_TFBS,
                                      epas1_mouse_msa_list$msa_func_TFBS, CRE_oi,
                                      'Mus_musculus',labels, 
                                      tfbs_size = 0.1) + ggtitle("Epas1:chr17_10063") + 
                                      coord_cartesian(xlim = c(0, 334))

epas1_seq_bar <- plot_seqID_trajectory(oCRE_mean_act, df_mouse_tiles_mpra, control_mean_act, CRE_oi, keep.tips) + 
                    coord_cartesian(clip = 'off') + labs(x = '% sequence similarity\nto mouse ortholog', y = '')

#### phyloP cum dist
library(ggbeeswarm)
xmax <- max(epas1_func_tfbs.phyloP$phyloP, na.rm = TRUE)
phyloP.plot <- ggplot() + geom_quasirandom_rast(data = epas1_func_tfbs.phyloP, 
                                                aes(x = phyloP, y = rank_TF, color = TF_name), 
            dodge.width=1, size = 1) +  # Adjust the width for separation 
            geom_vline(xintercept = 0 ,color = "black", linetype = 'dashed') + 
            geom_text(
                data = epas1_func_tfbs.phyloP_avg,
                aes(x = xmax + 0.5, y = rank_TF, label = round(avg_phyloP, 3)),
                inherit.aes = FALSE, hjust = 0, size = 3
              ) + 
                scale_color_manual(values = col_TFs, guide = 'none') + 
            theme_classic() + coord_cartesian(clip = 'off') +
                            labs(x = 'phyloP', y = 'rank TFBS') + 
                theme(legend.position = 'none')


### mpra and fc trajectory
fc.plot <- ggplot() + 
                    geom_hline(yintercept = 0, color = "black", linetype = 'dashed') +
                    geom_errorbar(data = stepwise_summary, aes(x = species2, ymin=avg_log2FC_steps-sd_log2FC_steps, 
                                                           ymax=avg_log2FC_steps+sd_log2FC_steps), 
                                  width=.2, position=position_dodge(0.05)) + 
                    geom_line(data = stepwise_summary, aes(x = species2, y = avg_log2FC_steps, group = group)) + 
                    geom_star(data = stepwise_summary, aes(x = species2, y = avg_log2FC_steps, 
                                                           starshape = marker, size = marker, fill = marker)) +  
                    scale_starshape_manual(values = c('Other' = 15, "Transitions" = 9, "Mouse" = 1, 'Rat' = 1), guide = guide_legend(override.aes = list(size = 2, fill = 'black')), name = "") +
                    scale_size_manual(values = c('Other' = 1, "Transitions" = 3, "Mouse" = 3, 'Rat' = 3), guide = 'none') +
                    scale_fill_manual(values = c('Other' = 'black', "Transitions" = 'red', "Mouse" = 'red', 'Rat' = 'dark gray'), guide = 'none') +
                    labs(y = "log2FC MPRA activity", x = "") + coord_cartesian(clip="off") +
                    theme_classic() + theme(legend.position = 'bottom', legend.box = "horizontal", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

mpra.plot <- ggplot() + 
                    geom_rect(data = epas1_avg_lineage_mpra %>% filter(species == 'Mus_musculus'), aes(ymin = mean_MPRA_act-sd_MPRA_act, ymax = mean_MPRA_act+sd_MPRA_act, color = NULL), xmin = -Inf, xmax = Inf,
                   alpha = .2, fill = "#008000") + 
                    geom_hline(data = epas1_avg_lineage_mpra %>% filter(species == 'Mus_musculus') , aes(yintercept = mean_MPRA_act), color = '#008000', linetype = 'dashed') +
                    geom_rect(data = control_mean_act %>% filter(CRE_id == 'minP'), aes(ymin = mean_MPRA_act-sd_MPRA_act, ymax = mean_MPRA_act+sd_MPRA_act, color = NULL), xmin = -Inf, xmax = Inf,
                   alpha = .2,fill = "#0057e7") +
                    geom_hline(data = control_mean_act %>% filter(CRE_id == 'minP'), aes(yintercept = mean_MPRA_act), color = '#0057e7', linetype = 'dashed') +
                    geom_line(data = epas1_avg_lineage_mpra, aes(x = species2, y = mean_MPRA_act, group = group)) + 
                    geom_errorbar(data = epas1_avg_lineage_mpra, aes(x = species2, ymin=mean_MPRA_act-sd_MPRA_act, 
                                                           ymax=mean_MPRA_act+sd_MPRA_act), 
                                  width=.2, position=position_dodge(0.05)) + 
                    geom_star(data = epas1_avg_lineage_mpra, aes(x = species2, y = mean_MPRA_act, starshape = marker, size = marker, fill = marker)) +  
                    scale_starshape_manual(values = c('Other' = 15, "Transitions" = 9, "Mouse" = 1, 'Rat' = 1), guide = guide_legend(override.aes = list(size = 2, fill = 'black')), name = "") +
                    scale_size_manual(values = c('Other' = 1, "Transitions" = 3, "Mouse" = 3, 'Rat' = 3), guide = 'none') +
                    scale_fill_manual(values = c('Other' = 'black', "Transitions" = 'red', "Mouse" = 'red', 'Rat' = 'dark gray'), guide = 'none') +
                    labs(y = "MPRA activity", x = "") + coord_cartesian(clip="off") + 
                    theme_classic() + theme(legend.position = 'none', legend.box = "horizontal", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



library(ggtext)
### fig. 2d
mouse_epas1_gg_tfbs <- plot_tfbs_trajectory(epas1_mouse_msa_list$msa_df_TFBS %>% filter(species == 'Mus_musculus'),
                                      epas1_mouse_msa_list$msa_func_TFBS, 
                                      CRE_oi, 'Mus_musculus',
                                      c('Mouse'), 
                                      tfbs_size = 0.1) + labs(x = "", y ="") + 
                 geom_rect(data = tfbs_map %>% head(n=1), aes(xmin = msa_start, xmax = msa_end+3, ymin = -Inf, ymax = Inf),
                          inherit.aes = FALSE, fill = NA, color = 'black', linetype = 'dashed', alpha = 0.7) + coord_cartesian(xlim = c(0, 334)) 



epas1_full_msa_plot <- ggplot(fixed_epas1_mouse_msa_df, aes(x = curpos, y = species2, fill = log_change)) +
  geom_tile_rast() +
  geom_rect(data = tfbs_map %>% head(1), aes(xmin = msa_start, xmax = msa_end+3, ymin = -Inf, ymax = Inf),
                          inherit.aes = FALSE, fill = NA, color = 'black', linetype = 'dashed', alpha = 0.7) +
  # Black borders for "First Match" only
  geom_tile_rast(
    data = subset(fixed_epas1_mouse_msa_df, match2 == "First Match"),
    fill = NA,
    color = "black",
    size = 0.1, alpha = 0.1
  ) +
  geom_segment(
    aes(x = curpos - 0.5, xend = curpos + 0.5, 
        y = as.numeric(species2) - 0.5, 
        yend = as.numeric(species2) - 0.5),
    inherit.aes = FALSE,
    color = "black",
    linewidth = 0.1
  ) +
  scale_fill_gradient2(low = "dodgerblue2",
                             mid = "white",
                             high = "firebrick1",
                             midpoint = 0,
                             space = "Lab",
                             na.value = "gainsboro",
                             limits=c(-1.5,1.5),
                             breaks = c(-1.5, 0, 1.5),
                             oob=squish, name = expression(log[2]~"FC (MPRA)")) +
  theme_minimal() + coord_cartesian(xlim = c(0, 334)) +
  labs(x = "msa bp position", y = "") +
  theme(
    panel.grid = element_blank(),
    axis.ticks.y = element_blank(), # Optionally remove ticks if necessary
      legend.position = 'bottom'
  )

epas1_snp_plot <- ggplot(epas1_snp_counts, aes(y=species2)) + #, color=mean_MPRA_act)) + 
            geom_point_rast(aes(x = total_nontfbs_snps), size = 1, color = 'black') +
            geom_point_rast(aes(x = total_tfbs_snps), size = 1, color = 'red') +
            labs(x = '# of mouse\nmutations gained') + theme_classic() +
            geom_vline(xintercept = 0, color = "black", linetype = 'dashed') +
            scale_x_continuous(breaks = c(0, 15 , 30)) + coord_cartesian(clip = 'off') +
            theme(legend.position = 'none',
                  axis.text.y = element_blank(),
                  axis.line.y = element_blank(),
                  axis.ticks.y = element_blank(),
                  axis.title.y=element_blank())

epas1_fc_plot <- ggplot(epas1_group_fc, aes(y=species2, x = avg_log_change, color = tf_call)) +
  geom_vline(xintercept = 0, color = "black", linetype = 'dashed') +
  
  # Plot black points first
  geom_point_rast(size = 1, alpha = 0.6) +
  scale_color_manual(values = c('non func. TFBS' = 'black', 'func. TFBS' = 'red')) +
  labs(x = expression(log[2] ~ "FC (MPRA)")) +
  theme_classic() +  coord_cartesian(clip = 'off') +
  theme(
    legend.position = 'none',
    axis.text.y = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank()
  )


## zoom in mutations in Anc238
epas1_triplet_map <- epas1_mouse_msa_list$msa_func_TFBS %>% arrange(msa_start) %>% dplyr::select(TF_name, rank_TF, msa_start, msa_end) %>% 
                            filter(msa_start >= 239 & msa_end <= 279)

triplet_plot <- ggplot(fixed_epas1_mouse_msa_df %>% filter(curpos >= 239) %>% filter(curpos <= 279), 
                       aes(x = curpos, y = species2, label = mut2, fill = log_change)) +
  geom_tile_rast() +
  geom_tile_rast(
    data = subset(fixed_epas1_mouse_msa_df, curpos >= 239 & curpos <= 279 & match2 == "First Match"),
    fill = NA,
    color = "black",
    linewidth = 0.15,   # use linewidth (not size) for borders
    width = 1, height = 1,   # match one grid cell exactly
  ) +
  scale_fill_gradient2(low = "dodgerblue2",
                             mid = "white",
                             high = "firebrick1",
                             midpoint = 0,
                             space = "Lab",
                             na.value = "gainsboro",
                             limits=c(-1.5,1.5),
                             breaks = c(-1.5, 0, 1.5),
                             oob=squish) +
  geom_text(aes(fontface = ifelse(match2 == "First Match", "bold",
                                 ifelse(match2 == "DEL", "bold", "plain")), 
                alpha = ifelse(match2 == "First Match", 1,
                             ifelse(match2 != "Match", 1, 0.9)),
               size = I(ifelse(mut2 == "*", 5, 3))  # treat as fixed, no scaling
               ),color = 'black',
           show.legend = FALSE) +
  geom_rect(data = epas1_triplet_map, aes(xmin = msa_start-0.5, xmax = msa_end+0.5, ymin = -Inf, ymax = Inf, color = TF_name),
                          inherit.aes = FALSE, fill = NA, linetype = 'dashed', alpha = 0.7) +
  scale_color_manual(values = col_TFs) +
  scale_x_continuous(breaks = seq(240, 279, 4)) + 
  theme_minimal() + 
  labs(x = "msa bp position", y = "") +
  theme(
    panel.grid = element_blank(),
    axis.ticks.y = element_blank(), # Optionally remove ticks if necessary
      legend.position = 'none'
  )

fixed_epas1_anc_msa_df <- fixed_epas1_mouse_msa_df %>% 
    filter(species %in%c('fullTreeAnc238','fullTreeAnc239','fullTreeAnc115','fullTreeAnc69','fullTreeAnc67','fullTreeAnc59'))
epas1_ap1_map <- epas1_mouse_msa_list$msa_func_TFBS %>% arrange(msa_start) %>% dplyr::select(TF_name, rank_TF, msa_start, msa_end) %>% 
                            filter(msa_start >= 173 & msa_end <= 200)

ap1_triplet_plot <- ggplot(fixed_epas1_anc_msa_df %>% filter(curpos >= 173) %>% filter(curpos <= 200), 
                       aes(x = curpos, y = species2, label = mut2, fill = log_change)) +
  geom_tile_rast() +
  geom_tile_rast(
    data = subset(fixed_epas1_anc_msa_df, curpos >= 173 & curpos <= 200 & match2 == "First Match"),
    fill = NA,
    color = "black",
    linewidth = 0.15,   # use linewidth (not size) for borders
    width = 1, height = 1,   # match one grid cell exactly
  ) +
  scale_fill_gradient2(low = "dodgerblue2",
                             mid = "white",
                             high = "firebrick1",
                             midpoint = 0,
                             space = "Lab",
                             na.value = "gainsboro",
                             limits=c(-1.5,1.5),
                             breaks = c(-1.5, 0, 1.5),
                             oob=squish) +
  geom_text(aes(fontface = ifelse(match2 == "First Match", "bold",
                                 ifelse(match2 == "DEL", "bold", "plain")), 
                alpha = ifelse(match2 == "First Match", 1,
                             ifelse(match2 != "Match", 1, 0.9)),
               size = I(ifelse(mut2 == "*", 5, 3))  # treat as fixed, no scaling
               ),color = 'black',
           show.legend = FALSE) +
  geom_rect(data = epas1_ap1_map, aes(xmin = msa_start-0.5, xmax = msa_end+0.5, ymin = -Inf, ymax = Inf, color = TF_name),
                          inherit.aes = FALSE, fill = NA, linetype = 'dashed', alpha = 0.7) +
  scale_color_manual(values = col_TFs) +
  scale_x_continuous(breaks = seq(174, 200, 4)) + 
  theme_minimal() + 
  labs(x = "msa bp position", y = "") +
  theme(
    panel.grid = element_blank(),
    axis.ticks.y = element_blank(), # Optionally remove ticks if necessary
      legend.position = 'none'
  )

### supp aff plot
triplet_aff.plot <- ggplot(new_aff.mat %>% filter(msa_start >= 239) %>% filter(msa_end <= 279), 
                                    aes(x = log2_fc_steps, y = species2, fill = TF_name)) +
  geom_bar(stat="identity") + 
  scale_fill_manual(values = col_TFs) +
  facet_wrap(~TF_name, ncol = 3) +
  labs(x = expression(log[2]~"FC predicted TFBS affinity"), y = "") +
  geom_vline(xintercept = 0, color = "black", linetype = 'dashed') +
#   scale_x_continuous(breaks = c(0, 0.5, 1)) + 
  theme_classic() + 
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
#     axis.ticks.y = element_blank(), # Optionally remove ticks if necessary
      legend.position = 'none'
  )


anc_aff.mat <- new_aff.mat %>% filter(species%in%c('fullTreeAnc238','fullTreeAnc239','fullTreeAnc115','fullTreeAnc69','fullTreeAnc68','fullTreeAnc67','fullTreeAnc59'))
ap1_anc_aff.mat <- anc_aff.mat %>% filter(msa_start >= 173) %>% filter(msa_end <= 200) %>% 
                mutate(TF_name2 = factor(TF_name2, levels = c('Jun_Atf3_2','Jun_Atf3_5')))

ap1_aff.plot <- ggplot(anc_aff.mat %>% filter(msa_start >= 173) %>% filter(msa_end <= 200), 
                                    aes(x = log2_fc_steps, y = species2, fill = TF_name)) +
  geom_bar(stat="identity") + 
  scale_fill_manual(values = col_TFs) +
  facet_wrap(~TF_name2, ncol = 2) +
  labs(x = expression(log[2]~"FC predicted TFBS affinity"), y = "") +
  geom_vline(xintercept = 0, color = "black", linetype = 'dashed') +
#   scale_x_continuous(breaks = c(0, 0.5, 1)) + 
  theme_classic() + 
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
      legend.position = 'none'
  )

layout <- 'AAAAAABCCCC
           AAAAAABCCCC
           AAAAAABCCCC
           AAAAAABCCCC
           AAAAAABCCCC
           DDDDDD#####
           EEEEEE#FFGG
           EEEEEE#FFGG
           EEEEEE#FFGG
           EEEEEE#FFGG
           HHHHHHHIIII
           HHHHHHHIIII
           HHHHHHHIIII
           HHHHHHHIIII'
p1 <-  epas1_gg_tfbs + epas1_seq_bar + phyloP.plot  + mouse_epas1_gg_tfbs + epas1_full_msa_plot + 
    epas1_snp_plot + epas1_fc_plot + 
    triplet_plot + triplet_aff.plot + plot_layout(design = layout) + plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(size = 12, face='bold'))
ggsave('sup_fig6.pdf', plot = p1, width = 210, height = 260, useDingbats = F, units = 'mm', device = 'pdf') 

