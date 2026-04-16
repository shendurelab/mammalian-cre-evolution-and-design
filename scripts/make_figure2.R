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
source("utils/trajectory_utils.R")
source("utils/figure1_func.R")
source("utils/make_tree_functions.R")
source("utils/analyze_figure2_data.R")

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

### plot activity changes
### fig 2b-c
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
#                     scale_y_log10(
#                        breaks = c(10^-1, 10^0, 10^1),
#                        labels = scales::trans_format("log10", scales::math_format(10^.x))
#                      ) + annotation_logticks(sides = 'l') +
                    labs(y = "log2FC MPRA activity", x = "") + coord_cartesian(clip="off") +
                    theme_classic() + theme(legend.position = 'bottom', legend.box = "horizontal", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

mpra.plot <- ggplot() + 
                    geom_rect(data = gata4_avg_lineage_mpra %>% filter(species == 'Mus_musculus'), aes(ymin = mean_MPRA_act-sd_MPRA_act, ymax = mean_MPRA_act+sd_MPRA_act, color = NULL), xmin = -Inf, xmax = Inf,
                   alpha = .2, fill = "#008000") + 
                    geom_hline(data = gata4_avg_lineage_mpra %>% filter(species == 'Mus_musculus') , aes(yintercept = mean_MPRA_act), color = '#008000', linetype = 'dashed') +
                    geom_rect(data = control_mean_act %>% filter(CRE_id == 'minP'), aes(ymin = mean_MPRA_act-sd_MPRA_act, ymax = mean_MPRA_act+sd_MPRA_act, color = NULL), xmin = -Inf, xmax = Inf,
                   alpha = .2,fill = "#0057e7") +
                    geom_hline(data = control_mean_act %>% filter(CRE_id == 'minP'), aes(yintercept = mean_MPRA_act), color = '#0057e7', linetype = 'dashed') +
                    geom_line(data = gata4_avg_lineage_mpra, aes(x = species2, y = mean_MPRA_act, group = group)) + 
                    geom_errorbar(data = gata4_avg_lineage_mpra, aes(x = species2, ymin=mean_MPRA_act-sd_MPRA_act, 
                                                           ymax=mean_MPRA_act+sd_MPRA_act), 
                                  width=.2, position=position_dodge(0.05)) + 
                    geom_star(data = gata4_avg_lineage_mpra, aes(x = species2, y = mean_MPRA_act, starshape = marker, size = marker, fill = marker)) +  
                    scale_starshape_manual(values = c('Other' = 15, "Transitions" = 9, "Mouse" = 1, 'Rat' = 1), guide = guide_legend(override.aes = list(size = 2, fill = 'black')), name = "") +
                    scale_size_manual(values = c('Other' = 1, "Transitions" = 3, "Mouse" = 3, 'Rat' = 3), guide = 'none') +
                    scale_fill_manual(values = c('Other' = 'black', "Transitions" = 'red', "Mouse" = 'red', 'Rat' = 'dark gray'), guide = 'none') +
#                     scale_y_log10(
#                        breaks = c(10^-1, 10^0, 10^1),
#                        labels = scales::trans_format("log10", scales::math_format(10^.x))
#                      ) + annotation_logticks(sides = 'l') +
                    labs(y = "MPRA activity", x = "") + coord_cartesian(clip="off") + 
                    theme_classic() + theme(legend.position = 'none', legend.box = "horizontal", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


### supplementary fig 5b
### plot gata4 tree as example
### plot subtree
### supplementary fig 5d
gata.subtree <- ggtree(sub.tree, aes(color = mean_MPRA_act), size = 0.8) + #geom_treescale(x=0, y=235) + 
          layout_dendrogram() +
          geom_text2(aes(label = label2), hjust = -0.15, vjust = -0.5, color = "black", size = 3) +   
          geom_treescale(x = 0, y = 1.5) +
          scale_color_gradient(breaks = c(0.1, 1, 10), limits = c(0.05, 13),  # Quantile limits from your range
                               labels = scales::trans_format("log10", scales::math_format(10^.x)),
                                oob = scales::squish,  # To handle values outside the specified limits
                                trans = "log10",
                               low = '#42047e',high = '#07f49e', na.value = 'gainsboro') +
        theme(legend.position = 'bottom')
gata.subtree <- gata.subtree %<+% sub.tips + geom_star(
                mapping=aes(fill=mean_MPRA_act, starshape = Shape), size = 2,
                position="identity",starstroke=1, color = 'NA', show_guide = FALSE) +
    scale_fill_gradient(breaks = c(0.1, 1, 10), limits = c(0.05, 13),  # Quantile limits from your range
                          labels = scales::trans_format("log10", scales::math_format(10^.x)),
                            oob = scales::squish,  # To handle values outside the specified limits
                            trans = "log10",
                           low = '#42047e',high = '#07f49e', na.value = 'gainsboro') +
    labs(color = 'MPRA activity', fill = 'MPRA activity') + 
    theme(legend.position = 'bottom',
            legend.spacing = unit(0, "cm"),        # Remove extra spacing between legend items
            legend.box.margin = margin(-10, 0, 0, 0), # Reduce margin below the legend
            plot.title.position = "plot",         # Align title closer to the plot
            plot.margin = margin(5, 5, 5, 5),     # Reduce overall plot margin
          panel.spacing = unit(0.1, "lines"))

#### supplementary fig 5b-c rat
rat_fc.plot <- ggplot() + 
                    geom_hline(yintercept = 0, color = "black", linetype = 'dashed') +
                    geom_errorbar(data = rat_stepwise_summary, aes(x = species2, ymin=avg_log2FC_steps-sd_log2FC_steps, 
                                                           ymax=avg_log2FC_steps+sd_log2FC_steps), 
                                  width=.2, position=position_dodge(0.05)) + 
                    geom_line(data = rat_stepwise_summary, aes(x = species2, y = avg_log2FC_steps, group = group)) + 
                    geom_star(data = rat_stepwise_summary, aes(x = species2, y = avg_log2FC_steps, 
                                                           starshape = marker, size = marker, fill = marker)) +  
                    scale_starshape_manual(values = c('Other' = 15, "Transitions" = 9, "Mouse" = 1, 'Rat' = 1), guide = guide_legend(override.aes = list(size = 2, fill = 'black')), name = "") +
                    scale_size_manual(values = c('Other' = 1, "Transitions" = 3, "Mouse" = 3, 'Rat' = 3), guide = 'none') +
                    scale_fill_manual(values = c('Other' = 'black', "Transitions" = 'red', "Mouse" = 'red', 'Rat' = 'dark gray'), guide = 'none') +
#                     scale_y_log10(
#                        breaks = c(10^-1, 10^0, 10^1),
#                        labels = scales::trans_format("log10", scales::math_format(10^.x))
#                      ) + annotation_logticks(sides = 'l') +
                    labs(y = "log2FC MPRA activity", x = "") + coord_cartesian(clip="off") +
                    theme_classic() + theme(legend.position = 'bottom', legend.box = "horizontal", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

rat_mpra.plot <- ggplot() + 
                    geom_rect(data = df_mouse_tiles_mpra %>% filter(full_CRE_id == CRE_oi), aes(ymin = mean_MPRA_act-sd_MPRA_act, ymax = mean_MPRA_act+sd_MPRA_act, color = NULL), xmin = -Inf, xmax = Inf,
                   alpha = .2, fill = "#008000") + 
                    geom_hline(data = df_mouse_tiles_mpra %>% filter(full_CRE_id == CRE_oi) , aes(yintercept = mean_MPRA_act), color = '#008000', linetype = 'dashed') +
                    geom_rect(data = control_mean_act %>% filter(CRE_id == 'minP'), aes(ymin = mean_MPRA_act-sd_MPRA_act, ymax = mean_MPRA_act+sd_MPRA_act, color = NULL), xmin = -Inf, xmax = Inf,
                   alpha = .2,fill = "#0057e7") +
                    geom_hline(data = control_mean_act %>% filter(CRE_id == 'minP'), aes(yintercept = mean_MPRA_act), color = '#0057e7', linetype = 'dashed') +
                    geom_line(data = gata4_avg_rat_lineage_mpra, aes(x = species2, y = mean_MPRA_act, group = group)) + 
                    geom_errorbar(data = gata4_avg_rat_lineage_mpra, aes(x = species2, ymin=mean_MPRA_act-sd_MPRA_act, 
                                                           ymax=mean_MPRA_act+sd_MPRA_act), 
                                  width=.2, position=position_dodge(0.05)) + 
                    geom_star(data = gata4_avg_rat_lineage_mpra, aes(x = species2, y = mean_MPRA_act, starshape = marker, size = marker, fill = marker)) +  
                    scale_starshape_manual(values = c('Other' = 15, "Transitions" = 9, "Mouse" = 1, 'Rat' = 1), guide = guide_legend(override.aes = list(size = 2, fill = 'black')), name = "") +
                    scale_size_manual(values = c('Other' = 1, "Transitions" = 3, "Mouse" = 3, 'Rat' = 3), guide = 'none') +
                    scale_fill_manual(values = c('Other' = 'black', "Transitions" = 'red', "Mouse" = 'red', 'Rat' = 'dark gray'), guide = 'none') +
#                     scale_y_log10(
#                        breaks = c(10^-1, 10^0, 10^1),
#                        labels = scales::trans_format("log10", scales::math_format(10^.x))
#                      ) + annotation_logticks(sides = 'l') +
                    labs(y = "MPRA activity", x = "") + coord_cartesian(clip="off") + 
                    theme_classic() + theme(legend.position = 'none', legend.box = "horizontal", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

library(ggtext)
### fig. 2d
gata4_full_msa_plot <- ggplot(final_gata4_msa_df, aes(x = curpos, y = species2, fill = log_change)) +
  geom_tile_rast() +
  geom_rect(data = tfbs_map, aes(xmin = msa_start, xmax = msa_end+3, ymin = -Inf, ymax = Inf),
                          inherit.aes = FALSE, fill = NA, color = 'black', linetype = 'dashed', alpha = 0.7) +
  # Black borders for "First Match" only
  geom_tile_rast(
    data = subset(final_gata4_msa_df, match2 == "First Match"),
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
  theme_minimal() + coord_cartesian(xlim = c(0, 310)) +
  labs(x = "msa bp position", y = "") +
  theme(
    panel.grid = element_blank(),
    axis.ticks.y = element_blank(), # Optionally remove ticks if necessary
      legend.position = 'bottom'
  )

gata4_snp_plot <- ggplot(final_gata4_snp_counts, aes(y=species2)) + #, color=mean_MPRA_act)) + 
            geom_point_rast(aes(x = total_nontfbs_snps), size = 1, color = 'black') +
            geom_point_rast(aes(x = total_tfbs_snps), size = 1, color = 'red') +
            labs(x = '# of mouse\nmutations\ngained') + theme_classic() +
            geom_vline(xintercept = 0, color = "black", linetype = 'dashed') +
            scale_x_continuous(breaks = c(0, 6 , 12)) + coord_cartesian(clip = 'off') +
            theme(legend.position = 'none',
                  axis.text.y = element_blank(),
                  axis.line.y = element_blank(),
                  axis.ticks.y = element_blank(),
                  axis.title.y=element_blank())

gata4_fc_plot <- ggplot(final_gata4_group_fc, aes(y=species2, x = avg_log_change, color = tf_call)) +
  geom_vline(xintercept = 0, color = "black", linetype = 'dashed') +
  
  # Plot black points first
  geom_point_rast(size = 1, alpha = 0.6) +
  scale_color_manual(values = c('non func. TFBS' = 'black', 'func. TFBS' = 'red')) +
  labs(x = expression(atop(log[2] ~ "FC", "(MPRA)"))) +
  theme_classic() +  coord_cartesian(clip = 'off') +
  theme(
    legend.position = 'none',
    axis.text.y = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank()
  )


## fig. 2e
gata4_klf_plot <- ggplot(final_gata4_msa_df %>% filter(curpos >= 26) %>% filter(curpos <= 37), aes(x = curpos, y = species2, label = mut2, fill = log_change)) +
  geom_tile_rast() +
  geom_tile_rast(
    data = subset(final_gata4_msa_df, curpos >= 26 & curpos <= 37 & match2 == "First Match"),
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
  scale_x_continuous(breaks = c(26, 30, 34)) + 
  theme_minimal() + 
  labs(x = "msa bp position", y = "") +
  theme(
    panel.grid = element_blank(),
    axis.ticks.y = element_blank(), # Optionally remove ticks if necessary
      legend.position = 'none'
  )

gata4_sox17_pair_plot <- ggplot(final_gata4_msa_df %>% filter(curpos >= 79) %>% filter(curpos <= 92), aes(x = curpos, y = species2, label = mut2, fill = log_change)) +
  geom_tile_rast() +
  geom_tile_rast(
    data = subset(final_gata4_msa_df, curpos >= 79 & curpos <= 92 & match2 == "First Match"),
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
  scale_x_continuous(breaks = c(80, 84, 88, 92)) + 
  theme_minimal() + 
  labs(x = "msa bp position", y = "") +
  theme(
    panel.grid = element_blank(),
    axis.ticks.y = element_blank(), # Optionally remove ticks if necessary
    axis.text.y = element_blank(), # Optionally remove ticks if necessary
      legend.position = 'none'
  )

new_sox17_plot <- ggplot(final_gata4_msa_df %>% filter(curpos >= 175) %>% filter(curpos <= 186), aes(x = curpos, y = species2, label = mut2, fill = log_change)) +
  geom_tile_rast() +
  geom_tile_rast(
    data = subset(final_gata4_msa_df, curpos >= 175 & curpos <= 186 & match2 == "First Match"),
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
  scale_x_continuous(breaks = c(176, 180, 184)) + 
  theme_minimal() + 
  labs(x = "msa bp position", y = "") +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(), # Optionally remove ticks if necessary
      legend.position = 'none'
  )

new_ap1_plot <- ggplot(final_gata4_msa_df %>% filter(curpos >= 241) %>% filter(curpos <= 252), 
                       aes(x = curpos, y = species2, label = mut2, fill = log_change)) +
  geom_tile_rast() +
  geom_tile_rast(
    data = subset(final_gata4_msa_df, curpos >= 241 & curpos <= 252 & match2 == "First Match"),
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
  scale_x_continuous(breaks = c(242, 246, 250)) + 
  theme_minimal() + 
  labs(x = "msa bp position", y = "") +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(), # Optionally remove ticks if necessary
      legend.position = 'none'
  )

### make mouse TFBS map
### goes with fig. 2d
mouse_gata4_gg_tfbs <- plot_tfbs_trajectory(gata4_msa_list$msa_df_TFBS %>% filter(species == 'Mus_musculus'),
                                      gata4_msa_list$msa_func_TFBS, 
                                      CRE_oi, 'Mus_musculus',
                                      c('Mouse'), 
                                      tfbs_size = 0.1) + labs(x = "", y ="") + 
                 geom_rect(data = tfbs_map, aes(xmin = msa_start, xmax = msa_end+3, ymin = -Inf, ymax = Inf),
                          inherit.aes = FALSE, fill = NA, color = 'black', linetype = 'dashed', alpha = 0.7) + coord_cartesian(xlim = c(0, 310)) 

### fig. 2f
gata4_klf_aff.plot <- ggplot(final_aff.mat %>% filter(msa_start >= 26) %>% filter(msa_end <= 37), 
                                    aes(x = log2_fc_steps, y = species2, fill = TF_name)) +
  geom_bar(stat="identity") + 
  scale_fill_manual(values = col_TFs) +
  facet_wrap(~TF_name, ncol = 2) +
  labs(x = "", y = "") +
  geom_vline(xintercept = 0, color = "black", linetype = 'dashed') +
#   scale_x_continuous(breaks = c(0, 0.5, 1)) + 
  theme_classic() + 
  theme(
    panel.grid = element_blank(),
    axis.ticks.y = element_blank(), # Optionally remove ticks if necessary
      legend.position = 'none'
  )

gata4_sox17_pair_aff.plot <- ggplot(final_aff.mat %>% filter(msa_start >= 79) %>% filter(msa_end <= 92), 
                                    aes(x = log2_fc_steps, y = species2, fill = TF_name)) +
  geom_bar(stat="identity") + 
  scale_fill_manual(values = col_TFs) +
  facet_wrap(~TF_name, ncol = 2) +
  labs(x = expression(log[2]~"FC predicted TFBS affinity"), y = "") +
  geom_vline(xintercept = 0, color = "black", linetype = 'dashed') +
#   scale_x_continuous(breaks = c(0, 0.5, 1)) + 
  theme_classic() + 
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(), # Optionally remove ticks if necessary
      legend.position = 'none'
  )

new_sox17_aff.plot <- ggplot(final_aff.mat %>% filter(msa_start >= 175) %>% filter(msa_end <= 186), 
                             aes(x = log2_fc_steps, y = species2, fill = TF_name)) +
  geom_bar(stat="identity") + 
  scale_fill_manual(values = col_TFs) +
  facet_wrap(~TF_name, ncol = 1) +
  labs(x = '', y = '') + 
#   labs(x = expression("Fractional " * Delta * " predicted TFBS affinity"), y = "") +
  geom_vline(xintercept = 0, color = "black", linetype = 'dashed') +
#   scale_x_continuous(breaks = c(0, 0.5, 1, 1.5)) + 
  theme_classic() + 
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(), # Optionally remove ticks if necessary
      legend.position = 'none'
  )

new_ap1_aff.plot <- ggplot(final_aff.mat %>% filter(msa_start >= 241) %>% filter(msa_end <= 252), 
                             aes(x = log2_fc_steps, y = species2, fill = TF_name)) +
  geom_bar(stat="identity") + 
  scale_fill_manual(values = col_TFs) +
  facet_wrap(~TF_name, ncol = 1) +
  geom_vline(xintercept = 0, color = "black", linetype = 'dashed') +
  scale_x_continuous(breaks = c(0, 0.5, 1)) + 
  theme_classic() + 
  labs(x = "", y = "") +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(), # Optionally remove ticks if necessary
      legend.position = 'none'
  )

### supplementary fig 5a
CRE_oi <- "Gata4_chr14_5729"
fasta_file <- '/net/shendure/vol10/projects/tli/nobackup/EB_dev_project/oCREs_v2/nobackup/reoriented_fasta/Gata4_chr14_5729_oCRE.fasta'

gata4_msa_list <- make_msa_tfbs_traj_df(fasta_file, mouse.nodes, 'Mus_musculus', CRE_oi, cactus.tfbs, mouse_func_tfbs)
labels <- collect_labels(gata4_msa_list$msa_df_TFBS, side = T)
gata4_gg_tfbs <- plot_tfbs_trajectory(gata4_msa_list$msa_df_TFBS, gata4_msa_list$msa_func_TFBS, CRE_oi, 'Mus_musculus',labels, 
                                      tfbs_size = 0.1) + ggtitle("Gata4:chr14_5729")
gata4_seq_bar <- plot_seqID_trajectory(oCRE_mean_act, df_mouse_tiles_mpra, control_mean_act, CRE_oi, mouse.nodes) + 
                    coord_cartesian(clip = 'off') + labs(x = '% sequence similarity\nto mouse ortholog', y = '')

#### make fig2a schematic
labels <- c()
for(species_oi in c('fullTreeAnc239','fullTreeAnc38','Mus_musculus')){
    name_oi <- str_replace(species_oi, "fullTree","")
    if(species_oi%in%select.images$label){
        img_path <- select.images %>% filter(label == species_oi) %>% pull(uid)
        name_oi <- select.images %>% filter(label == species_oi) %>% pull(name)
        labels <- c(labels, paste0(
                      "<img src='", img_path, "' width='25' style='vertical-align: middle;'/>"
                    ))
    } else{
        labels <- c(labels, name_oi)
    }
}

sub_gata4_avg_lineage_mpra <- gata4_avg_lineage_mpra %>% filter(species2%in%c('Anc239','Anc38','Mus_musculus'))
step.example <- ggplot() + 
                    geom_line(data = sub_gata4_avg_lineage_mpra, aes(x = species2, y = mean_MPRA_act, group = group),
                             color = 'cyan',linetype = 'dashed') + 
                    geom_star(data = sub_gata4_avg_lineage_mpra, aes(x = species2, y = mean_MPRA_act, starshape = marker, size = marker, fill = marker)) +  
                    scale_starshape_manual(values = c('Other' = 15, "Transitions" = 9, "Mouse" = 1, 'Rat' = 1), guide = guide_legend(override.aes = list(size = 2, fill = 'black')), name = "") +
                    scale_size_manual(values = c('Other' = 1, "Transitions" = 3, "Mouse" = 3, 'Rat' = 3), guide = 'none') +
                    scale_fill_manual(values = c('Other' = 'black', "Transitions" = 'red', "Mouse" = 'red', 'Rat' = 'dark gray'), guide = 'none') +
#                     scale_y_log10(
#                        breaks = c(10^-1, 10^0, 10^1),
#                        labels = scales::trans_format("log10", scales::math_format(10^.x))
#                      ) + annotation_logticks(sides = 'l') +
                    labs(y = "MPRA activity", x = "") + coord_cartesian(clip="off") + 
                    scale_x_discrete(labels = labels) + theme_classic() +
                            theme(legend.position = 'none',
                                  axis.text.x = element_markdown(color = "black", size = 9),
                                  panel.grid = element_blank(),
                                  axis.ticks.x = element_blank(),
                                  axis.text.y = element_blank(),
                                  axis.ticks.y = element_blank())

example.tfbs <- plot_tfbs_trajectory(gata4_msa_list$msa_df_TFBS %>% 
                                     filter(species%in%c('fullTreeAnc239','fullTreeAnc38','Mus_musculus')),
                                      gata4_msa_list$msa_func_TFBS, 
                                      CRE_oi, 'Mus_musculus',
                                      rev(labels), 
                                      tfbs_size = 0.5) + coord_cartesian(xlim = c(235, 255))
                        
subset_gata4_msa_df <- fixed_gata4_msa_df %>% filter(species%in%c('fullTreeAnc239','fullTreeAnc38','Mus_musculus'))
example_mut_plot <- ggplot(subset_gata4_msa_df %>% filter(curpos >= 241) %>% filter(curpos <= 252), 
                       aes(x = curpos, y = species2, label = mut2, fill = log_change)) +
  geom_tile() +
  geom_tile(
    data = subset(subset_gata4_msa_df, curpos >= 241 & curpos <= 252 & match2 == "First Match"),
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
                             ifelse(match2 == "DEL", 1, 0.9)),
               size = I(ifelse(mut2 == "*", 5, 3))  # treat as fixed, no scaling
               ),color = 'black',
           show.legend = FALSE) +
  scale_x_continuous(breaks = c(242, 246, 250)) + 
  theme_minimal() + 
  labs(x = "", y = "") +
  scale_y_discrete(labels = rev(labels)) + 
                            theme(legend.position = 'none',
                                  axis.text.y = element_markdown(color = "black", size = 9),
                                  panel.grid = element_blank(),
                                  axis.ticks.x = element_blank(),
                                  axis.text.x = element_blank(),
                                  axis.ticks.y = element_blank())

layout <- '###BBB
           AAACCC'
p <- step.example + example.tfbs + example_mut_plot + plot_layout(design = layout)
ggsave('fig2_schematic.pdf', plot = p, width = 120, height = 80, useDingbats = F, units = 'mm', device = 'pdf')

layout <- '#########AAAA
           #########AAAA
           #########AAAA
           #########BBBB
           #########BBBB
           #########BBBB
           CCCCCCCCCC###
           DDDDDDDDDDE#F
           DDDDDDDDDDE#F
           DDDDDDDDDDE#F
           DDDDDDDDDDE#F
           GGGHHHHIIIJJJ
           GGGHHHHIIIJJJ
           GGGHHHHIIIJJJ
           GGGHHHHIIIJJJ
           KKKLLLLMMMNNN
           KKKLLLLMMMNNN
           KKKLLLLMMMNNN
           KKKLLLLMMMNNN'
p1 <-  mpra.plot + fc.plot + mouse_gata4_gg_tfbs + gata4_full_msa_plot + gata4_snp_plot + gata4_fc_plot + gata4_klf_plot + gata4_sox17_pair_plot + new_sox17_plot + new_ap1_plot + gata4_klf_aff.plot + gata4_sox17_pair_aff.plot + new_sox17_aff.plot + new_ap1_aff.plot + plot_layout(design = layout) + plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(size = 12, face='bold'))
ggsave('fig2.pdf', plot = p1, width = 210, height = 295, useDingbats = F, units = 'mm', device = 'pdf') 



layout <- 'AAAAAABDDDD
           AAAAAABDDDD
           AAAAAABEEEE
           AAAAAABEEEE
           CCCCCCCC###
           CCCCCCCC###
           CCCCCCCC###
           CCCCCCCC###
           CCCCCCCC###
           CCCCCCCC###'
p1 <-  gata4_gg_tfbs + gata4_seq_bar + gata.subtree + rat_mpra.plot + rat_fc.plot + plot_layout(design = layout) + plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(size = 12, face='bold'))
ggsave('sup_fig5.pdf', plot = p1, width = 210, height = 210, useDingbats = F, units = 'mm', device = 'pdf') 