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
source("utils/analyze_sup_figure7_data.R")

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

### make species mpra/logfc trajectory
plot_epas1_lineage_trajectory_act <- function(epas1_avg_lineage_mpra, stepwise_summary){
    fc.plot <- ggplot() + 
                        geom_hline(yintercept = 0, color = "black", linetype = 'dashed') +
                        geom_errorbar(data = stepwise_summary, aes(x = species2, ymin=avg_log2FC_steps-sd_log2FC_steps, 
                                                               ymax=avg_log2FC_steps+sd_log2FC_steps), 
                                      width=.2, position=position_dodge(0.05)) + 
                        geom_line(data = stepwise_summary, aes(x = species2, y = avg_log2FC_steps, group = group)) + 
                        geom_star(data = stepwise_summary, aes(x = species2, y = avg_log2FC_steps, 
                                                               starshape = marker, size = marker, fill = marker)) +  
                        scale_starshape_manual(values = c('Other' = 15, "Transition" = 9, "Extant" = 1), guide = guide_legend(override.aes = list(size = 2, fill = 'black')), name = "") +
                        scale_size_manual(values = c('Other' = 1, "Transition" = 3, "Extant" = 3), guide = 'none') +
                        scale_fill_manual(values = c('Other' = 'black', "Transition" = 'red', "Extant" = 'red'), guide = 'none') +
    #                     scale_y_log10(
    #                        breaks = c(10^-1, 10^0, 10^1),
    #                        labels = scales::trans_format("log10", scales::math_format(10^.x))
    #                      ) + annotation_logticks(sides = 'l') +
                        labs(y = "log2FC MPRA activity\nrelative to prior ancestor", x = "") + coord_cartesian(clip="off") +
                        theme_classic() + theme(legend.position = 'none', legend.box = "horizontal", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

    mpra.plot <- ggplot() + 
                        geom_rect(data = df_mouse_tiles_mpra %>% filter(full_CRE_id == 'Epas1_chr17_10063'), aes(ymin = mean_MPRA_act-sd_MPRA_act, ymax = mean_MPRA_act+sd_MPRA_act, color = NULL), xmin = -Inf, xmax = Inf,
                       alpha = .2, fill = "#008000") + 
                        geom_hline(data = df_mouse_tiles_mpra %>% filter(full_CRE_id == 'Epas1_chr17_10063'), aes(yintercept = mean_MPRA_act), color = '#008000', linetype = 'dashed') +
                        geom_rect(data = control_mean_act %>% filter(CRE_id == 'minP'), aes(ymin = mean_MPRA_act-sd_MPRA_act, ymax = mean_MPRA_act+sd_MPRA_act, color = NULL), xmin = -Inf, xmax = Inf,
                       alpha = .2,fill = "#0057e7") +
                        geom_hline(data = control_mean_act %>% filter(CRE_id == 'minP'), aes(yintercept = mean_MPRA_act), color = '#0057e7', linetype = 'dashed') +
                        geom_line(data = epas1_avg_lineage_mpra, aes(x = species2, y = mean_MPRA_act, group = group)) + 
                        geom_errorbar(data = epas1_avg_lineage_mpra, aes(x = species2, ymin=mean_MPRA_act-sd_MPRA_act, 
                                                               ymax=mean_MPRA_act+sd_MPRA_act), 
                                      width=.2, position=position_dodge(0.05)) + 
                        geom_star(data = epas1_avg_lineage_mpra, aes(x = species2, y = mean_MPRA_act, starshape = marker, size = marker, fill = marker)) +  
                        scale_starshape_manual(values = c('Other' = 15, "Transition" = 9, "Extant" = 1), guide = guide_legend(override.aes = list(size = 2, fill = 'black')), name = "") +
                        scale_size_manual(values = c('Other' = 1, "Transition" = 3, "Extant" = 3), guide = 'none') +
                        scale_fill_manual(values = c('Other' = 'black', "Transition" = 'red', "Extant" = 'red'), guide = 'none') +
    #                     scale_y_log10(
    #                        breaks = c(10^-1, 10^0, 10^1),
    #                        labels = scales::trans_format("log10", scales::math_format(10^.x))
    #                      ) + annotation_logticks(sides = 'l') +
                        labs(y = "MPRA activity", x = "") + coord_cartesian(clip="off") + 
                        theme_classic() + theme(legend.position = 'none', legend.box = "horizontal", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    return(list(mpra_plot = mpra.plot, fc_plot = fc.plot))
}


### make mouse fc plot sup fig 7a
mouse_act_plot <- plot_epas1_lineage_trajectory_act(epas1_avg_lineage_mpra_list$Mus_musculus,
                                                       stepwise_summary_list$Mus_musculus)
### make mouse fc plot sup fig 7b
epas1_ap1_map <- epas1_mouse_msa_list$msa_func_TFBS %>% arrange(msa_start) %>% dplyr::select(TF_name, rank_TF, msa_start, msa_end) %>% 
                            filter(msa_start >= 158 & msa_end <= 200)

ap1_triplet_plot <- ggplot(fixed_epas1_mouse_msa_df %>% filter(curpos >= 158) %>% filter(curpos <= 200), 
                       aes(x = curpos, y = species2, label = mut2, fill = log_change)) +
  geom_tile_rast() +
  geom_tile_rast(
    data = subset(fixed_epas1_mouse_msa_df, curpos >= 158 & curpos <= 200 & match2 == "First Match"),
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
  scale_x_continuous(breaks = seq(158, 200, 4)) + 
  theme_minimal() + 
  labs(x = "msa bp position", y = "") +
  theme(
    panel.grid = element_blank(),
    axis.ticks.y = element_blank(), # Optionally remove ticks if necessary
      legend.position = 'none'
  )
ap1_aff.mat <- new_aff.mat %>% filter(msa_start >= 158) %>% filter(msa_end <= 200) %>% 
                mutate(TF_name2 = factor(TF_name2, levels = c('Jun_Atf3_8','Jun_Atf3_2','Jun_Atf3_5')))

ap1_aff.plot <- ggplot(ap1_aff.mat %>% filter(msa_start >= 158) %>% filter(msa_end <= 200), 
                                    aes(x = log2_fc_steps, y = species2, fill = TF_name)) +
  geom_bar(stat="identity") + 
  scale_fill_manual(values = col_TFs) +
  facet_wrap(
    ~ TF_name2, 
    ncol = 3,
    labeller = as_labeller(c(
      "Jun_Atf3_8" = "AP-1\n1st site",
      "Jun_Atf3_2" = "AP-1\n2nd site",
      "Jun_Atf3_5" = "AP-1\n3rd site"
    ))
  ) +
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


### sup fig7e
dog_act_plot <- plot_epas1_lineage_trajectory_act(epas1_avg_lineage_mpra_list$Lycaon_pictus,
                                                       stepwise_summary_list$Lycaon_pictus)

species.pair <- c("Lycaon_pictus",'fullTreeAnc239')
dog.nodes <- make.nodepath(tree, species.pair)
dog.nodes2 <- c(dog.nodes, 'Mus_musculus')
dog_msa_df_TFBS <- epas1_msa_list$msa_df_TFBS %>% filter(species%in%dog.nodes2) %>% 
                    mutate(species = factor(species, levels = dog.nodes2))
labels <- collect_labels(dog_msa_df_TFBS, side = T)
dog_epas1_gg_tfbs <- plot_tfbs_trajectory(dog_msa_df_TFBS,
                                      epas1_msa_list$msa_func_TFBS, CRE_oi,
                                      'Mus_musculus',labels, 
                                      tfbs_size = 0.1, rm_species = c('Mus_musculus')) + ggtitle("Epas1:chr17_10063")

dog_seq_bar <- plot_seqID_trajectory(oCRE_mean_act, df_mouse_tiles_mpra, control_mean_act, CRE_oi, dog.nodes) + 
                    coord_cartesian(clip = 'off') + labs(x = '% sequence similarity\nto mouse ortholog', y = '')


layout <- 'AAAAABBBB
           AAAAABBBB
           CCCCCCDDD
           CCCCCCDDD
           CCCCCCDDD
           EEEEEFFFF
           EEEEEFFFF'
p1 <-  mouse_act_plot$mpra_plot + mouse_act_plot$fc_plot + ap1_triplet_plot + ap1_aff.plot  + 
    dog_act_plot$mpra_plot + dog_act_plot$fc_plot +
    plot_layout(design = layout) + plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(size = 12, face='bold'))
ggsave('sup_fig7.pdf', plot = p1, width = 210, height = 200, useDingbats = F, units = 'mm', device = 'pdf') 

layout <- 'AAAAAAABB
           AAAAAAABB'

p1 <-  dog_epas1_gg_tfbs + dog_seq_bar + 
    plot_layout(design = layout) + plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(size = 12, face='bold'))
ggsave('sup_fig7_bottom.pdf', plot = p1, width = 210, height = 80, useDingbats = F, units = 'mm', device = 'pdf') 

### multiply all log_change numbers in klf, sox, and ap1 sites
# prod(c(1.7,3.4, 1.3, 2.1))
