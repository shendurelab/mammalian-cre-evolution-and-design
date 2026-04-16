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
source("utils/figure2_func.R")
source("utils/analyze_figure3_data.R")

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

#### make Epas1 tree Fig3a
epas.tree <- ggtree(sub.tree, aes(color = mean_MPRA_act, shape = capture), size = 0.4) + #geom_treescale(x=0, y=235) + 
#             layout_dendrogram() +
            geom_treescale(x = 0.4, y = 0) +
            scale_color_gradient(breaks = c(0.1, 1, 10), limits = c(0.05, 13),  # Quantile limits from your range
                                   labels = scales::trans_format("log10", scales::math_format(10^.x)),
                                    oob = scales::squish,  # To handle values outside the specified limits
                                    trans = "log10",
                                   low = '#42047e',high = '#07f49e', na.value = 'gainsboro') +
            theme(legend.position = 'bottom')
epas.tree <- epas.tree %<+% sub.tips + geom_star(
                        mapping=aes(fill=mean_MPRA_act, starshape = Shape), size = 3,
                        position="identity",starstroke=1, color = 'black', show_guide = FALSE) +
            scale_fill_gradient(breaks = c(0.1, 1, 10), limits = c(0.05, 13),  # Quantile limits from your range
                                  labels = scales::trans_format("log10", scales::math_format(10^.x)),
                                    oob = scales::squish,  # To handle values outside the specified limits
                                    trans = "log10",
                                   low = '#42047e',high = '#07f49e', na.value = 'gainsboro') +
            labs(color = 'MPRA activity', fill = 'MPRA activity') + ggtitle("Epas1:chr17_10063") + 
            theme(legend.position = 'bottom',
                    legend.spacing = unit(0, "cm"),        # Remove extra spacing between legend items
                    legend.box.margin = margin(-10, 0, 0, 0), # Reduce margin below the legend
                    plot.title.position = "plot",         # Align title closer to the plot
                    plot.margin = margin(5, 5, 5, 5),     # Reduce overall plot margin
                  panel.spacing = unit(0.1, "lines"))
d <- data.frame(node=c(293, 247,429,383), Order=c("RODENTIA", "PRIMATES",'CARNIVORA','CETARTIODACTYLA'),
               active = c('yes','no','yes','no'))
epas.tree <- epas.tree + new_scale_color() + 
    geom_hilight(data=d, aes(node=node, color=active),fill =NA,
                            type = "rect", linetype = 'dashed') +
  scale_color_manual(values = c("yes" = "red", "no" = "blue"), name = 'Active Clade?')


# heatmap for triplet pairs and 3' end AP1 aff
aff.hm <- ggplot(aff.mat.tree, aes(x = TF_name2, y = label, fill = hue_color)) +
  geom_tile() +
  # white-to-color gradient mapped by TF_name
    scale_fill_identity() +  # Use precomputed colors
  theme_classic() + labs(x = "predicted TFBS\naffinity", y = "") + 
  # add vertical borders between columns
  geom_vline(
    xintercept = seq(1.5, length(unique(aff.mat.tree$TF_name2)) - 0.5, 1),
    color = "black", linewidth = 0.5
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.y = element_blank(),
  axis.ticks.y = element_blank(),
  axis.title.y=element_blank(),
  axis.text.y=element_blank(),
      legend.position = 'bottom',
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)  # <-- add black border
  )

ap1.hm <- ggplot(ap1.mat.tree, aes(x = TF_name2, y = label, fill = hue_color)) +
  geom_tile() +
  # white-to-color gradient mapped by TF_name
    scale_fill_identity() +  # Use precomputed colors
  theme_classic() + labs(x = "mean predicted\n3' AP-1 affinity", y = "") + 
  theme(
    axis.text.x = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.y = element_blank(),
  axis.ticks.y = element_blank(),
  axis.title.y=element_blank(),
  axis.text.y=element_blank(),
      legend.position = 'bottom',
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)  # <-- add black border
  )

### template mouse TFBS map
mouse_epas1_gg_tfbs <- plot_tfbs_trajectory(epas1_msa_list$msa_df_TFBS %>% filter(species == 'Mus_musculus'),
                                      epas1_msa_list$msa_func_TFBS, 
                                      CRE_oi, 'Mus_musculus',
                                      c('Mouse'), 
                                      tfbs_size = 0.1) + labs(x = "", y ="") + 
                 geom_rect(data = tfbs_map, aes(xmin = msa_start, xmax = msa_end+9, ymin = -Inf, ymax = Inf),
                          inherit.aes = FALSE, fill = NA, color = 'black', linetype = 'dashed', alpha = 0.7) +
                theme(axis.line.y = element_blank(),
                  axis.ticks.y = element_blank(),
                  axis.title.y=element_blank(),
                      axis.text.x = element_blank(),
                     axis.line.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.title.x=element_blank())


library(aplot)
options("aplot_guides" = "keep")
# p <- aff.hm %>% insert_left(klf4.hm, width = 0.2) %>% insert_left(epas.tree, width = 1.7)  %>% 
#             insert_right(ap1.hm, width = 0.2)
p <- aff.hm  %>% insert_left(epas.tree, width = 1.7) %>% 
        insert_right(ap1.hm, width = 0.2)
ggsave("fig3_tree_map.pdf", plot = p, width = 130, height = 160,
       units = "mm", device = 'pdf', useDingbats = FALSE)
ggsave("fig3_mouse_tfbs_legend.pdf", plot = mouse_epas1_gg_tfbs, width = 160, height = 10,
       units = "mm", device = 'pdf', useDingbats = FALSE)

### plot primate lineage
gata_primate_plot <- ggplot(primate.epas1_msa_info$msa_df %>% filter(curpos >= 674) %>% filter(curpos <= 687), aes(x = curpos, y = species2, label = mut2, fill = log_change)) +
  geom_tile_rast() +
  geom_tile_rast(
    data = subset(primate.epas1_msa_info$msa_df, curpos >= 674 & curpos <= 687 & match2 == "Mutation"),
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
  geom_text(aes(fontface = ifelse(match2 == "Mutation", "bold",
                                 ifelse(match2 == "DEL", "bold", "plain")), 
                alpha = ifelse(match2 == "Mutation", 1,
                             ifelse(match2 != "Match", 1, 0.9)),
               size = I(ifelse(mut2 == "*", 5, 3))  # treat as fixed, no scaling
               ),color = 'black',
           show.legend = FALSE) +
  scale_x_continuous(breaks = seq(660, 686, 4)) + 
  theme_minimal() + 
  labs(x = "", y = "") + ggtitle('MONKEYS') +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid = element_blank(),
    axis.ticks.y = element_blank(), # Optionally remove ticks if necessary
      legend.position = 'none'
  )

mouse_y <- which(levels(primate.epas1_msa_info$msa_df$species2) == "Mouse")

gata_primate_plot <- gata_primate_plot +
  geom_hline(yintercept = mouse_y - 0.5,     # line just *below* Mouse row
             color = "black", linetype = "solid", linewidth = 0.2)


gata4_primate_aff.plot <- ggplot(primate.epas1_msa_info$msa_func_aff_TFBS %>% filter(msa_start >= 674) %>% filter(msa_end <= 687), aes(x = TF_name2, y = species2, fill = hue_color)) +
  geom_point(shape = 21, color = "black", size = 3, stroke = 0.5) +  # black border
  scale_fill_identity() +  # Use precomputed fill colors
  theme_classic() + labs(x = "", y = "") + 
  # add vertical borders between columns
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(), # Optionally remove ticks if necessary
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(), # Optionally remove ticks if necessary
      axis.line = element_blank(),
      legend.position = 'none'
  )

sox_lemur_plot <- ggplot(lemur.epas1_msa_info$msa_df %>% filter(curpos >= 660) %>% filter(curpos <= 672), aes(x = curpos, y = species2, label = mut2, fill = log_change)) +
  geom_tile_rast() +
  geom_tile_rast(
    data = subset(lemur.epas1_msa_info$msa_df, curpos >= 660 & curpos <= 672 & match2 == "Mutation"),
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
  geom_text(aes(fontface = ifelse(match2 == "Mutation", "bold",
                                 ifelse(match2 == "DEL", "bold", "plain")), 
                alpha = ifelse(match2 == "Mutation", 1,
                             ifelse(match2 != "Match", 1, 0.9)),
               size = I(ifelse(mut2 == "*", 5, 3))  # treat as fixed, no scaling
               ),color = 'black',
           show.legend = FALSE) +
  scale_x_continuous(breaks = seq(660, 672, 4)) + 
  theme_minimal() + 
  labs(x = "", y = "") + ggtitle('LEMURS') +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5),
    axis.ticks.y = element_blank(), # Optionally remove ticks if necessary
      legend.position = 'none'
  )

mouse_y <- which(levels(lemur.epas1_msa_info$msa_df$species2) == "Mouse")

sox_lemur_plot <- sox_lemur_plot +
  geom_hline(yintercept = mouse_y - 0.5,     # line just *below* Mouse row
             color = "black", linetype = "solid", linewidth = 0.2)

sox_lemur_aff.plot <- ggplot(lemur.epas1_msa_info$msa_func_aff_TFBS %>% filter(msa_start >= 660) %>% filter(msa_end <= 672), aes(x = TF_name2, y = species2, fill = hue_color)) +
  geom_point(shape = 21, color = "black", size = 3, stroke = 0.5) +  # black border
  scale_fill_identity() +  # Use precomputed fill colors
  theme_classic() + labs(x = "", y = "") + 
  # add vertical borders between columns
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(), # Optionally remove ticks if necessary
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(), # Optionally remove ticks if necessary
      axis.line = element_blank(),
      legend.position = 'none'
  )

### plot giraffe lineage
sox_giraffe_plot <- ggplot(giraffe.epas1_msa_info$msa_df %>% filter(curpos >= 660) %>% filter(curpos <= 672), aes(x = curpos, y = species2, label = mut2, fill = log_change)) +
  geom_tile_rast() +
  geom_tile_rast(
    data = subset(giraffe.epas1_msa_info$msa_df, curpos >= 660 & curpos <= 672 & match2 == "Mutation"),
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
  geom_text(aes(fontface = ifelse(match2 == "Mutation", "bold",
                                 ifelse(match2 == "DEL", "bold", "plain")), 
                alpha = ifelse(match2 == "Mutation", 1,
                             ifelse(match2 != "Match", 1, 0.9)),
               size = I(ifelse(mut2 == "*", 5, 3))  # treat as fixed, no scaling
               ),color = 'black',
           show.legend = FALSE) +
  scale_x_continuous(breaks = seq(660, 672, 4)) + 
  theme_minimal() + ggtitle('CETARTIODACTYLA') +
  labs(x = "msa bp position", y = "") +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5),
    axis.ticks.y = element_blank(), # Optionally remove ticks if necessary
      legend.position = 'none'
  )

mouse_y <- which(levels(giraffe.epas1_msa_info$msa_df$species2) == "Mouse")

sox_giraffe_plot <- sox_giraffe_plot +
  geom_hline(yintercept = mouse_y - 0.5,     # line just *below* Mouse row
             color = "black", linetype = "solid", linewidth = 0.2)

sox_giraffe_aff.plot <- ggplot(giraffe.epas1_msa_info$msa_func_aff_TFBS %>% filter(msa_start >= 660) %>% filter(msa_end <= 672), aes(x = TF_name2, y = species2, fill = hue_color)) +
  geom_point(shape = 21, color = "black", size = 3, stroke = 0.5) +  # black border
  scale_fill_identity() +  # Use precomputed fill colors
  theme_classic() + labs(x = "predicted\naffinity", y = "") + 
  # add vertical borders between columns
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(), # Optionally remove ticks if necessary
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(), # Optionally remove ticks if necessary
      axis.line = element_blank(),
      legend.position = 'none'
  )

ap1_end_sand_rat_plot <- ggplot(sand_rat.epas1_msa_info$msa_df %>% filter(curpos >= 753), aes(x = curpos, y = species2, label = mut2, fill = log_change)) +
  geom_tile_rast() +
  geom_tile_rast(
    data = subset(sand_rat.epas1_msa_info$msa_df, curpos >= 753 & match2 == "Mutation"),
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
  geom_text(aes(fontface = ifelse(match2 == "Mutation", "bold",
                                 ifelse(match2 == "DEL", "bold", "plain")), 
                alpha = ifelse(match2 == "Mutation", 1,
                             ifelse(match2 != "Match", 1, 0.9)),
               size = I(ifelse(mut2 == "*", 5, 3))  # treat as fixed, no scaling
               ),color = 'black',
           show.legend = FALSE) +
  scale_x_continuous(breaks = seq(754, 776, 4)) + 
  theme_minimal() + 
  labs(x = "", y = "") + 
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5),
    axis.ticks.y = element_blank(), # Optionally remove ticks if necessary
      legend.position = 'none'
  )

mouse_y <- which(levels(sand_rat.epas1_msa_info$msa_df$species2) == "Mouse") - 1

ap1_end_sand_rat_plot <- ap1_end_sand_rat_plot +
  geom_hline(yintercept = mouse_y - 0.5,     # line just *below* Mouse row
             color = "black", linetype = "solid", linewidth = 0.2)

ap1_sand_rat_aff.plot <- ggplot(new_sand_rat_ap1.mat, aes(x = log2_fc_steps, y = species2, fill = TF_name)) +
  geom_bar(stat="identity") + 
  scale_fill_manual(values = col_TFs) +
  labs(x = '', y = "") +
  geom_vline(xintercept = 0, color = "black", linetype = 'dashed') +
#   scale_x_continuous(breaks = c(0, 0.5, 1)) + 
  theme_classic() + 
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(), # Optionally remove ticks if necessary
      legend.position = 'none'
  )

high_ap1_site <- data.frame(msa_start = c(744), msa_end = c(764))
ap1_end_dog_plot <- ggplot(dog.epas1_msa_info$msa_df %>% filter(curpos >= 706 & curpos <=767), aes(x = curpos, y = species2, label = mut2, fill = log_change)) +
  geom_tile_rast() +
  geom_tile_rast(
    data = subset(dog.epas1_msa_info$msa_df, curpos >= 706 & curpos <= 767 & match2 == "Mutation"),
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
  geom_text(aes(fontface = ifelse(match2 == "Mutation", "bold",
                                 ifelse(match2 == "DEL", "bold", "plain")), 
                alpha = ifelse(match2 == "Mutation", 1,
                             ifelse(match2 != "Match", 1, 0.9)),
               size = I(ifelse(mut2 == "*", 5, 3))  # treat as fixed, no scaling
               ),color = 'black',
           show.legend = FALSE) +
  geom_rect(data = high_ap1_site, aes(xmin = msa_start-0.5, xmax = msa_end+0.5, ymin = -Inf, ymax = Inf),
                          inherit.aes = FALSE, fill = NA, color = 'firebrick2', linetype = 'dashed', alpha = 0.7) +
  scale_x_continuous(breaks = seq(706, 766, 4)) + 
  theme_minimal() + 
  labs(x = "msa bp position", y = "") + 
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5),
    axis.ticks.y = element_blank(), # Optionally remove ticks if necessary
      legend.position = 'none'
  )

mouse_y <- which(levels(dog.epas1_msa_info$msa_df$species2) == "Mouse")

ap1_end_dog_plot <- ap1_end_dog_plot +
  geom_hline(yintercept = mouse_y - 0.5,     # line just *below* Mouse row
             color = "black", linetype = "solid", linewidth = 0.2)


ap1_dog_aff.plot <- ggplot(new_dog_ap1.mat, aes(x = TF_name2, y = species2, fill = hue_color)) +
  geom_point(shape = 21, color = "black", size = 3, stroke = 0.5) +  # black border
  scale_fill_identity() +  # Use precomputed fill colors
  theme_classic() + labs(x = "predicted\naffinity", y = "") + 
  # add vertical borders between columns
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(), # Optionally remove ticks if necessary
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(), # Optionally remove ticks if necessary
      axis.line = element_blank(),
      legend.position = 'none'
  )

dog_epas1_gg_tfbs <- plot_tfbs_trajectory(epas1_msa_list$msa_df_TFBS %>% 
                                          filter(species%in%c('Lycaon_pictus','Mus_musculus')),
                                      epas1_msa_list$msa_func_TFBS, 
                                      CRE_oi, 'Mus_musculus',
                                      c('African_dog'), 
                                      tfbs_size = 0.1, rm_species = c('Mus_musculus')) + labs(x = "", y ="") + 
                 geom_rect(data = tfbs_map, aes(xmin = msa_start, xmax = msa_end+9, ymin = -Inf, ymax = Inf),
                          inherit.aes = FALSE, fill = NA, color = 'black', linetype = 'dashed', alpha = 0.7) +
                theme(axis.line.y = element_blank(),
                  axis.ticks.y = element_blank(),
                  axis.title.y=element_blank(),
                      axis.text.x = element_blank(),
                     axis.line.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.title.x=element_blank())


layout <- '#######AAAB
           #######AAAB
           #######CCCD
           #######CCCD
           #######EEEF
           #######EEEF
           #######EEEF
           GGGGGGGGGGG  
           HHHHHHHHHHI
           HHHHHHHHHHI
           HHHHHHHHHHI'
p <- gata_primate_plot + gata4_primate_aff.plot + sox_lemur_plot + sox_lemur_aff.plot + sox_giraffe_plot + sox_giraffe_aff.plot + dog_epas1_gg_tfbs + ap1_end_dog_plot + ap1_dog_aff.plot + 
plot_layout(design = layout) + plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(size = 12, face='bold'))
ggsave('fig3.pdf', width = 210, height = 250,
       units = "mm", device = 'pdf', useDingbats = FALSE)
