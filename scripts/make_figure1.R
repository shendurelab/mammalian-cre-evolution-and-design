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
library(ggtreeExtra)
library(ggstar)
source("utils/mpra_tablemaker.R")
source("utils/load_figure1_2_data.R")
source("utils/make_tree_functions.R")

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

p_count <-0 

### make basic tree
order.colors <- c("steelblue","#cc0000","#ffc425",'#905aa3','#ae5a41','#ff7f50','#1c9e48')
    names(order.colors) <- c('PRIMATES','RODENTIA','EULIPOTYPHLA','CHIROPTERA',
                             'PERISSODACTYLA','CARNIVORA','CETARTIODACTYLA')
big.tree <- as_tibble(tree)
big.tree <- groupClade(big.tree, c(293, 247,353,383,429,346,459))
species.tips <- left_join(big.tree, cactus.annotation, by = c('label')) %>% 
                mutate(Shape = case_when(
                    label == 'Mus_musculus' ~ 1,
                    label == 'Homo_sapiens' ~ 1,
                    TRUE ~ 2)) %>% 
                mutate(Size = case_when(
                    label == 'Mus_musculus' ~ 5,
                    label == 'Homo_sapiens' ~ 5,
                    TRUE ~ 0))
species.tips$Shape <- factor(species.tips$Shape)
species.tips <- as_tibble(species.tips)
big.tree <- as.treedata(big.tree)

select.tips <- species.tips %>% filter(label%in%c('Mus_musculus','Homo_sapiens'))

simple.tree <- ggtree(big.tree ,aes(color = group, shape = label), size = 0.5) %<+% filter.images +
         layout_dendrogram() +
        geom_tiplab(aes(image=uid), geom="image", size = 0.15 ,angle = 0, align = 1) +
        scale_color_manual(values = c("black", "firebrick", "steelblue","#905aa3",
                                      "#1c9e48","#ff7f50","#ffc425",'#ae5a41')) +
        theme(text = element_text(size=5), legend.position = 'none') +
        theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))
simple.tree <- simple.tree %<+% select.tips + geom_star(
                        mapping=aes(fill=group, starshape = Shape), size = 6,
                        position="identity",starstroke=1, color = 'black') +
        scale_fill_manual(values = c("black", "firebrick", "steelblue","#905aa3",
                                      "#1c9e48","#ff7f50","#ffc425",'#ae5a41'))
ggsave('example_tree.pdf', height = 4, width = 8, plot =simple.tree, device = 'pdf', dpi = 300)

length(unique(oCRE_acc_pred$tile_id)) - 1 # minus control tile
# [1] 1871
oCRE_acc_pred <- oCRE_acc_pred %>% separate(tile_id, c("full_CRE_id",'species'), sep = "__")
oCRE_acc_pred %>% filter(full_CRE_id != 'control') %>% group_by(full_CRE_id) %>% summarise(n_species = n_distinct(species)) %>% 
            mutate(prct_coverage = n_species/(240 + 239))
# full_CRE_id       n_species prct_coverage
#   <chr>                 <int>         <dbl>
# 1 Bend5_chr4_8201         277         0.578
# 2 Epas1_chr17_10063       439         0.916
# 3 Gata4_chr14_5729        479         1    
# 4 Lama1_chr17_7784        266         0.555
# 5 Sparc_chr11_7211        410         0.856
                                             
tibble.tree <- as_tibble(tree)
# heatmap.colours <- c("red","blue",'#008000','purple','orange',"pink","magenta",'black',"grey")
# names(heatmap.colours) <- c("EUARCHONTA","LAURASIATHERIA", "XENARTHRA", "AFROTHERIA", "GLIRES",
#        "MONOTREMES","GALLOANSERAE", "NEOAVES","OTHER")                                                

library(ggstar)

shape.values <- c(1,24,13,15)
names(shape.values) <- c('mouse','no ortholog','not sequenced', 'sequenced')

phyloP.labels <- phyloP.df %>% mutate(label = paste('avg. phyloP =', round(phyloP, 3)))

### check how many pass minP expression
pass <- oCRE_mean_act %>% filter(mean_MPRA_act > control_mean_act %>% filter(CRE_id == 'minP') %>% pull(mean_MPRA_act))
dim(pass)
# [1] 955   9
dim(oCRE_mean_act)
# [1] 1867   5
frac.pass <- dim(pass)[1]/dim(oCRE_mean_act)[1]
frac.pass
# [1] 0.5115158

pass.endogenous <- oCRE_mean_act %>%
  # Join with data_frame_1 on full_CRE_id
  inner_join(df_mouse_tiles_mpra, by = "full_CRE_id", suffix = c("_df2", "_df1")) %>%
  # Filter rows where mean_winsor_act from df2 is greater than in df1
  filter(mean_MPRA_act_df2 >= mean_MPRA_act_df1) %>%
  # Select only the relevant columns from df2
  dplyr::select(full_CRE_id, species_df2, mean_MPRA_act_df2, sd_MPRA_act_df2, common_name_df2)
dim(pass.endogenous)
# [1] 108   5
frac.pass <- dim(pass.endogenous)[1]/dim(oCRE_mean_act)[1]
frac.pass
# [1] 0.05784681

gata.tree <- plot_tree(tibble.tree, oCRE_mean_act %>% filter(full_CRE_id == 'Gata4_chr14_5729'),
                       phyloP.labels %>% filter(full_CRE_id == 'Gata4_chr14_5729') %>% pull(label), thickness = 0.6, size = 4, max = 2, layout = 'dendo') + ggtitle("Gata4:chr14_5729")
epas.tree <- plot_tree(tibble.tree, oCRE_mean_act %>% filter(full_CRE_id == 'Epas1_chr17_10063'),
                       phyloP.labels %>% filter(full_CRE_id == 'Epas1_chr17_10063') %>% pull(label), thickness = 0.6, size = 4,layout = 'dendo') + ggtitle("Epas1:chr17_10063")
lama.tree <- plot_tree(tibble.tree, oCRE_mean_act %>% filter(full_CRE_id == 'Lama1_chr17_7784'),
                       phyloP.labels %>% filter(full_CRE_id == 'Lama1_chr17_7784') %>% pull(label), thickness = 0.6, size = 4,layout = 'dendo') + ggtitle("Lama1:chr17_7784")
sparc.tree <- plot_tree(tibble.tree, oCRE_mean_act %>% filter(full_CRE_id == 'Sparc_chr11_7211'),
                       phyloP.labels %>% filter(full_CRE_id == 'Sparc_chr11_7211') %>% pull(label), thickness = 0.6, size = 4,layout = 'dendo') + ggtitle("Sparc:chr11_7211")
bend5.tree <- plot_tree(tibble.tree, oCRE_mean_act %>% filter(full_CRE_id == 'Bend5_chr4_8201'),
                       phyloP.labels %>% filter(full_CRE_id == 'Bend5_chr4_8201') %>% pull(label), thickness = 0.6, size = 4,layout = 'dendo') + ggtitle("Bend5:chr4_8201")

order.x <- c('RODENTIA','PRIMATES','PERISSODACTYLA','CHIROPTERA','CETARTIODACTYLA','CARNIVORA','EULIPOTYPHLA')
order.colors <- c("steelblue","#cc0000","#ffc425",'#905aa3','#ae5a41','#ff7f50','#1c9e48')
names(order.colors) <- c('PRIMATES','RODENTIA','EULIPOTYPHLA','CHIROPTERA',
                         'PERISSODACTYLA','CARNIVORA','CETARTIODACTYLA')

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
clade_oCRE_mean_act <- clade_oCRE_mean_act %>% filter(Order%in%order.x)
clade_oCRE_mean_act$Order <- factor(clade_oCRE_mean_act$Order, levels = order.x)
clade_oCRE_mean_act <- clade_oCRE_mean_act
species.count <- clade_oCRE_mean_act %>% group_by(full_CRE_id, genome) %>% summarise(n_species = n_distinct(species)) %>% 
                mutate(label = case_when(
                        genome == 'extant' ~ paste0('e = ',n_species),
                        TRUE ~ paste0('a = ',n_species))) %>% group_by(full_CRE_id) %>% 
                summarise(label = paste(label, collapse = "\n"))

### get msa distance                                                      
mpra.id <- left_join(max.dist, clade_oCRE_mean_act %>% dplyr::select(-seqid), by = c("species",'full_CRE_id'))
mpra.id <- mpra.id %>% filter(Order%in%order.x)
mpra.id$Order <- factor(mpra.id$Order, levels = order.x)

#### plot & save main figure 1

gata.dot <- make_clade_hist('Gata4_chr14_5729', clade_oCRE_mean_act, df_mouse_tiles_mpra, control_mean_act)
gata.id <- make_seq_id('Gata4_chr14_5729', mpra.id, df_mouse_tiles_mpra, control_mean_act)
gata.plot <- gata.tree + (gata.dot | gata.id) + 
  plot_layout(widths = c(1, 1)) + # Unified legend
  theme(plot.margin = margin(2, 2, 2, 2)) + # Reduce plot margins
  plot_annotation(tag_levels = NULL) # Remove annotations

epas.dot <- make_clade_hist('Epas1_chr17_10063',clade_oCRE_mean_act, df_mouse_tiles_mpra, control_mean_act)
epas.id <- make_seq_id('Epas1_chr17_10063',mpra.id, df_mouse_tiles_mpra, control_mean_act)
epas.plot <- epas.tree + (epas.dot | epas.id) + 
  plot_layout(widths = c(1, 1)) + # Unified legend
  theme(plot.margin = margin(2, 2, 2, 2)) + # Reduce plot margins
  plot_annotation(tag_levels = NULL) # Remove annotations

lama.dot <- make_clade_hist('Lama1_chr17_7784',clade_oCRE_mean_act, df_mouse_tiles_mpra, control_mean_act)
lama.id <- make_seq_id('Lama1_chr17_7784',mpra.id, df_mouse_tiles_mpra, control_mean_act)
lama.plot <- lama.tree + (lama.dot | lama.id)+ 
  plot_layout(widths = c(1, 1)) + # Unified legend
  theme(plot.margin = margin(2, 2, 2, 2)) + # Reduce plot margins
  plot_annotation(tag_levels = NULL) # Remove annotations

combined_plot <- gata.plot / epas.plot / lama.plot + 
  plot_layout(heights = c(1, 1, 1))  # Equal heights for all three plots
combined_plot <- combined_plot + 
  theme(plot.margin = margin(2, 2, 2, 2))  # Reduce outer margins


ggsave('fig1.pdf', plot = combined_plot, width = 210, height = 300, useDingbats = F, units = 'mm')


### plot supplementary figure 3
sparc.dot <- make_clade_hist('Sparc_chr11_7211',clade_oCRE_mean_act, df_mouse_tiles_mpra, control_mean_act)
sparc.id <- make_seq_id('Sparc_chr11_7211',mpra.id, df_mouse_tiles_mpra, control_mean_act)
sparc.plot <- sparc.tree + (sparc.dot | sparc.id)+ 
  plot_layout(widths = c(1, 1)) + # Unified legend
  theme(plot.margin = margin(2, 2, 2, 2)) + # Reduce plot margins
  plot_annotation(tag_levels = NULL) # Remove annotations

bend5.dot <- make_clade_hist('Bend5_chr4_8201',clade_oCRE_mean_act, df_mouse_tiles_mpra, control_mean_act)
bend5.id <- make_seq_id('Bend5_chr4_8201',mpra.id, df_mouse_tiles_mpra, control_mean_act)
bend5.plot <- bend5.tree + (bend5.dot | bend5.id)+ 
  plot_layout(widths = c(1, 1)) + # Unified legend
  theme(plot.margin = margin(2, 2, 2, 2)) + # Reduce plot margins
  plot_annotation(tag_levels = NULL) # Remove annotations


mpra.id$full_CRE_id <- sub("_", ":", mpra.id$full_CRE_id)
mpra.id$full_CRE_id <- factor(mpra.id$full_CRE_id , levels = c("Gata4:chr14_5729","Epas1:chr17_10063","Lama1:chr17_7784",
                                                              "Sparc:chr11_7211","Bend5:chr4_8201"))
seqId.plot <- ggplot() + geom_quasirandom_rast(data = mpra.id, aes(y = Order, x = seqid, color = Order, shape = genome, alpha = genome), 
        dodge.width=1, size = 1, stroke = 1.2) +  # Adjust the width for separation 
            facet_wrap(~full_CRE_id, ncol = 5) + scale_color_manual(values = order.colors, guide = 'none') + theme_classic() +
            scale_alpha_manual(values = c(0.3, 0.8)) +
#             scale_fill_manual(values = c("ancestor" = "white", "extant" = "black")) + 
            scale_shape_manual(values = c(23, 21), guide = guide_legend(override.aes = list(size = 2, color = 'black'))) +
                        labs(x = '% sequence similarity to mouse ortholog', y = '') + theme(
                                                                                             legend.position = 'bottom',
                                                                                            panel.spacing = unit(0.7, "lines"))

combined_plot <- sparc.plot / bend5.plot  + 
  plot_layout(heights = c(1, 1))  # Equal heights for all three plots
combined_plot <- combined_plot + 
  theme(plot.margin = margin(2, 2, 2, 2))  # Reduce outer margins

ggsave('supp_fig3_tree.pdf', plot = combined_plot, width = 210, height = 200, useDingbats = F, units = 'mm')

cor.df <- mpra.id %>% group_by(full_CRE_id, Order) %>% summarise(r = cor(mean_MPRA_act, seqid, method = 'spearman'), 
                                                                 mean_seqid = mean(seqid),sd_MPRA_act = sd(mean_MPRA_act), 
                                                                 mean_MPRA_act = mean(mean_MPRA_act), n = n_distinct(species))

corr.plot <- ggplot(cor.df, aes(x = full_CRE_id, y = r, color = Order)) +
  geom_point_rast(size = 2) +
    scale_color_manual(values = order.colors) +  
    theme_classic() + scale_y_continuous(breaks = c(-1, -0.5, 0, 0.5, 1), limits = c(-1,1)) +
  labs(x = "CRE", y = "Spearman correlation\n(MPRA-to-mouse sequence similarity)") +
  theme(panel.spacing = unit(1, "lines"),
       axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
                    

cor.matrix <- dcast(cor.df, full_CRE_id ~ Order, value.var = "r")
rownames(cor.matrix) <- cor.matrix$full_CRE_id
cor.matrix <- cor.matrix[,-1]
cor.matrix <- cor.matrix[,order.x]#^2

library(pheatmap)
library(viridis)

order_colour <- list(Order = c(PRIMATES = "steelblue", RODENTIA = "#cc0000", CARNIVORA = "#ff7f50", CETARTIODACTYLA = "#1c9e48",
                           EULIPOTYPHLA = "#ffc425", CHIROPTERA = "#905aa3", PERISSODACTYLA = "#ae5a41"))
anno_row = data.frame(Order = order.x)
rownames(anno_row) <- anno_row$Order

library(ggplotify)
# Define blue-white-red color palette
blue_white_red <- colorRampPalette(c("blue", "white", "red"))(256)
heat <- pheatmap(cor.matrix, annotation_colors = order_colour, cluster_cols = F, cluster_rows = F , angle_col = 90,
                 annotation_col = anno_row, color  = blue_white_red,
                  breaks = seq(-1, 1, length.out = 257),  # Range from -1 to 1, 
                 show_colnames = F,
                 cellheight = 20, cellwidth = 40, fontsize = 9)
heat <- as.ggplot(heat)

layout <- 'AAAAAAAAA
           #########
           BBBBBBBBB'
combined_plot <- seqId.plot + heat + plot_layout(design = layout)
ggsave('supp_fig3.pdf', plot = combined_plot, width = 210, height = 180, useDingbats = F, units = 'mm')

### plot supplementary figure 4
dms.avg <- mouse_DMS_MPRA %>% filter(log2FC_MPRA != 0) %>% ### filter WT nucleotide affect which is zero
            filter(mut != 'del') %>% group_by(pos, CRE) %>% 
            summarise(mean_log2FC_MPRA = mean(log2FC_MPRA), median_log2FC_MPRA = median(log2FC_MPRA))


dms.phyloP <- left_join(dms.avg, singleBase.phyloP, by = c("pos",'CRE'))
mouse.TFBS_map <- mouse_func_tfbs %>%
  rowwise() %>%
  mutate(pos = list(seq(TFBS_start, TFBS_end))) %>%  # Create a sequence of positions
  unnest(pos) %>%  # Expand the list column
  group_by(CRE) %>%
  summarize(pos = list(unique(pos)))  # Get unique positions per CRE

new.dms_phyloP <- data.frame()
for(CRE_oi in unique(dms.phyloP$CRE)){
    cre_tfbs_map <- mouse.TFBS_map %>% filter(CRE == CRE_oi) %>% pull(pos)
    cre_phyloP <- dms.phyloP %>% filter(CRE == CRE_oi) %>% mutate(is_tfbs = case_when(
                                        pos%in%cre_tfbs_map[[1]] ~ TRUE,
                                        TRUE ~ FALSE))
    new.dms_phyloP <- bind_rows(new.dms_phyloP, cre_phyloP)
}

dms.phyloP <- new.dms_phyloP

# write.table(dms.phyloP, 'PYS2_DMS_CRE_phyloP_linked.txt', sep = '\t', row.names = F, quote = F)

dms.filter <- dms.phyloP %>% filter(!is.na(phyloP)) %>% #filter(CRE != 'Bend5_chr4_8201') %>% 
                separate(CRE, c("CRE","chr",'id')) %>% mutate(CRE = paste0(CRE,":",chr,"_",id)) %>% dplyr::select(-chr,-id)
dms.filter$CRE <- factor(dms.filter$CRE, levels = c("Gata4:chr14_5729","Epas1:chr17_10063","Lama1:chr17_7784","Sparc:chr11_7211",
                                                   'Bend5:chr4_8201'))
dms.cor <- dms.filter %>% group_by(CRE) %>% summarise(r = cor(mean_log2FC_MPRA, phyloP, method = 'spearman'))
dms.cor$CRE <- factor(dms.cor$CRE, levels = c("Gata4:chr14_5729","Epas1:chr17_10063","Lama1:chr17_7784","Sparc:chr11_7211",
                                             'Bend5:chr4_8201'))
dms.tfbs_stats <- dms.filter %>%
              group_by(CRE) %>%
              summarise(
                p_value = wilcox.test(phyloP[is_tfbs == TRUE], phyloP[is_tfbs == FALSE])$p.value
              ) %>% mutate(p_label = paste0("p = ", round(p_value, 3)))
dms.tfbs_stats$CRE <- factor(dms.tfbs_stats$CRE, levels = c("Gata4:chr14_5729","Epas1:chr17_10063","Lama1:chr17_7784","Sparc:chr11_7211",'Bend5:chr4_8201'))
dms.tfbs_ci <- dms.filter %>% group_by(CRE, is_tfbs) %>% summarise(mean_val = mean(phyloP), sd_val = sd(phyloP), count = n()) %>%
                mutate(ci_lower = mean_val - 1.96 * sd_val / sqrt(count),
                      ci_upper = mean_val - 1.96 * sd_val / sqrt(count))

dms.plot <- ggplot() + 
    geom_hline(yintercept = 0, color = 'black', alpha = 0.4, linetype = 'dashed') + 
    geom_point_rast(data = dms.filter, aes(x = phyloP, y = mean_log2FC_MPRA), 
               pch = 21, color = 'black', alpha = 0.6) + 
    geom_text(data = dms.cor, aes(label = paste0("rho = ", round(r, 3))), color = 'black', x = Inf, y = Inf, 
              hjust = 1, vjust = 1.5, check_overlap = T) + labs(x = 'phyloP', y = 'avg. log2(MPRA FC)') +
    facet_wrap(~CRE, ncol = 5) + theme_classic()

dms_cdf <- ggplot() + 
        stat_ecdf(data = dms.filter, aes(x = phyloP, color = is_tfbs), geom = 'step') +
        geom_text(data = dms.tfbs_stats, aes(label = p_label), color = 'black', x = Inf, y = -Inf, 
              hjust = 1, vjust = -0.7, check_overlap = T) +
        scale_color_manual(values = c('black','red')) + facet_wrap(~CRE, ncol = 5) +
        theme_classic() +
        labs(x = 'phyloP', y = "Cumulative distribution") +
        theme(legend.position = 'bottom')

final_dms.plot <- dms.plot / dms_cdf

ggsave('supp_fig4.pdf', plot = final_dms.plot, width = 210, height = 120, useDingbats = F, units = 'mm', device = 'pdf')

### write species.count as wide table: full_CRE_id x Order
species.count.by.order <- clade_oCRE_mean_act %>%
  group_by(full_CRE_id, Order) %>%
  summarise(n_species = n_distinct(species), .groups = 'drop')
species.count.wide <- dcast(species.count.by.order, full_CRE_id ~ Order, value.var = "n_species", fill = 0)
write.table(species.count.wide, 'species_count_by_order.txt', sep = '\t', row.names = F, quote = F)