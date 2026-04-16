### module load pcre2/10.39; module load R/4.3.1

library(stringr)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(ggrastr)
library(ggstar)
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


plot_tree <- function(tibble.tree, cre_mean_act, phylop.label , thickness = 0.75, size = 3, layout = 'circular', max= 10){
    shape.values <- c(1,24,13,15)
    names(shape.values) <- c('mouse','no ortholog','not sequenced', 'sequenced')
    sub.seq <- cre_mean_act[,c('species','mean_MPRA_act','common_name')]
    sub.seq$label <- sub.seq$species
    
    max.value <- cre_mean_act %>% filter(species == 'Mus_musculus') %>% pull(mean_MPRA_act)
    
    sub.seq <- sub.seq %>% mutate(capture = case_when(
                            species == 'Mus_musculus' ~ 'mouse',
                            is.na(mean_MPRA_act) ~ 'not sequenced',
                            TRUE ~ 'sequenced'))
    
    sub.tree <- full_join(tibble.tree, sub.seq, by = 'label')
    sub.tree <- sub.tree %>% mutate(capture = case_when(
                            is.na(capture) ~ 'no ortholog',
                            TRUE ~ capture))
    sub.tips <- sub.tree %>% filter(label%in%c('Mus_musculus','Homo_sapiens')) %>%
                mutate(Shape = case_when(
                    label == 'Mus_musculus' ~ 1,
                    label == 'Homo_sapiens' ~ 1,
                    TRUE ~ 2)) %>% 
                mutate(Size = case_when(
                    label == 'Mus_musculus' ~ 5,
                    label == 'Homo_sapiens' ~ 5,
                    TRUE ~ 0))
    sub.tips$Shape <- factor(sub.tips$Shape)
    sub.tips <- as.tibble(sub.tips)
    sub.tree <- as.treedata(sub.tree)
    
    if(layout == 'circular'){
        p <- ggtree(sub.tree, layout = 'circular', aes(color = mean_MPRA_act), size = thickness) +
            geom_text(x = Inf, y = Inf, label = phylop.label, col = "black", hjust = 1, vjust = 1, check_overlap = TRUE) + 
            scale_color_gradient(breaks = c(0.1, 1, 10), limits = c(0.05, 13),  # Quantile limits from your range
                                   labels = scales::trans_format("log10", scales::math_format(10^.x)),
                                    oob = scales::squish,  # To handle values outside the specified limits
                                    trans = "log10",
                                   low = '#42047e',high = '#07f49e', na.value = 'gainsboro') 
        p <- p %<+% sub.tips + geom_star(
                        mapping=aes(fill=mean_MPRA_act, starshape = Shape), size = size,
                        position="identity",starstroke=1, color = 'black', show_guide = FALSE) +
            scale_fill_gradient(breaks = c(0.1, 1, 10), #limits = c(0.05, max.value),  # Quantile limits from your range
                                  labels = scales::trans_format("log10", scales::math_format(10^.x)),
                                    oob = scales::squish,  # To handle values outside the specified limits
                                    trans = "log10",
                                   low = '#42047e',high = '#07f49e', na.value = 'gainsboro') +
            labs(color = 'MPRA activity', fill = 'MPRA activity') + 
            theme(legend.position = 'bottom')
    }else{
    
        p <- ggtree(sub.tree, aes(color = mean_MPRA_act, shape = capture), size = thickness) + #geom_treescale(x=0, y=235) + 
            layout_dendrogram() +
            geom_treescale(x = 0, y = 30) +
        geom_text(x = Inf, y = Inf, label = phylop.label, col = "black", hjust = 1, vjust = 1, check_overlap = TRUE) + 
    #         geom_tiplab(aes(label = common_name), color = 'black') + 
    #         geom_tippoint(aes(colour=mean_winsor_act+p_count)) + 
    #         geom_nodepoint(aes(colour=mean_winsor_act+p_count)) +
#             geom_star(aes(fill=mean_MPRA_act+p_count ,colour=mean_MPRA_act+p_count, starshape=capture),
#                             position="identity") + 
#         geom_cladelab(node=293, label="", barcolor = order.colors['RODENTIA'],
#                      offset = .1, barsize=2, angle=0, hjust=0.5, vjust = -1, fontsize=5, align = T) + 
#         geom_cladelab(node=247, label="", barcolor = order.colors['PRIMATES'],
#                      offset = .1, barsize=2, angle=0, hjust=0.5, vjust = -1, fontsize=5, align = T) + 
#         geom_cladelab(node=429, label="", barcolor = order.colors['CARNIVORA'],
#                      offset = .1, barsize=2, angle=0, hjust=0.5, vjust = -1, fontsize=5, align = T) + 
#         geom_cladelab(node=346, label="", barcolor = order.colors['EULIPOTYPHLA'],angle = 45,
#                      offset = .1, barsize=2, hjust=0.5, vjust = -1.5, fontsize=5, align = T) + 
#         geom_cladelab(node=459, label="", barcolor = order.colors['PERISSODACTYLA'],angle = 45,
#                      offset = .1, barsize=2, hjust=0.5, vjust = -1.5, fontsize=5, align = T) + 
#         geom_cladelab(node=353, label="", barcolor = order.colors['CHIROPTERA'],
#                      offset = .1, barsize=2, angle=0, hjust=0.5, vjust = -1, fontsize=5, align = T) + 
#         geom_cladelab(node=383, label="", barcolor = order.colors['CETARTIODACTYLA'],
#                      offset = .1, barsize=2, angle=0, hjust=0.5, vjust = -1, fontsize=5, align = T) + 
              scale_color_gradient(breaks = c(0.1, 1, 10), limits = c(0.05, 13),  # Quantile limits from your range
                                   labels = scales::trans_format("log10", scales::math_format(10^.x)),
                                    oob = scales::squish,  # To handle values outside the specified limits
                                    trans = "log10",
                                   low = '#42047e',high = '#07f49e', na.value = 'gainsboro') +
#               scale_fill_gradient(breaks = c(0.1, 1, 10), limits = c(0.05, 10),  # Quantile limits from your range
#                                   labels = scales::trans_format("log10", scales::math_format(10^.x)),
#                                     oob = scales::squish,  # To handle values outside the specified limits
#                                     trans = "log10",
#                                    low = '#42047e',high = '#07f49e', na.value = 'gainsboro') +

#             labs(color = 'MPRA act.', fill = 'MPRA act.') + 
#             scale_starshape_manual(values = shape.values,name = '',
#                                   guide = guide_legend(override.aes = list(size = 5, fill = 'black'))) + 
            theme(legend.position = 'bottom')
        p <- p %<+% sub.tips + geom_star(
                        mapping=aes(fill=mean_MPRA_act, starshape = Shape), size = size,
                        position="identity",starstroke=1, color = 'black', show_guide = FALSE) +
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
    }
    return(p)
}

plot_subtrees <- function(nodes,cre){
    p <- plot_spacer()
    sub.seq <- oCRE_mean_act %>% filter(full_CRE_id == cre)
    sub.seq <- sub.seq[,c('species','mean_winsor_act','common_name')]
    sub.seq$label <- sub.seq$species

    sub.seq <- sub.seq %>% mutate(capture = case_when(
                            species == 'Mus_musculus' ~ 'mouse',
                            is.na(mean_winsor_act) ~ 'not sequenced',
                            TRUE ~ 'sequenced'))
    select.images <- Sys.glob('/net/shendure/vol10/projects/tli/nobackup/EB_dev_project/util_scripts/species_images/main_tree_images/*.svg')
    select.images <- data.frame(uid = select.images,
                           label = c("Hipposideros_armiger","Felis_catus","Erinaceus_europaeus",
                                    "Equus_caballus","Homo_sapiens","Mus_musculus","Delphinapterus_leucas"))
    rownames(select.images) <- select.images$label
    select.images <- left_join(big.tree, select.images, by = 'label')
    select.images <- as.tibble(select.images)
    
    
    for(n in nodes){
        subtree <- tree_subset(tree, n, levels_back = 4)
        tibble.sub_tree <- as_tibble(subtree)

        sub.tree <- left_join(tibble.sub_tree, sub.seq, by = 'label')
        sub.tree <- as.treedata(sub.tree)

        x <- ggtree(sub.tree, aes(color = mean_winsor_act), size = 0.75) %<+% select.images + 
#             layout_dendrogram() +
              scale_color_gradient(breaks = c(0.01, 0.1, 1), limits = c(0.01, 3),  # Quantile limits from your range
                                    oob = scales::squish,  # To handle values outside the specified limits
                                    trans = "log10",
                                   low = '#42047e',high = '#07f49e', na.value = 'gainsboro') +
            labs(color = 'MPRA activity') + 
            new_scale_color() +
            geom_tiplab(aes(image=uid), geom="image", size = 0.05) +             
            theme(legend.position = 'top')
        
        p <- p + x + plot_layout(guides = "collect") & theme(legend.position = 'bottom')
    }
    return(p)
}

order.x <- c('RODENTIA','PRIMATES','PERISSODACTYLA','CHIROPTERA','CETARTIODACTYLA','CARNIVORA','EULIPOTYPHLA')

make_clade_hist <- function(CRE_oi, oCRE_mean_act, mouse_wt_act, control_mean_act){
    library(ggbeeswarm)
    x_breaks = c(10^-1, 10^0, 10^1, 10^2)
    cre_mean_act <- oCRE_mean_act %>% filter(full_CRE_id == CRE_oi)
    mouse_wt_act <- mouse_wt_act %>% filter(full_CRE_id == CRE_oi)
    cre.count <- species.count %>% filter(full_CRE_id == CRE_oi)
    cre_mean_act$Order <- factor(cre_mean_act$Order, levels = order.x)

    dot.plot <- ggplot() + geom_quasirandom_rast(data = cre_mean_act, aes(x = Order, y = mean_MPRA_act, color = Order, 
                                                                          shape = genome, alpha = genome), 
            dodge.width=1, size = 1, stroke = 1.1) +  # Adjust the width for separation 
                geom_rect(data = mouse_wt_act, aes(ymin = mean_MPRA_act-sd_MPRA_act, ymax = mean_MPRA_act+sd_MPRA_act), xmin = -Inf, xmax = Inf,
               alpha = .2,fill = "#008000") +
                geom_rect(data = control_mean_act %>% filter(CRE_id == 'minP'), aes(ymin = mean_MPRA_act-sd_MPRA_act, ymax = mean_MPRA_act+sd_MPRA_act), xmin = -Inf, xmax = Inf,
               alpha = .1,fill = "blue")  +
            geom_hline(yintercept = control_mean_act %>% filter(CRE_id == 'minP') %>% pull(mean_MPRA_act) ,color = "blue", linetype = 'dashed') + 
                geom_hline(data = mouse_wt_act, aes(yintercept = mean_MPRA_act), color = "#008000", linetype = 'dashed') + 
                scale_y_log10(
                   breaks = c(10^-1, 10^0, 10^1), limits = c(0.1, 13),
                   labels = scales::trans_format("log10", scales::math_format(10^.x))
                 ) + annotation_logticks(sides="l") + 
                geom_text(data = cre.count, aes(label = label), x = Inf, y = Inf,
                col = "black", hjust = 1, vjust = 1.2, size = 3.5, check_overlap = TRUE) +
                scale_color_manual(values = order.colors, guide = 'none') + 
#                 scale_fill_manual(values = c("ancestor" = "white", "extant" = "black")) +
                scale_shape_manual(values = c(23, 21), guide = guide_legend(override.aes = list(size = 1, color = 'black'))) +
                scale_alpha_manual(values = c(0.3, 0.8)) + theme_classic() + scale_x_discrete(drop = FALSE)  +
                            labs(y = 'MPRA activity', x = '') + 
                theme(legend.position = 'bottom', 
                      axis.text.x = element_text(angle = 45, hjust=1, color = "black"),
                      axis.text.y = element_text(color = "black"))
    return(dot.plot)
}

make_seq_id <- function(CRE_oi, seq.id, mouse_wt_act, control_mean_act){
    cre_id <- seq.id %>% filter(full_CRE_id == CRE_oi)
    mouse_wt_act <- mouse_wt_act %>% filter(full_CRE_id == CRE_oi)
    
    seq.plot <- ggplot(cre_id) + geom_point_rast(aes(x = seqid, y = mean_MPRA_act, color = Order, shape = genome, alpha = genome), 
                                            size = 1, stroke = 1.1) +
        scale_alpha_manual(values = c(0.3, 0.8)) + 
        geom_rect(data = control_mean_act %>% filter(CRE_id == 'minP'), aes(ymin = mean_MPRA_act-sd_MPRA_act, ymax = mean_MPRA_act+sd_MPRA_act), xmin = -Inf, xmax = Inf,
               alpha = .1,fill = "blue")  +
        geom_hline(yintercept = control_mean_act %>% filter(CRE_id == 'minP') %>% pull(mean_MPRA_act) ,color = "blue", linetype = 'dashed') + 
        geom_rect(data = mouse_wt_act, aes(ymin = mean_MPRA_act-sd_MPRA_act, ymax = mean_MPRA_act+sd_MPRA_act), xmin = -Inf, xmax = Inf,
               alpha = .2,fill = "#008000") + 
        geom_hline(data = mouse_wt_act, aes(yintercept = mean_MPRA_act), color = "#008000", linetype = 'dashed') + 
        scale_x_continuous(limits = c(40, 100), breaks = c(40,60, 80, 100)) +
        scale_y_log10(
           breaks = c(10^-1, 10^0, 10^1), limits = c(0.1, 13),
           labels = scales::trans_format("log10", scales::math_format(10^.x))
         ) + annotation_logticks(sides="l") + scale_color_manual(values = order.colors) + theme_classic() + 
#         scale_fill_manual(values = c("ancestor" = "white", "extant" = "black")) +
        scale_shape_manual(values = c(23, 21), guide = guide_legend(override.aes = list(size = 2, color = 'black'))) +
        labs(y = 'MPRA activity', x = '% sequence similarity\nto mouse ortholog') + 
        theme(legend.position = 'none')
}
