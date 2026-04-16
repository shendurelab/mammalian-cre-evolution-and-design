library(tidyverse)       # Includes ggplot2, stringr, dplyr, readxl (replaces individual loading of these packages)
library(cowplot)
library(ggrastr)
library(castor)
library(reshape2)         # Only needs to be loaded once
library(ggridges)
library(DescTools)
library(patchwork)
library(gtools)
suppressMessages(suppressWarnings(library(GenomicFeatures)))
suppressMessages(suppressWarnings(library(GenomicRanges)))
suppressMessages(suppressWarnings(library(rtracklayer)))
library(ape)
library(phytools)
library(ggtree)
library(tidytree)
library(ggtreeExtra)
library(IRanges)
library(seqinr)
library(Biostrings)
library(ggstar)
library(ggtext)
require(grid)
library(ggnewscale)

source("utils/mafft_aligning_functions.R")
source("utils/mpra_tablemaker.R")
source("utils/load_figure4_data.R")
source("utils/figure4_func.R")

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

#### plot schema
traj_oi <- 'Gata4_chr14_5729__fullTreeAnc239__Mus_musculus'                                     
example.traj <- evo_model_mean_act %>% filter(traj == traj_oi) %>% filter((rank != 0 & !startsWith(Event, 'target'))) %>% arrange(rank) %>% 
                filter(rank <= 5) %>% select(rank, ref_start) %>% distinct()
                                     
set.seed(42) # For reproducibility
example.traj$ref_start <- sapply(1:nrow(example.traj), function(i) {
  if (i > 1 && abs(example.traj$ref_start[i] - example.traj$ref_start[i - 1]) < 10) {
    # Add or subtract 20 randomly
    example.traj$ref_start[i] + sample(c(-40, 40), 1)
  } else {
    # Keep the original value
    example.traj$ref_start[i]
  }
})

example.mod <- data.frame()
for(rank_oi in unique(example.traj$rank)){
    tmp <- example.traj %>% filter(rank == rank_oi)
    tmp <- bind_rows(example.mod, tmp)
    tmp$order <- rank_oi
    example.mod <- bind_rows(example.mod, tmp)
}
                                     
example.mod <- example.mod %>% distinct()
                               
model.scheme <- ggplot(example.mod, aes(y = order))+
                geom_hline(aes(yintercept = order), color = 'gainsboro', linewidth = 1) + ## add horizontal line
                geom_star(aes(x =ref_start + 1, y = order, fill = rank), 
                          size = 3, starshape = 12, alpha = 0.7, color = 'black', position = position_nudge(y = 0, x = 4)) + 
                scale_fill_gradient(low = "#D6E4FF",high = '#1F3A5F') + 
                scale_y_reverse() + coord_cartesian(xlim=c(0,300), clip="off") + theme_classic() +
                labs(x = 'Model\nreconstitution') + 
                theme(legend.position = 'none',
                          axis.text.y = element_blank(),
                          axis.line.y = element_blank(),
                          axis.ticks.y = element_blank(),
                          axis.title.y=element_blank(),
                          axis.text.x = element_blank(),
                          axis.line.x = element_blank(),
                          axis.ticks.x = element_blank()) 

set.seed(42)
random.traj <- example.traj
random.traj$ref_start <- sample(random.traj$ref_start)
random.mod <- data.frame()
for(rank_oi in unique(random.traj$rank)){
    tmp <- random.traj %>% filter(rank == rank_oi)
    tmp <- bind_rows(random.mod, tmp)
    tmp$order <- rank_oi
    random.mod <- bind_rows(random.mod, tmp)
} 
random.mod <- random.mod %>% distinct()
                                     
                                     
random.scheme <- ggplot(random.mod, aes(y = order))+
                geom_hline(aes(yintercept = order), color = 'gainsboro', linewidth = 1) + ## add horizontal line
                geom_star(aes(x =ref_start + 1, y = order), starshape = 12, size = 3, alpha =0.7,
                                color = 'black', fill = '#B0B0B0', position = position_nudge(y = 0, x = 4)) + theme_classic() + 
                scale_y_reverse() + coord_cartesian(xlim=c(0,300), clip="off") + labs(x = 'Randomized\nreconstitution') + 
                theme(legend.position = 'none',
                          axis.text.y = element_blank(),
                          axis.line.y = element_blank(),
                          axis.ticks.y = element_blank(),
                          axis.title.y=element_blank(),
                          axis.text.x = element_blank(),
                          axis.line.x = element_blank(),
                          axis.ticks.x = element_blank())

labels <- c()

for(species_oi in c('fullTreeAnc239',"Mus_musculus")){
    if(species_oi%in%select.images$label){
        img_path <- select.images %>% filter(label == species_oi) %>% pull(uid)
        name_oi <- select.images %>% filter(label == species_oi) %>% pull(name)
        labels <- c(labels, paste0(
                              "<span style='vertical-align: middle;'>*", name_oi, "*</span> ",
                              "<img src='", img_path, "' width='35' style='vertical-align: middle;'/>"
                            ))
    } else{
        labels <- c(labels, NA)
    }
} 
                                    
example.aln <- example.traj 
example.aln$rank <- 'Mus_musculus'
example.aln <- data.frame(bind_rows(example.aln, data.frame(rank = c("fullTreeAnc239"), ref_start = c(NA))))
example.aln$rank <- factor(example.aln$rank, levels = c('fullTreeAnc239',"Mus_musculus"))                               
aln.scheme <- ggplot(example.aln, aes(y = rank))+
                geom_hline(aes(yintercept = rank, color = rank), linewidth = 1) + ## add horizontal line
                scale_color_manual(values = c('gainsboro','red')) + 
                geom_star(aes(x =ref_start), starshape = 12, size = 3, 
                                color = 'white', fill = 'black', alpha = 0.7) + 
                # Add arrows pointing down from the stars
                  geom_segment(
                    aes(x = ref_start, xend = ref_start, y = as.numeric(rank) - 0.2, yend = as.numeric(rank) - 0.7),  # Adjust length with yend
                    arrow = arrow(type = "closed", angle = 30, length = unit(0.15, "cm")),
                    color = "black", linewidth = 0.6
                  ) +
                coord_cartesian(xlim=c(0,300), clip = 'off') + theme_classic() +
                scale_y_discrete(labels = labels) + 
                theme(legend.position = 'none',
                          axis.text.y = element_markdown(color = "black", size = 10),
                          axis.line.y = element_blank(),
                          axis.ticks.y = element_blank(),
                          axis.title.y=element_blank(),
                          axis.text.x = element_blank(),
                          axis.line.x = element_blank(),
                          axis.ticks.x = element_blank(),
                          axis.title.x=element_blank())

top.plot <-  model.scheme + random.scheme +
  plot_layout(design = "
              AABB
              AABB
              ") 
ggsave('fig5_aln_scheme.pdf', plot = top.plot, width = 60, height = 40, useDingbats = F, units = 'mm', device = 'pdf')

### get start, end, and best results from model
cre_levels <- c('Gata4:chr14_5729','Epas1:chr17_10063','Lama1:chr17_7784')
mpra.points <- data.frame()
random.points <- data.frame()
atac.points <- data.frame()

for(traj_oi in traj_oi_list){
    best_pid <- evo_model_mean_act %>% filter(traj == traj_oi) %>% 
                    filter(mean_MPRA_act >= ref.events %>% filter(traj == traj_oi) %>% pull(mean_MPRA_act)) %>%
                    arrange(PID) %>% head(n = 1) %>% mutate(point_type = "best")
    last_pid <- evo_model_mean_act %>% filter(traj == traj_oi) %>% filter(rank == max(rank)) %>% mutate(point_type = "end")
    target_pid <- target.events %>% filter(traj == traj_oi) %>% mutate(point_type = "start")
    mpra.points <- bind_rows(mpra.points, last_pid)
    mpra.points <- bind_rows(mpra.points, best_pid)
    mpra.points <- bind_rows(mpra.points, target_pid)
}

for(traj_oi in traj_oi_list){
    best_pid <- evo_random_mean_act %>% filter(traj == traj_oi) %>% 
                    filter(mean_MPRA_act >= ref.events %>% filter(traj == traj_oi) %>% pull(mean_MPRA_act)) %>%
                    arrange(PID) %>% head(n = 1) %>% mutate(point_type = "best")
    last_pid <- evo_random_mean_act %>% filter(traj == traj_oi) %>% filter(rank == max(rank)) %>% mutate(point_type = "end")
    target_pid <- target.events %>% filter(traj == traj_oi) %>% mutate(point_type = "start")
    random.points <- bind_rows(random.points, last_pid)
    random.points <- bind_rows(random.points, best_pid)
    random.points <- bind_rows(random.points, target_pid)
}

for(traj_oi in traj_oi_list){
    best_pid <- evo_model_mean_act %>% filter(traj == traj_oi) %>% 
                    filter(norm_footprint >= ref.events %>% filter(traj == traj_oi) %>% pull(norm_footprint)) %>%
                    arrange(PID) %>% head(n = 1) %>% mutate(point_type = "best")
    last_pid <- evo_model_mean_act %>% filter(traj == traj_oi) %>% filter(rank == max(rank)) %>% mutate(point_type = "end")
    target_pid <- target.events %>% filter(traj == traj_oi) %>% mutate(point_type = "start")
    atac.points <- bind_rows(atac.points, last_pid)
    atac.points <- bind_rows(atac.points, best_pid)
    atac.points <- bind_rows(atac.points, target_pid)
}
mpra.points <- mpra.points %>%
  mutate(
    CRE = sub("_", ":", CRE)  # Replace underscores with colons
  )
mpra.points$CRE <- factor(mpra.points$CRE, levels = cre_levels)

random.points <- random.points %>%
  mutate(
    CRE = sub("_", ":", CRE)  # Replace underscores with colons
  )
random.points$CRE <- factor(random.points$CRE, levels = cre_levels)

atac.points <- atac.points %>%
  mutate(
    CRE = sub("_", ":", CRE)  # Replace underscores with colons
  )
atac.points$CRE <- factor(atac.points$CRE, levels = cre_levels)

best_mpra.points <- bind_rows(mpra.points, random.points) %>% filter(point_type == 'best')
                                     
#### plot atac vs. mpra evo model correlations
all.select <- bind_rows(evo_model_mean_act.select %>% mutate(point_type = 'Model'), 
                        new_evo_random_mean_act.select %>% mutate(point_type = 'Random'))
all.select <- bind_rows(all.select, 
                       oCRE_mean_select %>% mutate(point_type = 'Mouse_Lineage'))
all.select$CRE <- factor(all.select$CRE, levels = cre_levels)
corr.df <- all.select %>% group_by(CRE) %>% summarise(rho = cor(log2(norm_footprint), mean_MPRA_act, method = 'spearman'))
all.select <- all.select %>% 
    left_join(corr.df, by = 'CRE') %>% 
  mutate(corr_label = paste0(CRE, "\nrho = ", round(rho, 3)))  # use for labeling only
all.select$point_type <- factor(all.select$point_type, levels = c('Model','Mouse_Lineage','Random'))

corr.plot <- ggplot() +
                geom_point_rast(data = all.select, aes(y = log2(norm_footprint), x = mean_MPRA_act, color = point_type, 
                                                       alpha = point_type), 
                                shape = 21, fill = 'white') +
                geom_hline(data = select.ref, aes(yintercept = log2(norm_footprint)), color = '#008000', linetype = 'dashed', show.legend = FALSE) +
                geom_vline(data = select.ref, aes(xintercept = mean_MPRA_act), color = '#008000', linetype = 'dashed', show.legend = FALSE) +
                geom_rect(data = select.ref, aes(xmin = mean_MPRA_act-sd_MPRA_act, xmax = mean_MPRA_act+sd_MPRA_act, color = NULL), 
                          ymin = -Inf, ymax = Inf, show.legend = FALSE, alpha = .2, fill = "#008000") + 
                scale_color_manual(values = c('Model' = '#1F3A5F', 'Random' = '#B0B0B0', 'Mouse_Lineage' = 'firebrick'), name = "") +
                scale_alpha_manual(values = c('Model' = 0.8, 'Random' = 0.5, 'Mouse_Lineage' = 0.8), name = "") +
                geom_star(data = atac.points, aes(x = mean_MPRA_act, y = log2(norm_footprint), fill = point_type, starshape = point_type), 
                              size = 3, show.legend = FALSE) +  
                scale_fill_manual(values = c('end' = '#1F3A5F', 'best' = '#ffd700', 'start' = '#D6E4FF')) + 
                scale_starshape_manual(values = c('start' = 15, 'best' = 1, 'end' = 15)) + 
                labs(y = expression(paste("predicted ", log[2], "(accessibility)")), 
                     x = 'MPRA activity') + 
#                 geom_label(data = corr.df, aes(label = paste0('rho = ', round(rho, 2))), x= Inf, y = -Inf,
#                           hjust = 1,  # Align to the right
#                         vjust = -0.7, size = 5) +   # Align to the bottom  
                facet_wrap(~CRE, ncol = 3, labeller = as_labeller(setNames(all.select$corr_label, all.select$CRE))) + 
                scale_x_log10(
                       breaks = c(10^-2,10^-1, 10^0, 10^1), limits = c(0.01, 10),
                       labels = scales::trans_format("log10", scales::math_format(10^.x))
                     ) + annotation_logticks(sides = 'b') + 
                scale_y_continuous(breaks = c(0,5,10)) +  # Format x-axis to show whole numbers
                theme_classic() + theme(#legend.position = 'bottom',
                                           panel.spacing = unit(1, "lines"),
                                       legend.position = c(0.08, 0.85),  # Adjust position inside the plot
                                      legend.background = element_rect(fill = NA, color = NA), # Semi-transparent background
                                      legend.key = element_blank()) # Remove legend key background
                                  
# corr2.plot <- ggplot() +
#                 geom_point_rast(data = all.select, aes(y = log2(norm_footprint), x = mean_MPRA_act, color = point_type), 
#                                 alpha = 0.6, shape = 21, fill = 'white') +
#                 geom_hline(data = select.ref, aes(yintercept = log2(norm_footprint)), color = '#008000', linetype = 'dashed', show.legend = FALSE) +
#                 geom_vline(data = select.ref, aes(xintercept = mean_MPRA_act), color = '#008000', linetype = 'dashed', show.legend = FALSE) +
#                 geom_rect(data = select.ref, aes(xmin = mean_MPRA_act-sd_MPRA_act, xmax = mean_MPRA_act+sd_MPRA_act, color = NULL), 
#                           ymin = -Inf, ymax = Inf, show.legend = FALSE, alpha = .2, fill = "#008000") + 
#                 scale_color_manual(values = c('Model' = '#1F3A5F', 'Random' = '#B0B0B0'), name = "") +
#                 labs(y = 'predicted log2(accessibility)', x = 'MPRA activity') + 
#                 geom_label(data = corr.df, aes(label = paste0('rho = ', round(rho, 2))), x= Inf, y = -Inf,
#                           hjust = 1,  # Align to the right
#                         vjust = -0.7, size = 5) +   # Align to the bottom  
#                 facet_wrap(~CRE, ncol = 3) + 
#                 scale_x_log10(
#                        breaks = c(10^-2,10^-1, 10^0, 10^1), limits = c(0.01, 10),
#                        labels = scales::trans_format("log10", scales::math_format(10^.x))
#                      ) + annotation_logticks(sides = 'b') + 
#                 scale_y_continuous(breaks = c(0,5,10)) +  # Format x-axis to show whole numbers
#                 theme_classic() + theme(#legend.position = 'bottom',
#                                            panel.spacing = unit(1, "lines"),
#                                        legend.position = c(0.05, 0.97),  # Adjust position inside the plot
#                                       legend.background = element_rect(fill = NA, color = NA), # Semi-transparent background
#                                       legend.key = element_blank()) # Remove legend key background
# ggsave("presentation_plots/corr_acc_act.pdf", plot = corr2.plot, height = 4, width = 12, dpi = 300, device = 'pdf')

## plot ATAC prediction plots
atac.pid <- ggplot() +   
                    geom_step(data = evo_model_mean_act.select, aes(x = 100-PID, y = log2(norm_footprint), color = 'Model')) + 
                    geom_step(data = new_evo_random_mean_act.select, aes(x = 100-PID, y = log2(norm_footprint), group = iter, color = 'Random'), 
                              alpha = 0.5) +
                    geom_step(data = oCRE_mean_select, aes(x = 100-PID, y = log2(norm_footprint), color = 'Mouse_Lineage')) +
                    scale_color_manual(values = c('Model' = '#1F3A5F', 'Random' = '#B0B0B0', 'Mouse_Lineage' = 'firebrick'), name = "") +  
#                     scale_color_gradient(limits = c(0,50), oob = scales::squish, name = "model-ranked\nmutations",
#                                         breaks = c(0, 25, 50), labels = c(0, 25, '50+')) +
                    geom_hline(data = select.ref, aes(yintercept = log2(norm_footprint)), color = '#008000', linetype = 'dashed') +
                    geom_hline(data = select.target, aes(yintercept = log2(norm_footprint)), color = "black", linetype = 'dashed') +
                    geom_star(data = atac.points, aes(x = 100-PID, y = log2(norm_footprint), fill = point_type, starshape = point_type), 
                              size = 3, show.legend = FALSE) +  
                    scale_fill_manual(values = c('end' = '#1F3A5F', 'best' = '#ffd700', 'start' = '#D6E4FF')) + 
                    scale_starshape_manual(values = c('start' = 15, 'best' = 1, 'end' = 15)) + 
                    labs(x = '% sequence difference to mouse ortholog', y = expression(paste("predicted ", log[2], "(accessibility)"))) +
                    scale_y_continuous(breaks = c(0, 5, 10)) + scale_x_reverse() +
                    facet_wrap(~CRE, ncol = 3, scales = "free_x") +  # Allow each facet to have its own x scale
                    theme_classic() + theme(legend.position = 'bottom', legend.box = "horizontal",
                                           panel.spacing = unit(1, "lines"))
# atac.rank <- ggplot() +   
#                     geom_step(data = evo_model_mean_act.select, aes(x = rank, y = log2(norm_footprint), color = 'Model')) + 
#                     geom_step(data = new_evo_random_mean_act.select, aes(x = rank, y = log2(norm_footprint), group = iter, color = 'Random'), 
#                               alpha = 0.5) +
#                     scale_color_manual(values = c('Model' = '#1F3A5F', 'Random' = '#B0B0B0')) + 
# #                     scale_color_gradient(limits = c(0,50), oob = scales::squish, name = "model-ranked\nmutations",
# #                                         breaks = c(0, 25, 50), labels = c(0, 25, '50+')) +
#                     geom_hline(data = select.ref, aes(yintercept = log2(norm_footprint)), color = '#008000', linetype = 'dashed') +
#                     geom_hline(data = select.target, aes(yintercept = log2(norm_footprint)), color = "black", linetype = 'dashed') +
#                     geom_star(data = atac.points, aes(x = rank, y = log2(norm_footprint), fill = point_type, starshape = point_type), 
#                               size = 3, show.legend = FALSE) +  
# #                     new_scale_fill() +
#                     scale_fill_manual(values = c('end' = '#1F3A5F', 'best' = '#ffd700', 'start' = '#D6E4FF')) + 
#                     scale_starshape_manual(values = c('start' = 15, 'best' = 1, 'end' = 15)) + 
#                     labs(x = 'mutation step', y = 'predicted log2(accessibility)') + scale_y_continuous(breaks = c(0, 5, 10)) +
#                     facet_wrap(~CRE, ncol = 1, scales = "free_x") +  # Allow each facet to have its own x scale
#                     theme_classic() + theme(legend.position = 'bottom', legend.box = "horizontal",
#                                            panel.spacing = unit(1, "lines"))
## plot MPRA plots 
mpra.pid <- ggplot() +   
                    geom_step(data = new_evo_random_mean_act.select, aes(x = 100-PID, y = mean_MPRA_act, group = iter, color = 'Random'), 
                              alpha = 0.5) +
                    geom_step(data = oCRE_mean_select, aes(x = 100-PID, y = mean_MPRA_act, color = 'Mouse_Lineage')) +
                    geom_step(data = evo_model_mean_act.select, aes(x = 100-PID, y = mean_MPRA_act, color = 'Model')) + 
                    scale_color_manual(values = c('Model' = '#1F3A5F', 'Random' = '#B0B0B0', 'Mouse_Lineage' = 'firebrick'), name = "") + 
#                     scale_color_gradient(limits = c(0,50), oob = scales::squish, name = "model-ranked\nmutations",
#                                         breaks = c(0, 25, 50), labels = c(0, 25, '50+')) +
                    geom_rect(data = select.ref, aes(ymin = mean_MPRA_act-sd_MPRA_act, ymax = mean_MPRA_act+sd_MPRA_act, color = NULL), 
                              xmin = -Inf, xmax = Inf, alpha = .2, fill = "#008000", show.legend = FALSE) + 
                    geom_hline(data = select.ref, aes(yintercept = mean_MPRA_act), color = '#008000', linetype = 'dashed', show.legend = FALSE) +
                    geom_rect(data = select.target, aes(ymin = mean_MPRA_act-sd_MPRA_act, ymax = mean_MPRA_act+sd_MPRA_act, color = NULL), 
                              xmin = -Inf, xmax = Inf, show.legend = FALSE, alpha = .2,fill = "black") +
                    geom_hline(data = select.target, aes(yintercept = mean_MPRA_act), color = "black", linetype = 'dashed',
                              show.legend = FALSE) +
                    geom_star(data = mpra.points, aes(x = 100-PID, y = mean_MPRA_act, fill = point_type, starshape = point_type), 
                              size = 3, show.legend = FALSE) +  
#                     new_scale_fill() +
                    scale_fill_manual(values = c('end' = '#1F3A5F', 'best' = '#ffd700', 'start' = '#D6E4FF')) + 
                    scale_starshape_manual(values = c('start' = 15, 'best' = 1, 'end' = 15)) + 
                    labs(x = '% sequence difference to mouse ortholog', y = 'MPRA activity') +
                    scale_y_log10(
                       breaks = c(10^-2,10^-1, 10^0, 10^1), limits = c(0.01, 10),
                       labels = scales::trans_format("log10", scales::math_format(10^.x))
                     ) + annotation_logticks(sides = 'l') + scale_x_reverse() +
                    facet_wrap(~CRE, ncol = 3, scales = "free_x") +  # Allow each facet to have its own x scale
                    theme_classic() + theme(legend.position = 'bottom', legend.box = "horizontal",
                                           panel.spacing = unit(1, "lines")) 
# mpra.rank <- ggplot() +   
#                     geom_step(data = evo_model_mean_act.select, aes(x = rank, y = mean_MPRA_act, color = 'Model')) + 
#                     geom_step(data = new_evo_random_mean_act.select, aes(x = rank, y = mean_MPRA_act, group = iter, color = 'Random'), 
#                               alpha = 0.5) +
#                     scale_color_manual(values = c('Model' = '#1F3A5F', 'Random' = '#B0B0B0')) + 
# #                     scale_color_gradient(limits = c(0,50), oob = scales::squish, name = "model-ranked\nmutations",
# #                                         breaks = c(0, 25, 50), labels = c(0, 25, '50+')) +
#                     geom_rect(data = select.ref, aes(ymin = mean_MPRA_act-sd_MPRA_act, ymax = mean_MPRA_act+sd_MPRA_act, color = NULL), xmin = -Inf, xmax = Inf,
#                    alpha = .2, fill = "#008000") + 
#                     geom_hline(data = select.ref, aes(yintercept = mean_MPRA_act), color = '#008000', linetype = 'dashed') +
#                     geom_rect(data = select.target, aes(ymin = mean_MPRA_act-sd_MPRA_act, ymax = mean_MPRA_act+sd_MPRA_act, color = NULL), xmin = -Inf, xmax = Inf,
#                    alpha = .2,fill = "black") +
#                     geom_hline(data = select.target, aes(yintercept = mean_MPRA_act), color = "black", linetype = 'dashed') +
#                     geom_star(data = mpra.points, aes(x = rank, y = mean_MPRA_act, fill = point_type, starshape = point_type), 
#                               size = 3, show.legend = FALSE) +  
# #                     new_scale_fill() +
#                     scale_fill_manual(values = c('end' = '#1F3A5F', 'best' = '#ffd700', 'start' = '#D6E4FF')) + 
#                     scale_starshape_manual(values = c('start' = 15, 'best' = 1, 'end' = 15)) + 
#                     labs(x = 'mutation step', y = 'MPRA activity') + 
#                     scale_y_log10(
#                        breaks = c(10^-2,10^-1, 10^0, 10^1), limits = c(0.01, 10),
#                        labels = scales::trans_format("log10", scales::math_format(10^.x))
#                      ) + annotation_logticks(sides = 'l') +
#                     facet_wrap(~CRE, ncol = 1, scales = "free_x") +  # Allow each facet to have its own x scale
#                     theme_classic() + theme(legend.position = 'bottom', legend.box = "horizontal",
#                                            panel.spacing = unit(1, "lines"))                                      

# mpra2.pid <- ggplot() +   
#                     geom_step(data = evo_model_mean_act.select, aes(x = 100-PID, y = mean_MPRA_act, color = 'Model')) + 
#                     geom_step(data = new_evo_random_mean_act.select, aes(x = 100-PID, y = mean_MPRA_act, group = iter, color = 'Random'), 
#                               alpha = 0.5) +
#                     geom_step(data = oCRE_mean_select, aes(x = 100-PID, y = mean_MPRA_act, color = 'Mouse_Lineage'), 
#                               linetype = 'dashed') +
#                     scale_color_manual(values = c('Model' = '#1F3A5F', 'Random' = '#B0B0B0', 'Mouse_Lineage' = '#2A9D8F'), name = "") + 
# #                     scale_color_gradient(limits = c(0,50), oob = scales::squish, name = "model-ranked\nmutations",
# #                                         breaks = c(0, 25, 50), labels = c(0, 25, '50+')) +
#                     geom_rect(data = select.ref, aes(ymin = mean_MPRA_act-sd_MPRA_act, ymax = mean_MPRA_act+sd_MPRA_act, color = NULL), 
#                               xmin = -Inf, xmax = Inf, alpha = .2, fill = "#008000", show.legend = FALSE) + 
#                     geom_hline(data = select.ref, aes(yintercept = mean_MPRA_act), color = '#008000', linetype = 'dashed', show.legend = FALSE) +
#                     geom_rect(data = select.target, aes(ymin = mean_MPRA_act-sd_MPRA_act, ymax = mean_MPRA_act+sd_MPRA_act, color = NULL), 
#                               xmin = -Inf, xmax = Inf, show.legend = FALSE, alpha = .2,fill = "black") +
#                     geom_hline(data = select.target, aes(yintercept = mean_MPRA_act), color = "black", linetype = 'dashed',
#                               show.legend = FALSE) +
#                     geom_star(data = mpra.points, aes(x = 100-PID, y = mean_MPRA_act, fill = point_type, starshape = point_type), 
#                               size = 5, show.legend = FALSE) +  
# #                     new_scale_fill() +
#                     scale_fill_manual(values = c('end' = '#1F3A5F', 'best' = '#ffd700', 'start' = '#D6E4FF')) + 
#                     scale_starshape_manual(values = c('start' = 15, 'best' = 1, 'end' = 15)) + 
#                     labs(x = '% sequence difference (mouse CRE)', y = 'MPRA activity') +
#                     scale_y_log10(
#                        breaks = c(10^-2,10^-1, 10^0, 10^1), limits = c(0.01, 10),
#                        labels = scales::trans_format("log10", scales::math_format(10^.x))
#                      ) + annotation_logticks(sides = 'l') + scale_x_reverse() +
#                     facet_wrap(~CRE, ncol = 3, scales = "free_x") +  # Allow each facet to have its own x scale
#                     theme_classic() + theme(legend.position = 'bottom', legend.box = "horizontal", text = element_text(family = "sans", size = 18),
#                                            panel.spacing = unit(1, "lines")) 

epas1_mpra.pid <- ggplot() +   
                    geom_step(data = new_evo_random_mean_act.select %>% filter(CRE == 'Epas1:chr17_10063'), 
                              aes(x = 100-PID, y = mean_MPRA_act, group = iter, color = 'Random'), 
                              alpha = 0.5) +
                    geom_step(data = oCRE_mean_select %>% filter(CRE == 'Epas1:chr17_10063'), 
                              aes(x = 100-PID, y = mean_MPRA_act, color = 'Mouse_Lineage')) +
                    geom_step(data = evo_model_mean_act.select %>% filter(CRE == 'Epas1:chr17_10063'), 
                              aes(x = 100-PID, y = mean_MPRA_act, color = 'Model')) + 
                    scale_color_manual(values = c('Model' = '#1F3A5F', 'Random' = '#B0B0B0', 'Mouse_Lineage' = 'firebrick'), name = "") + 
#                     scale_color_gradient(limits = c(0,50), oob = scales::squish, name = "model-ranked\nmutations",
#                                         breaks = c(0, 25, 50), labels = c(0, 25, '50+')) +
                    geom_rect(data = select.ref %>% filter(CRE == 'Epas1:chr17_10063'), 
                              aes(ymin = mean_MPRA_act-sd_MPRA_act, ymax = mean_MPRA_act+sd_MPRA_act, color = NULL), 
                              xmin = -Inf, xmax = Inf, alpha = .2, fill = "#008000", show.legend = FALSE) + 
                    geom_hline(data = select.ref %>% filter(CRE == 'Epas1:chr17_10063'), 
                               aes(yintercept = mean_MPRA_act), color = '#008000', linetype = 'dashed', show.legend = FALSE) +
                    geom_rect(data = select.target %>% filter(CRE == 'Epas1:chr17_10063'), 
                              aes(ymin = mean_MPRA_act-sd_MPRA_act, ymax = mean_MPRA_act+sd_MPRA_act, color = NULL), 
                              xmin = -Inf, xmax = Inf, show.legend = FALSE, alpha = .2,fill = "black") +
                    geom_hline(data = select.target %>% filter(CRE == 'Epas1:chr17_10063'),
                               aes(yintercept = mean_MPRA_act), color = "black", linetype = 'dashed',
                              show.legend = FALSE) +
                    geom_star(data = mpra.points %>% filter(CRE == 'Epas1:chr17_10063'), 
                              aes(x = 100-PID, y = mean_MPRA_act, fill = point_type, starshape = point_type), 
                              size = 3, show.legend = FALSE) +  
#                     new_scale_fill() +
                    scale_fill_manual(values = c('end' = '#1F3A5F', 'best' = '#ffd700', 'start' = '#D6E4FF')) + 
                    scale_starshape_manual(values = c('start' = 15, 'best' = 1, 'end' = 15)) + 
                    labs(x = '% sequence difference to mouse ortholog', y = 'MPRA activity') +
                    scale_y_log10(
                       breaks = c(10^-2,10^-1, 10^0, 10^1), limits = c(0.01, 10),
                       labels = scales::trans_format("log10", scales::math_format(10^.x))
                     ) + annotation_logticks(sides = 'l') + scale_x_reverse() +
                    facet_wrap(~CRE, ncol = 3, scales = "free_x") +  # Allow each facet to have its own x scale
                    theme_classic() + theme(legend.position = 'bottom', legend.box = "horizontal",
                                           panel.spacing = unit(1, "lines")) 

# ggsave("presentation_plots/mpra_traj_act.pdf", plot = mpra2.pid, height = 4, width = 12, dpi = 300, device = 'pdf')

### plot log2FC mutations to mapped functional TFBS
load("../extdata/evo_delta_info.RData")
### plot functional tfbs with respect to MPRA for trajetories
load("../extdata/evo_func_tfbs_info.RData")
###
load("../extdata/all_evo_delta_info.RData")

model.all <- left_join(rank_func_tfbs.count, all.events %>% dplyr::select(tile_name,frac_mean_delta,func_tfbs), by = c("tile_name"))
model.all <- model.all %>%
  mutate(
    CRE = sub("_", ":", CRE)  # Replace underscores with colons
  )
model.all$CRE <- factor(model.all$CRE, levels = cre_levels)

model.maturation_stats <- model.all %>% filter(rank == 0) %>% filter(tile_name%in%target.events$tile_name) %>% distinct() %>% 
                    mutate(status = 'no TFBS')

for(CRE_oi in unique(model.maturation_stats$CRE)){
    cre_func_tfbs.count <- model.all %>% filter(rank > 0) %>% filter(CRE == CRE_oi)
    current.count <- model.maturation_stats %>% filter(CRE == CRE_oi)
    current.func_tfbs <- current.count %>% pull(n_func_tfbs)
    current.mature_func_tfbs <- current.count %>% pull(n_func_aff_tfbs)
    
    cre_update <- data.frame()
    for(rank_oi in seq(1:max(cre_func_tfbs.count$rank))){
        rank.count <- cre_func_tfbs.count %>% filter(rank == rank_oi) %>% distinct() %>% mutate(status = case_when(
                                    n_func_aff_tfbs > current.mature_func_tfbs ~ 'optimizing TFBS',
                                    n_func_tfbs > current.func_tfbs ~ 'recreating TFBS',
                                    func_tfbs == 'func' ~ 'recreating TFBS',
                                    TRUE ~ 'no TFBS'))
        if(nrow(rank.count) == 0){
        }else{
            ## reset counts to the next rank
            current.count <- rank.count 
            current.func_tfbs <- current.count %>% pull(n_func_tfbs)
            current.mature_func_tfbs <- current.count %>% pull(n_func_aff_tfbs)
            cre_update <- bind_rows(cre_update, rank.count)
        }
    }
    model.maturation_stats <- bind_rows(model.maturation_stats, cre_update)
}

random.all <- left_join(random_func_tfbs.count, all.events %>% dplyr::select(tile_name,frac_mean_delta,func_tfbs), by = c("tile_name"))
random.all <- random.all %>%
  mutate(
    CRE = sub("_", ":", CRE)  # Replace underscores with colons
  )
random.all$CRE <- factor(random.all$CRE, levels = cre_levels)

random.maturation_stats <- random.all %>%
  filter(rank == 0) %>% filter(iter == 1) %>%
  distinct() %>%
  filter(!grepl("ref", tile_name)) %>%  # Corrected grep to grepl
  mutate(status = 'no TFBS')

for(CRE_oi in unique(random.maturation_stats$CRE)){
    cre_func_tfbs.count <- random.all %>% filter(rank > 0) %>% filter(CRE == CRE_oi)
    current.count <- random.maturation_stats %>% filter(CRE == CRE_oi)
    current.func_tfbs <- current.count %>% pull(n_func_tfbs)
    current.mature_func_tfbs <- current.count %>% pull(n_func_aff_tfbs)
    cre_update <- data.frame()
    for(iter_oi in seq(1:10)){
        iter_func_tfbs.count <- cre_func_tfbs.count %>% filter(iter == iter_oi)
        iter_update <- data.frame()
        
        iter.count <- current.count
        iter.func_tfbs <- current.func_tfbs
        iter.mature_func_tfbs <- current.mature_func_tfbs
        
        for(rank_oi in seq(1:max(iter_func_tfbs.count$rank))){
            rank.count <- iter_func_tfbs.count %>% filter(rank == rank_oi) %>% distinct() %>% mutate(status = case_when(
                                        n_func_aff_tfbs > iter.mature_func_tfbs ~ 'optimizing TFBS',
                                        n_func_tfbs > iter.func_tfbs ~ 'recreating TFBS',
                                        func_tfbs == 'func' ~ 'recreating TFBS',
                                        TRUE ~ 'no TFBS'))
            if(nrow(rank.count) == 0){
            }else{
                ## reset counts to the next rank
                iter.count <- rank.count 
                iter.func_tfbs <- iter.count %>% pull(n_func_tfbs)
                iter.mature_func_tfbs <- iter.count %>% pull(n_func_aff_tfbs)
                iter_update <- bind_rows(iter_update, rank.count)
            }
        }
        cre_update <- bind_rows(cre_update, iter_update)
    }
    random.maturation_stats <- bind_rows(random.maturation_stats, cre_update)
}

maturation.stats <- bind_rows(model.maturation_stats,random.maturation_stats)
all.events <- all.events %>%
  mutate(
    CRE = sub("_", ":", CRE)  # Replace underscores with colons
  )
all.events$CRE <- factor(all.events$CRE, levels = cre_levels)

all.events2 <- left_join(all.events, maturation.stats %>% dplyr::select(tile_name, status)) %>% 
            dplyr::select(tile_name, rank, step, center_pos, mean_MPRA_act, sd_MPRA_act, mean_delta,sd_delta, frac_mean_delta, CRE, status, func_tfbs, point_type, iter) %>% distinct()
all.func <- all.func %>%
  mutate(
    CRE = sub("_", ":", CRE)  # Replace underscores with colons
  )
all.func$CRE <- factor(all.func$CRE, levels = cre_levels)

epas1_lollipop <- ggplot() +
#         geom_rect(data= gata_func_tfbs ,aes(xmin=TFBS_start,xmax=TFBS_end,fill=TF_name), ymin = -Inf, ymax = Inf,,alpha=0.3,color=NA)+
        geom_segment(data=all.func %>% filter(CRE == 'Epas1:chr17_10063'),
                     aes(x=TFBS_start,xend=TFBS_end,
                         y=TFBS_prop_WT,yend=TFBS_prop_WT,color=TF_name),
                     linewidth=1.5)+    
        #         scale_alpha_manual(values=c(0.1,0.3))+
        #         scale_fill_manual(values=col_TFs)+
        scale_colour_manual(values=col_TFs, name = "")+
        #         geom_point(data=mpra_result,
        #                    aes(x=pos,y=delta_MPRA), color = 'black', alpha = 0.2)+
#         guides(alpha="none")+
        #         new_scale_fill()+
        new_scale_color()+
        geom_hline(yintercept = 0, color = "black", linetype = 'dashed') +
        geom_point_rast(data=all.events2 %>% filter(CRE == 'Epas1:chr17_10063'),
                   aes(x=center_pos,y=frac_mean_delta, color = status, alpha = func_tfbs), shape = 21)+
        scale_color_manual(values = c("optimizing TFBS" = "#E69F00", "recreating TFBS" = "#56B4E9", 'no TFBS' = "black"), name = "") + 
#         scale_shape_manual(values = c("deletion" = 7, "mismatch" = 1, "insertion" = 10), name = "") + 
        scale_alpha_manual(values = c('not func' = 0.1, 'func' = 0.8), guide = 'none') + 
        coord_cartesian(xlim= c(0,300)) + 
        labs(y = expression("Fractional " * Delta * " MPRA activity"), x = "bp position") + theme_classic() + 
        theme(legend.position = 'bottom', legend.box="vertical", legend.margin=margin())

layout <- 'AABBBB'
epas1.plot <- epas1_mpra.pid + epas1_lollipop +
  plot_layout(design = layout)
ggsave("presentation_plots/epas1_example.pdf", plot = epas1.plot, height = 4, width = 12, dpi = 300, device = 'pdf')

# plt_lollipop <- ggplot() +
# #         geom_rect(data= gata_func_tfbs ,aes(xmin=TFBS_start,xmax=TFBS_end,fill=TF_name), ymin = -Inf, ymax = Inf,,alpha=0.3,color=NA)+
#         geom_segment(data=all.func,
#                      aes(x=TFBS_start,xend=TFBS_end,
#                          y=TFBS_prop_WT,yend=TFBS_prop_WT,color=TF_name),
#                      linewidth=1.5)+    
#         #         scale_alpha_manual(values=c(0.1,0.3))+
#         #         scale_fill_manual(values=col_TFs)+
#         scale_colour_manual(values=col_TFs, name = "")+
#         #         geom_point(data=mpra_result,
#         #                    aes(x=pos,y=delta_MPRA), color = 'black', alpha = 0.2)+
# #         guides(alpha="none")+
#         #         new_scale_fill()+
#         new_scale_color()+
#         geom_hline(yintercept = 0, color = "black", linetype = 'dashed') +
#         geom_point_rast(data=all.events2,
#                    aes(x=center_pos,y=frac_mean_delta, color = status, alpha = func_tfbs), shape = 21)+
#         scale_color_manual(values = c("optimizing TFBS" = "#E69F00", "recreating TFBS" = "#56B4E9", 'no TFBS' = "black"), name = "") + 
# #         scale_shape_manual(values = c("deletion" = 7, "mismatch" = 1, "insertion" = 10), name = "") + 
#         scale_alpha_manual(values = c('not func' = 0.1, 'func' = 0.8), guide = 'none') + 
#         coord_cartesian(xlim= c(0,300)) + 
#         facet_wrap(~CRE, ncol = 1) + 
#         labs(y = expression("Fractional " * Delta * " MPRA activity"), x = "bp position") + theme_classic() + 
#         theme(legend.position = 'bottom', legend.box="vertical", legend.margin=margin())

func_fc_tbl <- all.events2 %>% filter(!is.na(func_tfbs)) %>%  
  group_by(CRE, func_tfbs) %>%
  summarise(mean_delta = mean(frac_mean_delta), .groups = "drop") %>%
  group_by(CRE) %>%
  mutate(
    baseline = mean_delta[func_tfbs == "not func"],
    fold_change = (mean_delta) / baseline
  ) %>% filter(func_tfbs == 'func')

func_pval_tbl <- all.events2 %>% filter(!is.na(func_tfbs)) %>% 
  group_by(CRE) %>%
  summarise(
    p_val = wilcox.test(
      frac_mean_delta[func_tfbs == "func"],
      frac_mean_delta[func_tfbs == "not func"],
      exact = FALSE
    )$p.value
  ) %>%
  mutate(
    m_label = paste0(
      "p = ", signif(p_val, 3)
    )
  )


func_tfbs_delta <- ggplot() + 
        geom_quasirandom_rast(data = all.events2 %>% filter(!is.na(func_tfbs)), aes(x = func_tfbs, y = frac_mean_delta, color = func_tfbs), 
                              shape = 21) +
        scale_color_manual(values = c('not func' = '#414487', 'func' = '#22A884'), name = "",
                          labels = c("Func. TFBS", "Non-func. TFBS")) +  
        scale_x_discrete(drop = FALSE, labels = c("Func. TFBS","Non-func. TFBS"))  +
        geom_hline(yintercept = 0, color = "black", linetype = 'dashed') + 
        geom_text(data = func_pval_tbl, aes(label = m_label), color = 'black', x = Inf, y = -Inf, 
              hjust = 1, vjust = -1, check_overlap = T) +
        theme_classic() + facet_wrap(~CRE, ncol = 3) + 
        labs(y = expression(atop("Fractional ", Delta * " MPRA activity")), x = "") + 
        theme(legend.position = 'none', legend.box="vertical", legend.margin=margin())

### stacked bar plot
sum.events <- all.events2 %>%
              filter(!is.na(func_tfbs)) %>%
              mutate(iter = if_else(is.na(iter), "Model", paste0("Random ", iter))) %>%
              group_by(CRE, iter, func_tfbs) %>%
              summarise(sum_delta = sum(frac_mean_delta), .groups = "drop_last") %>%
              mutate(
                prct_delta = abs(sum_delta) / sum(abs(sum_delta))
              ) %>%
              ungroup()
events.total <- sum.events %>%
                group_by(CRE, func_tfbs) %>% summarise(mean_prct_delta = mean(prct_delta), sd_prct_delta = sd(prct_delta),
                                                      max_prct_delta = max(prct_delta), min_prct_delta = min(prct_delta)) #%>%
                #group_by(CRE) %>% summarise(sd_delta_events =sd(sum_delta))
sum.events$iter <- factor(sum.events$iter, levels = rev(c("Model", paste0('Random ', seq(1:10)))))
delta_barstack <- ggplot(data=sum.events, aes(x=prct_delta*100, y=iter, fill=func_tfbs)) +
          geom_bar(stat="identity") + facet_wrap(~CRE, ncol = 3) +
          labs(x = expression("% Fractional " * Delta * " MPRA contribution"), y = "") + 
          scale_fill_manual(values = c('not func' = '#414487', 'func' = '#22A884'), name = "",
                          labels = c("Func. TFBS", "Non-func. TFBS")) + theme_classic() +
          theme(legend.position = 'bottom', legend.box="vertical", legend.margin=margin())
### TFBS length
interval_union_length <- function(starts, ends) {
  ord <- order(starts)
  starts <- starts[ord]
  ends   <- ends[ord]

  total <- 0
  cur_start <- starts[1]
  cur_end   <- ends[1]

  for (i in 2:length(starts)) {
    if (starts[i] <= cur_end + 1) {
      # overlap or adjacency → extend
      cur_end <- max(cur_end, ends[i])
    } else {
      # no overlap → close previous interval
      total <- total + (cur_end - cur_start + 1)
      cur_start <- starts[i]
      cur_end   <- ends[i]
    }
  }

  # add final interval
  total + (cur_end - cur_start + 1)
}


tfbs_coverage <- mouse_func_tfbs %>%  
    mutate(
        CRE = sub("_", ":", CRE)  # Replace underscores with colons
      ) %>% filter(CRE%in%cre_levels) %>% 
  group_by(CRE) %>%
  summarise(
    func_length = interval_union_length(
      TFBS_start,
      TFBS_end
    ),
    .groups = "drop"
  ) %>% mutate(nonfunc_length = 300 - func_length)
tfbs_coverage <- melt(tfbs_coverage) 
colnames(tfbs_coverage) <- c("CRE","func_tfbs","length")
tfbs_coverage <- tfbs_coverage %>% group_by(CRE) %>% mutate(prct_length = length/300)
tfbs_coverage$CRE <- factor(tfbs_coverage$CRE, levels = cre_levels)

length_barstack <- ggplot(tfbs_coverage, aes(x = 1, y = prct_length*100, fill = func_tfbs)) +
  geom_col(width = 0.7) +
  facet_wrap(~CRE, ncol = 3) +
  coord_flip() +
  scale_y_continuous(limits = c(0, 100), breaks = c(0, 50, 100)) +
  scale_x_continuous(breaks = NULL) +
  scale_fill_manual(
    values = c(nonfunc_length = "#414487", func_length = "#22A884"),
    name = "",
    labels = c("Func. TFBS", "Non-func. TFBS")
  ) +
  labs(y = "% of CRE length", x = NULL) +
  theme_classic() +
  theme(legend.position = "none")
# layout <- "AAAAAAA
#            AAAAAAA
#            AAAAAAA
#            BBBBBBB"
# stackbars <- delta_barstack + length_barstack +
#   plot_layout(design = layout) 

# ggsave('test_bar.pdf', plot = stackbars, width = 210, height = 100,  useDingbats = F, units = 'mm', device = 'pdf')




# all.range <- all.range %>%
#   mutate(
#     CRE = sub("_", ":", CRE)  # Replace underscores with colons
#   )
# all.range$CRE <- factor(all.range$CRE, levels = cre_levels)

# range.stats <- all.range %>%
#   group_by(CRE) %>%
#   summarise(
#     p_value = wilcox.test(frac_mean_delta_range, frac_max_sigma)$p.value,
#     mean_range = mean(frac_mean_delta_range, na.rm = TRUE),
#     mean_sigma = mean(frac_max_sigma, na.rm = TRUE),
#     fold_change = mean_range / mean_sigma,
#     over = mean(frac_max_sigma > frac_mean_delta_range)
#   ) %>%
#   mutate(
#     m_label = paste0(
#       "p = ", signif(p_value, 3)
#     )
#   )


# gg_cdf <- ggplot() + 
#         stat_ecdf(data = all.range, aes(x = frac_mean_delta_range, color = "observed"), geom = 'step', linewidth = 1) +
#         stat_ecdf(data = all.range, aes(x = frac_max_sigma, color = "expected"), geom = 'step', linetype = "dashed", linewidth = 1) +
#         geom_text(data = range.stats, aes(label = m_label), color = 'black', x = Inf, y = -Inf, 
#               hjust = 1, vjust = -1.5, check_overlap = T) +
#         scale_color_manual(values = c('observed' = 'black','expected' = 'grey'), name = "") + facet_wrap(~CRE, ncol = 3) +
#         scale_x_log10(breaks = c(10^-2,10^-1, 10^0),
#                        labels = scales::trans_format("log10", scales::math_format(10^.x))) + annotation_logticks(side = 'b') +
#         theme_classic() +
#         labs(x = expression("Fractional " * Delta * " MPRA activity range"), y = "Cumulative distribution") +
#         guides(color = guide_legend(override.aes = list(linewidth = 0.5))) +  # Reduces linewidth in legend
#         theme(legend.position = c(0.07, 0.90),  # Adjust position inside the plot
#                                       legend.background = element_rect(fill = NA, color = NA), # Semi-transparent background
#                                       legend.key = element_blank()) 


### example plot
cre.events <- all_delta %>% filter(traj%in%traj_oi_list) %>%
                mutate(iter = case_when(
                        is.na(iter) ~ 'Model',
                        TRUE ~ paste0('Random ', iter)))
cre.tech_var <- cre.events %>% filter(target_mut_id != "") %>% group_by(CRE, iter, target_mut_id) %>% summarise(sd_delta = sd(delta)) %>%
                    mutate(variation = 'Variation across\ntransfection\nreplicates')
cre.traj_var <- cre.events %>% filter(target_mut_id != "") %>% group_by(CRE, target_mut_id) %>% summarise(sd_delta = sd(delta)) %>% 
            mutate(iter = 'Order-dependence', variation = 'Variation across\nmutational\norders')
cre.var <- bind_rows(cre.traj_var, cre.tech_var) 
cre.var$iter <- factor(cre.var$iter, levels = c("Order-dependence", "Model",paste0("Random ", seq(1:10))))
cre.var <- cre.var %>%
  mutate(
    CRE = sub("_", ":", CRE)  # Replace underscores with colons
  )
cre.var$CRE <- factor(cre.var$CRE, levels = cre_levels)

spread.stats <- cre.var %>%   # replace df with your data.frame name
  group_by(CRE) %>%
  summarise(
    mean_obs = mean(sd_delta[variation == "Variation across\nmutational\norders"], na.rm = TRUE),
    mean_exp = mean(sd_delta[variation == "Variation across\ntransfection\nreplicates"], na.rm = TRUE),
    sd_obs = mean(sd_delta[variation == "Variation across\nmutational\norders"], na.rm = TRUE),
    sd_exp = mean(sd_delta[variation == "Variation across\ntransfection\nreplicates"], na.rm = TRUE),
    fold_change = mean_obs / mean_exp,
    log2_fc = log2(fold_change),

    p_value = wilcox.test(
      sd_delta[variation == "Variation across\nmutational\norders"],
      sd_delta[variation == "Variation across\ntransfection\nreplicates"],
      exact = FALSE
    )$p.value,

    .groups = "drop"
  ) %>%
  mutate(
    m_label = paste0(
      "p = ", signif(p_value, 3)
    )
  )

cre_delta_spread <- ggplot() + 
        geom_quasirandom_rast(data = cre.var, aes(y = sd_delta, x = variation, color = variation), 
                              shape = 21, alpha = 0.7) +
        geom_boxplot(data = cre.var, aes(y = sd_delta, x = variation), fill = NA, color = 'black',
                     outlier.shape = NA) + ### remove outlier points since they will be plotted below
        scale_color_manual(values = c('Variation across\nmutational\norders' = 'black', 'Variation across\ntransfection\nreplicates' = 'grey'), name = "", labels = c("observed", "expected")) + 
        scale_y_log10(breaks = c(10^-2,10^-1, 10^0),
                       labels = scales::trans_format("log10", scales::math_format(10^.x))) + annotation_logticks(side = 'l') +
        geom_text(data = spread.stats, aes(label = m_label), color = 'black', x = -Inf, y = -Inf, 
              hjust = -0.2, vjust = -1.5, check_overlap = T) +
#         geom_hline(yintercept = 0, color = "black", linetype = 'dashed') + 
        theme_classic() + facet_wrap(~CRE, ncol = 3) + 
        labs(
              y = expression(atop("Standard deviation of", Delta * " MPRA activity")),
              x = ""
            ) + 
        theme(legend.position = 'none', legend.box="vertical", legend.margin=margin())

epas1.func <- all.func %>% filter(CRE == 'Epas1:chr17_10063')
epas1_lollipop <- ggplot() +
#         geom_rect(data= gata_func_tfbs ,aes(xmin=TFBS_start,xmax=TFBS_end,fill=TF_name), ymin = -Inf, ymax = Inf,,alpha=0.3,color=NA)+
        geom_segment(data=all.func %>% filter(CRE == 'Epas1:chr17_10063'),
                     aes(x=TFBS_start,xend=TFBS_end,
                         y=TFBS_prop_WT,yend=TFBS_prop_WT,color=TF_name),
                     linewidth=1.5)+    
        #         scale_alpha_manual(values=c(0.1,0.3))+
        #         scale_fill_manual(values=col_TFs)+
        scale_colour_manual(values=col_TFs, name = "")+
        #         geom_point(data=mpra_result,
        #                    aes(x=pos,y=delta_MPRA), color = 'black', alpha = 0.2)+
#         guides(alpha="none")+
        #         new_scale_fill()+
        new_scale_color()+
        geom_hline(yintercept = 0, color = "black", linetype = 'dashed') +
        geom_point_rast(data=epas1.range,
                   aes(x=center_pos,y=frac_mean_delta_range, color = 'observed'))+
        geom_point_rast(data=epas1.range,
                   aes(x=center_pos,y=frac_max_sigma, color = 'expected'))+
        scale_color_manual(values = c('observed' = 'black','expected' = 'grey'), name = "") + 
#         scale_shape_manual(values = c("deletion" = 7, "mismatch" = 1, "insertion" = 10), name = "") + 
        coord_cartesian(xlim= c(0,300)) + 
        labs(y = expression("Fractional " * Delta * " MPRA activity range"), x = "bp position") + theme_classic() + 
        theme(legend.position = 'bottom', legend.box="vertical", legend.margin=margin())


# gg_cdf2 <- ggplot() + 
#         stat_ecdf(data = all.range, aes(x = frac_mean_delta_range, color = "observed"), geom = 'step', linewidth = 1) +
#         stat_ecdf(data = all.range, aes(x = frac_max_sigma, color = "expected"), geom = 'step', linetype = "dashed", linewidth = 1) +
#         geom_text(data = range.stats, aes(label = m_label), color = 'black', x = Inf, y = -Inf, 
#               hjust = 1, vjust = -1, size = 5, check_overlap = T) +
#         scale_color_manual(values = c('observed' = 'black','expected' = 'grey'), name = "") + facet_wrap(~CRE, ncol = 3) +
#         scale_x_log10(breaks = c(10^-2,10^-1, 10^0),
#                        labels = scales::trans_format("log10", scales::math_format(10^.x))) + annotation_logticks(side = 'b') +
#         theme_classic() +
#         labs(x = expression("Fractional " * Delta * " MPRA activity range"), y = "Cumulative distribution") +
#         guides(color = guide_legend(override.aes = list(linewidth = 0.5))) +  # Reduces linewidth in legend
#         theme(text = element_text(family = "sans", size = 18), legend.position = c(0.07, 0.98),  # Adjust position inside the plot
#                                       legend.background = element_rect(fill = NA, color = NA), # Semi-transparent background
#                                       legend.key = element_blank()) 
# ggsave("presentation_plots/signal_noise_traj.pdf", plot = gg_cdf2, height = 4, width = 12, dpi = 300, device = 'pdf')


point.x <- left_join(mpra.points, rank_func_tfbs.count %>% select(tile_name, n_map_tfbs, n_map_aff_tfbs, n_func_tfbs, n_func_aff_tfbs), 
               by = c('tile_name')) %>% distinct()
update.points <- data.frame()
inital_random_func_tfbs.count <- data.frame()
rank_func_tfbs.count <- rank_func_tfbs.count %>%
  mutate(
    CRE = sub("_", ":", CRE)  # Replace underscores with colons
  )
rank_func_tfbs.count$CRE <- factor(rank_func_tfbs.count$CRE, levels = cre_levels)
random_func_tfbs.count <- random_func_tfbs.count %>%
  mutate(
    CRE = sub("_", ":", CRE)  # Replace underscores with colons
  )
random_func_tfbs.count$CRE <- factor(random_func_tfbs.count$CRE, levels = cre_levels)


for(cre_oi in unique(mpra.points$CRE)){
    cre.points <- point.x %>% filter(CRE == cre_oi)
    random.tmp <- random_func_tfbs.count %>% filter(CRE == cre_oi) %>% filter(rank == 0)
    end.points <- cre.points %>% filter(point_type == 'end')
    max.rank <- rank_func_tfbs.count %>% filter(CRE == cre_oi) %>% filter(rank == max(rank)) %>% pull(rank)
    max.tfbs <- rank_func_tfbs.count %>% filter(CRE == cre_oi) %>% filter(rank == max(rank)) %>% pull(n_func_tfbs)
    end.points$rank <- max.rank
    end.points$n_func_tfbs <- max.tfbs
    other.points <- cre.points %>% filter(point_type != 'end')
    inital.mpra <- cre.points %>% filter(point_type == 'start') %>% pull(mean_MPRA_act)
    update.points <- bind_rows(update.points, other.points)
    update.points <- bind_rows(update.points, end.points)
    random.tmp$mean_MPRA_act <- inital.mpra
    inital_random_func_tfbs.count <- bind_rows(inital_random_func_tfbs.count, random.tmp)
}

random_func_tfbs.count <- bind_rows(random_func_tfbs.count %>% filter(rank != 0),
                                   inital_random_func_tfbs.count)
update.points$CRE <- factor(update.points$CRE, levels = cre_levels)
update.points <- update.points %>% 
  arrange(factor(point_type, levels = c("start", "end", "best")))  # Ensure "best" is last

library(ggbeeswarm)

func_tfbs.count <- bind_rows(random_func_tfbs.count %>% mutate(point_type = 'Random'),
                            rank_func_tfbs.count %>% mutate(point_type = 'Model'))
func_tfbs.count <- func_tfbs.count %>% 
          arrange(factor(point_type, levels = c("Random", "Model")))  # Ensure "Model" is last

# Normalize n_aff_tfbs per CRE
func_tfbs.count <- func_tfbs.count %>% group_by(CRE) %>%
  mutate(n_func_aff_tfbs_scaled = pmin(n_func_aff_tfbs, 10),
        n_func_tfbs_scaled = pmin(n_func_tfbs, 10)) # Cap at 95th percentile
simple_func_tfbs.count <- func_tfbs.count %>% dplyr::select(CRE, tile_name, n_func_tfbs ,n_func_aff_tfbs, n_func_aff_tfbs_scaled, 
                                                            mean_MPRA_act, sd_MPRA_act, n_func_tfbs_scaled) %>% distinct()

func_tfbs_points <- ggplot() + 
        geom_quasirandom_rast(data = simple_func_tfbs.count, aes(x = n_func_aff_tfbs, y = mean_MPRA_act), 
                              color = 'grey', shape = 21, alpha = 0.5) +
#         geom_point(data = random_func_tfbs.count, aes(x = n_func_tfbs, y = mean_MPRA_act, color = 'Random'), alpha = 0.5, shape =21) + 
#         geom_point(data = rank_func_tfbs.count, aes(x = n_func_tfbs, y = mean_MPRA_act, color = 'Model'), shape =21) +
#         scale_color_manual(values = c('Model' = '#222f5b', 'Random' = '#FFA07A'), name = "") +  
#         scale_alpha_manual(values = c('Model' = 1, 'Random' = 0.5), name = "", guide = "none") + 
        geom_rect(data = select.ref, aes(ymin = mean_MPRA_act-sd_MPRA_act, ymax = mean_MPRA_act+sd_MPRA_act, color = NULL), xmin = -Inf, xmax = Inf,
                   alpha = .2, fill = "#008000", show.legend = FALSE) + 
        geom_hline(data = select.ref, aes(yintercept = mean_MPRA_act), color = '#008000', linetype = 'dashed', show.legend = FALSE) +
        geom_rect(data = select.target, aes(ymin = mean_MPRA_act-sd_MPRA_act, ymax = mean_MPRA_act+sd_MPRA_act, color = NULL), 
                  xmin = -Inf, xmax = Inf, show.legend = FALSE,alpha = .2,fill = "black") +
        geom_hline(data = select.target, aes(yintercept = mean_MPRA_act), color = "black", linetype = 'dashed', show.legend = FALSE) +
        geom_star(data = update.points, aes(x = n_func_aff_tfbs, y = mean_MPRA_act, fill = point_type, starshape = point_type), 
                  size = 3, show.legend = FALSE) +  
        scale_fill_manual(values = c('end' = '#222f5b', 'best' = '#ffd700', 'start' = '#D6E4FF')) + 
        scale_starshape_manual(values = c('start' = 15, 'best' = 1, 'end' = 15)) + 
        theme_classic() + facet_wrap(~ CRE, ncol = 3, scales = "free_x") + 
        scale_y_log10(breaks = c(10^-2,10^-1, 10^0, 10^1),
                   labels = scales::trans_format("log10", scales::math_format(10^.x)),
                          limits = c(0.01,10)) + annotation_logticks(sides = 'l') + 
        scale_x_continuous(breaks = scales::breaks_pretty(5)) +
        labs(x = "# of optimized functional TFBS", y = "MPRA activity") + 
        theme(panel.spacing = unit(1, "lines"),
             legend.position = 'bottom')

layout <- "
  AAAAAAAAAAA
  BBBBBBBBBBB
"
                                      
top.plot <-  mpra.pid + func_tfbs_points +
  plot_layout(design = layout) + plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(size = 12, face='bold'))
ggsave('fig4.pdf', plot = top.plot, width = 210, height = 160, useDingbats = F, units = 'mm', device = 'pdf')


fc_tbl <- maturation.stats %>% filter(!is.na(frac_mean_delta)) %>% 
  group_by(CRE, status) %>%
  summarise(mean_delta = mean(frac_mean_delta), .groups = "drop") %>%
  group_by(CRE) %>%
  mutate(
    baseline = mean_delta[status == "no TFBS"],
    fold_change = (mean_delta) / baseline
  )

pval_tbl <- maturation.stats %>% filter(!is.na(frac_mean_delta)) %>% 
  group_by(CRE) %>%
  summarise(
    p_recreating = wilcox.test(
      frac_mean_delta[status == "recreating TFBS"],
      frac_mean_delta[status == "no TFBS"],
      exact = FALSE
    )$p.value,
    p_optimizing = wilcox.test(
      frac_mean_delta[status == "optimizing TFBS"],
      frac_mean_delta[status == "no TFBS"],
      exact = FALSE
    )$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    m_label = paste0(
      "p optimizing = ", signif(p_optimizing, 3),
      "\np recreating = ", signif(p_recreating, 3)
    )
  )

maturation_cdf <- ggplot() + 
        stat_ecdf(data = maturation.stats, aes(x = frac_mean_delta, color = status, linetype = status), geom = 'step',
                 linewidth = 0.8) +
        scale_color_manual(values = c('optimizing TFBS' = '#E69F00',
                                     'recreating TFBS' = '#56B4E9',
                                     'no TFBS' = 'black'), name = "") +
        scale_linetype_manual(values = c('optimizing TFBS' = 'solid',
                                     'recreating TFBS' = 'dashed',
                                     'no TFBS' = 'dotted'), name = "") +
        geom_text(data = pval_tbl, aes(label = m_label), color = 'black', x = Inf, y = -Inf, 
              hjust = 1, vjust = -0.5, check_overlap = T) +
        theme_classic() + 
        facet_wrap(~CRE, ncol = 3) +
#         coord_cartesian(xlim=c(-0.3,1)) +
        scale_x_continuous(limit = c(-0.3, 1), breaks = c(0, 0.5, 1)) +
        labs(x = expression("Fractional " * Delta * " MPRA activity"), y = "Cumulative distribution") +
        theme(legend.position = 'bottom',  # Adjust position inside the plot
                                      legend.background = element_rect(fill = NA, color = NA), # Semi-transparent background
                                      legend.key = element_blank())

#### plot supp fig 8
layout <- 'AAAAAA
           BBBBBB'
sup8.plot <- atac.pid + corr.plot +
  plot_layout(design = layout) + plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(size = 12, face='bold'))
ggsave('supp_fig8.pdf', plot = sup8.plot, width = 210, height = 160, useDingbats = F, units = 'mm', device = 'pdf')

#### plot supp fig 9
layout <- 'AAAAAA
           AAAAAA
           BBBBBB
           CCCCCC
           CCCCCC
           CCCCCC
           DDDDDD
           DDDDDD
           EEEEEE
           EEEEEE'
sup9.plot <- func_tfbs_delta + length_barstack + delta_barstack + maturation_cdf + cre_delta_spread + 
  plot_layout(
    design = layout,
    heights = c(
      1, 1,        # A
      0.3,         # B  ← smaller
      1, 1, 1,     # C
      1, 1,        # D
      1, 1         # E
    )
  ) + plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(size = 12, face='bold'))
ggsave('supp_fig9.pdf', plot = sup9.plot, width = 210, height = 290, useDingbats = F, units = 'mm', device = 'pdf')