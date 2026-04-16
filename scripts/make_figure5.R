### module load pcre2/10.39; module load R/4.3.1

library(tidyverse)
library(ggplot2)
library(cowplot)
library(ggrastr)
library(castor)
library(ape)
library(phytools)
library(reshape2)
library(DescTools)
library(readxl)
library(gtools)
library(patchwork)
library(ggnewscale)
library(viridis)
library(ggtext)
library(ggbeeswarm)
library(ggstar)
source("utils/mpra_tablemaker.R")
source("utils/load_figure5_data.R")
source("utils/mafft_aligning_functions.R")
source("utils/figure5_func.R")
load("../extdata/CRE_traj_func_tfbs_info.RData")
load('../extdata/CRE_traj_delta_info.RData')


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

traj_mean_act <- traj_mean_act %>% distinct()

##### make figure scheme #####
CRE_oi <- 'Epas1_chr17_10063'                                     
example.max_traj <- traj_mean_act %>% filter(CRE == CRE_oi) %>% filter(gen > 0) %>% filter(selection == 'Maximizing') %>% distinct() %>% 
            arrange(gen) %>% head(n = 10)
example.min_traj <- traj_mean_act %>% filter(CRE == CRE_oi) %>% filter(gen > 0) %>% filter(selection == 'Minimizing') %>% distinct() %>% 
            arrange(gen) %>% head(n = 10)
                                     
example.max_mod <- data.frame()
example.min_mod <- data.frame()

for(gen_oi in seq(1:5)){
    tmp.max <- example.max_traj %>% filter(gen == gen_oi) %>% dplyr::select(gen, mapping)
    tmp.min <- example.min_traj %>% filter(gen == gen_oi) %>% dplyr::select(gen, mapping)
    
    tmp.max <- bind_rows(example.max_mod, tmp.max)
    tmp.min <- bind_rows(example.min_mod, tmp.min)
    
    tmp.max$order <- gen_oi
    example.max_mod <- bind_rows(example.max_mod, tmp.max)
    tmp.min$order <- gen_oi
    example.min_mod <- bind_rows(example.min_mod, tmp.min)
}
                
example.max_mod <- example.max_mod %>% distinct()
example.min_mod <- example.min_mod %>% distinct()

max.scheme <- ggplot(example.max_mod, aes(y = order))+
                geom_hline(aes(yintercept = order), color = 'black', linewidth = 1) + ## add horizontal line
                geom_star(aes(x = mapping + 1, y = order, fill = gen), 
                          size = 3, starshape = 12, alpha = 0.7, color = 'black', position = position_nudge(y = 0, x = 4)) + 
                scale_fill_gradient(low = "white",high = 'red') + 
                scale_y_reverse() + coord_cartesian(xlim=c(0,300), clip="off") + theme_classic() +
                labs(x = 'Enhancing\nObjective') + 
                theme(legend.position = 'none',
                          axis.text.y = element_blank(),
                          axis.line.y = element_blank(),
                          axis.ticks.y = element_blank(),
                          axis.title.y=element_blank(),
                          axis.text.x = element_blank(),
                          axis.line.x = element_blank(),
                          axis.ticks.x = element_blank()) 

min.scheme <- ggplot(example.min_mod, aes(y = order))+
                geom_hline(aes(yintercept = order), color = 'black', linewidth = 1) + ## add horizontal line
                geom_star(aes(x = mapping + 1, y = order, fill = gen), 
                          size = 3, starshape = 12, alpha = 0.7, color = 'black', position = position_nudge(y = 0, x = 4)) + 
                scale_fill_gradient(low = "white",high = 'blue') + 
                scale_y_reverse() + coord_cartesian(xlim=c(0,300), clip="off") + theme_classic() +
                labs(x = 'Ablation\nObjective') + 
                theme(legend.position = 'none',
                          axis.text.y = element_blank(),
                          axis.line.y = element_blank(),
                          axis.ticks.y = element_blank(),
                          axis.title.y=element_blank(),
                          axis.text.x = element_blank(),
                          axis.line.x = element_blank(),
                          axis.ticks.x = element_blank()) 

labels <- c()

for(species_oi in c("Mus_musculus")){
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
                                    
example.seq <- example.max_mod %>% head(n =1) 
example.seq$gen <- 'Mus_musculus'
             

seq.scheme <- ggplot(example.seq, aes(y = gen))+
                geom_hline(aes(yintercept = gen), color = 'black',linewidth = 1) + ## add horizontal line
                geom_star(aes(x =mapping), starshape = 12, size = 3, 
                                color = 'white', fill = 'black', alpha = 0.7) + 
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


max.point <- traj_mean_act %>% group_by(CRE) %>% filter(mean_MPRA_act == max(mean_MPRA_act))
min.point <- traj_mean_act %>%
  group_by(CRE) %>%
  filter(
    mean_MPRA_act <=
      control_mean_act %>%
        filter(CRE_id == "minP") %>%
        pull(mean_MPRA_act)
  ) %>%
  slice_min(gen, n = 1, with_ties = FALSE) %>%
  ungroup()

max.wt <- data.frame()
min.wt <- data.frame()

for(CRE_oi in unique(traj_mean_act$CRE)){
    wt_act <- wt_mean_act %>% filter(CRE == CRE_oi) %>% pull(mean_MPRA_act)
    max_cre_act <- traj_mean_act %>% filter(CRE == CRE_oi) %>% filter(mean_MPRA_act > wt_act) %>% filter(gen == min(gen))
    min_cre_act <- traj_mean_act %>% filter(CRE == CRE_oi) %>% 
            filter(mean_MPRA_act < control_mean_act %>% filter(CRE_id == 'minP') %>% pull(mean_MPRA_act)) %>% filter(gen == min(gen))
    max.wt <- bind_rows(max.wt, max_cre_act)
    min.wt <- bind_rows(min.wt, min_cre_act)
}

wt_mean_act <- traj_mean_act %>% filter(gen == 0) %>% dplyr::select(CRE, mean_MPRA_act, sd_MPRA_act, norm_footprint, pluri_pred, ecto_pred, meso_pred) %>% distinct()
imp.points <- bind_rows(max.point %>% mutate(point_type = 'Max'),
                       min.point %>%  mutate(point_type = 'Min'))
imp.points <- left_join(imp.points, wt_mean_act %>% dplyr::select(CRE, WT_MPRA_act = mean_MPRA_act), by = 'CRE') %>% 
                mutate(FC_WT = mean_MPRA_act/WT_MPRA_act, neg_FC_WT = WT_MPRA_act/mean_MPRA_act)

corr.df <- traj_mean_act[complete.cases(traj_mean_act$mean_MPRA_act),] %>% distinct() %>% 
    group_by(CRE, selection) %>% 
    summarise(R = cor(norm_footprint, mean_MPRA_act, method = 'spearman'))

max_corr.df <- traj_mean_act %>% filter(selection == 'Maximizing')
min_corr.df <- traj_mean_act %>% filter(selection == 'Minimizing')

link.max <- max_corr.df %>% filter(gen == 0 | gen == max(gen))
link.max <- bind_rows(link.max, max.point) %>% mutate(point_type = case_when(
                        gen == 0 ~ 'start',
                        gen == 50 ~ 'end',
                        TRUE ~ 'best'))
link.min <- min_corr.df %>% filter(gen == 0 | gen == max(gen))
link.min <- bind_rows(link.min, min.point) %>% mutate(point_type = case_when(
                        gen == 0 ~ 'start',
                        gen == 50 ~ 'end',
                        TRUE ~ 'best'))

link.min %>% filter(point_type == 'best') %>% summarise(m = median(gen))
# median 10 mutations needed for strongest ablation
link.max %>% filter(point_type == 'best') %>% summarise(m = median(gen))   
# median 7 mutations needed for strongest enhancement 

max_corr.df <- max_corr.df %>% 
    left_join(corr.df %>% filter(selection == 'Maximizing'), by = 'CRE') %>% 
    mutate(corr_label = paste0(CRE, "\nrho = ", round(R, 3)))  # use for labeling only

min_corr.df <- min_corr.df %>% 
    left_join(corr.df %>% filter(selection == 'Minimizing'), by = 'CRE') %>% 
    mutate(corr_label = paste0(CRE, "\nrho = ", round(R, 3)))  # use for labeling only

max_corr.plot <- ggplot() +
                # Arrows replacing lines for directionality
                geom_point_rast(data = max_corr.df, aes(y = mean_MPRA_act, x = norm_footprint, fill = gen), 
                                shape = 21, color = 'white', alpha = 0.8) +
                geom_line(data = link.max %>% filter(gen == 0 | gen != 50), aes(y = mean_MPRA_act, x = norm_footprint), 
                          linetype = 'dashed', linewidth = 0.8, alpha = 0.6, color = '#ee4035') +  # Connect points by gen order
                geom_line(data = link.max %>% filter(gen != 0 | gen == 50), aes(y = mean_MPRA_act, x = norm_footprint), 
                          linetype = 'dashed', linewidth = 0.8, alpha = 0.8, color = '#ee4035') +  # Connect points by gen order
                geom_vline(data = wt_mean_act, aes(xintercept = norm_footprint), color = '#008000', linetype = 'dashed', show.legend = FALSE) +
                geom_rect(data = wt_mean_act, aes(ymin = mean_MPRA_act-sd_MPRA_act, ymax = mean_MPRA_act+sd_MPRA_act, color = NULL), xmin = -Inf, xmax = Inf,
                   alpha = .4, fill = "#008000", show.legend = FALSE) + 
        geom_hline(data = wt_mean_act, aes(yintercept = mean_MPRA_act), color = '#008000', linetype = 'dashed', show.legend = FALSE) +
                geom_rect(data = control_mean_act %>% filter(CRE_id == 'minP'), 
                  aes(ymin = mean_MPRA_act-sd_MPRA_act, ymax = mean_MPRA_act+sd_MPRA_act), xmin = -Inf, xmax = Inf,
               alpha = .4,fill = "blue")  +
        geom_hline(yintercept = control_mean_act %>% filter(CRE_id == 'minP') %>% pull(mean_MPRA_act), color = 'blue', linetype = 'dashed') +
                geom_star(data = link.max, aes(x = norm_footprint, y = mean_MPRA_act, fill = gen, starshape = point_type), 
                  size = 3, show.legend = FALSE) +  
                scale_starshape_manual(values = c('start' = 15, 'best' = 1, 'end' = 15)) +
                scale_fill_gradient(breaks = c(0, 10, 30, 50),name = "# enhancing mutations", low = 'black', high = "#ee4035") + 
                labs(x = 'predicted log2(accessibility)', y = 'MPRA activity') +  
                facet_wrap(~CRE, ncol = 1, , labeller = as_labeller(setNames(max_corr.df$corr_label, max_corr.df$CRE))) + 
                scale_x_log10(breaks = c(10^0, 10^2, 10^4),
                     labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
                scale_y_log10(breaks = c(10^-1, 10^0, 10^1),
                   labels = scales::trans_format("log10", scales::math_format(10^.x)),
                          limits = c(0.05,15)) +
                annotation_logticks() + 
                theme_classic() + theme(legend.position = 'bottom',
                                           panel.spacing = unit(1, "lines"))

min_corr.plot <- ggplot() +
                # Arrows replacing lines for directionality
                geom_point_rast(data = min_corr.df, aes(y = mean_MPRA_act, x = norm_footprint, fill = gen), 
                                shape = 21, color = 'white', alpha = 0.8) +
                geom_line(data = link.min %>% filter(gen == 0 | gen != 50), aes(y = mean_MPRA_act, x = norm_footprint), 
                          linetype = 'dashed', linewidth = 0.8, alpha = 0.6, color = '#0392cf') +  # Connect points by gen order
                geom_line(data = link.min %>% filter(gen != 0 | gen == 50), aes(y = mean_MPRA_act, x = norm_footprint), 
                          linetype = 'dashed', linewidth = 0.8, alpha = 0.8, color = '#0392cf') +  # Connect points by gen order
                geom_vline(data = wt_mean_act, aes(xintercept = norm_footprint), color = '#008000', linetype = 'dashed', show.legend = FALSE) +
                geom_rect(data = wt_mean_act, aes(ymin = mean_MPRA_act-sd_MPRA_act, ymax = mean_MPRA_act+sd_MPRA_act, color = NULL), xmin = -Inf, xmax = Inf,
                   alpha = .4, fill = "#008000", show.legend = FALSE) + 
        geom_hline(data = wt_mean_act, aes(yintercept = mean_MPRA_act), color = '#008000', linetype = 'dashed', show.legend = FALSE) +
            geom_rect(data = control_mean_act %>% filter(CRE_id == 'minP'), 
                  aes(ymin = mean_MPRA_act-sd_MPRA_act, ymax = mean_MPRA_act+sd_MPRA_act), xmin = -Inf, xmax = Inf,
               alpha = .4,fill = "blue")  +
        geom_hline(yintercept = control_mean_act %>% filter(CRE_id == 'minP') %>% pull(mean_MPRA_act), color = 'blue', linetype = 'dashed') +
                geom_star(data = link.min, aes(x = norm_footprint, y = mean_MPRA_act, fill = gen, starshape = point_type), 
                  size = 3, show.legend = FALSE) +  
                scale_starshape_manual(values = c('start' = 15, 'best' = 1, 'end' = 15)) +
                scale_fill_gradient(breaks = c(0, 10, 30, 50),name = "# ablating mutations", low = 'black', high = "#0392cf") +           
                labs(x = 'predicted log2(accessibility)', y = 'MPRA activity') + 
                facet_wrap(~CRE, ncol = 1, labeller = as_labeller(setNames(min_corr.df$corr_label, min_corr.df$CRE))) + 
                scale_x_log10(breaks = c(10^-1,10^0, 10^1),
                     labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
                scale_y_log10(breaks = c(10^-1, 10^0, 10^1),
                   labels = scales::trans_format("log10", scales::math_format(10^.x)),
                          limits = c(0.05,15)) +
                annotation_logticks() + 
                theme_classic() + theme(legend.position = 'bottom',
                                           panel.spacing = unit(1, "lines"))



acc.plot <- ggplot() +
        geom_step(data = traj_mean_act, aes(x = gen, y = norm_footprint, color = selection)) +
        geom_step(data = traj_mean_act, aes(x = gen, y = pluri_pred, color = selection), linetype = 'dashed') +
        geom_hline(data = wt_mean_act, aes(yintercept = norm_footprint), color = '#008000', linetype = 'dashed', show.legend = FALSE) +
        geom_point(data = wt_mean_act, aes(y = norm_footprint), x = 0, color = 'black', size = 4) +
        scale_y_log10(breaks = c(10^0, 10^2, 10^4),
                     labels = scales::trans_format("log10", scales::math_format(10^.x))) + annotation_logticks(sides = 'l') + 
        coord_cartesian(xlim = c(0,50)) +
        facet_wrap(~CRE,ncol = 1) + theme_classic() + 
        scale_color_manual(
            values = selection.color,
            labels = c(
              "Maximizing" = "Enhancing",
              "Minimizing" = "Ablating"
            ),
            name = ""
          ) +  
        labs(x = '# of mutations', y = 'predicted log2(accessibility)') +
        theme(legend.position = 'bottom', panel.spacing = unit(1, "lines"))

act.plot <- ggplot() +
        geom_step(data = traj_mean_act, aes(x = gen, y = mean_MPRA_act, color = selection)) +
        geom_rect(data = wt_mean_act, aes(ymin = mean_MPRA_act-sd_MPRA_act, ymax = mean_MPRA_act+sd_MPRA_act, color = NULL), xmin = -Inf, xmax = Inf,
                   alpha = .4, fill = "#008000", show.legend = FALSE) + 
        geom_hline(data = wt_mean_act, aes(yintercept = mean_MPRA_act), color = '#008000', linetype = 'dashed', show.legend = FALSE) +
        geom_rect(data = control_mean_act %>% filter(CRE_id == 'minP'), 
                  aes(ymin = mean_MPRA_act-sd_MPRA_act, ymax = mean_MPRA_act+sd_MPRA_act), xmin = -Inf, xmax = Inf,
               alpha = .4,fill = "blue")  +
        geom_hline(yintercept = control_mean_act %>% filter(CRE_id == 'minP') %>% pull(mean_MPRA_act), color = 'blue', linetype = 'dashed') +
        geom_point(data = wt_mean_act, aes(y = mean_MPRA_act), x = 0, color = 'black', size = 4) +
        geom_errorbar(data = wt_mean_act, aes(ymin=mean_MPRA_act-sd_MPRA_act, ymax=mean_MPRA_act+sd_MPRA_act, x = 0), width=.2,
                         position=position_dodge(0.05)) +
        geom_star(data = imp.points, aes(x = gen, y = mean_MPRA_act, fill = point_type, starshape = point_type), 
                  size = 3, show.legend = FALSE, starshape = 1) +  
        scale_fill_manual(values = c('Max' = '#ee4035','Min' = '#0392cf')) + 
        scale_y_log10(breaks = c(10^-1, 10^0, 10^1),
                   labels = scales::trans_format("log10", scales::math_format(10^.x)),
                          limits = c(0.05,15)) + annotation_logticks(sides = 'l') + 
        coord_cartesian(xlim = c(0,50)) +
        facet_wrap(~CRE,ncol = 5) + theme_classic() + 
        scale_color_manual(
            values = selection.color,
            labels = c(
              "Maximizing" = "Enhancing",
              "Minimizing" = "Ablating"
            ),
            name = ""
          ) +  
        labs(x = '# of mutations', y = 'MPRA activity') +
        theme(legend.position = 'bottom', panel.spacing = unit(1, "lines"),)


match_traj_dms.df <- traj_mean_act
wt_DMS_MPRA <- read.delim("../data/mouse_WT_DMS_tile_CRE_act.txt", sep = '\t')
wt_DMS_MPRA <- wt_DMS_MPRA %>% 
    mutate(
    CRE = sub("_", ":", CRE)  # Replace underscores with colons
  )
mouse_DMS_MPRA <- mouse_DMS_MPRA %>% 
    mutate(
    CRE = sub("_", ":", CRE)  # Replace underscores with colons
  )

dms_MPRA <- left_join(mouse_DMS_MPRA, wt_DMS_MPRA %>% dplyr::select(CRE , WT_MPRA_act = mean_WT_act), by = c("CRE")) %>%
            mutate(DMS_delta = MPRA_act - WT_MPRA_act, DMS_prod = MPRA_act/WT_MPRA_act , pos = pos - 1) %>% 
            dplyr::select(CRE, mapping = pos, mut_nuc = mut, DMS_log2FC_MPRA = log2FC_MPRA , 
                          mean_DMS_MPRA_act = MPRA_act, sd_DMS_MPRA_act = sd_MPRA_act, DMS_delta, DMS_prod)
dms_MPRA[is.na(dms_MPRA)] <- 0

dms_MPRA$mut_nuc <- toupper(dms_MPRA$mut_nuc)

match_traj_dms.df <- left_join(match_traj_dms.df, dms_MPRA, by = c("CRE","mapping","mut_nuc"))
match_traj_gen0 <- left_join(traj_mean_act %>% filter(gen == 0),
                            wt_mean_act %>% dplyr::select(CRE, mean_DMS_MPRA_act = mean_MPRA_act, sd_DMS_MPRA_act = sd_MPRA_act),
                            by = c("CRE")) %>% mutate(DMS_delta = 0)
match_traj_dms.df <- bind_rows(match_traj_dms.df %>% filter(gen != 0), match_traj_gen0)

combinatorial_model.df <- left_join(match_traj_dms.df, wt_mean_act %>% dplyr::select(CRE , WT_MPRA_act = mean_MPRA_act), by = c("CRE")) %>%
                distinct() %>% arrange(gen) %>%
                filter(gen > 0) %>% 
                group_by(CRE, selection) %>% 
                mutate(delta_additive_model = cumsum(DMS_delta) + WT_MPRA_act, prod_multiplicative_model = cumprod(DMS_prod) * WT_MPRA_act) %>%
                ungroup() %>% dplyr::select(-WT_MPRA_act)

combinatorial_model_filtered.df <- combinatorial_model.df %>% filter(gen <= 10)
## add 0 point
wt_model <- data.frame(selection = c(rep('Maximizing', 5), rep('Minimizing', 5)),
                                                      wt_delta = rep(0,10), traj_delta = rep(0,10), DMS_delta = rep(0, 10), gen = rep(0,10),
                                                       func_tfbs = rep("not func", 10),
                                                      CRE = rep(cre_levels, 2))
wt_model <- left_join(wt_model, wt_mean_act %>% dplyr::select(CRE , delta_additive_model = mean_MPRA_act, mean_MPRA_act= mean_MPRA_act,
                                                              prod_multiplicative_model = mean_MPRA_act), by = c("CRE"))
combinatorial_model_filtered.df <- bind_rows(combinatorial_model_filtered.df,
                                            wt_model)


combinatorial_model_filtered.df$CRE <- factor(combinatorial_model_filtered.df$CRE, levels = cre_levels)

model.corr <- combinatorial_model_filtered.df %>% filter(selection == 'Maximizing') %>% group_by(CRE) %>% 
                summarise(rho_add = cor(mean_MPRA_act, delta_additive_model, method = 'spearman'),
                         rho_mult = cor(mean_MPRA_act, prod_multiplicative_model, method = 'spearman')) %>% ungroup() %>%
                mutate(label = paste0("rho[add] = ", round(rho_add, 2), "\n",
                        "rho[mult.] = ", round(rho_mult, 2)))
model.corr$CRE <- factor(model.corr$CRE, levels = cre_levels)

cum_additive.plot <- ggplot() + 
        geom_point_rast(data =combinatorial_model_filtered.df %>% filter(selection == 'Maximizing'), 
                        aes(x = gen, y = prod_multiplicative_model, color = "multiplicative", shape = func_tfbs), alpha = 0.7, size = 2) +
        geom_line(data =combinatorial_model_filtered.df %>% filter(selection == 'Maximizing'), 
                        aes(x = gen, y = prod_multiplicative_model, color = "multiplicative", linetype = "expected"), alpha = 0.7) +
        geom_point_rast(data =combinatorial_model_filtered.df %>% filter(selection == 'Maximizing'), 
                        aes(x = gen, y = delta_additive_model, color = "additive", shape = func_tfbs), alpha = 0.7, size = 2) +
        geom_line(data =combinatorial_model_filtered.df %>% filter(selection == 'Maximizing'), 
                        aes(x = gen, y = delta_additive_model, color = "additive", linetype = "expected"), alpha = 0.7) +
        geom_point_rast(data =combinatorial_model_filtered.df %>% filter(selection == 'Maximizing'), 
                        aes(x = gen, y = mean_MPRA_act, color = "enhancing", shape = func_tfbs), size = 2) +
        geom_errorbar(data =combinatorial_model_filtered.df %>% filter(selection == 'Maximizing'),
                      aes(x = gen, ymin=mean_MPRA_act-sd_MPRA_act, ymax=mean_MPRA_act+sd_MPRA_act, color = "enhancing"), width=.2,
                         position=position_dodge(0.05)) +
        geom_line(data =combinatorial_model_filtered.df %>% filter(selection == 'Maximizing'), 
                        aes(x = gen, y = mean_MPRA_act, color = "enhancing", linetype = "observed")) +
        scale_color_manual(values = c("additive" = '#8c510a', 'multiplicative' = '#666666', "enhancing" = "#ee4035"), 
                           name = "", guide = guide_legend(override.aes = list(shape = NA, size = 1))) +
        scale_shape_manual(values = c("func"=18 ,"not func" = 4), name = "",
                          guide = guide_legend(override.aes = list(size = 4))) +
        scale_linetype_manual(values = c("expected" = 'dashed', "observed" = 'solid'), guide = 'none') +
        geom_point_rast(data = wt_mean_act, x = 0, aes(y = mean_MPRA_act), color = 'black', size = 3, show.legend = FALSE) +
        geom_errorbar(data =wt_mean_act, color = 'black',
                      aes(x = 0, ymin=mean_MPRA_act-sd_MPRA_act, ymax=mean_MPRA_act+sd_MPRA_act), width=.2,
                         position=position_dodge(0.05)) +
        geom_text(data = model.corr, 
                  aes(label = label), color = 'black', x = Inf, y = -Inf, 
              hjust = 1,  # Align to the right
                        vjust = -0.2, size = 3, check_overlap = T) +
        scale_y_log10(breaks = c(10^-1, 10^0, 10^1, 10^2),
                   labels = scales::trans_format("log10", scales::math_format(10^.x)),
                          limits = c(0.05,100)) + annotation_logticks(sides = 'l') + 
        facet_wrap(~CRE, ncol = 5) +
        scale_x_continuous(breaks=c(0,2,4,6,8,10)) +
        labs(x = '# of mutations', 
             y = 'Cumulative activity') + theme_classic() +
        theme(legend.position = 'bottom', panel.spacing = unit(0.5, "lines"))


library(ggrepel)

#### calculate delta number of tfbs within a CRE

delta_tfbs_df <- traj_tfbs_info %>% distinct() %>%
  group_by(CRE, gen, selection) %>% summarise(n_map_tfbs = sum(n_map_tfbs), n_func_tfbs = sum(n_func_tfbs)) %>%
  arrange(gen, .by_group = TRUE) %>%
  group_by(CRE, selection) %>%
  mutate(
    delta_n_map_tfbs  = n_map_tfbs  - dplyr::lag(n_map_tfbs,  default = dplyr::first(as.numeric(n_map_tfbs))),
    delta_n_func_tfbs = n_func_tfbs - dplyr::lag(n_func_tfbs, default = dplyr::first(as.numeric(n_func_tfbs)))
  )

mouse_func_tfbs <- mouse_func_tfbs %>% 
    mutate(
    CRE = sub("_", ":", CRE)  # Replace underscores with colons
  )
mouse_func_tfbs$CRE <- factor(mouse_func_tfbs$CRE, levels = cre_levels)
new_delta_traj_model_act <- distinct(new_delta_traj_model_act)

mut.catalog <- new_delta_traj_model_act %>% filter(gen > 10) %>% group_by(CRE, selection, func_tfbs) %>%
            summarise(n = n()) %>% group_by(CRE, selection) %>% mutate(prct = n/sum(n))
mut.catalog$CRE <- factor(mut.catalog$CRE, levels = rev(cre_levels))

mut.stack <- ggplot(mut.catalog, 
                              aes(y = CRE, x = prct*100, fill = func_tfbs)) +
            geom_bar(stat="identity") + 
            facet_wrap(
                ~ selection,
                ncol = 1,
                labeller = as_labeller(c(
                  "Maximizing" = "Enhancing",
                  "Minimizing" = "Ablating"
                ))
          ) +
          labs(x = "% of early mutations overlapping TFBSs", y = "") + 
          scale_fill_manual(values = c('not func' = '#414487', 'func' = '#22A884'), name = "",
                          labels = c("Func. TFBS", "Non-func. TFBS")) + theme_classic() +
          theme(legend.position = 'bottom', legend.box="vertical", legend.margin=margin())


func_pval_tbl <- new_delta_traj_model_act %>% filter(gen <= 10) %>% 
  group_by(selection) %>%
  summarise(
    mean_func    = mean(log2FC[func_tfbs == "func"], na.rm = TRUE),
    mean_notfunc = mean(log2FC[func_tfbs == "not func"], na.rm = TRUE),
    fold_change  = mean_func/mean_notfunc,
    p_val = wilcox.test(
      log2FC[func_tfbs == "func"],
      log2FC[func_tfbs == "not func"],
      exact = FALSE
    )$p.value
  ) %>%
  mutate(
    m_label = paste0(
      "p = ", signif(p_val, 3)
    )
  )

func_pval_by_CRE_tbl <- new_delta_traj_model_act %>% filter(gen <= 10) %>% 
  group_by(CRE, selection) %>%
  summarise(
    mean_func    = mean(log2FC[func_tfbs == "func"], na.rm = TRUE),
    mean_notfunc = mean(log2FC[func_tfbs == "not func"], na.rm = TRUE),
    fold_change  = mean_func/mean_notfunc,
    p_val = wilcox.test(
      log2FC[func_tfbs == "func"],
      log2FC[func_tfbs == "not func"],
      exact = FALSE
    )$p.value
  ) %>%
  mutate(
    m_label = paste0(
      "p = ", signif(p_val, 3)
    )
  )

delta_mpra_dot.plot <- ggplot() + 
        geom_quasirandom_rast(data = new_delta_traj_model_act %>% filter(gen <= 10), 
                                             aes(y = log2FC, x = func_tfbs, color = func_tfbs), 
        dodge.width=1, shape = 21) +  # Adjust the width for separation 
            geom_hline(yintercept = 0, color = "black", linetype = 'dashed') +
        scale_color_manual(values = c('not func' = '#414487', 'func' = '#22A884'), name = "") +  
        geom_text(data = func_pval_tbl, aes(label = m_label), color = 'black', x = Inf, y = -Inf, 
              hjust = 1, vjust = -1, check_overlap = T) +
            theme_classic() + scale_x_discrete(drop = FALSE, labels = c("Func. TFBS","Non-func. TFBS"))  +
                        labs(y = expression(atop(log[2] ~ "FC (MPRA)","relative to prior mutation")), x = '') + 
            facet_wrap(
            ~ selection,
            ncol = 1,
            labeller = as_labeller(c(
              "Maximizing" = "Enhancing",
              "Minimizing" = "Ablating"
            ))
              ) +
            theme(legend.position = 'none')

layout <- "
  AAAAAAAA
"
                                      
fig5.plot <-  act.plot + 
  plot_layout(design = layout) + plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(size = 12, face='bold'))
ggsave('fig5.pdf', plot = fig5.plot, width = 210, height = 70, useDingbats = F, units = 'mm', device = 'pdf') 

layout <- 'AABBCC
           AABBCC
           AABBCC
           AABBCC
           '

sup10.plot <- acc.plot + max_corr.plot + min_corr.plot + plot_layout(design = layout) + plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(size = 12, face='bold'))
ggsave('supp_fig10.pdf', plot = sup10.plot, width = 210, height = 230, useDingbats = F, units = 'mm', device = 'pdf')


layout <- 'AAAABBBB
           AAAABBBB
           AAAABBBB
           AAAABBBB
           CCCCCCCC
           CCCCCCCC'

sup11.plot <- mut.stack + delta_mpra_dot.plot + cum_additive.plot + plot_layout(design = layout) + plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(size = 12, face='bold'))
ggsave('supp_fig11.pdf', plot = sup11.plot, width = 210, height = 160, useDingbats = F, units = 'mm', device = 'pdf')
                        
