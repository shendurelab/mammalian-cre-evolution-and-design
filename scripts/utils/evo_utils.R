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
library(viridis)
source("/net/shendure/vol8/projects/tli/PYS2_oCRE_BC_libs_20240904/mafft_aligning_functions.R")
set.seed(42)

shape.events <- c(3,4,16,15)
names(shape.events) <- c('insertion','deletion','target','mismatch')

calc_step_changes <- function(data, target.events){
    fc_df <- data.frame()
    for(traj_oi in unique(data$traj)){
        mean_fc <- c()
        max_fc <- c()
        min_fc <- c()
        mean_delta <- c()
        max_delta <- c()
        min_delta <- c()
        sigma_delta <- c()
        sd_delta <- c()
        if('rank'%in%colnames(data)){
            traj_df <- data %>% filter(traj == traj_oi) %>% filter((rank != 0 & !startsWith(Event, 'target'))) %>% arrange(rank)
            mean_traj_act <- traj_df %>% pull(mean_MPRA_act)
            max_traj_act <- traj_df %>% pull(max_MPRA_act)
            min_traj_act <- traj_df %>% pull(min_MPRA_act)
            sigma_traj_act <- traj_df %>% pull(sigma_MPRA_act)
            sd_traj_act <- traj_df %>% pull(sd_MPRA_act)
            CRE_oi <- unique(traj_df$CRE)[1]
            target_oi <- unique(traj_df$target)[1]
            start_mean_act <- target.events %>% filter(CRE == CRE_oi) %>% filter(target == target_oi) %>% head(n=1) %>% pull(mean_MPRA_act)
            start_max_act <- target.events %>% filter(CRE == CRE_oi) %>% filter(target == target_oi) %>% head(n=1) %>% pull(max_MPRA_act)
            start_min_act <- target.events %>% filter(CRE == CRE_oi) %>% filter(target == target_oi) %>% head(n=1) %>% pull(min_MPRA_act)
            start_sigma_act <- target.events %>% filter(CRE == CRE_oi) %>% filter(target == target_oi) %>% head(n=1) %>% pull(sigma_MPRA_act)
            start_sd_act <- target.events %>% filter(CRE == CRE_oi) %>% filter(target == target_oi) %>% head(n=1) %>% pull(sd_MPRA_act)
            
            
            for(i in 1:length(mean_traj_act)){
                if(i == 1){
                    mean_fc <- c(mean_fc, mean_traj_act[i]/start_mean_act)
                    max_fc <- c(max_fc, max_traj_act[i]/start_max_act)
                    min_fc <- c(min_fc, min_traj_act[i]/start_min_act)
                    
                    mean_delta <- c(mean_delta, mean_traj_act[i] - start_mean_act)
                    max_delta <- c(max_delta, max_traj_act[i] - start_max_act)
                    min_delta <- c(min_delta, min_traj_act[i] - start_min_act)
                    sigma_delta <- c(sigma_delta, sigma_traj_act[i] - start_sigma_act)
                    sd_delta <- c(sd_delta, sd_traj_act[i] - start_sd_act)
                }else{
                    mean_fc <- c(mean_fc, mean_traj_act[i]/mean_traj_act[i-1])
                    max_fc <- c(max_fc, max_traj_act[i]/max_traj_act[i-1])
                    min_fc <- c(min_fc, min_traj_act[i]/min_traj_act[i-1])
                    
                    mean_delta <- c(mean_delta, mean_traj_act[i] - mean_traj_act[i-1])
                    max_delta <- c(max_delta, max_traj_act[i] - max_traj_act[i-1])
                    min_delta <- c(min_delta, min_traj_act[i] - min_traj_act[i-1])
                    sigma_delta <- c(sigma_delta, sigma_traj_act[i] - sigma_traj_act[i-1])
                    sd_delta <- c(sd_delta, sd_traj_act[i] - sd_traj_act[i-1])
                }
            }
            traj_df$mean_fc <- mean_fc
            traj_df$max_fc <- max_fc
            traj_df$min_fc <- min_fc
            
            traj_df$mean_delta <- mean_delta
            traj_df$max_delta <- max_delta
            traj_df$min_delta <- min_delta
            traj_df$sigma_delta <- sigma_delta
            traj_df$sd_delta <- sd_delta
            fc_df <- bind_rows(fc_df, traj_df)
        }else{
        }
    }
    fc_df$log2_mean_fc <- log2(fc_df$mean_fc)
    fc_df$log2_mean_delta <- log2(fc_df$mean_delta)
    return(fc_df)
}

plot_group_mpra_traj <- function(evo_model_df, evo_random_df, traj_oi, target_oi, ref.events, target.events, ncol = 1, type = 'model'){
    evo_model_df <- evo_model_df %>% filter(traj%in% traj_oi) %>% filter((rank != 0 & !startsWith(Event, 'target')))
    evo_random_df <- evo_random_df %>% filter(traj%in% traj_oi) %>% filter((rank != 0 & !startsWith(Event, 'target'))) #%>% 
#                             group_by(rank, traj, CRE) %>% summarise(sd_MPRA_act = sd(mean_MPRA_act),
#                                                                                 mean_MPRA_act = mean(mean_MPRA_act))
    reference_oi <- unique(evo_model_df$reference)
    CRE_oi <- unique(evo_model_df$CRE)
    select.ref <- ref.events %>% filter(traj%in% traj_oi)
    select.target <- target.events %>% filter(traj%in% traj_oi) 
    
    mpra.plot <- ggplot() + 
                    geom_line(data = evo_random_df, aes(x = rank, y = mean_MPRA_act, group = iter), color = 'tomato', alpha = 0.3) + 
                    geom_point(data = evo_random_df, aes(x = rank, y = mean_MPRA_act, group = iter), shape = 16, color = 'tomato', alpha = 0.3) +  
                    geom_line(data = evo_model_df, aes(x = rank, y = mean_MPRA_act, color = rank)) + 
                    geom_point(data = evo_model_df, aes(x = rank, y = mean_MPRA_act, color = rank, shape = Event)) +
                    geom_errorbar(data = evo_model_df, aes(x = rank, ymin=mean_MPRA_act-sd_MPRA_act, 
                                                           ymax=mean_MPRA_act+sd_MPRA_act, color = rank), width=.2, position=position_dodge(0.05)) + 
                    scale_shape_manual(values = shape.events, guide = guide_legend(override.aes = list(size = 3, fill = 'black'))) +
                    scale_color_gradient(limits = c(0,50), oob = scales::squish, name = "model-ranked mutations",
                                        breaks = c(0, 25, 50), labels = c(0, 25, '50+')) +
                    geom_rect(data = select.ref, aes(ymin = mean_MPRA_act-sd_MPRA_act, ymax = mean_MPRA_act+sd_MPRA_act, color = NULL), xmin = -Inf, xmax = Inf,
                   alpha = .2, fill = "#008000") + 
                    geom_hline(data = select.ref, aes(yintercept = mean_MPRA_act), color = '#008000', linetype = 'dashed') +
                    geom_rect(data = select.target, aes(ymin = mean_MPRA_act-sd_MPRA_act, ymax = mean_MPRA_act+sd_MPRA_act, color = NULL), xmin = -Inf, xmax = Inf,
                   alpha = .2,fill = "black") + scale_x_continuous(breaks = c(0, 30, 60, 90, 120), limits = c(0,120)) + 
                    guides(shape = guide_legend(title.position = "top", title.hjust=0.5)) + 
                    geom_hline(data = select.target, aes(yintercept = mean_MPRA_act), color = "black", linetype = 'dashed') +
                    scale_y_log10(
                       breaks = c(10^-1, 10^0, 10^1),
                       labels = scales::trans_format("log10", scales::math_format(10^.x))
                     ) + annotation_logticks(sides = 'l') + labs(x = '# sequential mutations', y = 'MPRA act.') +
                    facet_wrap(~CRE, ncol = ncol) +
                    theme_cowplot() + theme(legend.position = 'bottom', legend.box = "horizontal", text = element_text(family = "sans", size = 16))
#     if(type == 'random_avg'){
#         select_model <- select_model %>% group_by(rank, traj, CRE) %>% summarise(sd_MPRA_act = sd(mean_MPRA_act),
#                                                                                 mean_MPRA_act = mean(mean_MPRA_act))
#         mpra.plot <- ggplot(select_model, aes(x = rank, y = mean_MPRA_act)) + 
#                     geom_line(color = 'tomato') + geom_point(color = 'tomato') + 
#                     #geom_errorbar(aes(ymin=mean_MPRA_act-sd_MPRA_act, ymax=mean_MPRA_act+sd_MPRA_act), width=.2,
#                     #     position=position_dodge(0.05), color = 'tomato') + 
#                     geom_rect(data = select.ref, aes(ymin = mean_MPRA_act-sd_MPRA_act, ymax = mean_MPRA_act+sd_MPRA_act, color = NULL), xmin = -Inf, xmax = Inf,
#                    alpha = .2, fill = "#008000") + 
#                     geom_hline(data = select.ref, aes(yintercept = mean_MPRA_act), color = '#008000', linetype = 'dashed') +
#                     geom_rect(data = select.target, aes(ymin = mean_MPRA_act-sd_MPRA_act, ymax = mean_MPRA_act+sd_MPRA_act, color = NULL), xmin = -Inf, xmax = Inf,
#                    alpha = .2,fill = "black") + scale_x_continuous(breaks = c(0, 30, 60, 90, 120), limits = c(0,120)) + 
#                     geom_hline(data = select.target, aes(yintercept = mean_MPRA_act), color = "black", linetype = 'dashed') +
#                     scale_y_log10(
#                        breaks = c(10^-1, 10^0, 10^1),
#                        labels = scales::trans_format("log10", scales::math_format(10^.x))
#                      ) + annotation_logticks(sides = 'l') + labs(x = '# random iterative mutations', y = 'MPRA act.') +
#                     facet_wrap(~CRE, ncol = ncol) +
#                     theme_cowplot() + theme(legend.position = 'bottom')
        
#     } 
#     else{
    
#         mpra.plot <- ggplot(select_model, aes(x = rank, y = mean_MPRA_act, shape = Event, group = traj, color = rank)) + 
#                     geom_line() + geom_point() + geom_errorbar(aes(ymin=mean_MPRA_act-sd_MPRA_act, ymax=mean_MPRA_act+sd_MPRA_act), width=.2,
#                          position=position_dodge(0.05)) + 
#                     geom_rect(data = select.ref, aes(ymin = mean_MPRA_act-sd_MPRA_act, ymax = mean_MPRA_act+sd_MPRA_act, color = NULL), xmin = -Inf, xmax = Inf,
#                    alpha = .2, fill = "#008000") + 
#                     geom_hline(data = select.ref, aes(yintercept = mean_MPRA_act), color = '#008000', linetype = 'dashed') +
#                     geom_rect(data = select.target, aes(ymin = mean_MPRA_act-sd_MPRA_act, ymax = mean_MPRA_act+sd_MPRA_act, color = NULL), xmin = -Inf, xmax = Inf,
#                    alpha = .2,fill = "black") + 
#                     geom_hline(data = select.target, aes(yintercept = mean_MPRA_act), color = "black", linetype = 'dashed') +
#                     ylim(c(0,10)) + 
#                     scale_y_log10(
#                        breaks = c(10^-1, 10^0, 10^1),
#                        labels = scales::trans_format("log10", scales::math_format(10^.x))
#                      ) + annotation_logticks(sides = 'l') + labs(x = '# model-selected iterative mutations', y = 'MPRA act.') +
#                     facet_wrap(~CRE, ncol = ncol) + scale_x_continuous(breaks = c(0, 30, 60, 90, 120), limits = c(0,120)) + 
#                     scale_color_gradient(limits = c(0,50), oob = scales::squish, name = "",
#                                         breaks = c(0, 25, 50),
#                        labels = c(0, 25, '50+')) +
#                     scale_shape_manual(values = shape.events, guide = guide_legend(override.aes = list(size = 3, fill = 'black'))) + 
#                     theme_cowplot() + theme(legend.position = 'bottom')
#     }
    return(mpra.plot)
}

plot_group_acc_traj <- function(evo_df, traj_oi, target_oi, ref.events, target.events){
    select_model <- evo_df %>% filter(traj%in% traj_oi) %>% filter((rank != 0 & !startsWith(Event, 'target')))
    reference_oi <- unique(select_model$reference)
    CRE_oi <- unique(select_model$CRE)
    select.ref <- ref.events %>% filter(reference%in%reference_oi) %>% filter(CRE%in%CRE_oi) %>% filter(Event == 'ref1')
    select.target <- target.events %>% filter(target == target_oi) %>% filter(CRE == CRE_oi) %>% filter(Event == 'target1')

    pred.plot <- ggplot(select_model, aes(x = rank, y = log2_footprint, shape = Event, group = traj)) + 
                geom_line() + geom_point() + 
                geom_hline(data = select.ref, aes(yintercept = log2_footprint), color = '#008000', linetype = 'dashed') +
                geom_hline(data = select.target, aes(yintercept = log2_footprint), color = "black", linetype = 'dashed') +
                labs(x = '# mutations', y = 'log2(chrombpnet pred.)') + 
                ggtitle('Predicted chromatin accessibility') + facet_wrap(~CRE, ncol = 1) + 
                scale_shape_manual(values = shape.events) + theme_cowplot() + theme(legend.position = 'bottom', 
                                                                                    text = element_text(family = "sans",size = 16))
    
    return(pred.plot)
}

get_besk_rank <- function(evo_act, act_cutoff){
    best_rank <- evo_act %>% filter(mean_MPRA_act >= act_cutoff) %>%
                arrange(rank) %>% head(n = 1) %>% pull(rank)
    return(best_rank)
}

check_tfbs_addition_order <- function(evo_act, traj_oi, act_cutoff, tfbs_map, CRE_oi, diff = 1){
    traj_act <- evo_act %>% filter(traj == traj_oi)
    best_rank <- get_besk_rank(traj_act, act_cutoff)
    print(paste0("Rank threshold above act. cutoff = ", best_rank))
    best_delta_rank <- traj_act %>% filter(rank <= best_rank) %>% arrange(-rank_delta) %>% head(n = 1) %>% pull(rank)
    print(paste0("Rank w. best delta = ", best_delta_rank))
    range_check <- traj_act %>% filter(rank >= best_delta_rank - diff & rank <= best_rank + diff) %>%
                select(tile_name, rank, mean_MPRA_act, sd_MPRA_act, rank_delta) %>% arrange(rank) %>% 
                mutate(rank_label = case_when(rank %in% c(11,12,13) ~ paste0(rank, "th m."),
                        rank %% 10 == 1 ~ paste0(rank, 'st m.'),
                        rank %% 10 == 2 ~ paste0(rank, 'nd m.'),
                        rank %% 10 == 3 ~ paste0(rank, 'rd m.'),
                        TRUE ~ paste0(rank, "th m.")))
    range_label_order <- range_check %>% pull(rank_label)
    range_check$species <- range_check$rank_label
    range_check$species <- factor(range_check$species, levels = range_label_order)
    mpra_bar <- ggplot(range_check, aes(x=mean_MPRA_act, y=species, color=mean_MPRA_act)) + 
        geom_point(size = 2) + 
        geom_errorbar(aes(xmin=mean_MPRA_act-sd_MPRA_act, xmax=mean_MPRA_act+sd_MPRA_act), width=.2,
                     position=position_dodge(0.05)) +
#         geom_rect(xmin = cre_oi_wt_act-cre_oi_wt_sd, xmax = cre_oi_wt_act+cre_oi_wt_sd, ymin = -Inf, ymax = Inf,
#            alpha = .2,fill = "#008000") + 
#             geom_rect(xmin = minP_oi_act - minP_oi_sd, xmax = minP_oi_act + minP_oi_sd, ymin = -Inf, ymax = Inf,
#            alpha = .1,fill = "blue")  +
        geom_vline(xintercept = act_cutoff ,color = "#008000", linetype = 'dashed') +
        scale_color_gradient(breaks = c(10^-1, 10^0, 10^1),
               labels = scales::trans_format("log10", scales::math_format(10^.x)), limits = c(0.05, 15),  # Quantile limits from your range
                                    oob = scales::squish,  # To handle values outside the specified limits
                                    trans = "log10",
                                   low = '#42047e',high = '#07f49e', na.value = 'gainsboro') +
        scale_x_log10(breaks = c(10^-1, 10^0, 10^1),
               labels = scales::trans_format("log10", scales::math_format(10^.x)),
                      limits = c(0.05,15))+ annotation_logticks(sides = 'b')+ theme_classic() + 
            labs(x = 'MPRA act.', color = 'MPRA act.') + 
            theme(text = element_text(size = 16),
                          axis.text.y = element_blank(),
                          axis.line.y = element_blank(),
                          axis.ticks.y = element_blank(),
                          axis.title.y=element_blank()) 
    tfbs_map <- tfbs_map %>% filter(tile_name%in%range_check$tile_name) %>% 
                            mutate(center_pos = (start_pos + end_pos)/2)
    mm <- match(tfbs_map$tile_name, range_check$tile_name)
    tfbs_map$species <- range_check$rank_label[mm]
    tfbs_map$rank <- range_check$rank[mm]
    tfbs_map$species <- factor(tfbs_map$species, levels = range_label_order)
    max_width = max(tfbs_map$width)
    best_tfbs <- tfbs_map %>% filter(rank == best_delta_rank) %>% select(center_pos, TF_name, species, norm_affinity, width)
    best_tfbs_minus1 <- tfbs_map %>% filter(rank == best_delta_rank - 1) %>% select(center_pos, TF_name, species, norm_affinity, width)
    
    tolerance_width <- abs(max(best_tfbs_minus1 %>% pull(width)) - max(best_tfbs %>% pull(width)))

    non_overlapping_rows <- best_tfbs %>%
                      rowwise() %>% 
                      filter(!any(abs(center_pos - best_tfbs_minus1$center_pos) <= tolerance_width &
                                    TF_name == best_tfbs_minus1$TF_name &
                                  abs(norm_affinity - best_tfbs_minus1$norm_affinity) <= 0.1))
    print(non_overlapping_rows)
    
    tfbs_order_check <- make_tfbs_plot(tfbs_map, max_width, marker = non_overlapping_rows) + ggtitle(CRE_oi)
    layout <- '
        AAAAAAABB
        AAAAAAABB'
    p <- tfbs_order_check + mpra_bar + plot_layout(design = layout)    
    return(p)
    
}

fix_tfbs_overlaps <- function(df, buffer = 2) {
  select.df <- df %>% select(center_pos, species, TF_name, norm_affinity)
  for (i in 1:nrow(select.df)) {
      for (j in 1:nrow(select.df)) {
        if (i != j && abs(select.df$center_pos[i] - select.df$center_pos[j]) <= 1) {
          # If there is an overlap, modify the center_pos by adding/subtracting 1 to avoid overlap
          select.df$center_pos[i] <- select.df$center_pos[i] + buffer
        }
      }
    }
  return(select.df)
}


make_tfbs_plot <- function(tfbs_map, max_width, labels = NULL, marker = NULL){
    
    final_TFs <- c("Jun:atf3","Foxa2","Gata6","Sox17","Klf4")
    col_TFs <- c("#DC3220","#009E73","#FFC20A","#56B4E9","#5D3A9B")
    names(col_TFs) <- final_TFs 
    shaps_TFs <- c(29,13,15,28,11)
    names(shaps_TFs) <- final_TFs
    
#     fix_tfbs_map <- fix_tfbs_overlaps(tfbs_map)
    
    gg_tfbs <- ggplot(tfbs_map, aes(y=species)) +
                  geom_hline(aes(yintercept = species), color = "black", alpha = 0.5) +  # Add horizontal line
                  geom_star(aes(x =center_pos, starshape=TF_name, fill = TF_name, alpha = norm_affinity), size = 3,
                                    color = 'black', position=position_jitter(h=0.2)) + 
                  labs(x = 'bp position', fill = "TFBS", shape = "TFBS", size = "affinity binding") +
                  scale_fill_manual(values = col_TFs) + 
                            theme_cowplot() + 
                            scale_starshape_manual(values = shaps_TFs, name = "TFBS",
                                                  guide = guide_legend(override.aes = list(size = 3))) + # coord_cartesian(xlim = c(0, 300)) + 
                  scale_alpha_continuous(name = "TFBS affinity", limits = c(0.1, 1), breaks = c(0.25, 0.5, 0.75, 1), oob=scales::squish) + # Affinity binding size scaling
                    guides(size = guide_legend(override.aes = list(shape = 21, fill = "black"))) + 
                    coord_cartesian(xlim = c(0, max_width)) + 
                    labs(x = 'bp position') + 
                        theme(legend.position = 'bottom',legend.box="vertical", legend.margin=margin(),
                              axis.text.y = element_markdown(color = "black", size = 10),
                              axis.line.y = element_blank(),
                              axis.ticks.y = element_blank(),
                              axis.title.y=element_blank(), 
                              plot.title = element_text(hjust = 0.5, face = 'plain'),
                             text = element_text(family = "sans", size = 16))
    if(!is.null(labels)){
        gg_tfbs <- gg_tfbs + scale_y_discrete(labels = labels) 

    }
    if(!is.null(marker)){
        gg_tfbs <- gg_tfbs + geom_point(data = marker, aes(y = as.numeric(species) + 0.5, x = center_pos), shape = 25, 
                                        size = 5, color = 'black', fill = 'red', position=position_jitter(h=0.1))
    }
    return(gg_tfbs)
}


map_TFBS_to_msa_pos <- function(tfbs.map, msa.map){
    mm.start <- match(tfbs.map$start_pos, msa.map$query_pos)
    mm.end <- match(tfbs.map$end_pos, msa.map$query_pos)
    
    tfbs.map$msa_start <- msa.map$curpos[mm.start]
    tfbs.map$msa_end <- msa.map$curpos[mm.end]
    return(tfbs.map)
}

plot_traj_tfbs_lfc <- function(evo_df, msa_map, tfbs_map, ref_tile, target_tile, final_tile, 
                               ref_name, target_name, traj_oi, CRE_oi){
    select_model <- evo_df %>% filter(traj == traj_oi) %>% filter((rank != 0 & !startsWith(Event, 'target')))
    pos.mut <- select_model %>% filter(rank_delta > 0) %>% select(ref_start, ref_end)
    pos.mut <- IRanges(start = pos.mut$ref_start, end = pos.mut$ref_end)
    
    ref_start <- select_model %>% select(ref_start,Event, rank, log2_rank_fc, rank_delta, phyloP)
    colnames(ref_start) <- c("curpos",'Event','rank','fc','delta','phyloP')
    ref_end <- select_model %>% select(ref_end, Event, rank, log2_rank_fc, rank_delta, phyloP)
    colnames(ref_end) <- c("curpos",'Event','rank','fc','delta','phyloP')
    lfc_pos <- full_join(ref_start, ref_end, by = c("curpos",'Event','rank','phyloP'))
    lfc_pos$fc_combined <- ifelse(
      is.na(lfc_pos$fc.x), lfc_pos$fc.y, # If log2_rank_fc.x is NA, take log2_rank_fc.y
      ifelse(
        is.na(lfc_pos$fc.y), lfc_pos$fc.x, # If log2_rank_fc.y is NA, take log2_rank_fc.x
        ifelse(
          abs(lfc_pos$fc.x) >= abs(lfc_pos$fc.y), lfc_pos$fc.x, lfc_pos$fc.y # Otherwise, take the max abs value
        )
      )
    )
    lfc_pos$delta_combined <- ifelse(
      is.na(lfc_pos$delta.x), lfc_pos$delta.y, # If log2_rank_fc.x is NA, take log2_rank_fc.y
      ifelse(
        is.na(lfc_pos$delta.y), lfc_pos$delta.x, # If log2_rank_fc.y is NA, take log2_rank_fc.x
        ifelse(
          abs(lfc_pos$delta.x) >= abs(lfc_pos$delta.y), lfc_pos$delta.x, lfc_pos$delta.y # Otherwise, take the max abs value
        )
      )
    )

    lfc_pos$species <- ref_name
    lfc_pos$curpos <- lfc_pos$curpos + 1
    ref_msa <- msa_map %>% filter(species == ref_name) %>% select(msa_pos = curpos, curpos = ref_pos)
    lfc_pos <- left_join(lfc_pos, ref_msa, by = c("curpos"))
    max_width = max(msa_map$curpos)
    shape.events <- c(3,4)
    names(shape.events) <- c('insertion','deletion')
    
#     event_pos <- lfc_pos %>% mutate(Event = case_when(
#                     Event == 'mismatch' ~ NA,
#                     TRUE ~ Event))
#     gg_event <- ggplot(event_pos, aes(y=species)) +
#                   geom_point(aes(x =curpos, shape=Event),
#                                     color = 'black') + 
#                   labs(x = 'bp position', shape = "Event", color = "Event") +
#                   scale_fill_manual(values = col_TFs) + 
#                             theme_cowplot() + 
#                             scale_shape_manual(values = shape.events, name = "Event",
#                                                   guide = guide_legend(override.aes = list(size = 5))) + coord_cartesian(xlim = c(0, 300)) + 
#                 ggplot2::theme(legend.key = ggplot2::element_rect(colour = "black", size = 1.0),
#                     axis.line.y = element_blank(),axis.text.y=element_blank(),
#                               axis.ticks.y = element_blank(),legend.position = 'top',
#                               axis.title.y=element_blank(),
#                               axis.text.x=element_blank(),
#                               axis.line.x = element_blank(),
#                               axis.ticks.x = element_blank(),
#                               axis.title.x=element_blank(), text = element_text(family = "sans", size = 16)) 
    
    gg_phyloP <- ggplot(lfc_pos) +
              geom_tile_rast(aes(x = msa_pos, y = species, fill = phyloP), width = 1.5, height = 1, alpha = 1) +
              scale_fill_gradient2(low = '#6624c0', mid = "#a67be0",
                                      high = '#58b580', midpoint = 0, space = "Lab",
                                      na.value = "gainsboro", 
                                     limits=c(-2,2),  breaks = c(-2,0,2)) + 
              scale_y_discrete(labels = c("phyloP")) + 
              labs(fill = "phyloP") +
              ggplot2::xlab("bp position") +  theme_cowplot() + 
              coord_cartesian(xlim = c(0, max_width)) + 
              ggplot2::theme(legend.key = ggplot2::element_rect(colour = "black", size = 1.0),
                    axis.line.y = element_blank(),
                              axis.ticks.y = element_blank(),legend.position = 'top',
                              axis.title.y=element_blank(),
                              axis.text.x=element_blank(),
                              axis.line.x = element_blank(),
                              axis.ticks.x = element_blank(),
                              axis.title.x=element_blank(), text = element_text(family = "sans", size = 16)) 
#     best_delta <-  lfc_pos %>% filter(delta_combined  == max(delta_combined)) %>% head(n = 1)
#     print(best_delta)
    gg_delta <- ggplot() +            
#               geom_point(data = best_delta, aes(x = msa_pos, y = as.numeric(species)+0.1), color = 'black', fill = 'red', shape = 25, size = 4) +
              geom_tile_rast(data = lfc_pos, aes(x = msa_pos, y = species, fill = delta_combined), width = 1.5, height = 1, alpha = 1) +
              scale_fill_gradient2(low = "dodgerblue2",  
                                     mid = "white",  
                                     high = "firebrick1",  
                                     midpoint = 0,  
                                     space = "Lab",  
                                     na.value = "gainsboro",  
                                     limits=c(-0.5,0.5),  breaks = c(-0.5,0,0.5), 
                                     oob=scales::squish) + 
              scale_y_discrete(labels = c(expression(Delta ~ (m[i] - m[i-1])))) + 
              labs(fill = expression(Delta ~ (m[i] - m[i-1]))) +
              ggplot2::xlab("bp position") +  theme_cowplot() + 
              coord_cartesian(xlim = c(0, max_width)) + 
              ggplot2::theme(legend.key = ggplot2::element_rect(colour = "black", size = 1.0),
                    axis.line.y = element_blank(),
                              axis.ticks.y = element_blank(),legend.position = 'top',
                              axis.title.y=element_blank(), text = element_text(family = "sans", size = 16))
    lfc_filter_pos <-  lfc_pos %>% filter(delta_combined  > 0)
    gg_density <- ggplot(lfc_filter_pos) +
              stat_density(aes(x = msa_pos, y = species,fill = stat(density)), geom = "raster", position = "identity") +
              scale_fill_viridis(name = "mutation density", 
                                 breaks = c(0,0.002,0.003), ) +
                scale_y_discrete(labels = c("beneficial mutation\ndensity")) + 
              ggplot2::xlab("bp position") +  theme_cowplot() + 
              coord_cartesian(xlim = c(0, max_width)) +
              ggplot2::theme(legend.position = 'none',
                    axis.line.y = element_blank(),
                              axis.ticks.y = element_blank(),
                              axis.title.y=element_blank(),
                              axis.text.x=element_blank(),
                              axis.line.x = element_blank(),
                              axis.ticks.x = element_blank(),
                              axis.title.x=element_blank(), text = element_text(family = "sans", size = 16)) 
    
    ref_TFBS <- tfbs_map %>% filter(tile_name == ref_tile)
    ref_msa <- msa_map %>% filter(species == ref_name)
    ref_TFBS <- map_TFBS_to_msa_pos(ref_TFBS, ref_msa) %>% mutate(species = ref_name)
    
    target_TFBS <- tfbs_map %>% filter(tile_name == target_tile)
    target_msa <- msa_map %>% filter(species == target_name)
    target_TFBS <- map_TFBS_to_msa_pos(target_TFBS, target_msa) %>% mutate(species = target_name)
    final_TFBS <- bind_rows(ref_TFBS, target_TFBS) %>% 
                        mutate(center_pos = (msa_start + msa_end)/2)
    final_TFBS$species <- factor(final_TFBS$species, levels = c(ref_name, target_name))
    
    select.images <- sort(Sys.glob('../PYS2_oCRE_BC_libs_20240904/imgs/*.png'))
    select.images <- data.frame(uid = select.images,
                               label = c("fullTreeAnc223","Felis_catus","Canis_lupus_familiaris","fullTreeAnc115",'fullTreeAnc239',
                                        "Equus_caballus","Homo_sapiens",'fullTreeAnc237',"Mus_musculus","fullTreeAnc110","Rattus_norvegicus","fullTreeAnc68"),
                               name = c("Carnivora","Felis_catus","Canis_lupus_familiaris","Euarchontoglires",
                                       "First_Ancestor","Equus_caballus","Homo_sapiens",'Laurasiatheria',
                                        "Mus_musculus","Primates","Rattus_norvegicus","Rodentia"))
    
    labels <- c()

    for(species_oi in c(ref_name, target_name)){
        img_path <- select.images %>% filter(label == species_oi) %>% pull(uid)
        labels <- c(labels, paste0("<img src='", img_path,  "' width='25' /><br>*", species_oi,"*"))
    }
    gg_tfbs <- make_tfbs_plot(final_TFBS, max_width, labels) + ggtitle(CRE_oi) +
                        theme(axis.text.x=element_blank(),
                                  axis.line.x = element_blank(),
                              axis.title.x=element_blank(),
                              axis.ticks.x = element_blank())

    ref_img_path <- select.images %>% filter(label == ref_name) %>% pull(uid)
    target_img_path <- select.images %>% filter(label == target_name) %>% pull(uid)
    labels <- c(paste0("<img src='", target_img_path,  "' width='25' />", 
                       "<span class='plus-sign'>+</span>",
                       "<img src='", ref_img_path, "' width='25' /><br>*",
                      paste0(max(lfc_pos$rank), " mutations"),"*"))
    
    final_TFBS <- tfbs_map %>% filter(tile_name == final_tile)
    final_msa <- msa_map %>% filter(species == final_tile)
    final_TFBS <- map_TFBS_to_msa_pos(final_TFBS, final_msa) %>% 
                        mutate(center_pos = (msa_start + msa_end)/2) %>% mutate(species = 'final')
    final_tfbs <- make_tfbs_plot(final_TFBS, max_width, labels) + 
                        theme(axis.text.x=element_blank(),
                                  axis.line.x = element_blank(),
                              axis.title.x=element_blank(),
                              axis.ticks.x = element_blank())

    p <- gg_tfbs + gg_density + gg_phyloP + gg_delta + final_tfbs +
                plot_layout(heights = c(0.3, 0.1 ,0.1, 0.1, 0.1), ncol = 1,guides = "collect") & theme(legend.position = 'bottom',
                                                                                             legend.box="vertical", legend.margin=margin())
    return(p)
}