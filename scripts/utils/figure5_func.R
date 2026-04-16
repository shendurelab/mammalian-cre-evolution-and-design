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
library(ggstar)
source("utils/trajectory_utils.R")

#### load TFBS mapping
final_TFs <- c("Jun_Atf3","Foxa2","Gata4/6","Sox17","Klf4","Hnf1b")
col_TFs <- c("firebrick2","forestgreen","darkorange2","dodgerblue1","mediumorchid3",'goldenrod1')
names(col_TFs) <- final_TFs 

select.images <- Sys.glob('../imgs/*.png')
select.images <- data.frame(uid = select.images,
                           label = c("Psammomys_obesus","fullTreeAnc237","fullTreeAnc223","Felis_catus","Canis_lupus_familiaris","fullTreeAnc115",'fullTreeAnc239',
                                    "Equus_caballus","Homo_sapiens","Mus_musculus","fullTreeAnc110","Rattus_norvegicus","fullTreeAnc68"),
                           name = c('Psammomys_obesus',"Laurasiatheria","Carnivora","Felis_catus","Canis_lupus_familiaris","Euarchontoglires",
                                   "First_Ancestor","Equus_caballus","Homo_sapiens","Mus_musculus","Primates","Rattus_norvegicus","Rodentia"))


### load functions for plotting
make_msa_tfbs_traj_df <- function(fasta.file, tile.list, ref_name, CRE_oi, tfbs_df, func_tfbs_df){
    colnames(tfbs_df)[8] <- 'species'
    colnames(tfbs_df)[17] <- 'ori_CRE_id'
    msa_df <- make_alignment_dataframe(fasta.file, tile.list, ref_name, CRE_oi)
    msa_df_TFBS <- make_cre_tfbs_alignment_dataframe(msa_df, tfbs_df, CRE_oi)
#     msa_df_TFBS$species <- factor(msa_df_TFBS$species, levels = tile.list[tile.list%in%unique(msa_df_TFBS$species)])
    
    ref_msa_df <- msa_df %>% filter(species == ref_name)
    msa_ref_df_TFBS <- msa_df_TFBS %>% filter(species == ref_name)
    msa_ref_df_TFBS_gr <- GRanges(msa_ref_df_TFBS %>% dplyr::select(chr = TF_name, start = msa_start, end = msa_end))
    
    functional_tfbs_cre <- func_tfbs_df %>% filter(CRE == CRE_oi)
    functional_tfbs_cre$rank_TF = rank(-functional_tfbs_cre$TFBS_activity, ties.method = "random")
    mm.start <- match(functional_tfbs_cre$TFBS_start, ref_msa_df$ref_pos)
    mm.end <- match(functional_tfbs_cre$TFBS_end, ref_msa_df$ref_pos)
    functional_tfbs_cre$msa_start <- ref_msa_df$curpos[mm.start]
    functional_tfbs_cre$msa_end <- ref_msa_df$curpos[mm.end]
    msa_ref_df_TFBS <- left_join(msa_ref_df_TFBS, functional_tfbs_cre %>% dplyr::select(TF_name, msa_start, msa_end, rank_TF),
                                by = c("TF_name",'msa_start','msa_end'))
    
    
    final_df_TFBS <- data.frame()
    for(species_oi in unique(msa_df_TFBS$species)){
        species_tfbs_df <- msa_df_TFBS %>% filter(species == species_oi)
        species_tfbs_gr <- GRanges(species_tfbs_df %>% dplyr::select(chr = TF_name, start = msa_start, end = msa_end))
        
        overlaps <- findOverlaps(msa_ref_df_TFBS_gr, species_tfbs_gr, type = 'start')
        species_tfbs_df_w_overlaps <- species_tfbs_df[unique(subjectHits(overlaps)),]
        species_tfbs_df_w_overlaps$norm_affinity_ref <- msa_ref_df_TFBS$norm_affinity[unique(queryHits(overlaps))]
        species_tfbs_df_w_overlaps$rank_TF <- msa_ref_df_TFBS$rank_TF[unique(queryHits(overlaps))]
        species_tfbs_df <- left_join(species_tfbs_df, species_tfbs_df_w_overlaps)
        
        final_df_TFBS <- bind_rows(final_df_TFBS, species_tfbs_df)
        
    }
    
    final_df_TFBS <- final_df_TFBS %>% 
                mutate(binding_div = abs(norm_affinity_ref - norm_affinity)/(norm_affinity_ref + norm_affinity))
    final_df_TFBS <- final_df_TFBS %>% mutate(TF_name2 = case_when(
                            is.na(rank_TF) ~ NA,
                            TRUE ~ paste0(TF_name,"_",rank_TF)))
    return(list(msa_df = msa_df, msa_df_TFBS = final_df_TFBS, msa_func_TFBS = functional_tfbs_cre))
}

plot_tfbs_trajectory <- function(msa_tfbs_df, func_tfbs_df, CRE_oi, ref_name, labels, 
                                 aff_cutoff = 0.05, tfbs_size = 0.5 ,rm_species = NULL, order.tiles = NULL){
    
    msa_filter_tfbs_df <- msa_tfbs_df %>% filter(norm_affinity > aff_cutoff)
    
    arrow_forward <- arrow(type = "closed", angle = 30, length = unit(0.2, "cm"))
    arrow_reverse <- arrow(type = "closed", angle = 30, length = unit(0.2, "cm"))
    
    cre_func_tfbs <- data.frame()
    for(species_oi in unique(msa_filter_tfbs_df$species)){
        species_func_tfbs <- make_aff_corr_df(msa_filter_tfbs_df, func_tfbs_df, species_oi, ref_name) %>%
                            dplyr::select(TF_name,msa_start,msa_end, norm_affinity_query, TFBS_orientation_query, functional_TF) %>% 
                            dplyr::select(TF_name,msa_start,msa_end,norm_affinity = norm_affinity_query, 
                                          TFBS_orientation = TFBS_orientation_query, functional_TF)
        species_func_tfbs$species = species_oi
        cre_func_tfbs <- bind_rows(cre_func_tfbs, species_func_tfbs)
    }
    
    cre_func_tfbs <- cre_func_tfbs %>% filter(functional_TF) %>% mutate(center_pos = (msa_start + msa_end)/2)
    
    # Calculate global min and max for each TF_name across both datasets
    tf_min_max <- bind_rows(
      msa_filter_tfbs_df %>% dplyr::select(TF_name, norm_affinity),
      cre_func_tfbs %>% dplyr::select(TF_name, norm_affinity)
    ) %>%
      group_by(TF_name) %>%
      summarise(
        global_min = min(norm_affinity, na.rm = TRUE),
        global_max = max(norm_affinity, na.rm = TRUE)
      )
    
    final_msa_filter_tfbs_df <- data.frame()
    final_cre_func_tfbs <- data.frame()
    for(TF_oi in unique(msa_filter_tfbs_df$TF_name)){        
        tf_global_min <- tf_min_max %>% filter(TF_name == TF_oi) %>% pull(global_min)
        tf_global_max <- tf_min_max %>% filter(TF_name == TF_oi) %>% pull(global_max)
        
        tf_df <- msa_filter_tfbs_df %>% filter(TF_name == TF_oi) %>% 
                mutate(adjusted_affinity = scales::rescale(norm_affinity, to = c(0.2, 1), from = c(tf_global_min, tf_global_max)))
        func_tf_df <- cre_func_tfbs %>% filter(TF_name == TF_oi) %>% 
                mutate(adjusted_affinity = scales::rescale(norm_affinity, to = c(0.2, 1), from = c(tf_global_min, tf_global_max)))
        final_msa_filter_tfbs_df <- bind_rows(final_msa_filter_tfbs_df, tf_df)
        final_cre_func_tfbs <- bind_rows(final_cre_func_tfbs, func_tf_df)
    }
    
    
    # Create a new column for row-specific colors
    final_msa_filter_tfbs_df$hue_color <- mapply(
      function(tf, affinity) {
        tf_color <- col_TFs[tf]  # Get the base color for the TF
        colorRampPalette(c("white", tf_color))(100)[round(affinity * 99) + 1]  # Interpolate color based on norm_affinity
      },
      final_msa_filter_tfbs_df$TF_name,  # Transcription factor name
      final_msa_filter_tfbs_df$adjusted_affinity  # Normalized affinity
    )
    
    # Repeat for `cre_func_tfbs` if necessary
    final_cre_func_tfbs$hue_color <- mapply(
      function(tf, affinity) {
        tf_color <- col_TFs[tf]  # Get the base color for the TF
        colorRampPalette(c("white", tf_color))(100)[round(affinity * 99) + 1]  # Interpolate color based on norm_affinity
      },
      final_cre_func_tfbs$TF_name,  # Transcription factor name
      final_cre_func_tfbs$adjusted_affinity  # Normalized affinity
    )
    
    if(is.null(rm_species)){
    }else{
        final_msa_filter_tfbs_df <- final_msa_filter_tfbs_df %>% filter(!species%in%rm_species)
        final_cre_func_tfbs <- final_cre_func_tfbs %>% filter(!species%in%rm_species)
    }
    
    final_msa_filter_tfbs_df$species <- factor(final_msa_filter_tfbs_df$species, levels = rev(order.tiles))
    final_cre_func_tfbs$species <- factor(final_cre_func_tfbs$species, levels = rev(order.tiles))
    
    library(ggtext)

    gg_tfbs <- ggplot() +
                  geom_hline(data = final_msa_filter_tfbs_df, aes(yintercept = species), color = "black", alpha = 0.5) +  # Add horizontal line
                  geom_segment(
                    data = subset(final_msa_filter_tfbs_df, TFBS_orientation == "+"),
                    aes(x = center_pos, xend = center_pos + 2, y = species, yend = species, color = hue_color),
                    arrow = arrow_forward, size = tfbs_size, position = position_nudge(y = 0, x = 3)
                  ) +
                  geom_segment(
                    data = subset(final_msa_filter_tfbs_df, TFBS_orientation == "-"),
                    aes(x = center_pos, xend = center_pos - 2, y = species, yend = species, color = hue_color),
                    arrow = arrow_reverse, size = tfbs_size, position = position_nudge(y = 0, x = 3)
                  ) +
                  geom_segment(
                    data = subset(final_cre_func_tfbs, TFBS_orientation == "+"),
                    aes(x = center_pos, xend = center_pos + 2, y = species, yend = species),
                    arrow = arrow_forward, size = (tfbs_size * 3) + 1, color = "black", position = position_nudge(y = 0, x = 3)
                  ) +
                  geom_segment(
                    data = subset(final_cre_func_tfbs, TFBS_orientation == "+"),
                    aes(x = center_pos, xend = center_pos + 2, y = species, yend = species, color = hue_color),
                    arrow = arrow_forward, size = (tfbs_size * 3), position = position_nudge(y = 0, x = 3)
                  ) +
                  geom_segment(
                    data = subset(final_cre_func_tfbs, TFBS_orientation == "-"),
                    aes(x = center_pos, xend = center_pos - 2, y = species, yend = species),
                    arrow = arrow_reverse, size = (tfbs_size * 3) + 1, color = "black", position = position_nudge(y = 0, x = 3)
                  ) +
                  geom_segment(
                    data = subset(final_cre_func_tfbs, TFBS_orientation == "-"),
                    aes(x = center_pos, xend = center_pos - 2, y = species, yend = species, color = hue_color),
                    arrow = arrow_reverse, size = (tfbs_size * 3), position = position_nudge(y = 0, x = 3)
                  ) + labs(x= 'msa bp position') +
                  scale_color_identity() +  # Use precomputed colors
                  scale_x_continuous(limits = c(0, 300)) + 
#                  scale_color_manual(values = col_TFs) + 
                theme_cowplot() + scale_y_discrete(labels = labels) + 
                            theme(#legend.position = 'bottom',legend.box="horizontal", legend.margin=margin(),
                                  axis.text.y = element_markdown(color = "black", size = 8),
                                  axis.line.y = element_blank(),
                                  axis.ticks.y = element_blank(),
                                  axis.title.y=element_blank())
    return(gg_tfbs)
}

plot_mpra_act_trajectory <- function(mpra_act_df, ref_mpra_df, control_mpra_df, CRE_oi, species.list, labels = NULL, pic_labels = F, flip = F){
    colnames(mpra_act_df)[1] <- 'species'
    species_act <- mpra_act_df %>% filter(CRE == CRE_oi) %>% filter(species%in%species.list)
    species_act$species <- factor(species_act$species, levels = species.list[species.list%in%species_act$species])
    cre_oi_wt_act <- ref_mpra_df %>% filter(CRE == CRE_oi) %>% pull(mean_MPRA_act)
    cre_oi_wt_sd <- ref_mpra_df %>% filter(CRE == CRE_oi) %>% pull(sd_MPRA_act)
    minP_oi_act <- control_mpra_df %>% filter(CRE_id == 'minP') %>% pull(mean_MPRA_act)
    minP_oi_sd <- control_mpra_df %>% filter(CRE_id == 'minP') %>% pull(sd_MPRA_act)
    
    species_act <- species_act %>% mutate(point_type = case_when(
                                    mean_MPRA_act == max(mean_MPRA_act) ~ 'Max',
                                    mean_MPRA_act == min(mean_MPRA_act) ~ 'Min',
                                    TRUE ~ 'none'))
    if(flip){
        mpra_bar <- ggplot(species_act, aes(y=mean_MPRA_act, x=species)) + #, color=mean_MPRA_act)) + 
            geom_star(aes(starshape = point_type, fill = point_type, size = point_type)) + 
            geom_errorbar(aes(ymin=mean_MPRA_act-sd_MPRA_act, ymax=mean_MPRA_act+sd_MPRA_act, color = point_type), width=.2,
                         position=position_dodge(0.05)) +
            scale_size_manual(values = c('Max' = 3,'Min' = 3,'none' = 1.5)) +
            scale_starshape_manual(values = c('Max' = 1,'Min' = 1,'none' = 15)) + 
            scale_fill_manual(values = c('Max' = '#ee4035','Min' = '#0392cf','none' = 'gainsboro')) + 
            scale_color_manual(values = c('Max' = '#ee4035','Min' = '#0392cf','none' = 'gainsboro')) + 
        #         geom_rect(xmin = cre_oi_wt_act-cre_oi_wt_sd, xmax = cre_oi_wt_act+cre_oi_wt_sd, ymin = -Inf, ymax = Inf,
        #            alpha = .2,fill = "#008000") + 
        #             geom_rect(xmin = minP_oi_act - minP_oi_sd, xmax = minP_oi_act + minP_oi_sd, ymin = -Inf, ymax = Inf,
        #            alpha = .1,fill = "blue")  +
                geom_hline(yintercept = minP_oi_act ,color = "blue", linetype = 'dashed')  +
            geom_hline(yintercept = cre_oi_wt_act ,color = "#008000", linetype = 'dashed') +
    #         scale_color_gradient(breaks = c(10^-1, 10^0, 10^1),
    #                labels = scales::trans_format("log10", scales::math_format(10^.x)), limits = c(0.05, 15),  # Quantile limits from your range
    #                                     oob = scales::squish,  # To handle values outside the specified limits
    #                                     trans = "log10",
    #                                    low = '#42047e',high = '#07f49e', na.value = 'gainsboro') +
            scale_y_log10(breaks = c(10^-1, 10^0, 10^1),
                   labels = scales::trans_format("log10", scales::math_format(10^.x)),
                          limits = c(0.05,15))+ annotation_logticks(sides = 'l')+ theme_classic() + 
                labs(y = 'MPRA act.', color = 'MPRA act.')
        if(is.null(labels)){
            mpra_bar <- mpra_bar + 
                        theme(legend.position = 'none',
                                      axis.text.x = element_blank(),
                                      axis.line.x = element_blank(),
                                      axis.ticks.x = element_blank(),
                                      axis.title.x=element_blank())
        }else{
            mpra_bar <- mpra_bar + scale_x_discrete(labels = labels)
            if(pic_labels){
                mpra_bar <- mpra_bar + theme(legend.position = 'none',legend.box="horizontal", legend.margin=margin(),
                                  axis.text.x = element_markdown(color = "black", size = 8),
                                  axis.line.x = element_blank(),
                                  axis.ticks.x = element_blank(),
                                  axis.title.x=element_blank())
            }else{
                mpra_bar <- mpra_bar + theme(axis.text.x = element_text(angle = 45, hjust = 1))
            }
        }
    }else{
        mpra_bar <- ggplot(species_act, aes(x=mean_MPRA_act, y=species)) + #, color=mean_MPRA_act)) + 
            geom_star(aes(starshape = point_type, fill = point_type, size = point_type)) + 
            geom_errorbar(aes(xmin=mean_MPRA_act-sd_MPRA_act, xmax=mean_MPRA_act+sd_MPRA_act, color = point_type), width=.2,
                         position=position_dodge(0.05)) +
            scale_size_manual(values = c('Max' = 3,'Min' = 3,'none' = 1.5)) +
            scale_starshape_manual(values = c('Max' = 1,'Min' = 1,'none' = 15)) + 
            scale_fill_manual(values = c('Max' = '#ee4035','Min' = '#0392cf','none' = 'gainsboro')) + 
            scale_color_manual(values = c('Max' = '#ee4035','Min' = '#0392cf','none' = 'gainsboro')) + 
        #         geom_rect(xmin = cre_oi_wt_act-cre_oi_wt_sd, xmax = cre_oi_wt_act+cre_oi_wt_sd, ymin = -Inf, ymax = Inf,
        #            alpha = .2,fill = "#008000") + 
        #             geom_rect(xmin = minP_oi_act - minP_oi_sd, xmax = minP_oi_act + minP_oi_sd, ymin = -Inf, ymax = Inf,
        #            alpha = .1,fill = "blue")  +
                geom_vline(xintercept = minP_oi_act ,color = "blue", linetype = 'dashed')  +
            geom_vline(xintercept = cre_oi_wt_act ,color = "#008000", linetype = 'dashed') +
    #         scale_color_gradient(breaks = c(10^-1, 10^0, 10^1),
    #                labels = scales::trans_format("log10", scales::math_format(10^.x)), limits = c(0.05, 15),  # Quantile limits from your range
    #                                     oob = scales::squish,  # To handle values outside the specified limits
    #                                     trans = "log10",
    #                                    low = '#42047e',high = '#07f49e', na.value = 'gainsboro') +
            scale_x_log10(breaks = c(10^-1, 10^0, 10^1),
                   labels = scales::trans_format("log10", scales::math_format(10^.x)),
                          limits = c(0.05,15))+ annotation_logticks(sides = 'b')+ theme_classic() + 
                labs(x = 'MPRA act.', color = 'MPRA act.')
        if(is.null(labels)){
            mpra_bar <- mpra_bar + 
                        theme(legend.position = 'none',
                              axis.text.y = element_blank(),
                              axis.line.y = element_blank(),
                              axis.ticks.y = element_blank(),
                              axis.title.y=element_blank())
        }else{
            mpra_bar <- mpra_bar + scale_y_discrete(labels = labels)
            if(pic_labels){
                mpra_bar <- mpra_bar + theme(legend.position = 'none',legend.box="horizontal", legend.margin=margin(),
                                  axis.text.y = element_markdown(color = "black", size = 8),
                                  axis.line.y = element_blank(),
                                  axis.ticks.y = element_blank(),
                                  axis.title.y=element_blank())
            }else{
                mpra_bar <- mpra_bar
            }
        }
    }
    return(mpra_bar)
}

### load functions for plotting
make_msa_tfbs_tile_df <- function(fasta.file, tile.list, ref_name, CRE_oi, tfbs_df, func_tfbs_df){
    print("make MSA dataframe")
    msa_df <- make_alignment_dataframe(fasta.file, tile.list, ref_name,CRE_oi)
    ref_msa_df <- msa_df %>% filter(species == ref_name)
    
    print("make MSA TFBS dataframe")
    msa_df_TFBS <- make_cre_tfbs_alignment_dataframe(msa_df, tfbs_df %>% select(species = tile_name, ori_CRE_id = CRE,start = start_pos, end = end_pos, TF_name, TFBS_orientation, norm_affinity), CRE_oi)
    msa_ref_df_TFBS <- msa_df_TFBS %>% filter(species == ref_name)
    msa_ref_df_TFBS_gr <- GRanges(msa_ref_df_TFBS %>% dplyr::select(chr = TF_name, start = msa_start, end = msa_end))
    
    print("make functional TFBS dataframe")
    functional_tfbs_cre <- func_tfbs_df %>% filter(CRE == CRE_oi)
    functional_tfbs_cre <- functional_tfbs_cre %>%
              arrange(desc(functional_TF), desc(TFBS_activity)) %>%  # Prioritize TRUE first, then sort by TFBS_activity
              mutate(rank_TF = rank(-TFBS_activity, ties.method = "random"))
    #functional_tfbs_cre$rank_TF = rank(-functional_tfbs_cre$TFBS_activity, ties.method = "random")

    mm.start <- match(functional_tfbs_cre$TFBS_start, ref_msa_df$ref_pos)
    mm.end <- match(functional_tfbs_cre$TFBS_end, ref_msa_df$ref_pos)
    functional_tfbs_cre$msa_start <- ref_msa_df$curpos[mm.start]
    functional_tfbs_cre$msa_end <- ref_msa_df$curpos[mm.end]
    msa_ref_df_TFBS <- left_join(msa_ref_df_TFBS, functional_tfbs_cre %>% dplyr::select(TF_name, msa_start, msa_end, rank_TF),
                                by = c("TF_name",'msa_start','msa_end'))
    
    print("add functional TFBS to MSA TFBS dataframe")
    final_df_TFBS <- data.frame()
    for(species_oi in unique(msa_df_TFBS$species)){
        species_tfbs_df <- msa_df_TFBS %>% filter(species == species_oi)
        species_tfbs_gr <- GRanges(species_tfbs_df %>% dplyr::select(chr = TF_name, start = msa_start, end = msa_end))
        
        overlaps <- findOverlaps(msa_ref_df_TFBS_gr, species_tfbs_gr, type = 'start')
        species_tfbs_df_w_overlaps <- species_tfbs_df[unique(subjectHits(overlaps)),]
        species_tfbs_df_w_overlaps$norm_affinity_ref <- msa_ref_df_TFBS$norm_affinity[unique(queryHits(overlaps))]
        species_tfbs_df_w_overlaps$rank_TF <- msa_ref_df_TFBS$rank_TF[unique(queryHits(overlaps))]
        species_tfbs_df <- left_join(species_tfbs_df, species_tfbs_df_w_overlaps)
        
        final_df_TFBS <- bind_rows(final_df_TFBS, species_tfbs_df)
        
    }
    
    final_df_TFBS <- final_df_TFBS %>% mutate(TF_name2 = case_when(
                            is.na(rank_TF) ~ NA,
                            TRUE ~ paste0(TF_name,"_",rank_TF)))
    
    
    return(list(msa_df_TFBS = final_df_TFBS, msa_func_TFBS = functional_tfbs_cre))
}

mark_func_overlap_df <- function(events_df, func_tfbs_df){
    events_df <- events_df %>% mutate(tag = row_number())
    query_gr <- GRanges(events_df %>% dplyr::select(start = ref_start, end = ref_end) %>% mutate(chr = 'chr1'))

    func_tfbs_gr <- GRanges(func_tfbs_df %>% dplyr::select(start = TFBS_start, end = TFBS_end) %>% mutate(chr = 'chr1'))
    overlaps <- findOverlaps(query_gr, func_tfbs_gr, type = 'any')
    
    func_events <- events_df[queryHits(overlaps),] %>% mutate(func_tfbs = 'func')
    non_func_events <- events_df %>% filter(!tag%in%unique(func_events$tag)) %>% mutate(func_tfbs = 'not func')
    # if any overlap with func tfbs mark as so
    events_df <- bind_rows(func_events, non_func_events) %>% dplyr::select(-tag)
    return(events_df)
}


make_msa_tfbs_cactus_df <- function(fasta.file, species.list, ref_name, CRE_oi, tfbs_df, func_tfbs_df){
    msa_df <- make_alignment_dataframe(fasta.file, species.list, ref_name,CRE_oi)
    ref_msa_df <- msa_df %>% filter(species == ref_name)
    
    msa_df_TFBS <- make_cre_tfbs_alignment_dataframe(msa_df, tfbs_df %>% select(species, ori_CRE_id = CRE,start = start_pos, end = end_pos, TF_name, TFBS_orientation, norm_affinity), CRE_oi)
    msa_ref_df_TFBS <- msa_df_TFBS %>% filter(species == ref_name)
    msa_ref_df_TFBS_gr <- GRanges(msa_ref_df_TFBS %>% dplyr::select(chr = TF_name, start = msa_start, end = msa_end))
    
    functional_tfbs_cre <- func_tfbs_df %>% filter(CRE == CRE_oi)
    functional_tfbs_cre$rank_TF = rank(-functional_tfbs_cre$TFBS_activity, ties.method = "random")
    mm.start <- match(functional_tfbs_cre$TFBS_start, ref_msa_df$ref_pos)
    mm.end <- match(functional_tfbs_cre$TFBS_end, ref_msa_df$ref_pos)
    functional_tfbs_cre$msa_start <- ref_msa_df$curpos[mm.start]
    functional_tfbs_cre$msa_end <- ref_msa_df$curpos[mm.end]
    msa_ref_df_TFBS <- left_join(msa_ref_df_TFBS, functional_tfbs_cre %>% dplyr::select(TF_name, msa_start, msa_end, rank_TF),
                                by = c("TF_name",'msa_start','msa_end'))
    
    final_df_TFBS <- data.frame()
    for(species_oi in unique(msa_df_TFBS$species)){
        species_tfbs_df <- msa_df_TFBS %>% filter(species == species_oi)
        species_tfbs_gr <- GRanges(species_tfbs_df %>% dplyr::select(chr = TF_name, start = msa_start, end = msa_end))
        
        overlaps <- findOverlaps(msa_ref_df_TFBS_gr, species_tfbs_gr, type = 'start')
        species_tfbs_df_w_overlaps <- species_tfbs_df[unique(subjectHits(overlaps)),]
        species_tfbs_df_w_overlaps$norm_affinity_ref <- msa_ref_df_TFBS$norm_affinity[unique(queryHits(overlaps))]
        species_tfbs_df_w_overlaps$rank_TF <- msa_ref_df_TFBS$rank_TF[unique(queryHits(overlaps))]
        species_tfbs_df <- left_join(species_tfbs_df, species_tfbs_df_w_overlaps)
        
        final_df_TFBS <- bind_rows(final_df_TFBS, species_tfbs_df)
        
    }
    
    final_df_TFBS <- final_df_TFBS %>% mutate(TF_name2 = case_when(
                            is.na(rank_TF) ~ NA,
                            TRUE ~ paste0(TF_name,"_",rank_TF)))
    
    
    return(list(msa_df_TFBS = final_df_TFBS, msa_func_TFBS = functional_tfbs_cre))
}

make_aff_corr_df <- function(msa_tfbs_df, func_tfbs_df, query_name, ref_name){
    query_tfbs_df <- msa_tfbs_df %>% filter(species == query_name) %>% 
            dplyr::select(TF_name, msa_start, msa_end, norm_affinity, TFBS_orientation)
    query_tfbs_gr <- GRanges(query_tfbs_df %>% dplyr::select(chr = TF_name, start = msa_start, end = msa_end))

    ref_tfbs_df <- msa_tfbs_df %>% filter(species == ref_name) %>% dplyr::select(TF_name, msa_start, msa_end, norm_affinity, TFBS_orientation)
    ref_tfbs_df <- left_join(ref_tfbs_df, func_tfbs_df %>% dplyr::select(TF_name, msa_start, msa_end, rank_TF, functional_TF),
                                by = c("TF_name",'msa_start','msa_end'))
    ref_tfbs_gr <- GRanges(ref_tfbs_df %>% dplyr::select(chr = TF_name, start = msa_start, end = msa_end))
    overlaps <- findOverlaps(ref_tfbs_gr, query_tfbs_gr, type = 'any')
    # Compute intersection widths
    intersection_widths <- width(pintersect(ref_tfbs_gr[queryHits(overlaps)], query_tfbs_gr[subjectHits(overlaps)]))
    # Create a data frame of overlaps with widths
    overlap_df <- data.frame(
      queryHit = queryHits(overlaps),
      subjectHit = subjectHits(overlaps),
      width = intersection_widths
    )

    # Identify the subjectHit with the maximum overlap width for each queryHit
    max_overlap_df <- overlap_df %>%
      group_by(queryHit) %>%
      slice_max(width, n = 1) %>%  # Select the subjectHit with the maximum overlap
      ungroup()

    # Subset reference and query TFBS data
    comp_tfbs_df <- ref_tfbs_df[max_overlap_df$queryHit, ]
    comp_tfbs_df$norm_affinity_query <- query_tfbs_df$norm_affinity[max_overlap_df$subjectHit]
    comp_tfbs_df$TFBS_orientation_query <- query_tfbs_df$TFBS_orientation[max_overlap_df$subjectHit]
#     comp_tfbs_df <- comp_tfbs_df %>% mutate(func_tfbs = case_when(
#                             is.na(rank_TF) ~ FALSE,
#                             TRUE ~ TRUE))
    comp_tfbs_df$TF_name <- factor(comp_tfbs_df$TF_name, levels = final_TFs)
    return(comp_tfbs_df)
}

make_aff_corr_bar <- function(msa_tfbs_df, func_tfbs_df, species.list, ref_name, labels = NULL, pic_labels = F, ncol =1 , nrow = 5){
    corr.df <- data.frame()
    for(species_oi in species.list){
        tmp.df <- make_aff_corr_df(msa_tfbs_df, func_tfbs_df, species_oi, 'Mus_musculus')
        comp_aff_corr <- tmp.df %>% group_by(TF_name) %>% summarise(R = cor(norm_affinity_query, norm_affinity, method = 'spearman'))
        comp_aff_corr$species <- species_oi
        corr.df <- bind_rows(corr.df, comp_aff_corr)
    }
    
    corr.df$species <- factor(corr.df$species, levels = species.list)

    corr_bar <- ggplot(corr.df, aes(x=R, y=species)) + 
                geom_point(size = 2) + facet_wrap(~TF_name, ncol = ncol, nrow = nrow) + 
                scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0",0.5,'1'),
                              limits = c(0,1)) + theme_classic() + 
                    labs(x = 'Spearman Correlation Coefficient (R)')
    if(is.null(labels)){
            corr_bar <- corr_bar + 
                        theme(legend.position = 'bottom',
                              axis.text.y = element_blank(),
                              axis.line.y = element_blank(),
                              axis.ticks.y = element_blank(),
                              axis.title.y=element_blank())
    }else{
        corr_bar <- corr_bar + scale_y_discrete(labels = labels)
        if(pic_labels){
            corr_bar <- corr_bar + theme(legend.position = 'bottom',legend.box="horizontal", legend.margin=margin(),
                              axis.text.y = element_markdown(color = "black", size = 8),
                              axis.line.y = element_blank(),
                              axis.ticks.y = element_blank(),
                              axis.title.y=element_blank())
        }else{
            corr_bar <- corr_bar
        }
    }
    return(corr_bar)
}

plot_aff_corr_plot <- function(comp_tfbs_df, ncol = 1, nrow = 5){
    comp_aff_corr <- comp_tfbs_df %>% group_by(TF_name) %>% summarise(R = cor(norm_affinity_query, norm_affinity, method = 'spearman'))
    print(paste0("Avg. corr for TFBS affinities = ", mean(comp_aff_corr$R)))

    cor.aff <- ggplot() +
                geom_point_rast(data = comp_tfbs_df, aes(x = norm_affinity, y = norm_affinity_query, 
                                                                color = func_tfbs, alpha = func_tfbs, size = func_tfbs)) +             
                scale_color_manual(values = c('gainsboro','red')) + 
                scale_alpha_manual(values = c(0.3,0.7)) +
                scale_size_manual(values = c(1, 2)) +
#                 labs(x = 'Mouse binding affinity', y = 'First Ancestor\nbinding affinity') + 
                geom_vline(xintercept = 0.05, color = 'grey', linetype = 'dashed') + 
                geom_hline(yintercept = 0.05, color = 'grey', linetype = 'dashed') + 
                geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") + # Add diagonal dashed line
                geom_text(data = comp_aff_corr, aes(label = paste0('R = ', round(R, 2))), x= Inf, y = Inf,
                          hjust = 1,  # Align to the right
                        vjust = 1) +   # Align to the bottom  
                facet_wrap(~TF_name, ncol = ncol, nrow = nrow) + 
                scale_x_log10(breaks = c(10^-7, 10^-4, 10^-1),
                       labels = scales::trans_format("log10", scales::math_format(10^.x))
                     ) + 
                scale_y_log10(breaks = c(10^-7, 10^-4, 10^-1),
                       labels = scales::trans_format("log10", scales::math_format(10^.x))
                     ) + annotation_logticks() + 
                theme_classic() + theme(strip.text = element_text(size = 10))
    return(cor.aff)
}

collect_labels <- function(msa_df_TFBS, written = F){
    labels <- c()

    for(species_oi in levels(msa_df_TFBS$species)){
        if(species_oi%in%select.images$label){
            img_path <- select.images %>% filter(label == species_oi) %>% pull(uid)
            name_oi <- select.images %>% filter(label == species_oi) %>% pull(name)
            if(written){
                labels <- c(labels, name_oi)
            }else{
                labels <- c(labels, paste0("<img src='", img_path,  "' width='", 25, "' />"))            }
        } else{
            if(written){
                labels <- c(labels, "")
            }else{
                labels <- c(labels, NA)
            }
        }
    }
    return(labels)
}

plot_div <- function(msa_tfbs_df, species.list, labels){
    div.data <- msa_tfbs_df[complete.cases(msa_tfbs_df$rank_TF),]
    div.data$species <- factor(div.data$species, levels = rev(species.list[species.list%in%unique(div.data$species)]))
    div.data$TF_name <- factor(div.data$TF_name, levels = final_TFs)

    div.plot <- ggplot(div.data, 
                       aes(x = species, y = binding_div, color = rank_TF, group = TF_name2)) + 
                geom_step() + facet_wrap(~TF_name, ncol = 1) + 
                scale_color_gradient(breaks = c(2, 4, 6, 8, 10), 
                    labels=c(2,4,6,8,"+10"),
                    low = "black",  # Light gray for low values
                    high = "gray90",   # Dark gray for high values
                    name = 'TF rank'
                  ) +
                theme_classic() + labs(y = 'fractional binding divergence') +
                scale_x_discrete(labels = rev(labels)) + scale_y_continuous(breaks = c(0,0.5,1), limits = c(0,1), labels = c('0', '0.5', '1')) +
              theme(axis.text.x = element_text(angle = 45, hjust = 1),
                                  axis.line.x = element_blank(),
                                  axis.ticks.x = element_blank(),
                                  axis.title.x =element_blank())
    return(div.plot)
}

plot_aff_trace <- function(msa_tfbs_df, func_tfbs_df, query_name, ref_name){
    aff_trace_query <- msa_tfbs_df %>% filter(species == query_name) %>% 
              rowwise() %>%
                  mutate(pos = list(msa_start:msa_end)) %>% # Expand positions into a list column
                  unnest(pos) %>%                          # Flatten the list column
                  dplyr::select(species, pos, norm_affinity, TF_name) %>%
                  group_by(species, pos, TF_name) %>%                        # Group by position
                  summarise(
                    rolling_avg_norm_aff = mean(norm_affinity, na.rm = TRUE) # Rolling average
              )
    aff_trace_ref <- msa_tfbs_df %>% filter(species == ref_name) %>%
              rowwise() %>%
                  mutate(pos = list(msa_start:msa_end)) %>% # Expand positions into a list column
                  unnest(pos) %>%                          # Flatten the list column
                  dplyr::select(species, pos, norm_affinity, TF_name) %>%
                  group_by(species, pos, TF_name) %>%                        # Group by position
                  summarise(
                    rolling_avg_norm_aff = mean(norm_affinity, na.rm = TRUE) # Rolling average
              )
    aff_trace_query$TF_name <- factor(aff_trace_query$TF_name, levels = final_TFs)
    aff_trace_ref$TF_name <- factor(aff_trace_ref$TF_name, levels = final_TFs)
    func_tfbs_df$TF_name <- factor(func_tfbs_df$TF_name, levels = final_TFs)
    func_tfbs_df$species <- ref_name
    
    aff_plot <- ggplot() + 
                  geom_step(data = aff_trace_query, aes(x = pos, y = rolling_avg_norm_aff, color = species)) + 
                  geom_step(data = aff_trace_ref, aes(x = pos, y = rolling_avg_norm_aff, color = species)) + 
                  geom_rect(data = func_tfbs_df, 
                        aes(
                        xmin = msa_start,
                        xmax = msa_end), ymin = -Inf, ymax = Inf,
                        colour = "black", fill = NA, linetype = "dashed", alpha = 0.6) + 
                  geom_hline(yintercept = 0.05, color = 'black', linetype = 'dashed') + 
                  facet_wrap(~TF_name, ncol = 1) + 
                  theme_classic() +
                  scale_color_manual(
                    values = setNames(c("gainsboro", "#cc0000"), c(query_name, ref_name)),  # Define colors for the legend
                    name = "",  # Legend title
                    guide = guide_legend(override.aes = list(size = 5))
                  ) + scale_y_log10(breaks = c(10^-5, 10^-3, 10^-1),
                           labels = scales::trans_format("log10", scales::math_format(10^.x))
                         ) + annotation_logticks(sides = 'l') + 
                  coord_cartesian(xlim=c(1,max(msa_tfbs_df$msa_end))) + 
                  labs(x = 'msa bp position', y = 'norm. binding affinity') #+ theme(legend.position = 'bottom')
}

pairwise_functional_TFBS <- function(msa_tfbs_df, species.list, labels){
    library(pheatmap)
    library(viridis)
    
    # Generate a viridis color palette
    color_palette <- viridis(100)

    # Generate breaks to span the range 0 to 1
    breaks <- seq(0, 1, length.out = 101)
    species.list <- species.list[species.list%in%unique(msa_tfbs_df$species)]
    func_tfbs_mat <- msa_tfbs_df[complete.cases(msa_tfbs_df$rank_TF),]
    func_tfbs_cast <- dcast(func_tfbs_mat, TF_name2 ~ species, value.var = "norm_affinity")
    # Replace NA values with 0
    func_tfbs_cast[is.na(func_tfbs_cast)] <- 0
    rownames(func_tfbs_cast) <- func_tfbs_cast$TF_name2
    func_tfbs_cast <- func_tfbs_cast[,-1]
    func_tfbs_cor <- cor(func_tfbs_cast, method = 'spearman')
    func_tfbs_cor <- func_tfbs_cor[rev(species.list), rev(species.list)]
    # Generate the pheatmap silently
    pheatmap_obj <- pheatmap(
        func_tfbs_cor, 
        cluster_rows = FALSE, 
        cluster_cols = FALSE, 
        labels_row = rev(labels),
        labels_col = rev(labels),
        color = color_palette, 
        breaks = breaks,
        angle_col = 45,
        silent = TRUE
    )
    
    # Convert the pheatmap object into a ggplot object
    gg_heatmap <- as.ggplot(pheatmap_obj$gtable)
    
    return(gg_heatmap)
}


