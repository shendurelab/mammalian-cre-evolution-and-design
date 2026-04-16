library(tidyverse)
library(castor)
suppressMessages(suppressWarnings(library(dplyr)))
library(reshape2)
library(ape)
library(phytools)
require(grid)
library(ggtree)
library(tidytree)
library(ggtreeExtra)
library(stringr)
library(ggplot2)
library(patchwork)
library(ggnewscale)
source("/net/shendure/vol8/projects/tli/PYS2_oCRE_BC_libs_20240904/mafft_aligning_functions.R")

select.images <- Sys.glob('/net/shendure/vol8/projects/tli/PYS2_oCRE_BC_libs_20240904/imgs/*.png')
select.images <- data.frame(uid = select.images,
                           label = c("Psammomys_obesus","fullTreeAnc237","fullTreeAnc223","Felis_catus","Canis_lupus_familiaris","fullTreeAnc115",'fullTreeAnc239',
                                    "Equus_caballus","Homo_sapiens","Mus_musculus","fullTreeAnc110","Rattus_norvegicus","fullTreeAnc68"),
                           name = c('Psammomys_obesus',"Laurasiatheria","Carnivora","Felis_catus","Canis_lupus_familiaris","Euarchontoglires",
                                   "First_Ancestor","Equus_caballus","Homo_sapiens","Mus_musculus","Primates","Rattus_norvegicus","Rodentia"))

order.x <- c('PRIMATES','RODENTIA','EULIPOTYPHLA','CHIROPTERA',
                             'PERISSODACTYLA','CARNIVORA','CETARTIODACTYLA')

### make basic tree
order.colors <- c("steelblue","#cc0000","#ffc425",'#905aa3','#ae5a41','#ff7f50','#1c9e48')
names(order.colors) <- order.x


make.nodepath <- function(tree, pair){
    all.labels <- c(tree$tip.label, tree$node.label)
    tip.numbers <- match(pair, all.labels)
    pathway <- nodepath(tree, tip.numbers[1], tip.numbers[2])
    return(all.labels[pathway])
}

plot_dms_tfbs_act <- function(msa_df_dms_TFBS, msa_dms_df, species_list,
                          CRE_oi, shaps_TFs, labels = NULL, tile_height = 0.3, tile_width = 0.4){
    
    msa_df_dms_TFBS$species <- factor(msa_df_dms_TFBS$species, levels = species_list)
    if('rollmean_median'%in%colnames(msa_dms_df)){
        msa_dms_df <- msa_dms_df %>% select(center_pos = curpos, species, log2FC_MPRA = rollmean_median)
    } else{
        msa_dms_df <- msa_dms_df %>% select(center_pos = curpos, species, log2FC_MPRA = log2FC_MPRA)
    }
    msa_dms_df$log2FC_MPRA <- as.numeric(msa_dms_df$log2FC_MPRA)
    msa_dms_df$species <- factor(msa_dms_df$species, levels = species_list)
    
    
    gg_tfbs <- ggplot() +
              geom_hline(data = msa_df_dms_TFBS, aes(yintercept = species), color = "black", alpha = 0.5) +  # Add horizontal line
              geom_tile(data = msa_dms_df, aes(x= center_pos, y = as.numeric(species) -0.5, fill = log2FC_MPRA), width = tile_width, height = tile_height) + 
              labs(x = 'MSA position', fill = "log2(variant/WT)", shape = "TFBS", size = "affinity binding") +
              scale_fill_gradient2(low = "dodgerblue2",  
                                 mid = "white",  
                                 high = "firebrick1",  
                                 midpoint = 0,  
                                 space = "Lab",  
                                 na.value = "gainsboro",  
                                 limits=c(-1,1),  breaks = c(-1,0,1), 
                                 oob=scales::squish) + 
                new_scale("fill") +  # Needed to add a new color scale for a different layer
            geom_star(data = msa_df_dms_TFBS, aes(x =center_pos, y = species, starshape=TF_name, fill = avg_log2FC_MPRA, size = norm_affinity),
                                color = 'black', position=position_jitter(h=0.1)) + 
                scale_fill_gradient2(low = "#008080", mid = "white", high = "orange", 
                       midpoint = 0, space = "Lab", na.value = "gainsboro", name = "TFBS log2(variant/WT)",
                       limits = c(-0.1, 0.1), breaks = c(-0.1, 0, 0.1), oob = scales::squish) +
                        theme_cowplot() + scale_starshape_manual(values = shaps_TFs, name = "TFBS", 
                                                                 guide = guide_legend(override.aes = list(size = 5, fill = 'black'))) +
                        scale_size_continuous(limits = c(0, 5), breaks = c(0.5, 1, 3, 5)) + # Affinity binding size scaling
                        guides(size = guide_legend(override.aes = list(shape = 21, fill = "black"))) 
    if(is.null(labels)){
        gg_tfbs <- gg_tfbs + theme(legend.position = 'bottom',legend.box="vertical", legend.margin=margin(),
                              text = element_text(size = 16),
                              axis.line.y = element_blank(),
                              axis.ticks.y = element_blank(),
                              axis.title.y=element_blank())
    }else{
        gg_tfbs <- gg_tfbs + scale_y_discrete(labels = labels) + 
                        theme(legend.position = 'bottom',legend.box="horizontal", legend.margin=margin(),
                              text = element_text(size = 16),
                              axis.text.y = element_markdown(color = "black", size = 10),
                              axis.line.y = element_blank(),
                              axis.ticks.y = element_blank(),
                              axis.title.y=element_blank()) 
    }
    
    species_act <- oCRE_mean_act %>% filter(full_CRE_id == CRE_oi) %>% filter(species%in%species_list)
    species_act$species <- factor(species_act$species, levels = species_list)
    cre_oi_wt_act <- df_mouse_tiles_mpra %>% filter(full_CRE_id == CRE_oi) %>% pull(mean_MPRA_act)
    cre_oi_wt_sd <- df_mouse_tiles_mpra %>% filter(full_CRE_id == CRE_oi) %>% pull(sd_MPRA_act)
    minP_oi_act <- control_mean_act %>% filter(CRE_id == 'minP') %>% pull(mean_MPRA_act)
    minP_oi_sd <- control_mean_act %>% filter(CRE_id == 'minP') %>% pull(sd_MPRA_act)
    
    mpra_bar <- ggplot(species_act, aes(x=mean_MPRA_act, y=species, color=mean_MPRA_act)) + 
        geom_point(size = 2) + 
        geom_errorbar(aes(xmin=mean_MPRA_act-sd_MPRA_act, xmax=mean_MPRA_act+sd_MPRA_act), width=.2,
                     position=position_dodge(0.05)) +
#         geom_rect(xmin = cre_oi_wt_act-cre_oi_wt_sd, xmax = cre_oi_wt_act+cre_oi_wt_sd, ymin = -Inf, ymax = Inf,
#            alpha = .2,fill = "#008000") + 
#             geom_rect(xmin = minP_oi_act - minP_oi_sd, xmax = minP_oi_act + minP_oi_sd, ymin = -Inf, ymax = Inf,
#            alpha = .1,fill = "blue")  +
            geom_vline(xintercept = minP_oi_act ,color = "blue", linetype = 'dashed')  +
        geom_vline(xintercept = cre_oi_wt_act ,color = "#008000", linetype = 'dashed') +
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
    layout <- '
        AAAAAAABB
        AAAAAAABB'
    p <- gg_tfbs + mpra_bar + plot_layout(design = layout)
    
#     plot_grid(gg_tfbs, mpra_bar, align = 'h', rel_widths = c(1, 0.3))
    return(p)
}


plot_traj <- function(traj_act.df, df_mouse_tiles_mpra, minP, CRE_oi){
    shape_ori <- c(18, 16)
    names(shape_ori) <- c('ancestral','extant')
    
    order.x <- c('PRIMATES','RODENTIA','EULIPOTYPHLA','CHIROPTERA',
                             'PERISSODACTYLA','CARNIVORA','CETARTIODACTYLA')
    
    select_traj.act <- traj_act.df %>% filter(Order%in%order.x) %>% filter(full_CRE_id == CRE_oi)
    select_traj.act$Order <- factor(select_traj.act$Order, levels = order.x)
    
    cre_control <- df_mouse_tiles_mpra %>% filter(full_CRE_id == CRE_oi) %>% pull(mean_MPRA_act)
    
    dist.traj <- ggplot(select_traj.act[!duplicated(select_traj.act),], aes(x = distance, y = mean_MPRA_act, color = Order, group = Order, shape = origin)) + 
#         geom_errorbar(aes(ymin=mean_MPRA_act-sd_MPRA_act, ymax=mean_MPRA_act+sd_MPRA_act, color = Order), alpha = 0.) + 
        geom_point_rast(alpha = 0.6) + geom_smooth() + 
        scale_color_manual(values = order.colors, guide = 'none') + 
        scale_shape_manual(values = shape_ori, guide = guide_legend(override.aes = list(size = 5, fill = 'black'))) + 
        geom_hline(yintercept = minP  ,color = "black", linetype = 'dashed') + 
            geom_hline(yintercept = cre_control, color = "#008000", linetype = 'dashed') + 
        scale_x_continuous(breaks = c(0, 0.25, 0.5)) + 
        scale_y_log10(breaks = c(10^-1, 10^0, 10^1), limits = c(0.1, 15),
               labels = scales::trans_format("log10", scales::math_format(10^.x))) + annotation_logticks(sides = 'l') + facet_wrap(~Order, nrow = 1) +coord_cartesian() +  
        theme_classic() + labs(x = 'evolutionary distance', y = 'MPRA act.') + theme(panel.spacing = unit(1, "lines"))

    id.traj <- ggplot(select_traj.act[!duplicated(select_traj.act),], aes(x = 100-id, y = mean_MPRA_act, color = Order, group = Order, shape = origin)) + 
#         geom_errorbar(aes(ymin=mean_MPRA_act-sd_MPRA_act, ymax=mean_MPRA_act+sd_MPRA_act, color = Order), alpha = 0.) + 
        geom_point_rast(alpha = 0.6) + geom_smooth() + 
        scale_color_manual(values = order.colors, guide = 'none') + 
        scale_shape_manual(values = shape_ori, guide = guide_legend(override.aes = list(size = 5, fill = 'black'))) + 
        geom_hline(yintercept = minP  ,color = "black", linetype = 'dashed') + 
            geom_hline(yintercept = cre_control, color = "#008000", linetype = 'dashed') + 
            scale_x_continuous(breaks = c(0, 20, 40), limits = c(0,45)) + 
        scale_y_log10(breaks = c(10^-1, 10^0, 10^1), limits = c(0.1, 15),
               labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
        annotation_logticks(sides = 'l') + facet_wrap(~Order, nrow = 1) + coord_cartesian() + 
        theme_classic() + labs(x = 'Sequence difference (%, vs. First Mammal)', y = 'MPRA act.') + theme(panel.spacing = unit(1, "lines"))

    order.traj <- dist.traj + id.traj +
      plot_layout(ncol = 1, guides = 'collect') & theme(legend.position = 'top')

    return(order.traj)
}


generate_trajectory_df_2_plot <- function(species_list, ref_oi, fasta_file, CRE_oi, pic_size=25){
    shaps_TFs <- c(29,13,15,28,11)
    names(shaps_TFs) <- final_TFs
    if(ref_oi == 'Mus_musculus'){
        msa_df <- make_alignment_dataframe_from_fasta(fasta_file, species_list, "Mus_musculus",CRE_oi)
        msa_dms_df <- left_join(msa_df, df_DMS_all_CREs, by =c("mut",'ori_CRE_id','ref_pos'))
        msa_dms_df <- msa_dms_df %>% mutate(log2FC_MPRA = case_when( ### replace NA in log2FC_MPRA to 0 since its WT
                                mut_type == 'Match' ~ 0,
                                TRUE ~ log2FC_MPRA))
        msa_dms_df$log2FC_MPRA[msa_dms_df$log2FC_MPRA == -Inf] <- NA ### replace missing DEL expression with NA 

        msa_df_TFBS <- make_cre_tfbs_alignment_dataframe(msa_df, final_cactus_TFBS, CRE_oi)
        msa_df_dms_TFBS <- data.frame()

        for(species_oi in unique(msa_df_TFBS$species)){
            tmp_msa_df_TFBS <- msa_df_TFBS %>% filter(species == species_oi)
            tmp_msa_dms_df <- msa_dms_df %>% filter(species == species_oi)
            tmp_msa_df_dms_TFBS <- tmp_msa_df_TFBS %>%
              rowwise() %>%
              mutate(avg_log2FC_MPRA = mean(
                tmp_msa_dms_df %>% 
                  filter(curpos >= msa_start & curpos <= msa_end) %>%
                  pull(log2FC_MPRA), na.rm = TRUE)) %>%
              ungroup()  # ungroup to ensure no row-wise grouping remains
            msa_df_dms_TFBS <- bind_rows(msa_df_dms_TFBS, tmp_msa_df_dms_TFBS)
        }
        labels <- c()

        for(species_oi in species_list){
            if(species_oi%in%select.images$label){
                img_path <- select.images %>% filter(label == species_oi) %>% pull(uid)
                name_oi <- select.images %>% filter(label == species_oi) %>% pull(name)
                labels <- c(labels, paste0("<img src='", img_path,  "' width='",pic_size,"' /><br>*", name_oi,"*"))
            } else{
                labels <- c(labels, NA)
            }
        }

        p <- plot_dms_tfbs_act(msa_df_dms_TFBS, msa_dms_df, species_list, CRE_oi, shaps_TFs, labels)
    } else{
        msa_df <- make_alignment_dataframe_from_fasta(fasta_file, c(species_list,"Mus_musculus"), "Mus_musculus",CRE_oi)
        msa_df_TFBS <- make_cre_tfbs_alignment_dataframe(msa_df, final_cactus_TFBS, CRE_oi)

        msa_dms_df <- left_join(msa_df, df_DMS_all_CREs, by =c("mut",'ori_CRE_id','ref_pos'))
        msa_dms_df <- msa_dms_df %>% mutate(log2FC_MPRA = case_when( ### replace NA in log2FC_MPRA to 0 since its WT
                                mut_type == 'Match' ~ 0,
                                TRUE ~ log2FC_MPRA))
        msa_dms_df$log2FC_MPRA[msa_dms_df$log2FC_MPRA == -Inf] <- NA ### replace missing DEL expression with NA 
        ref_df_dms_TFBS <- data.frame()

        for(species_oi in unique(msa_df_TFBS$species)){
            tmp_msa_df_TFBS <- msa_df_TFBS %>% filter(species == species_oi)
            tmp_msa_dms_df <- msa_dms_df %>% filter(species == species_oi)
            tmp_msa_df_dms_TFBS <- tmp_msa_df_TFBS %>%
              rowwise() %>%
              mutate(avg_log2FC_MPRA = mean(
                tmp_msa_dms_df %>% 
                  filter(curpos >= msa_start & curpos <= msa_end) %>%
                  pull(log2FC_MPRA), na.rm = TRUE)) %>%
              ungroup()  # ungroup to ensure no row-wise grouping remains
            ref_df_dms_TFBS <- bind_rows(ref_df_dms_TFBS, tmp_msa_df_dms_TFBS)
        }

        msa_df <- make_alignment_dataframe_from_fasta(fasta_file, species_list, ref_oi,CRE_oi)
        msa_df_TFBS <- make_cre_tfbs_alignment_dataframe(msa_df, final_cactus_TFBS, CRE_oi)
        mm <-match(msa_df_TFBS$species, names(nodes.order))
        msa_df_TFBS$order <- nodes.order[mm]
        mm <-match(msa_df_TFBS$species, names(dist.mat))
        msa_df_TFBS$distance <- dist.mat[mm]
        msa_df_TFBS$id <- id.mat[mm]
        ref_dms_df <- msa_dms_df %>% select(mut, species, ori_CRE_id, query_pos,
                                   log2FC_MPRA, sd_log2FC_MPRA, n_BC, MPRA_act, sd_MPRA_act, min_act,
                                   max_act, wilcox_BH)
        ref_dms_df <- left_join(msa_df, ref_dms_df, by =c("mut",'species','ori_CRE_id','query_pos'))

        labels <- c()

        for(species_oi in species_list){
            if(species_oi%in%select.images$label){
                img_path <- select.images %>% filter(label == species_oi) %>% pull(uid)
                name_oi <- select.images %>% filter(label == species_oi) %>% pull(name)
                labels <- c(labels, paste0("<img src='", img_path,  "' width='",pic_size,"' /><br>*", name_oi,"*"))
            } else{
                labels <- c(labels, NA)
            }
        }

        keep.cols <- c("start_pos", "end_pos","kmer_seq", 
                       "affinity","TF_name","substring_l","ori_CRE_id", "species", "norm_affinity","avg_log2FC_MPRA")

        msa_df_dms_TFBS <- left_join(msa_df_TFBS, ref_df_dms_TFBS[,keep.cols], by =c("start_pos", "end_pos","kmer_seq",
                                                                               "affinity","TF_name","substring_l",
                                                                              "ori_CRE_id", "species", "norm_affinity")) 
        p <- plot_dms_tfbs_act(msa_df_dms_TFBS, ref_dms_df, nodes.keep, CRE_oi, shaps_TFs, labels)
    }
    return(p)
}