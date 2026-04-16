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
source("utils/load_figure1_2_data.R")

# Helper function to safely pull a value or return 0
pull_or_zero <- function(df, col) {
  if (nrow(df) == 0) return(0)
  val <- df[[col]]
  if (length(val) == 0) return(0)
  return(val)
}

CRE_oi <- "Epas1_chr17_10063"
big.tree <- as_tibble(tree)
big.tree <- groupClade(big.tree, c(293, 247,353,383,429,346,459))
anc.tips <- big.tree %>% filter(node%in%c(293, 247,353,383,429,346,459)) %>% pull(label)
anc.labels <- data.frame(species = anc.tips, anc_name = c('PRIMATES','RODENTIA','PERISSODACTYLA','CHIROPTERA',
                                                         'CETARTIODACTYLA','CARNIVORA','EULIPOTYPHLA'), group = c(2,1,3,4,5,6,7))

common_name.df <- data.frame(species = c('Mus_musculus','Giraffa_tippelskirchi','Homo_sapiens','Eulemur_fulvus','Lycaon_pictus','Psammomys_obesus'),
                            common_name = c('Mouse',"Giraffe",'Human','Brown_lemur','African_dog','Sand_rat'))

group_no <- 1
epas1_avg_lineage_mpra_list <- list()
stepwise_summary_list <- list()
for(species_oi in c('Mus_musculus','Lycaon_pictus','Psammomys_obesus')){
    ### plot Epas1 transitions
    order_oi <- cactus.annotation %>% filter(Species == species_oi) %>% pull(Order)
    anc_oi <- anc.labels %>% filter(anc_name == order_oi) %>% pull(species)
    common_name2 <- common_name.df %>% filter(species == species_oi) %>% pull(common_name)
    species.pair <- c(species_oi,'fullTreeAnc239')
    keep.nodes <- make.nodepath(tree, species.pair)
    ### filter nodes if in activity
    keep.nodes <- keep.nodes[keep.nodes%in%unique(oCRE_mean_act$species)]

    start_node <- "fullTreeAnc239"
    order_species <- rev(keep.nodes)  # traverse in your chosen order

    ### get replicate
    oCRE_replicate_act <- df_mpra_act %>% filter(class == 'oCRE')
    oCRE_replicate_act <- oCRE_replicate_act %>% separate(CRE_id, c('tile_id','something'), sep = "_20240416") %>% dplyr::select(-something)
    oCRE_replicate_act <- oCRE_replicate_act %>% separate(tile_id, c('full_CRE_id','species'), sep = "__")
    mm <- match(oCRE_replicate_act$species, cactus.annotation$Species)
    oCRE_replicate_act$common_name <- cactus.annotation$`Common Name`[mm]
    oCRE_replicate_act <- left_join(oCRE_replicate_act, max.dist, by =c("full_CRE_id","species"))

    # base table (ordered once)
    temp_avg_lineage_mpra <- oCRE_mean_act %>%
      filter(full_CRE_id == CRE_oi, species %in% keep.nodes) %>%
      mutate(
        species2  = case_when(
                species == species_oi ~ common_name2, 
                TRUE ~ str_replace(species, "fullTree","")
        ),
        marker = case_when(
          species%in%c('fullTreeAnc68','fullTreeAnc59','fullTreeAnc222')  ~ "Transition",
          species == species_oi                 ~ 'Extant',
          TRUE                                       ~ "Other"
        ),
        species2 = factor(species2, 
                          levels = c(
            ifelse(order_species == species_oi, common_name2,str_replace(order_species, "fullTree",""))
          )),
          lineage = species_oi,
            group = group_no
      ) %>%
      arrange(species)
    epas1_avg_lineage_mpra_list[[species_oi]] <- temp_avg_lineage_mpra

    epas1_lineage_mpra <- oCRE_replicate_act %>%
      filter(full_CRE_id == CRE_oi, species %in% keep.nodes) %>%
      mutate(
        species2  = case_when(
                species == species_oi ~ common_name2, 
                TRUE ~ str_replace(species, "fullTree","")
        ),
        marker = case_when(
          species%in%c('fullTreeAnc68','fullTreeAnc59','fullTreeAnc222')  ~ "Transition",
          species2 == common_name2                 ~ 'Extant',
          TRUE                                       ~ "Other"
        ),
        species2 = factor(species2, levels = c(
        ifelse(order_species == species_oi, common_name2, str_replace(order_species, "fullTree",""))
      )),
          lineage = species_oi
      ) %>%
      arrange(species)
    
    # constants
    start_act_df <- epas1_lineage_mpra %>%
      filter(species == start_node) %>%
      dplyr::select(biol_rep, start_act = MPRA_act)

    # Get mouse wild-type activity per replicate
    wt_act_df <- epas1_lineage_mpra %>%
      filter(species == species_oi) %>%
      dplyr::select(biol_rep, wt_act = MPRA_act)

    # stepwise deltas & fold-changes
    epas1_lineage_mpra_steps <- epas1_lineage_mpra %>%
      arrange(biol_rep, match(species, order_species)) %>%
      group_by(biol_rep) %>%
      mutate(
        prev_act     = lag(MPRA_act, default = NA),
        delta_steps  = MPRA_act - prev_act,
        fc_steps     = MPRA_act / prev_act,
        log2FC_steps = log2(fc_steps)
      ) %>%
      ungroup() %>%
      left_join(start_act_df, by = "biol_rep") %>%
      left_join(wt_act_df, by = "biol_rep") %>%
      mutate(
        delta_steps  = if_else(species == start_node, 0, delta_steps),
        fc_steps     = if_else(species == start_node, 1, fc_steps),
        log2FC_steps = if_else(species == start_node, 0, log2FC_steps),
        frac_delta   = delta_steps / (wt_act - start_act),
        species2     = factor(species2, levels = c(
            ifelse(order_species == species_oi, common_name2,str_replace(order_species, "fullTree",""))
          ))
      )

    stepwise_tmp <- epas1_lineage_mpra_steps %>%
      group_by(species) %>%
      summarise(
        avg_delta_steps     = mean(delta_steps, na.rm = TRUE),
        sd_delta_steps      = sd(delta_steps, na.rm = TRUE),
        avg_frac_delta      = mean(frac_delta, na.rm = TRUE),
        sd_frac_delta       = sd(frac_delta, na.rm = TRUE),
        avg_fc_steps        = mean(fc_steps, na.rm = TRUE),
        sd_fc_steps         = sd(fc_steps, na.rm = TRUE),
        avg_log2FC_steps    = mean(log2FC_steps, na.rm = TRUE),
        sd_log2FC_steps     = sd(log2FC_steps, na.rm = TRUE),
        p_delta_steps       = ifelse(n() > 1, t.test(delta_steps)$p.value, NA_real_),
        p_frac_delta        = ifelse(n() > 1, t.test(frac_delta)$p.value, NA_real_),
        p_log2FC_steps      = ifelse(n() > 1, t.test(log2FC_steps)$p.value, NA_real_)
      ) %>%
      ungroup()
    stepwise_tmp$p_adj_log2FC_steps <- p.adjust(stepwise_tmp$p_log2FC_steps, method = 'fdr')
    stepwise_tmp <- stepwise_tmp %>% 
        mutate(
        species2  = ifelse(species == species_oi, common_name2, str_replace(species, "fullTree","")),
        marker = case_when(
          species%in%c('fullTreeAnc68','fullTreeAnc59','fullTreeAnc222')  ~ "Transition",
          species == species_oi                 ~ 'Extant',
          TRUE                                       ~ "Other"
        ),
        species2 = factor(species2, 
                          levels = c(
            ifelse(order_species == species_oi, common_name2,str_replace(order_species, "fullTree",""))
          )),
          lineage = species_oi,
            group = group_no
      )
    
    stepwise_summary_list[[species_oi]] <- stepwise_tmp
    group_no <- group_no + 1
}

### plot Epas1 transitions
species.pair <- c("Mus_musculus",'fullTreeAnc239')
mouse.nodes <- make.nodepath(tree, species.pair)
### filter nodes if in activity
mouse.nodes <- mouse.nodes[mouse.nodes%in%unique(oCRE_mean_act$species)]

### Epas1 fasta file ###
CRE_oi <- "Epas1_chr17_10063"
fasta_file <- '/net/shendure/vol10/projects/tli/nobackup/EB_dev_project/oCREs_v2/nobackup/reoriented_fasta/Epas1_chr17_10063_oCRE.fasta'


### plot Epas1 MSA for mouse lineage
mouse_DMS_MPRA$mut <- toupper(mouse_DMS_MPRA$mut)
species.pair <- c("Mus_musculus",'fullTreeAnc239')
keep.tips <- make.nodepath(tree, species.pair)
epas1_mouse_msa_df <- make_alignment_dataframe(fasta_file, keep.tips, 'Mus_musculus',CRE_oi)
mouse_msa_df <- epas1_mouse_msa_df %>% filter(species == 'Mus_musculus')

epas1_mouse_msa_list <- make_msa_tfbs_traj_df(fasta_file, rev(keep.tips), 'Mus_musculus', CRE_oi, cactus.tfbs, mouse_func_tfbs)
labels <- collect_labels(epas1_mouse_msa_list$msa_df_TFBS, side = T)
epas1_mouse_func_tfbs <- mouse_func_tfbs %>% filter(CRE == CRE_oi)
mm.start <- match(epas1_mouse_func_tfbs$TFBS_start, mouse_msa_df$ref_pos)
mm.end <- match(epas1_mouse_func_tfbs$TFBS_end, mouse_msa_df$ref_pos)
epas1_mouse_func_tfbs$msa_start <- mouse_msa_df$curpos[mm.start]
epas1_mouse_func_tfbs$msa_end <- mouse_msa_df$curpos[mm.end]
epas1_mouse_msa_df <- left_join(epas1_mouse_msa_df, mouse_DMS_MPRA %>% filter(CRE == CRE_oi) %>%
                          dplyr::select(ref_pos = pos, mut, log2FC_MPRA, sd_log2FC_MPRA, MPRA_act, sd_MPRA_act, min_act,max_act, wilcox_BH),
                         by = c("ref_pos","mut"))

epas1_mouse_msa_df <- epas1_mouse_msa_df %>% mutate(mut2 = case_when(
                        mut == 'DEL' ~ '-',
                        TRUE ~ mut))
epas1_mouse_msa_df <- epas1_mouse_msa_df %>% mutate(species2 = ifelse(species == 'Mus_musculus', 'Mouse', str_replace(species, "fullTree","")))
epas1_mouse_msa_df$species2 <- factor(epas1_mouse_msa_df$species2, levels = c(
    ifelse(keep.tips == "Mus_musculus", "Mouse", str_replace(keep.tips, "fullTree",""))
  ))

epas1_mouse_msa_df <- epas1_mouse_msa_df %>%
  mutate(fill_type = case_when(
    mut2 == "-" & is.na(log2FC_MPRA) ~ "black",
    mut == ref_nuc ~ "grey",
    TRUE ~ "gradient"
  ))

fixed_epas1_mouse_msa_df <- data.frame()
### correct for missing data
keep.tips <- keep.tips[keep.tips%in%unique(epas1_mouse_msa_df$species)]

for(pos in seq(1:max(epas1_mouse_msa_df$curpos))){
    tmp_msa_df <- epas1_mouse_msa_df %>% filter(curpos == pos)
    all_mut_types <- unique(tmp_msa_df$mut_type)
    if(any(all_mut_types %in% c('Mismatch', 'DEL', 'INS'))){
        current_log2fc <- 0
        have_not_found_match <- TRUE
        log_change <- c()
        match2 <- c()
        for(species_oi in rev(keep.tips)){
            tmp_species <- tmp_msa_df %>% filter(species == species_oi)
            if(have_not_found_match){
                check_match <- tmp_species %>% filter(mut2 == ref_nuc)
                current_muttype <- tmp_species %>% pull(mut_type)
                if(nrow(check_match) > 0){
                    log_change <- c(log_change, -1*current_log2fc)
                    have_not_found_match <- FALSE
                    match2 <- c(match2, 'First Match')
                }else{
                    current_log2fc <- tmp_species %>% pull(log2FC_MPRA)
                    log_change <- c(log_change, 0)
                    match2 <- c(match2, current_muttype)
                }
            }else{
                log_change <- c(log_change, 0)
                match2 <- c(match2, 'Match')
            }
        }
        names(log_change) <- rev(keep.tips)
        names(match2) <- rev(keep.tips)
        mm <- match(tmp_msa_df$species, names(log_change))
        tmp_msa_df$log_change <- log_change[mm]
        tmp_msa_df$match2 <- match2[mm]
        fixed_epas1_mouse_msa_df <- bind_rows(fixed_epas1_mouse_msa_df, tmp_msa_df)
    }else{
        tmp_msa_df$log_change <- 0
        tmp_msa_df$match2 <- 'Match'
        fixed_epas1_mouse_msa_df <- bind_rows(fixed_epas1_mouse_msa_df, tmp_msa_df)
    }
}

fixed_epas1_mouse_msa_df <- fixed_epas1_mouse_msa_df %>% mutate(log_change = case_when(
                                        mut_type%in%c('DEL','INS') ~ NA,
                                        TRUE ~ log_change))
fixed_epas1_mouse_msa_df <- fixed_epas1_mouse_msa_df %>% 
            mutate(species2 = ifelse(species == 'Mus_musculus', 'Mouse', str_replace(species, "fullTree","")))
fixed_epas1_mouse_msa_df$species2 <- factor(fixed_epas1_mouse_msa_df$species2, levels = c(
    ifelse(keep.tips == "Mus_musculus", "Mouse", str_replace(keep.tips, "fullTree",""))
  ))


#### calculate number of SNPs that become the original mouse sequence at each step and mark which original mouse SNP confers a TFBS
### calculate avg log2FC for each mutation gained
all_mouse_tfbs_positions <- unlist(mapply(seq, epas1_mouse_func_tfbs$msa_start, epas1_mouse_func_tfbs$msa_end))
epas1_current_snps <- epas1_mouse_msa_df %>% filter(species == 'fullTreeAnc239') %>% filter(mut2 != ref_nuc) %>% pull(curpos)
gain_mouse_snps <- fixed_epas1_mouse_msa_df %>% filter(species == 'fullTreeAnc239') %>% filter(mut2 == ref_nuc) %>%
                            arrange(curpos) %>% mutate(group_id = cumsum(c(TRUE, diff(curpos) > 1)),
                                                      in_TFBS = case_when(
                                                          curpos%in%all_mouse_tfbs_positions ~ TRUE,
                                                          TRUE ~ FALSE))
fc_calc <- gain_mouse_snps %>% group_by(in_TFBS) %>% 
                        summarise(avg_log_change = mean(log_change, na.rm = TRUE), sd_log_change = sd(log_change, na.rm = TRUE))
epas1_all_fc <- gain_mouse_snps

epas1_snp_counts <- gain_mouse_snps %>%
                  group_by(species, group_id) %>%
                  summarise(
                    num_snps = 1,                     # each contiguous run counts as 1
                    tfbs_count = any(in_TFBS),        # TRUE if any in the run is in TFBS
                    .groups = "drop_last"
                  ) %>%
                  summarise(
                    n_snp_count = sum(num_snps),
                    total_tfbs_snps = sum(tfbs_count),# same as above unless you want total TRUEs
                    .groups = "drop"
                  ) %>% mutate(total_nontfbs_snps = n_snp_count - total_tfbs_snps)
epas1_avg_fc <- data.frame(
                                species         = 'fullTreeAnc239',
                                avg_tfbs_FC     = pull_or_zero(fc_calc %>% filter(in_TFBS), "avg_log_change"),
                                avg_nontfbs_FC  = pull_or_zero(fc_calc %>% filter(!in_TFBS), "avg_log_change"),
                                sd_tfbs_FC      = pull_or_zero(fc_calc %>% filter(in_TFBS), "sd_log_change"),
                                sd_nontfbs_FC   = pull_or_zero(fc_calc %>% filter(!in_TFBS), "sd_log_change")
                              )
for(species_oi in rev(keep.tips)){
    if(species_oi == 'fullTreeAnc239'){
    } else{
        temp_msa_df <- fixed_epas1_mouse_msa_df %>% filter(species == species_oi)
        current_snps <- temp_msa_df %>% filter(mut2 != ref_nuc) %>% pull(curpos)
        if(length(epas1_current_snps) > 0){
            gain_mouse_snps <- temp_msa_df %>% filter(curpos%in%epas1_current_snps) %>% filter(mut2 == ref_nuc) %>%
                            arrange(curpos) %>% mutate(group_id = cumsum(c(TRUE, diff(curpos) > 1)),
                                                      in_TFBS = case_when(
                                                          curpos%in%all_mouse_tfbs_positions ~ TRUE,
                                                          TRUE ~ FALSE))
            if (nrow(gain_mouse_snps) > 0) {
              fc_calc <- gain_mouse_snps %>% group_by(in_TFBS) %>% 
                        summarise(avg_log_change = mean(log_change, na.rm = TRUE), sd_log_change = sd(log_change, na.rm = TRUE))
              epas1_all_fc <- bind_rows(epas1_all_fc, gain_mouse_snps)
              snp_counts <- gain_mouse_snps %>%
                  group_by(species, group_id) %>%
                  summarise(
                    num_snps = 1,                     # each contiguous run counts as 1
                    tfbs_count = any(in_TFBS),        # TRUE if any in the run is in TFBS
                    .groups = "drop_last"
                  ) %>%
                  summarise(
                    num_snps = sum(num_snps),
                    total_tfbs_snps = sum(tfbs_count),# same as above unless you want total TRUEs
                    .groups = "drop"
                  ) %>% mutate(total_nontfbs_snps = num_snps - total_tfbs_snps)
                epas1_snp_counts <- bind_rows(epas1_snp_counts,
                                     data.frame(species = c(species_oi), 
                                    n_snp_count = c(snp_counts %>% pull(num_snps)),
                                    total_tfbs_snps = c(snp_counts %>% pull(total_tfbs_snps)),
                                     total_nontfbs_snps = c(snp_counts %>% pull(total_nontfbs_snps))))
                epas1_avg_fc <- bind_rows(epas1_avg_fc,
                              data.frame(
                                species         = species_oi,
                                avg_tfbs_FC     = pull_or_zero(fc_calc %>% filter(in_TFBS), "avg_log_change"),
                                avg_nontfbs_FC  = pull_or_zero(fc_calc %>% filter(!in_TFBS), "avg_log_change"),
                                sd_tfbs_FC      = pull_or_zero(fc_calc %>% filter(in_TFBS), "sd_log_change"),
                                sd_nontfbs_FC   = pull_or_zero(fc_calc %>% filter(!in_TFBS), "sd_log_change")
                              )
                            )
            }else{
                epas1_snp_counts <- bind_rows(epas1_snp_counts,
                                     data.frame(species = c(species_oi), 
                                    n_snp_count = 0,
                                    total_tfbs_snps = 0,
                                     total_nontfbs_snps = 0))
                epas1_avg_fc <- bind_rows(epas1_avg_fc, 
                                    data.frame(species = species_oi, 
                                     avg_tfbs_FC = 0, 
                                     avg_nontfbs_FC = 0, 
                                     sd_tfbs_FC = 0,
                                     sd_nontfbs_FC = 0))
                epas1_all_fc <- bind_rows(epas1_all_fc, data.frame(species = species_oi))
            }
        }else{
            epas1_snp_counts <- bind_rows(epas1_snp_counts,
                                     data.frame(species = c(species_oi), 
                                    n_snp_count = 0,
                                    total_tfbs_snps = 0,
                                     total_nontfbs_snps = 0))
            epas1_avg_fc <- bind_rows(epas1_avg_fc, 
                                    data.frame(species = species_oi, 
                                     avg_tfbs_FC = 0, 
                                     avg_nontfbs_FC = 0, 
                                     sd_tfbs_FC = 0,
                                     sd_nontfbs_FC = 0))
            epas1_all_fc <- bind_rows(epas1_all_fc, data.frame(species = species_oi))
        }
        epas1_current_snps <- current_snps
    }
}

epas1_avg_fc <- epas1_avg_fc %>%
  mutate(across(
    .cols = where(is.numeric),
    .fns = ~ replace_na(.x, 0))
  )


epas1_snp_counts <- epas1_snp_counts %>% mutate(species2 = ifelse(species == 'Mus_musculus', 'Mouse', str_replace(species, "fullTree","")))
epas1_snp_counts$species2 <- factor(epas1_snp_counts$species2, levels = c(
    ifelse(keep.tips == "Mus_musculus", "Mouse", str_replace(keep.tips, "fullTree",""))
  ))
epas1_avg_fc <- epas1_avg_fc %>% mutate(species2 = ifelse(species == 'Mus_musculus', 'Mouse', str_replace(species, "fullTree","")))
epas1_avg_fc$species2 <- factor(epas1_avg_fc$species2, levels = c(
    ifelse(keep.tips == "Mus_musculus", "Mouse", str_replace(keep.tips, "fullTree",""))))
epas1_all_fc <- epas1_all_fc %>% mutate(species2 = ifelse(species == 'Mus_musculus', 'Mouse', str_replace(species, "fullTree","")))
epas1_all_fc$species2 <- factor(epas1_all_fc$species2, levels = c(
    ifelse(keep.tips == "Mus_musculus", "Mouse", str_replace(keep.tips, "fullTree",""))))

tfbs_map <- data.frame(msa_start = c(239,173),
                      msa_end = c(279,200))

epas1_all_fc <- epas1_all_fc %>%
                    mutate(tf_call = case_when(
                            in_TFBS ~ 'func. TFBS',
                            TRUE ~ 'non func. TFBS'))
epas1_group_fc <- epas1_all_fc %>% group_by(species2, tf_call, group_id) %>% 
                        summarise(avg_log_change = mean(log_change, na.rm = TRUE), sd_log_change = sd(log_change, na.rm = TRUE))


aff.mat <- epas1_mouse_msa_list$msa_df_TFBS %>% filter(!is.na(TF_name2)) %>% 
            mutate(species2 = ifelse(species == 'Mus_musculus', 'Mouse', str_replace(species, "fullTree",""))) %>% 
                mutate(marker = case_when(
                        species2%in%c("Anc52",'Anc38','Anc37') ~ 'Transitions',
                        species2 == 'Mus_musculus' ~ 'Mouse',
                        TRUE ~ 'Other'))
aff.mat$species2 <- factor(aff.mat$species2, levels = c(
    ifelse(keep.tips == "Mus_musculus", "Mouse", str_replace(keep.tips, "fullTree",""))
  ))
new_aff.mat <- data.frame()

for(tf2 in unique(aff.mat$TF_name2)){
    tmp.mat <- aff.mat %>% filter(TF_name2 == tf2)
    fc_steps <- c(0)
    delta_steps <- c(0)
    initial_step <- tmp.mat %>% filter(species == 'fullTreeAnc239') %>% pull(norm_affinity)
    for(species_oi in rev(keep.tips)){
        if(species_oi == 'fullTreeAnc239'){
        }else{
            current_step <- tmp.mat %>% filter(species == species_oi) %>% pull(norm_affinity)
            fc = current_step/initial_step
            delta <- current_step - initial_step
            fc_steps <- c(fc_steps, fc)
            delta_steps <- c(delta_steps, delta)
            initial_step <- current_step
        }
    }
    names(delta_steps) <- rev(keep.tips)
    names(fc_steps) <- rev(keep.tips)
    mm <- match(tmp.mat$species, names(delta_steps))
    tmp.mat$delta_steps <- delta_steps[mm]
    tmp.mat$fc_steps <- fc_steps[mm]
    start_aff <- tmp.mat %>% filter(species == 'fullTreeAnc239') %>% pull(norm_affinity)
    wt_aff <- tmp.mat %>% filter(species == 'Mus_musculus') %>% pull(norm_affinity)
    tmp.mat <- tmp.mat %>% mutate(frac_delta = delta_steps/(wt_aff - start_aff))
    new_aff.mat <- bind_rows(new_aff.mat, tmp.mat)
}
new_aff.mat <- new_aff.mat %>% mutate(log2_fc_steps = ifelse(fc_steps == 0, 0, log2(fc_steps)))
new_aff.mat <- new_aff.mat %>% 
            mutate(species2 = ifelse(species == 'Mus_musculus', 'Mouse', str_replace(species, "fullTree","")))
new_aff.mat$species2 <- factor(new_aff.mat$species2, levels = c(
    ifelse(keep.tips == "Mus_musculus", "Mouse", str_replace(keep.tips, "fullTree",""))
  ))
new_aff.mat$TF_name <- factor(new_aff.mat$TF_name, levels = c("Foxa2","Sox17","Gata4/6","Jun_Atf3","Klf4"))


##### all extant species
fasta_file <- '/net/shendure/vol10/projects/tli/nobackup/EB_dev_project/oCREs_v2/nobackup/reoriented_fasta/Epas1_chr17_10063_oCRE.fasta'
all_species <- oCRE_mean_act %>% filter(full_CRE_id == CRE_oi) %>% pull(species)
epas1_msa_df <- make_alignment_dataframe(fasta_file, all_species, 'Mus_musculus',CRE_oi)
mouse_msa_df <- epas1_msa_df %>% filter(species == 'Mus_musculus')

epas1_msa_list <- make_msa_tfbs_traj_df(fasta_file, all_species, 'Mus_musculus', CRE_oi, cactus.tfbs, mouse_func_tfbs)
labels <- collect_labels(epas1_msa_list$msa_df_TFBS, side = T)
epas1_func_tfbs <- mouse_func_tfbs %>% filter(CRE == CRE_oi)
mm.start <- match(epas1_func_tfbs$TFBS_start, mouse_msa_df$ref_pos)
mm.end <- match(epas1_func_tfbs$TFBS_end, mouse_msa_df$ref_pos)
epas1_func_tfbs$msa_start <- mouse_msa_df$curpos[mm.start]
epas1_func_tfbs$msa_end <- mouse_msa_df$curpos[mm.end]
epas1_msa_df <- left_join(epas1_msa_df, mouse_DMS_MPRA %>% filter(CRE == CRE_oi) %>%
                          mutate(log2FC_MPRA = ifelse(is.infinite(log2FC_MPRA) & log2FC_MPRA < 0, NA, log2FC_MPRA)) %>% ### replace -Inf with NA
                          dplyr::select(ref_pos = pos, mut, log2FC_MPRA, sd_log2FC_MPRA, MPRA_act, sd_MPRA_act, min_act,max_act, wilcox_BH),
                         by = c("ref_pos","mut"))

epas1_msa_df <- epas1_msa_df %>% mutate(mut2 = case_when(
                        mut == 'DEL' ~ '-',
                        TRUE ~ mut))

