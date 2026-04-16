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

CRE_oi <- "Gata4_chr14_5729"
### plot Gata4 transitions
species.pair <- c("Mus_musculus",'fullTreeAnc239')
mouse.nodes <- make.nodepath(tree, species.pair)
### filter nodes if in activity
mouse.nodes <- mouse.nodes[mouse.nodes%in%unique(oCRE_mean_act$species)]

start_node <- "fullTreeAnc239"
order_species <- rev(mouse.nodes)  # traverse in your chosen order

### get replicate
oCRE_replicate_act <- df_mpra_act %>% filter(class == 'oCRE')
oCRE_replicate_act <- oCRE_replicate_act %>% separate(CRE_id, c('tile_id','something'), sep = "_20240416") %>% dplyr::select(-something)
oCRE_replicate_act <- oCRE_replicate_act %>% separate(tile_id, c('full_CRE_id','species'), sep = "__")
mm <- match(oCRE_replicate_act$species, cactus.annotation$Species)
oCRE_replicate_act$common_name <- cactus.annotation$`Common Name`[mm]
oCRE_replicate_act <- left_join(oCRE_replicate_act, max.dist, by =c("full_CRE_id","species"))


# base table (ordered once)
gata4_avg_lineage_mpra <- oCRE_mean_act %>%
  filter(full_CRE_id == CRE_oi, species %in% mouse.nodes) %>%
  mutate(
    species2  = ifelse(species == 'Mus_musculus', 'Mouse', str_replace(species, "fullTree","")),
    marker = case_when(
      species2 %in% c("Anc52", "Anc38", "Anc37") ~ "Transitions",
      species2 == "Mouse"                 ~ "Mouse",
      TRUE                                       ~ "Other"
    ),
    species2 = factor(species2, levels = c(
    ifelse(order_species == "Mus_musculus", "Mouse", str_replace(order_species, "fullTree",""))
  ))
  ) %>%
  arrange(species)

gata4_lineage_mpra <- oCRE_replicate_act %>%
  filter(full_CRE_id == CRE_oi, species %in% mouse.nodes) %>%
  mutate(
    species2 = ifelse(species == 'Mus_musculus', 'Mouse', str_replace(species, "fullTree","")),
    marker = case_when(
      species2 %in% c("Anc52", "Anc38", "Anc37") ~ "Transitions",
      species2 == "Mouse"                 ~ "Mouse",
      TRUE                                       ~ "Other"
    ),
    species = factor(species, levels = order_species)
  ) %>%
  arrange(species)

# constants
start_act_df <- gata4_lineage_mpra %>%
  filter(species == start_node) %>%
  dplyr::select(biol_rep, start_act = MPRA_act)

# Get mouse wild-type activity per replicate
wt_act_df <- gata4_lineage_mpra %>%
  filter(species == 'Mus_musculus') %>%
  dplyr::select(biol_rep, wt_act = MPRA_act)

# stepwise deltas & fold-changes
gata4_lineage_mpra_steps <- gata4_lineage_mpra %>%
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
    ifelse(order_species == "Mus_musculus", "Mouse", str_replace(order_species, "fullTree",""))
  ))
  )

stepwise_summary <- gata4_lineage_mpra_steps %>%
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
stepwise_summary$p_adj_log2FC_steps <- p.adjust(stepwise_summary$p_log2FC_steps, method = 'fdr')
stepwise_summary <- stepwise_summary %>% 
    mutate(
        species2 = ifelse(species == 'Mus_musculus', 'Mouse', str_replace(species, "fullTree","")),
        marker = case_when(
          species2 %in% c("Anc52", "Anc38", "Anc37") ~ "Transitions",
          species2 == "Mouse"                 ~ "Mouse",
          TRUE                                       ~ "Other"
        ),
        species2  = factor(species2, levels = c(
    ifelse(order_species == "Mus_musculus", "Mouse", str_replace(order_species, "fullTree",""))
  ))
    )


species.pair <- c("Rattus_norvegicus",'fullTreeAnc239')
rat.nodes <- make.nodepath(tree, species.pair)
### filter nodes if in activity
rat.nodes <- rat.nodes[rat.nodes%in%unique(oCRE_mean_act$species)]

order_species <- rev(rat.nodes)  # traverse in your chosen order

# base table (ordered once)
gata4_avg_rat_lineage_mpra <- oCRE_mean_act %>%
  filter(full_CRE_id == CRE_oi, species %in% rat.nodes) %>%
  mutate(
    species2 = ifelse(species == 'Rattus_norvegicus', 'Rat', str_replace(species, "fullTree","")),
    marker = case_when(
      species2 %in% c("Anc52", "Anc38") ~ "Transitions",
      species2 == "Rat"            ~ "Rat",
      TRUE                                       ~ "Other"
    ),
    species2 = factor(species2, levels = c(
    ifelse(order_species == "Rattus_norvegicus", "Rat", str_replace(order_species, "fullTree",""))
  ))
  ) %>%
  arrange(species)

gata4_rat_lineage_mpra <- oCRE_replicate_act %>%
  filter(full_CRE_id == CRE_oi, species %in% rat.nodes) %>%
  mutate(
    species2 = ifelse(species == 'Rattus_norvegicus', 'Rat', str_replace(species, "fullTree","")),
    marker = case_when(
      species2 %in% c("Anc52", "Anc38") ~ "Transitions",
      species2 == "Rat"            ~ "Rat",
      TRUE                                       ~ "Other"
    ),
    species2 = factor(species2, levels = c(
    ifelse(order_species == "Rattus_norvegicus", "Rat", str_replace(order_species, "fullTree",""))
  ))
  ) %>%
  arrange(species)

# constants
start_act_df <- gata4_rat_lineage_mpra %>%
  filter(species == start_node) %>%
  dplyr::select(biol_rep, start_act = MPRA_act)

# Get mouse wild-type activity per replicate
wt_act_df <- gata4_rat_lineage_mpra %>%
  filter(species == 'Rattus_norvegicus') %>%
  dplyr::select(biol_rep, wt_act = MPRA_act)

# stepwise deltas & fold-changes
gata4_rat_lineage_mpra_steps <- gata4_rat_lineage_mpra %>%
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
    ifelse(order_species == "Rattus_norvegicus", "Rat", str_replace(order_species, "fullTree",""))
  ))
  )

rat_stepwise_summary <- gata4_rat_lineage_mpra_steps %>%
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
rat_stepwise_summary$p_adj_log2FC_steps <- p.adjust(rat_stepwise_summary$p_log2FC_steps, method = 'fdr')
rat_stepwise_summary <- rat_stepwise_summary %>% 
    mutate(
    species2 = ifelse(species == 'Rattus_norvegicus', 'Rat', str_replace(species, "fullTree","")),
    marker = case_when(
      species2 %in% c("Anc52", "Anc38") ~ "Transitions",
      species2 == "Rat"            ~ "Rat",
      TRUE                                       ~ "Other"
    ),
    species2 = factor(species2, levels = c(
    ifelse(order_species == "Rattus_norvegicus", "Rat", str_replace(order_species, "fullTree",""))
  ))
  )

gata4_avg_lineage_mpra$group <- 1
stepwise_summary$group <- 1

gata4_avg_rat_lineage_mpra$group <- 2
rat_stepwise_summary$group <- 2

### generate subtree
rodent.tree <- tree_subset(tree, node = 'fullTreeAnc55', levels_back = 0)
rodent.tree <- as_tibble(rodent.tree)
cre_mean_act <- oCRE_mean_act %>% filter(full_CRE_id == 'Gata4_chr14_5729')
shape.values <- c(1,24,13,15)
names(shape.values) <- c('mouse','no ortholog','not sequenced', 'sequenced')
sub.seq <- cre_mean_act[,c('species','mean_MPRA_act','common_name')]
sub.seq$label <- sub.seq$species
sub.seq <- sub.seq %>% mutate(capture = case_when(
                        species == 'Mus_musculus' ~ 'mouse',
                        species == 'Rattus_norvegicus' ~ 'rat',
                        is.na(mean_MPRA_act) ~ 'not sequenced',
                        TRUE ~ 'sequenced'))

sub.tree <- left_join(rodent.tree, sub.seq, by = 'label')
sub.tree <- sub.tree %>% mutate(capture = case_when(
                        is.na(capture) ~ 'no ortholog',
                        TRUE ~ capture)) %>% 
            mutate(label2 = case_when(
                label == 'Rattus_norvegicus' ~ 'Rat',
                label == 'Mus_musculus' ~ 'Mouse',
                label%in%c(mouse.nodes,'Rattus_norvegicus') ~ str_replace(label, "fullTree",""),
                TRUE ~ NA))
sub.tips <- sub.tree %>% filter(label%in%c(mouse.nodes,'Rattus_norvegicus')) %>%
            mutate(Shape = case_when(
                label == 'Mus_musculus' ~ 1,
                label == 'Rattus_norvegicus' ~ 1,
                label%in%mouse.nodes ~ 9,
                TRUE ~ 15)) %>% 
            mutate(Size = case_when(
                label == 'Mus_musculus' ~ 5,
                label == 'Rattus_norvegicus' ~ 5,
                label%in%mouse.nodes ~ 5,
                TRUE ~ 0))
sub.tips$Shape <- factor(sub.tips$Shape)
sub.tips <- as.tibble(sub.tips)
sub.tree <- as.treedata(sub.tree)

### gata4 fasta file ###
CRE_oi <- "Gata4_chr14_5729"
fasta_file <- '../data/oCRE_fasta/Gata4_chr14_5729_oCRE.fasta'
species.pair <- c("Mus_musculus",'fullTreeAnc55')
anc.cutoff <- make.nodepath(tree, species.pair)

### plot Gata4 MSA for five rodents
mouse_DMS_MPRA$mut <- toupper(mouse_DMS_MPRA$mut)
species.pair <- c("Mus_musculus",'fullTreeAnc239')
keep.tips <- make.nodepath(tree, species.pair)
gata4_msa_df <- make_alignment_dataframe(fasta_file, keep.tips, 'Mus_musculus',CRE_oi)
mouse_msa_df <- gata4_msa_df %>% filter(species == 'Mus_musculus')

gata4_msa_list <- make_msa_tfbs_traj_df(fasta_file, rev(keep.tips), 'Mus_musculus', CRE_oi, cactus.tfbs, mouse_func_tfbs)
labels <- collect_labels(gata4_msa_list$msa_df_TFBS, side = T)
gata4_func_tfbs <- mouse_func_tfbs %>% filter(CRE == CRE_oi)
mm.start <- match(gata4_func_tfbs$TFBS_start, mouse_msa_df$ref_pos)
mm.end <- match(gata4_func_tfbs$TFBS_end, mouse_msa_df$ref_pos)
gata4_func_tfbs$msa_start <- mouse_msa_df$curpos[mm.start]
gata4_func_tfbs$msa_end <- mouse_msa_df$curpos[mm.end]
gata4_msa_df <- left_join(gata4_msa_df, mouse_DMS_MPRA %>% filter(CRE == CRE_oi) %>%
                          dplyr::select(ref_pos = pos, mut, log2FC_MPRA, sd_log2FC_MPRA, MPRA_act, sd_MPRA_act, min_act,max_act, wilcox_BH),
                         by = c("ref_pos","mut"))

gata4_msa_df <- gata4_msa_df %>% mutate(mut2 = case_when(
                        mut == 'DEL' ~ '-',
                        TRUE ~ mut))
gata4_msa_df <- gata4_msa_df %>% mutate(species2 = ifelse(species == 'Mus_musculus', 'Mouse', str_replace(species, "fullTree","")))
gata4_msa_df$species2 <- factor(gata4_msa_df$species2, levels = c(
    ifelse(keep.tips == "Mus_musculus", "Mouse", str_replace(keep.tips, "fullTree",""))
  ))

gata4_msa_df <- gata4_msa_df %>%
  mutate(fill_type = case_when(
    mut2 == "-" & is.na(log2FC_MPRA) ~ "black",
    mut == ref_nuc ~ "grey",
    TRUE ~ "gradient"
  ))

fixed_gata4_msa_df <- data.frame()

for(pos in seq(1:max(gata4_msa_df$curpos))){
    tmp_msa_df <- gata4_msa_df %>% filter(curpos == pos)
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
#                 current_muttype <- tmp_species %>% pull(mut_type)
#                 if(current_muttype != 'Match'){
#                     current_log2fc <- tmp_species %>% pull(log2FC_MPRA)
#                     log_change <- c(log_change, 0)
#                     match2 <- c(match2, current_muttype)
#                 }else{
#                     log_change <- c(log_change, -1*current_log2fc)
#                     have_not_found_match <- FALSE
#                     match2 <- c(match2, 'First Match')
#                 }
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
        fixed_gata4_msa_df <- bind_rows(fixed_gata4_msa_df, tmp_msa_df)
    }else{
        tmp_msa_df$log_change <- 0
        tmp_msa_df$match2 <- 'Match'
        fixed_gata4_msa_df <- bind_rows(fixed_gata4_msa_df, tmp_msa_df)
    }
}

fixed_gata4_msa_df <- fixed_gata4_msa_df %>% mutate(log_change = case_when(
                                        mut_type%in%c('DEL','INS') ~ NA,
                                        TRUE ~ log_change))
fixed_gata4_msa_df <- fixed_gata4_msa_df %>% 
            mutate(species2 = ifelse(species == 'Mus_musculus', 'Mouse', str_replace(species, "fullTree","")))
fixed_gata4_msa_df$species2 <- factor(fixed_gata4_msa_df$species2, levels = c(
    ifelse(keep.tips == "Mus_musculus", "Mouse", str_replace(keep.tips, "fullTree",""))
  ))


#### calculate number of SNPs that become the original mouse sequence at each step and mark which original mouse SNP confers a TFBS
### calculate avg log2FC for each mutation gained
all_mouse_tfbs_positions <- unlist(mapply(seq, gata4_func_tfbs$msa_start, gata4_func_tfbs$msa_end))
gata4_current_snps <- gata4_msa_df %>% filter(species == 'fullTreeAnc239') %>% filter(mut2 != ref_nuc) %>% pull(curpos)
gain_mouse_snps <- fixed_gata4_msa_df %>% filter(species == 'fullTreeAnc239') %>% filter(mut2 == ref_nuc) %>%
                            arrange(curpos) %>% mutate(group_id = cumsum(c(TRUE, diff(curpos) > 1)),
                                                      in_TFBS = case_when(
                                                          curpos%in%all_mouse_tfbs_positions ~ TRUE,
                                                          TRUE ~ FALSE))
fc_calc <- gain_mouse_snps %>% group_by(in_TFBS) %>% 
                        summarise(avg_log_change = mean(log_change, na.rm = TRUE), sd_log_change = sd(log_change, na.rm = TRUE))
gata4_all_fc <- gain_mouse_snps

gata4_snp_counts <- gain_mouse_snps %>%
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
gata4_avg_fc <- data.frame(
                                species         = 'fullTreeAnc239',
                                avg_tfbs_FC     = pull_or_zero(fc_calc %>% filter(in_TFBS), "avg_log_change"),
                                avg_nontfbs_FC  = pull_or_zero(fc_calc %>% filter(!in_TFBS), "avg_log_change"),
                                sd_tfbs_FC      = pull_or_zero(fc_calc %>% filter(in_TFBS), "sd_log_change"),
                                sd_nontfbs_FC   = pull_or_zero(fc_calc %>% filter(!in_TFBS), "sd_log_change")
                              )
for(species_oi in rev(keep.tips)){
    if(species_oi == 'fullTreeAnc239'){
    } else{
        temp_msa_df <- fixed_gata4_msa_df %>% filter(species == species_oi)
        current_snps <- temp_msa_df %>% filter(mut2 != ref_nuc) %>% pull(curpos)
        if(length(gata4_current_snps) > 0){
            gain_mouse_snps <- temp_msa_df %>% filter(curpos%in%gata4_current_snps) %>% filter(mut2 == ref_nuc) %>%
                            arrange(curpos) %>% mutate(group_id = cumsum(c(TRUE, diff(curpos) > 1)),
                                                      in_TFBS = case_when(
                                                          curpos%in%all_mouse_tfbs_positions ~ TRUE,
                                                          TRUE ~ FALSE))
            if (nrow(gain_mouse_snps) > 0) {
              fc_calc <- gain_mouse_snps %>% group_by(in_TFBS) %>% 
                        summarise(avg_log_change = mean(log_change, na.rm = TRUE), sd_log_change = sd(log_change, na.rm = TRUE))
              gata4_all_fc <- bind_rows(gata4_all_fc, gain_mouse_snps)
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
                gata4_snp_counts <- bind_rows(gata4_snp_counts,
                                     data.frame(species = c(species_oi), 
                                    n_snp_count = c(snp_counts %>% pull(num_snps)),
                                    total_tfbs_snps = c(snp_counts %>% pull(total_tfbs_snps)),
                                     total_nontfbs_snps = c(snp_counts %>% pull(total_nontfbs_snps))))
                gata4_avg_fc <- bind_rows(gata4_avg_fc,
                              data.frame(
                                species         = species_oi,
                                avg_tfbs_FC     = pull_or_zero(fc_calc %>% filter(in_TFBS), "avg_log_change"),
                                avg_nontfbs_FC  = pull_or_zero(fc_calc %>% filter(!in_TFBS), "avg_log_change"),
                                sd_tfbs_FC      = pull_or_zero(fc_calc %>% filter(in_TFBS), "sd_log_change"),
                                sd_nontfbs_FC   = pull_or_zero(fc_calc %>% filter(!in_TFBS), "sd_log_change")
                              )
                            )
            }else{
                gata4_snp_counts <- bind_rows(gata4_snp_counts,
                                     data.frame(species = c(species_oi), 
                                    n_snp_count = 0,
                                    total_tfbs_snps = 0,
                                     total_nontfbs_snps = 0))
                gata4_avg_fc <- bind_rows(gata4_avg_fc, 
                                    data.frame(species = species_oi, 
                                     avg_tfbs_FC = 0, 
                                     avg_nontfbs_FC = 0, 
                                     sd_tfbs_FC = 0,
                                     sd_nontfbs_FC = 0))
                gata4_all_fc <- bind_rows(gata4_all_fc, data.frame(species = species_oi))
            }
        }else{
            gata4_snp_counts <- bind_rows(gata4_snp_counts,
                                     data.frame(species = c(species_oi), 
                                    n_snp_count = 0,
                                    total_tfbs_snps = 0,
                                     total_nontfbs_snps = 0))
            gata4_avg_fc <- bind_rows(gata4_avg_fc, 
                                    data.frame(species = species_oi, 
                                     avg_tfbs_FC = 0, 
                                     avg_nontfbs_FC = 0, 
                                     sd_tfbs_FC = 0,
                                     sd_nontfbs_FC = 0))
            gata4_all_fc <- bind_rows(gata4_all_fc, data.frame(species = species_oi))
        }
        gata4_current_snps <- current_snps
    }
}

gata4_avg_fc <- gata4_avg_fc %>%
  mutate(across(
    .cols = where(is.numeric),
    .fns = ~ if_else(species == "fullTreeAnc51", .x, replace_na(.x, 0))
  ))


gata4_snp_counts <- gata4_snp_counts %>% mutate(species2 = ifelse(species == 'Mus_musculus', 'Mouse', str_replace(species, "fullTree","")))
gata4_snp_counts$species2 <- factor(gata4_snp_counts$species2, levels = c(
    ifelse(keep.tips == "Mus_musculus", "Mouse", str_replace(keep.tips, "fullTree",""))
  ))
gata4_avg_fc <- gata4_avg_fc %>% mutate(species2 = ifelse(species == 'Mus_musculus', 'Mouse', str_replace(species, "fullTree","")))
gata4_avg_fc$species2 <- factor(gata4_avg_fc$species2, levels = c(
    ifelse(keep.tips == "Mus_musculus", "Mouse", str_replace(keep.tips, "fullTree",""))))
gata4_all_fc <- gata4_all_fc %>% mutate(species2 = ifelse(species == 'Mus_musculus', 'Mouse', str_replace(species, "fullTree","")))
gata4_all_fc$species2 <- factor(gata4_all_fc$species2, levels = c(
    ifelse(keep.tips == "Mus_musculus", "Mouse", str_replace(keep.tips, "fullTree",""))))

tfbs_map <- data.frame(msa_start = c(26,79,175,241),
                      msa_end = c(37,92,186,252))


final_gata4_msa_df <- fixed_gata4_msa_df %>% filter(species%in%anc.cutoff)
final_gata4_snp_counts <- gata4_snp_counts %>% filter(species%in%anc.cutoff)
final_gata4_avg_fc <- gata4_avg_fc %>% filter(species%in%anc.cutoff)
final_gata4_all_fc <- gata4_all_fc %>% filter(species%in%anc.cutoff) %>% 
                    mutate(tf_call = case_when(
                            in_TFBS ~ 'func. TFBS',
                            TRUE ~ 'non func. TFBS'))
final_gata4_group_fc <- final_gata4_all_fc %>% group_by(species2, tf_call, group_id) %>% 
                        summarise(avg_log_change = mean(log_change, na.rm = TRUE), sd_log_change = sd(log_change, na.rm = TRUE))


aff.mat <- gata4_msa_list$msa_df_TFBS %>% filter(!is.na(TF_name2)) %>% 
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
    initial_step <- tmp.mat %>% filter(species == 'fullTreeAnc55') %>% pull(norm_affinity)
    for(species_oi in rev(keep.tips)){
        if(species_oi == 'fullTreeAnc55'){
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
    start_aff <- tmp.mat %>% filter(species == 'fullTreeAnc55') %>% pull(norm_affinity)
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
new_aff.mat$TF_name <- factor(new_aff.mat$TF_name, levels = c("Sox17","Gata4/6","Jun_Atf3","Klf4"))
final_aff.mat <- new_aff.mat %>% filter(species%in%anc.cutoff)


### for rat msa
species.pair <- c("Rattus_norvegicus",'fullTreeAnc239')
keep.tips <- make.nodepath(tree, species.pair)
rat_gata4_msa_df <- make_alignment_dataframe(fasta_file, c(keep.tips,'Mus_musculus'), 'Mus_musculus',CRE_oi)
mouse_msa_df <- rat_gata4_msa_df %>% filter(species == 'Mus_musculus')

rat_gata4_msa_list <- make_msa_tfbs_traj_df(fasta_file, rev(c(keep.tips,'Mus_musculus')), 'Mus_musculus', 
                                            CRE_oi, cactus.tfbs, mouse_func_tfbs)