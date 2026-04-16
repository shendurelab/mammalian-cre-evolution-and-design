### module load pcre2/10.39; module load R/4.3.1
library(tidyverse)
library(ggplot2)
library(cowplot)
library(ggrastr)
library(castor)
library(reshape2)
library(DescTools)
library(readxl)
library(gtools)
source("utils/mpra_tablemaker.R")
source("utils/figure4_func.R")

calc_gen_changes <- function(data){
    wt.event <- data %>% filter(gen == 0) %>% dplyr::select(CRE, wt_MPRA_act = mean_MPRA_act) %>% distinct()
    data <- left_join(data, wt.event, by = c("CRE"))
    data <- data %>% mutate(log2FC_WT_MPRA = log2(mean_MPRA_act/wt_MPRA_act))    
    
    final.df <- data %>% filter(gen == 0) ### keep WT
    
    for(CRE_oi in unique(data$CRE)){
        gen_df <- data %>% filter(CRE == CRE_oi) %>% filter(gen != 0) %>% arrange(gen)
        start_mean_act <- wt.event %>% filter(CRE == CRE_oi) %>% pull(wt_MPRA_act)
        mean_gen_act <- gen_df %>% pull(mean_MPRA_act)
        
        log2FC <- c()
        traj.delta <- c()
        wt.delta <- c()
        for(i in 1:length(mean_gen_act)){
            if(i == 1){
                log2FC <- c(log2FC, log2(mean_gen_act[i]/start_mean_act))
                traj.delta <- c(traj.delta, mean_gen_act[i] - start_mean_act)
                wt.delta <- c(wt.delta, mean_gen_act[i] - start_mean_act)
            }else{
                log2FC <- c(log2FC, log2(mean_gen_act[i]/mean_gen_act[i-1]))
                traj.delta <- c(traj.delta, mean_gen_act[i] - mean_gen_act[i-1])
                wt.delta <- c(wt.delta, mean_gen_act[i] - start_mean_act)
            }
        }
        gen_df$log2FC_TRAJ_MPRA <- log2FC
        gen_df$traj_delta <- traj.delta
        gen_df$wt_delta <- wt.delta
        gen_df <- gen_df %>% mutate(frac_traj_delta = traj_delta/start_mean_act,
                                   frac_wt_delta = wt_delta/start_mean_act) 
        final.df <- bind_rows(final.df, gen_df)
    }
    return(final.df)
}

add_func_mapping <- function(data){
    final.df <- data %>% filter(gen == 0) ### keep WT
    
    for(CRE_oi in unique(data$CRE)){
        cre_events <- data %>% filter(CRE == CRE_oi) %>% filter(gen != 0) %>% mutate(ref_start = mapping, ref_end = mapping)
        cre_func_tfbs <- mouse_func_tfbs %>% filter(CRE == CRE_oi)
        cre_events <- mark_func_overlap_df(cre_events, cre_func_tfbs) %>% dplyr::select(-ref_start, -ref_end)
        final.df <- bind_rows(final.df, cre_events)
    }
    return(final.df)
}

print("Loading model trajectory dataset...") 
df_counts_w_CREs2 <- read.table("../data/oCRE_and_CRE_optimization_MPRA_count_table.txt.gz",
                                header=TRUE)


# BC per CRE recovered?
coverage_lib <- calculate_MPRA_coverage(df_counts_w_CREs2)
df_wide <- generate_MPRA_wide_table(df_counts_w_CREs2, coverage_lib)
df_mpra_act <- generate_winsorize_MPRA_table(df_wide)

control_mean_act <- df_mpra_act %>% filter(class%in%c('neg_control','pos_control')) %>% group_by(CRE_id) %>% 
        summarise(mean_MPRA_act = mean(MPRA_act), sd_MPRA_act = sd(MPRA_act)) %>% dplyr::select(CRE_id, mean_MPRA_act, sd_MPRA_act)

### load mouse WT
cactus.annotation <- read_xlsx('../data/cactus_ucsc_species_list.xlsx')
cactus.annotation$Species <- gsub(' ','_',cactus.annotation$Species)
cactus.annotation$label <- cactus.annotation$Species
oCRE_mean_act <- df_mpra_act %>% filter(class == 'oCRE') %>% group_by(CRE_id) %>% 
        summarise(mean_MPRA_act = mean(MPRA_act), sd_MPRA_act = sd(MPRA_act), max_act = max(MPRA_act), min_act = min(MPRA_act))
oCRE_acc_pred <- read.delim("../data/oCRE_marginal_footprint.txt",
                           sep = '\t')

oCRE_mean_act <- oCRE_mean_act %>% separate(CRE_id, c('tile_id','something'), sep = "_20240416") %>% dplyr::select(-something)
oCRE_mean_act <- left_join(oCRE_mean_act, oCRE_acc_pred, by = c("tile_id"))
oCRE_mean_act <- oCRE_mean_act %>% separate(tile_id, c('full_CRE_id','species'), sep = "__")
mm <- match(oCRE_mean_act$species, cactus.annotation$Species)
oCRE_mean_act$common_name <- cactus.annotation$`Common Name`[mm]
oCRE_mean_act$log2_norm_score <- log2(oCRE_mean_act$norm_score)

### get mouse WT controls
df_mouse_tiles_mpra <- oCRE_mean_act %>% filter(species == 'Mus_musculus') 

### get trajectory analysis
traj_mean_act <- df_mpra_act %>% filter(class == 'CRE_traj') %>% group_by(CRE_id) %>% summarise(mean_MPRA_act = mean(MPRA_act),
                                                                           sd_MPRA_act = sd(MPRA_act), 
                                                                        max_act = max(MPRA_act), min_act = min(MPRA_act))
colnames(traj_mean_act) <- c("tile_name",'mean_MPRA_act','sd_MPRA_act','max_act','min_act')
traj_info_df <- read.delim("../data/oligo_info/CRE_optimization_order.txt", sep = '\t') %>% dplyr::select(-oligo_seq)

traj_mean_act <- left_join(traj_mean_act, traj_info_df, by = c('tile_name'))
wt_mean_act <- traj_mean_act %>% filter(gen == 0) %>% dplyr::select(CRE, mean_MPRA_act, sd_MPRA_act)
colnames(wt_mean_act) <- c('CRE','mean_MPRA_act','sd_MPRA_act')
traj_mean_act <- bind_rows(traj_mean_act, 
                          traj_mean_act %>% filter(gen == 0) %>% mutate(selection = 'Minimizing'))
traj_mean_act <- traj_mean_act %>% mutate(mapping = case_when(
                                mapping == 'WT_1' ~ NA,
                                TRUE ~ as.numeric(mapping)))
wt_seqs <- traj_mean_act %>% filter(is.na(mapping)) %>% filter(selection == 'Maximizing') %>% dplyr::select(CRE, wt_seq = seq)
traj_mean_act <- left_join(traj_mean_act, wt_seqs, by = c("CRE"))
traj_mean_act <- traj_mean_act %>% mutate(mut_nuc = case_when(
                                is.na(mapping) ~ NA,
                                TRUE ~ substr(seq, mapping + 1, mapping + 1)), 
                                         wt_nuc = case_when(
                                is.na(mapping) ~ NA,
                                TRUE ~ substr(wt_seq, mapping + 1, mapping + 1))) %>% 
                mutate(is_mutation = case_when(
                    mut_nuc != wt_nuc ~ TRUE,
                    TRUE ~ FALSE))
delta_mpra_act <- traj_mean_act %>%
  group_by(CRE, selection) %>%
  arrange(gen) %>%  # Ensure gen is in order
  mutate(traj_delta = mean_MPRA_act - lag(mean_MPRA_act)) %>%
  ungroup() %>% 
  mutate(traj_delta = case_when(
          is.na(traj_delta) ~ 0,
          TRUE ~ traj_delta))
delta_mpra_act <- left_join(delta_mpra_act, wt_mean_act %>% dplyr::select(CRE, wt_mean_MPRA_act = mean_MPRA_act), by = c("CRE"))
delta_mpra_act <- delta_mpra_act %>% mutate(wt_delta = mean_MPRA_act - wt_mean_MPRA_act,
                                           frac_traj_delta = traj_delta/(wt_mean_MPRA_act - control_mean_act %>% filter(CRE_id == 'minP') %>% pull(mean_MPRA_act)),
                                           frac_wt_delta = wt_delta/(wt_mean_MPRA_act - control_mean_act %>% filter(CRE_id == 'minP') %>% pull(mean_MPRA_act))) %>% dplyr::select(-wt_mean_MPRA_act)


traj_mean_act <- delta_mpra_act
traj_mean_act$CRE <- factor(traj_mean_act$CRE, levels = c("Gata4_chr14_5729","Epas1_chr17_10063","Lama1_chr17_7784","Sparc_chr11_7211","Bend5_chr4_8201"))


selection.color <- c("#ee4035","#0392cf")
names(selection.color) <- c("Maximizing","Minimizing")

print("Finish loading model trajectory dataset.") 


print("Loading model trajectory TFBS dataset...")

final_TFs <- c("Jun_Atf3","Foxa2","Gata4/6","Sox17","Klf4","Hnf1b")
col_TFs <- c("firebrick2","forestgreen","darkorange2","dodgerblue1","mediumorchid3",'goldenrod1')
names(col_TFs) <- final_TFs 

df_TF_max_aff <- read.table("../data/ParEndo_markerTF_cuttoffs.txt",
                            header=TRUE)
traj.tfbs <- read.delim(gzfile("../data/CRE_optimization_ProBound_ParEndo_TFBS_predictions.txt.gz"),
                       sep = '\t')
# joining Gata4/6 calls to avoid redundancy
gata.tfbs <- traj.tfbs %>% filter(TF_name %in% c("Gata4", "Gata6")) 
gata.tfbs <- gata.tfbs %>% mutate(TF_name = "Gata4/6") %>%
  group_by(start_pos, end_pos, CRE) %>%
  dplyr::slice_max(order_by = norm_affinity, with_ties = FALSE) %>%
  ungroup()
traj.tfbs <- bind_rows(traj.tfbs %>% filter(!TF_name %in% c("Gata4", "Gata6")),
                         gata.tfbs) 
traj.tfbs <- traj.tfbs %>% filter(TF_name%in%final_TFs)
traj.tfbs$start_pos <- traj.tfbs$start_pos + 1 ### fix 0-index position from probound
colnames(traj.tfbs)[8] <- 'tile_name'
traj.tfbs <- left_join(traj.tfbs, traj_info_df, 
                                   by =c('tile_name')) %>% dplyr::select(-seq)     
### load functional mouse TFBS data
mouse_DMS_tfbs <- read.delim("../data/mouse_TFBS_impact_and_affinity_CRM.txt", sep = '\t') 
gata.tfbs <- mouse_DMS_tfbs %>% filter(TF_name %in% c("Gata4", "Gata6")) 
gata.tfbs <- gata.tfbs %>% mutate(TF_name = "Gata4/6") %>%
  group_by(TFBS_start, TFBS_end, CRE) %>%
  dplyr::slice_max(order_by = norm_affinity, with_ties = FALSE) %>%
  ungroup()
mouse_DMS_tfbs <- bind_rows(mouse_DMS_tfbs %>% filter(!TF_name %in% c("Gata4", "Gata6")) , gata.tfbs) %>% 
                  filter(TF_name%in%final_TFs)
mouse_func_tfbs <- mouse_DMS_tfbs %>% filter(functional_TF) 
traj_mean_act <- add_func_mapping(traj_mean_act)
print("Finish loading model trajectory TFBS dataset")


print("Loading mouse DMS data...")
mouse_DMS_MPRA <- read.delim("../data/mouse_DMS_MPRA_data.txt", sep = '\t')  
print("Finished loading mouse DMS data.")


### add mEB predictions
pluri.predictions <- read.delim('../data/predictions/CRE_optimization_pluri_predictions.txt.gz', sep = '\t')
ecto.predictions <- read.delim('../data/predictions/CRE_optimization_ecto_predictions.txt.gz', sep = '\t')
meso.predictions <- read.delim('../data/predictions/CRE_optimization_meso_predictions.txt.gz', sep = '\t')

traj_mean_act <- traj_mean_act %>% mutate(MOTIF_NAME = case_when(
                                is.na(mapping) ~ paste(gen, selection, CRE, 'WT', sep = "_"),
                                TRUE ~ paste(gen, selection, CRE, mapping, sep = "_")))

traj_mean_act <- left_join(traj_mean_act, pluri.predictions %>% dplyr::select(MOTIF_NAME, pluri_pred = norm_footprint), by = 'MOTIF_NAME')
traj_mean_act <- left_join(traj_mean_act, ecto.predictions %>% dplyr::select(MOTIF_NAME, ecto_pred = norm_footprint), by = 'MOTIF_NAME')
traj_mean_act <- left_join(traj_mean_act, meso.predictions %>% dplyr::select(MOTIF_NAME, meso_pred = norm_footprint), by = 'MOTIF_NAME')

cre_levels <- c('Gata4:chr14_5729','Epas1:chr17_10063','Lama1:chr17_7784','Sparc:chr11_7211','Bend5:chr4_8201')
traj_mean_act <- traj_mean_act %>%
  mutate(
    CRE = sub("_", ":", CRE)  # Replace underscores with colons
  )
traj_mean_act$CRE <- factor(traj_mean_act$CRE, levels = cre_levels)

wt_mean_act <- wt_mean_act %>%
  mutate(
    CRE = sub("_", ":", CRE)  # Replace underscores with colons
  )
wt_mean_act$CRE <- factor(wt_mean_act$CRE, levels = cre_levels)

