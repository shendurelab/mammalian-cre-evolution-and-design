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
library(readxl)

source("utils/mafft_aligning_functions.R")
source("utils/mpra_tablemaker.R")

cre_levels <- c('Gata4:chr14_5729','Epas1:chr17_10063','Lama1:chr17_7784')

## load tree
print("Loading cactus tree...")
tree <- read.newick("../data/241-mammalian-2020v2.phast-242.nh")
cactus.annotation <- read_xlsx('../data/cactus_ucsc_species_list.xlsx')
cactus.annotation$Species <- gsub(' ','_',cactus.annotation$Species)
cactus.annotation$label <- cactus.annotation$Species
print("Finish loading cactus tree")

#### evo model driven activity
print("Loading evo model dataset...")
evo_model_mean_act <- read.delim("../data/PYS2_evo_model_MPRA_mean_act_20250214.txt.gz", sep = '\t')
colnames(evo_model_mean_act)[1] <- 'tile_name'
evo_model_twist_order <- read.delim('../data/twist_in_silico_evolution_model_driven_order_wpos_phyloP.txt', sep = '\t') %>% select(-oligo_seq, -oligo_seq_width, -step)

evo_model_mean_act <- left_join(evo_model_mean_act, evo_model_twist_order, 
                                   by =c('tile_name'))
evo_model_mean_act <- evo_model_mean_act
evo_model_mean_act$CRE <- factor(evo_model_mean_act$CRE, levels = c("Gata4_chr14_5729","Epas1_chr17_10063",
                                                                   "Lama1_chr17_7784","Bend5_chr4_8201","Sparc_chr11_7211"))
evo_model_mean_act$ref_start <- sapply(strsplit(gsub("[()]", "", evo_model_mean_act$ref_start), ", "), function(x) min(as.numeric(x)))
evo_model_mean_act$ref_end <- sapply(strsplit(gsub("[()]", "", evo_model_mean_act$ref_end), ", "), function(x) max(as.numeric(x)))
evo_model_mean_act$sigma_MPRA_act <- evo_model_mean_act$max_MPRA_act - evo_model_mean_act$min_MPRA_act
print("Finish loading evo model dataset")
                                     
                                     
select.images <- Sys.glob('../imgs/*.png')
select.images <- data.frame(uid = select.images,
                           label = c("Psammomys_obesus","fullTreeAnc237","fullTreeAnc223","Felis_catus","Canis_lupus_familiaris","fullTreeAnc115",'fullTreeAnc239',
                                    "Equus_caballus","Homo_sapiens","Mus_musculus","fullTreeAnc110","Rattus_norvegicus","fullTreeAnc68"),
                           name = c('Psammomys_obesus',"Laurasiatheria","Carnivora","Felis_catus","Canis_lupus_familiaris","Euarchontoglires",
                                   "First_Ancestor","Equus_caballus","Homo_sapiens","Mus_musculus","Primates","Rattus_norvegicus","Rodentia"))
                                     
#### plot mpra act
source("utils/evo_utils.R")
ref.events <- evo_model_mean_act %>% filter(Event == 'ref1') 

target.events <- evo_model_mean_act %>% filter(Event == 'target2')
                                     
# evo_model_mean_act <- calc_step_changes(evo_model_mean_act, target.events)

traj_oi_list <- c("Gata4_chr14_5729__fullTreeAnc239__Mus_musculus", "Epas1_chr17_10063__fullTreeAnc239__Mus_musculus",
                 "Lama1_chr17_7784__fullTreeAnc239__Mus_musculus")
# model_mpra_plot <- plot_group_mpra_traj(evo_model_mean_act, traj_oi_list,"fullTreeAnc239", ref.events, target.events, ncol = 3)
print("Loading evo random model dataset...")                       
evo_random_mean_act <- read.delim("../data/PYS2_evo_random_MPRA_mean_act_20250214.txt.gz", sep = '\t')
colnames(evo_random_mean_act)[1] <- 'tile_name'
evo_random_twist_order <- read.delim('../data/twist_in_silico_evolution_random_walk_order_wpos_phyloP_map_model_ranks.txt', sep = '\t') %>% select( -oligo_seq, -oligo_seq_width, -width)

evo_random_mean_act <- left_join(evo_random_mean_act, evo_random_twist_order, 
                                   by =c('tile_name'))
colnames(evo_random_mean_act)[29] <- 'model_rank'
colnames(evo_random_mean_act)[14] <- 'rank'
                                     

evo_random_mean_act$CRE <- factor(evo_random_mean_act$CRE, levels = c("Gata4_chr14_5729","Epas1_chr17_10063",
                                                                   "Lama1_chr17_7784","Bend5_chr4_8201","Sparc_chr11_7211"))
evo_random_mean_act$ref_start <- sapply(strsplit(gsub("[()]", "", evo_random_mean_act$ref_start), ", "), function(x) min(as.numeric(x)))
evo_random_mean_act$ref_end <- sapply(strsplit(gsub("[()]", "", evo_random_mean_act$ref_end), ", "), function(x) max(as.numeric(x)))
evo_random_mean_act$sigma_MPRA_act <- evo_random_mean_act$max_MPRA_act - evo_random_mean_act$min_MPRA_act

                                      
new_evo_random_mean_act <- evo_random_mean_act %>% filter(rank == 0)
for(i in seq(1:max(evo_random_mean_act[complete.cases(evo_random_mean_act),]$iter))){
    tmp <- evo_random_mean_act %>% filter(iter == i)
#     tmp <- calc_step_changes(tmp, target.events)
    new_evo_random_mean_act <- bind_rows(new_evo_random_mean_act, tmp)
}

# random_mpra_plot <- plot_group_mpra_traj(new_evo_random_mean_act, traj_oi_list,"fullTreeAnc239", ref.events, target.events, ncol =3,type = 'random_avg') 
evo_model_mean_act.select <- evo_model_mean_act %>% filter(traj%in%traj_oi_list) %>% filter((rank != 0 & !startsWith(Event, 'target')))
evo_model_mean_act.select <- evo_model_mean_act.select %>% 
            mutate(
    CRE = sub("_", ":", CRE)  # Replace underscores with colons
  )
evo_model_mean_act.select$CRE <- factor(evo_model_mean_act.select$CRE, levels = cre_levels)
                                      
new_evo_random_mean_act.select <- new_evo_random_mean_act %>% filter(traj%in%traj_oi_list) %>% filter((rank != 0 & !startsWith(Event, 'target')))
new_evo_random_mean_act.select <- new_evo_random_mean_act.select %>% 
            mutate(
    CRE = sub("_", ":", CRE)  # Replace underscores with colons
  )
new_evo_random_mean_act.select$CRE <- factor(new_evo_random_mean_act.select$CRE, levels = cre_levels)
                                      
reference_oi <- unique(evo_model_mean_act$reference)
CRE_oi <- unique(evo_model_mean_act$CRE)
select.ref <- ref.events %>% filter(traj%in% traj_oi_list)
select.target <- target.events %>% filter(traj%in% traj_oi_list) 
select.ref <- select.ref %>% 
            mutate(
    CRE = sub("_", ":", CRE)  # Replace underscores with colons
  )
select.ref$CRE <- factor(select.ref$CRE, levels = cre_levels)
select.target <- select.target %>% 
            mutate(
    CRE = sub("_", ":", CRE)  # Replace underscores with colons
  )
select.target$CRE <- factor(select.target$CRE, levels = cre_levels)

print("Finish loading evo random model dataset")
                                      
#### load/process oCRE MPRA data
print("Loading oCRE model dataset...")
# Prefer the pre-computed winsorized per-CRE activity table (ships in the
# repo; ~500 KB). Fall back to recomputing from the raw BC-level count
# table (~62 MB, gitignored) if the pre-computed file is not present.
if (file.exists("../data/oCRE_and_CRE_optimization_MPRA_act.txt.gz")) {
    df_mpra_act <- read.table("../data/oCRE_and_CRE_optimization_MPRA_act.txt.gz",
                              header = TRUE, sep = "\t")
} else {
    df_counts_w_CREs2 <- read.table("../data/oCRE_and_CRE_optimization_MPRA_count_table.txt.gz",
                                    header = TRUE)
    coverage_lib <- calculate_MPRA_coverage(df_counts_w_CREs2)
    df_wide <- generate_MPRA_wide_table(df_counts_w_CREs2, coverage_lib)
    df_mpra_act <- generate_winsorize_MPRA_table(df_wide)
}

#### extract control MPRA data
control_mean_act <- df_mpra_act %>% filter(class%in%c('neg_control','pos_control')) %>% group_by(CRE_id) %>% 
        summarise(mean_MPRA_act = mean(MPRA_act), sd_MPRA_act = sd(MPRA_act))

#### extract oCRE MPRA data
oCRE_mean_act <- df_mpra_act %>% filter(class == 'oCRE') %>% group_by(CRE_id) %>% 
        summarise(mean_MPRA_act = mean(MPRA_act), sd_MPRA_act = sd(MPRA_act))
oCRE_acc_pred <- read.delim("../data/oCRE_marginal_footprint.txt",
                           sep = '\t')

oCRE_mean_act <- oCRE_mean_act %>% separate(CRE_id, c('tile_id','something'), sep = "_20240416") %>% dplyr::select(-something)
oCRE_mean_act <- left_join(oCRE_mean_act, oCRE_acc_pred, by = c("tile_id"))
oCRE_mean_act <- oCRE_mean_act %>% separate(tile_id, c('full_CRE_id','species'), sep = "__")

### get mouse WT controls
df_mouse_tiles_mpra <- oCRE_mean_act %>% filter(species == 'Mus_musculus') 

dist.mat <- read.delim("../data/JB_ortholog_CREs_corrected_pairwise_pid.txt.gz",
                      sep = '\t')
                                      

dist.mat <- dist.mat %>% filter(species_x == 'Mus_musculus') %>% select(species = species_y, PID, full_CRE_id = CRE)
                                                             
oCRE_mean_act.x <- left_join(oCRE_mean_act, dist.mat, by = c("species",'full_CRE_id'))
print("Finish loading oCRE model dataset")                 
### get ancestral mouse lineage order
species.pair <- c("Mus_musculus",'fullTreeAnc239')
nodes.keep <- make.nodepath(tree, species.pair)
### filter nodes if in activity
nodes.keep <- nodes.keep[nodes.keep%in%unique(oCRE_mean_act$species)]
nodes.order <- seq(1:length(nodes.keep))
names(nodes.order) <- rev(nodes.keep)


oCRE_mean_select <- oCRE_mean_act.x %>% filter(species%in%nodes.keep) %>% 
                                      filter(full_CRE_id%in%c("Epas1_chr17_10063","Gata4_chr14_5729","Lama1_chr17_7784")) %>% 
                                      select(CRE = full_CRE_id, species, norm_footprint = norm_score, PID, mean_MPRA_act, sd_MPRA_act)
oCRE_mean_select <- oCRE_mean_select %>% 
            mutate(
    CRE = sub("_", ":", CRE)  # Replace underscores with colons
  )
oCRE_mean_select$CRE <- factor(oCRE_mean_select$CRE, levels = cre_levels)
    
final_TFs <- c("Jun_Atf3","Foxa2","Gata4/6","Sox17","Klf4","Hnf1b")
col_TFs <- c("firebrick2","forestgreen","darkorange2","dodgerblue1","mediumorchid3",'goldenrod1')
names(col_TFs) <- final_TFs 

df_TF_max_aff <- read.table("../data/ParEndo_markerTF_cuttoffs.txt",
                            header=TRUE)
print("Loading evo model TFBS dataset...")
final_evo.tfbs <- read.delim("../data/in_silico_evolution_model_driven_order.txt.gz",
                      sep = '\t')
# joining Gata4/6 calls to avoid redundancy
gata.tfbs <- final_evo.tfbs %>% filter(TF_name %in% c("Gata4", "Gata6")) 
gata.tfbs <- gata.tfbs %>% mutate(TF_name = "Gata4/6") %>%
  group_by(start_pos, end_pos, CRE) %>%
  dplyr::slice_max(order_by = norm_affinity, with_ties = FALSE) %>%
  ungroup()
final_evo.tfbs <- bind_rows(final_evo.tfbs %>% filter(!TF_name %in% c("Gata4", "Gata6")),
                         gata.tfbs) 
final_evo.tfbs <- final_evo.tfbs %>% filter(TF_name%in%final_TFs)
final_evo.tfbs$start_pos <- final_evo.tfbs$start_pos + 1 ### fix 0-index position from probound
colnames(final_evo.tfbs)[8] <- 'tile_name'
final_evo.tfbs <- left_join(final_evo.tfbs, evo_model_twist_order, 
                                   by =c('tile_name')) %>% select(-seq) 
final_evo.tfbs <-  final_evo.tfbs %>% 
            mutate(
    CRE = sub("_", ":", CRE)  # Replace underscores with colons
  )
print("Finish loading evo model TFBS dataset")
                                                                      
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

mouse_DMS_MPRA <- read.delim("../data/mouse_DMS_MPRA_data.txt", sep = '\t')
### need to back calculate the wild-type
mouse_DMS_MPRA2 <- mouse_DMS_MPRA %>% mutate(WT_MPRA_act = MPRA_act/(2^log2FC_MPRA))
wt_DMS_MPRA <- data.frame()
                                      
for(CRE_oi in unique(mouse_DMS_MPRA2$CRE)){
    tmp <- mouse_DMS_MPRA2 %>% filter(CRE == CRE_oi)
    tmp <- tmp[complete.cases(tmp$WT_MPRA_act),]
    tmp <- tmp[!is.infinite(tmp$WT_MPRA_act),]
    wt_mean <- mean(tmp$WT_MPRA_act)
    wt_sd <- sd(tmp$WT_MPRA_act)
    wt_DMS_MPRA <- bind_rows(wt_DMS_MPRA, data.frame(CRE = CRE_oi, WT_mean_MPRA_act = wt_mean, WT_sd_MPRA_act = wt_sd))
}
mouse_DMS_MPRA2 <- left_join(mouse_DMS_MPRA2, wt_DMS_MPRA, by = c("CRE"))                                      
mouse_DMS_MPRA2 <-  mouse_DMS_MPRA2 %>% 
            mutate(
    CRE = sub("_", ":", CRE)  # Replace underscores with colons
  )                                    
                                      
                                      
