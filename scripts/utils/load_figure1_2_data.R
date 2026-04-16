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
source("utils/mpra_tablemaker.R")
source("utils/make_tree_functions.R")
source("utils/trajectory_utils.R")
source("utils/figure1_func.R")

print("Loading cactus tree...")
tree <- read.newick("../data/241-mammalian-2020v2.phast-242.nh")
cactus.annotation <- read_xlsx('../data/cactus_ucsc_species_list.xlsx')
cactus.annotation$Species <- gsub(' ','_',cactus.annotation$Species)
cactus.annotation$label <- cactus.annotation$Species

big.tree <- as_tibble(tree)
big.tree <- groupClade(big.tree, c(293, 247,353,383,429,346,459))
big.tree <- as.treedata(big.tree)

select.images <- Sys.glob('../imgs/*.png')
select.images <- data.frame(uid = select.images,
                           label = c("Psammomys_obesus","fullTreeAnc237","fullTreeAnc223","Felis_catus","Lycaon_pictus","fullTreeAnc115",'fullTreeAnc239',
                                    "Equus_caballus","Homo_sapiens","Mus_musculus","fullTreeAnc110","Rattus_norvegicus","fullTreeAnc68"),
                           name = c('Sand_rat',"Laurasiatheria","Carnivora","Cat","African_dog","Euarchontoglires",
                                   "First_Ancestor","Horse","Human","Mouse","Primates","Rat","Rodentia"))
filter.images <- select.images %>% filter(label%in%c("fullTreeAnc239","Mus_musculus"))
rownames(filter.images) <- filter.images$label
filter.images <- left_join(big.tree, filter.images, by = 'label')
filter.images <- as_tibble(filter.images)
print("Finish loading cactus tree")

print("Loading phyloP data...")
phyloP.df <- read.delim("../data/DMS_300bp_tiles_phyloP.bed",
                       sep = '\t', header = F)
colnames(phyloP.df) <- c("chr",'start','end','full_CRE_id','phyloP')
singleBase.phyloP <- read.delim('../data/DMS_300bp_tiles_1bp_size_phyloP.bed',
                       sep = '\t', header = F) %>% separate(V4, c("CRE_id",'chr','id','pos'), "_") %>% 
                        mutate(CRE = paste0(CRE_id,"_",chr,"_",id)) %>% dplyr::select(-CRE_id,-chr,-id)
colnames(singleBase.phyloP) <- c("chr",'start','end','pos','phyloP','CRE')
singleBase.phyloP$pos <- as.numeric(singleBase.phyloP$pos)

adjusted.phyloP <- data.frame()
for(cre_oi in unique(singleBase.phyloP$CRE)){
    tmp.cre <- singleBase.phyloP %>% filter(CRE == cre_oi)
    if(cre_oi%in%c("Sparc_chr11_7211","Gata4_chr14_5729")){ ### reverse complement in mpra 
        tmp.cre <- tmp.cre %>%
              mutate(pos = max(pos) + 1 - pos)
        adjusted.phyloP <- bind_rows(adjusted.phyloP, tmp.cre)
    }else{
        adjusted.phyloP <- bind_rows(adjusted.phyloP, tmp.cre)
    }
    
}

singleBase.phyloP <- adjusted.phyloP

print("Finish loading phyloP dataframe")

print("Loading MSA dataframe...")
### get msa distance
max.dist <- data.frame()
dist.matrices <- list.files('../data/msa_dist_matrix', full.names = T)

for(d in dist.matrices){
    cre <- strsplit(basename(d),"_oCRE")[[1]][1]
    temp <- read.delim(d, sep = '\t', row.names = 1)
    temp <- temp[grepl("Mus_musculus", rownames(temp)),]
    temp <- melt(temp)
    temp$cre <- rep(cre, nrow(temp))
    max.dist <- bind_rows(max.dist, temp)
}

colnames(max.dist) <- c('species','seqid', 'full_CRE_id')                                                             
print("Finish loading MSA dataframe.")



print("Loading oCRE model dataset...") 
#### load/process MPRA data
df_counts_w_CREs2 <- read.table("../data/oCRE_and_CRE_optimization_MPRA_count_table.txt.gz",
                                header=TRUE)

coverage_lib <- calculate_MPRA_coverage(df_counts_w_CREs2)
df_wide <- generate_MPRA_wide_table(df_counts_w_CREs2, coverage_lib)
df_mpra_act <- generate_winsorize_MPRA_table(df_wide)

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
mm <- match(oCRE_mean_act$species, cactus.annotation$Species)
oCRE_mean_act$common_name <- cactus.annotation$`Common Name`[mm]
oCRE_mean_act$log2_norm_score <- log2(oCRE_mean_act$norm_score)
oCRE_mean_act <- left_join(oCRE_mean_act, max.dist, by =c("full_CRE_id","species"))

### get mouse WT controls
df_mouse_tiles_mpra <- oCRE_mean_act %>% filter(species == 'Mus_musculus') 
print("Finish loading oCRE model dataset") 

#### load TFBS mapping
final_TFs <- c("Jun_Atf3","Foxa2","Gata4/6","Sox17","Klf4","Hnf1b")
col_TFs <- c("firebrick2","forestgreen","darkorange2","dodgerblue1","mediumorchid3",'goldenrod1')
names(col_TFs) <- final_TFs 

print("Loading cactus TFBS dataset...")
df_TF_max_aff <- read.table("../data/ParEndo_markerTF_cuttoffs.txt",
                            header=TRUE)
cactus.tfbs <- read.delim(gzfile("../data/oCRE_ProBound_ParEndo_TFBS_predictions.txt.gz"),
                       sep = '\t')
CRE_oi <- cactus.tfbs %>% pull(CRE) %>% str_split("__") %>% lapply("[[",1) %>% unlist()
species_oi <- cactus.tfbs %>% pull(CRE) %>% str_split("__") %>% lapply("[[",2) %>% unlist() %>%
                str_split("_20240416") %>% lapply("[[",1) %>% unlist()
cactus.tfbs$ori_CRE_id <- CRE_oi
cactus.tfbs$species <- species_oi
cactus.tfbs <- cactus.tfbs %>% mutate(TF_ori_pos = case_when(
                                TFBS_orientation == '+' ~ 1,
                                TRUE ~ -1)) 
# joining Gata4/6 calls to avoid redundancy
gata.tfbs <- cactus.tfbs %>% filter(TF_name %in% c("Gata4", "Gata6")) %>% 
                    mutate(TF_ori_pos = case_when( ### correct for Gata4 and Gata6 orientation issue
                            TF_name == 'Gata6' ~ TF_ori_pos * -1,
                            TRUE ~ TF_ori_pos)) %>% 
                    mutate(TFBS_orientation = case_when(
                            TF_ori_pos == 1 ~ '+',
                            TRUE ~ '-'))
gata.tfbs <- gata.tfbs %>% mutate(TF_name = "Gata4/6") %>%
  group_by(start_pos, end_pos, CRE, ori_CRE_id, species) %>%
  dplyr::slice_max(order_by = norm_affinity, with_ties = FALSE) %>%
  ungroup()
cactus.tfbs <- bind_rows(cactus.tfbs %>% filter(!TF_name %in% c("Gata4", "Gata6")),
                         gata.tfbs) 
cactus.tfbs <- cactus.tfbs %>% filter(TF_name%in%final_TFs)
cactus.tfbs$start_pos <- cactus.tfbs$start_pos + 1 ### fix 0-index position from probound
### load functional mouse TFBS data
mouse_DMS_tfbs <- read.delim("../data/mouse_TFBS_impact_and_affinity_CRM.txt", sep = '\t')
mouse_DMS_tfbs <- mouse_DMS_tfbs %>% mutate(TF_ori_pos = case_when(
                                TFBS_orientation == '+' ~ 1,
                                TRUE ~ -1)) 
# joining Gata4/6 calls to avoid redundancy
gata.tfbs <- mouse_DMS_tfbs %>% filter(TF_name %in% c("Gata4", "Gata6")) %>% 
                    mutate(TF_ori_pos = case_when( ### correct for Gata4 and Gata6 orientation issue
                            TF_name == 'Gata6' ~ TF_ori_pos * -1,
                            TRUE ~ TF_ori_pos)) %>% 
                    mutate(TFBS_orientation = case_when(
                            TF_ori_pos == 1 ~ '+',
                            TRUE ~ '-'))
gata.tfbs <- gata.tfbs %>% mutate(TF_name = "Gata4/6") %>%
  group_by(TFBS_start, TFBS_end, CRE) %>%
  dplyr::slice_max(order_by = norm_affinity, with_ties = FALSE) %>%
  ungroup()
mouse_DMS_tfbs <- bind_rows(mouse_DMS_tfbs %>% filter(!TF_name %in% c("Gata4", "Gata6")) , gata.tfbs) %>% 
                  filter(TF_name%in%final_TFs)
mouse_func_tfbs <- mouse_DMS_tfbs %>% filter(functional_TF)
mouse_DMS_MPRA <- read.delim("../data/mouse_DMS_MPRA_data.txt", sep = '\t') %>% mutate(log2FC_MPRA = case_when(
            is.na(log2FC_MPRA) ~ 0,
            TRUE ~ log2FC_MPRA))
print("Finish loading cactus TFBS dataset")     