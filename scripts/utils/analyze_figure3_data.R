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
source("/net/shendure/vol8/projects/tli/PYS2_oCRE_BC_libs_20240904/mpra_tablemaker.R")
source("/net/shendure/vol8/projects/tli/PYS2_oCRE_BC_libs_20240904/make_tree_functions.R")
source("/net/shendure/vol8/projects/tli/PYS2_oCRE_BC_libs_20240904/trajectory_utils.R")
source("/net/shendure/vol8/projects/tli/PYS2_oCRE_BC_libs_20240904/figure2_func.R")
source("/net/shendure/vol8/projects/tli/PYS2_oCRE_BC_libs_20240904/make_tree_functions.R")
source("/net/shendure/vol8/projects/tli/PYS2_oCRE_BC_libs_20240904/load_figure1_2_data.R")

# Helper function to safely pull a value or return 0
pull_or_zero <- function(df, col) {
  if (nrow(df) == 0) return(0)
  val <- df[[col]]
  if (length(val) == 0) return(0)
  return(val)
}

### make Epas1 Tree
### select best representative species for each clade
CRE_oi <- "Epas1_chr17_10063"
epas1_oCRE_mean_act <- oCRE_mean_act %>%
  filter(full_CRE_id == CRE_oi) %>% left_join(cactus.annotation %>% dplyr::select(species = Species, Order, Clade))
top_order_epas1_act <- epas1_oCRE_mean_act %>% group_by(Order) %>% filter(mean_MPRA_act == max(mean_MPRA_act)) %>% 
                dplyr::select(Order, species, mean_MPRA_act) %>% arrange(-mean_MPRA_act)

tibble.tree <- as_tibble(tree)
shape.values <- c(1,24,13,15)
names(shape.values) <- c('mouse','no ortholog','not sequenced', 'sequenced')
sub.seq <- oCRE_mean_act %>% filter(full_CRE_id == 'Epas1_chr17_10063') %>%
        dplyr::select(species,mean_MPRA_act,common_name)
sub.seq$label <- sub.seq$species
sub.seq <- sub.seq %>% mutate(capture = case_when(
                        species == 'Mus_musculus' ~ 'mouse',
                        is.na(mean_MPRA_act) ~ 'not sequenced',
                        TRUE ~ 'sequenced'))
sub.tree <- full_join(tibble.tree, sub.seq, by = 'label')
sub.tree <- sub.tree %>% mutate(capture = case_when(
                        is.na(capture) ~ 'no ortholog',
                        TRUE ~ capture))
anc.keep <- c('fullTreeAnc110point5','fullTreeAnc238','fullTreeAnc68','fullTreeAnc194','fullTreeAnc223','fullTreeAnc79')
sub.tips <- sub.tree %>% filter(label%in%c('Mus_musculus','Homo_sapiens','Lycaon_pictus','Giraffa_tippelskirchi','Eulemur_fulvus',
                                          anc.keep)) %>%
            mutate(Shape = case_when(
                label == 'Mus_musculus' ~ 1,
                label == 'Homo_sapiens' ~ 1,
                label == 'Lycaon_pictus' ~1,
                label == 'Eulemur_fulvus' ~1,
                label == 'Giraffa_tippelskirchi'~1,
                label%in%anc.keep ~ 9,
                TRUE ~ 2)) %>% 
            mutate(Size = case_when(
                label == 'Mus_musculus' ~ 3,
                label == 'Homo_sapiens' ~ 3,
                label == 'Lycaon_pictus' ~3,
                label == 'Eulemur_fulvus' ~3,
                label == 'Giraffa_tippelskirchi' ~ 3,
                label%in%anc.keep ~ 3,
                TRUE ~ 0))
sub.tips$Shape <- factor(sub.tips$Shape)
sub.tips <- as.tibble(sub.tips)
sub.tree <- as.treedata(sub.tree)

### make affinity heatmap
mouse_DMS_MPRA$mut <- toupper(mouse_DMS_MPRA$mut)
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

### format aff mat for gheatmap table
tree_species <- tibble.tree %>%
  filter(!startsWith(label, "full")) %>%
  dplyr::select(species = label)

aff.mat <- epas1_msa_list$msa_df_TFBS %>% filter(!is.na(TF_name2)) %>% 
            dplyr::select(species, TF_name,TF_name2, norm_affinity) %>%
            group_by(TF_name2) %>% ### normalize by functional TFBS
            mutate(norm_affinity_scaled = (norm_affinity - min(norm_affinity, na.rm = TRUE)) /
                                            (max(norm_affinity, na.rm = TRUE) - min(norm_affinity, na.rm = TRUE))) %>%
            ungroup()
aff.mat$hue_color <- mapply(
  function(tf, affinity) {
    if (is.na(affinity)) {
      # if no affinity, return grey
      return("grey80")
    } else {
      tf_color <- col_TFs[tf]   # lookup TF-specific base color
      # interpolate: scale affinity (0–1) to 1–100 index
      pal <- colorRampPalette(c("white", tf_color))(100)
      pal[round(affinity * 99) + 1]
    }
  },
  aff.mat$TF_name,  # Transcription factor name
  aff.mat$norm_affinity_scaled  # Normalized affinity
)

### calculate 3' end AP-1 activity
ap1.mat <- epas1_msa_list$msa_df_TFBS %>% filter(TF_name == 'Jun_Atf3') %>% 
            filter(msa_start > 700) %>% group_by(species, TF_name) %>% 
            summarise(avg_norm_affinity = mean(norm_affinity), sd_norm_affinity = sd(norm_affinity)) %>% 
            mutate(TF_name2 = "3' AP-1") %>% 
            dplyr::select(species, TF_name,TF_name2, norm_affinity = avg_norm_affinity) %>%
            group_by(TF_name2) %>%  ### normalize by functional TFBS
            mutate(norm_affinity_scaled = (norm_affinity - min(norm_affinity, na.rm = TRUE)) /
                                            (max(norm_affinity, na.rm = TRUE) - min(norm_affinity, na.rm = TRUE))) %>%
            ungroup()

ap1.mat$hue_color <- mapply(
  function(tf, affinity) {
    if (is.na(affinity)) {
      # if no affinity, return grey
      return("grey80")
    } else {
      tf_color <- col_TFs[tf]   # lookup TF-specific base color
      # interpolate: scale affinity (0–1) to 1–100 index
      pal <- colorRampPalette(c("white", tf_color))(100)
      pal[round(affinity * 99) + 1]
    }
  },
  ap1.mat$TF_name,  # Transcription factor name
  ap1.mat$norm_affinity_scaled  # Normalized affinity
)

### calculate 5' end Klf4 activity
klf.mat <- epas1_msa_list$msa_df_TFBS %>% filter(TF_name == 'Klf4') %>% 
            filter(msa_start < 100) %>% group_by(species, TF_name) %>% 
            summarise(avg_norm_affinity = mean(norm_affinity), sd_norm_affinity = sd(norm_affinity)) %>% 
            mutate(TF_name2 = "5' Klf4") %>% 
            dplyr::select(species, TF_name,TF_name2, norm_affinity = avg_norm_affinity) %>%
            group_by(TF_name2) %>%  ### normalize by functional TFBS
            mutate(norm_affinity_scaled = (norm_affinity - min(norm_affinity, na.rm = TRUE)) /
                                            (max(norm_affinity, na.rm = TRUE) - min(norm_affinity, na.rm = TRUE))) %>%
            ungroup()

extant_aff.mat <- aff.mat %>% mutate(species = as.character(species)) %>% filter(!startsWith(species, "full"))
extant_ap1.mat <- ap1.mat %>% mutate(species = as.character(species)) %>% filter(!startsWith(species, "full"))
extant_klf4.mat <- klf.mat %>% mutate(species = as.character(species)) %>% filter(!startsWith(species, "full"))

# join so that every species × TF_name2 combination exists
aff.mat.tree <- tree_species %>%
  crossing(TF_name2 = unique(extant_aff.mat$TF_name2)) %>%
  left_join(aff.mat, by = c("species", "TF_name2")) %>% filter(TF_name2%in%c('Foxa2_3','Sox17_1','Gata4/6_4',
                                                                            'Jun_Atf3_8','Jun_Atf3_2','Jun_Atf3_5',"3' AP-1")) %>%
  mutate(TF_name2 = factor(TF_name2, levels = c('Jun_Atf3_8','Jun_Atf3_2','Jun_Atf3_5',
                                               'Foxa2_3','Sox17_1','Gata4/6_4',"3' AP-1"))) %>%
  dplyr::select(label = species, TF_name2, TF_name, norm_affinity_scaled)

aff.mat.tree$hue_color <- mapply(
      function(tf, affinity) {
        if (is.na(affinity)) {
          # if no affinity, return grey
          return("grey80")
        } else {
          tf_color <- col_TFs[tf]   # lookup TF-specific base color
          # interpolate: scale affinity (0–1) to 1–100 index
          pal <- colorRampPalette(c("white", tf_color))(100)
          pal[round(affinity * 99) + 1]
        }
      },
      aff.mat.tree$TF_name,  # Transcription factor name
      aff.mat.tree$norm_affinity_scaled  # Normalized affinity
    )

ap1.mat.tree <- tree_species %>%
  crossing(TF_name2 = unique(extant_ap1.mat$TF_name2)) %>%
  left_join(extant_ap1.mat, by = c("species", "TF_name2")) %>%
  dplyr::select(label = species, TF_name2, TF_name, norm_affinity_scaled)

ap1.mat.tree$hue_color <- mapply(
      function(tf, affinity) {
        if (is.na(affinity)) {
          # if no affinity, return grey
          return("grey80")
        } else {
          tf_color <- col_TFs[tf]   # lookup TF-specific base color
          # interpolate: scale affinity (0–1) to 1–100 index
          pal <- colorRampPalette(c("white", tf_color))(100)
          pal[round(affinity * 99) + 1]
        }
      },
      ap1.mat.tree$TF_name,  # Transcription factor name
      ap1.mat.tree$norm_affinity_scaled  # Normalized affinity
    )

klf4.mat.tree <- tree_species %>%
  crossing(TF_name2 = unique(extant_klf4.mat$TF_name2)) %>%
  left_join(extant_klf4.mat, by = c("species", "TF_name2")) %>%
  dplyr::select(label = species, TF_name2, TF_name, norm_affinity_scaled)

klf4.mat.tree$hue_color <- mapply(
      function(tf, affinity) {
        if (is.na(affinity)) {
          # if no affinity, return grey
          return("grey80")
        } else {
          tf_color <- col_TFs[tf]   # lookup TF-specific base color
          # interpolate: scale affinity (0–1) to 1–100 index
          pal <- colorRampPalette(c("white", tf_color))(100)
          pal[round(affinity * 99) + 1]
        }
      },
      klf4.mat.tree$TF_name,  # Transcription factor name
      klf4.mat.tree$norm_affinity_scaled  # Normalized affinity
    )

tfbs_map <- data.frame(msa_start = c(622,506, 700),
                      msa_end = c(687,565, max(epas1_msa_list$msa_df_TFBS$msa_end)))

common_name.df <- data.frame(species = c('Mus_musculus','Giraffa_tippelskirchi','Homo_sapiens','Eulemur_fulvus','Lycaon_pictus','Psammomys_obesus'),
                            common_name = c('Mouse',"Giraffe",'Human','Brown_lemur','African_dog','Sand_rat'))

get_msa_aff_lineage <- function(species_input, anc_root, add_mouse = F){
    species.pair <- c(species_input, anc_root)
    species.nodes <- make.nodepath(tree, species.pair)
    lineage.epas1_msa_df <- epas1_msa_df %>% filter(species%in%species.nodes)
    root.anc <- lineage.epas1_msa_df %>% filter(species == anc_root) %>% dplyr::select(curpos, anc_mut = mut2)
    lineage.epas1_msa_df <- lineage.epas1_msa_df %>% left_join(root.anc)
    ### update species.nodes in case some are missing
    species.nodes <- species.nodes[species.nodes%in%unique(epas1_msa_df$species)]
    
    fixed_lineage.epas1_msa_df <- data.frame()
    tmp.lineage_msa <- lineage.epas1_msa_df
    # store previous species' mutation set (curpos + mutated base)
    prev_mut_tbl <- NULL
    
    for(species_oi in rev(species.nodes)){
        if(species_oi == anc_root){ ### skip first ancestor
            tmp_species <- tmp.lineage_msa %>% filter(species == species_oi) %>%
                    mutate(match2 = 'Match', log_change = ifelse(mut_type%in%c("DEL","INS"), NA, 0))
            fixed_lineage.epas1_msa_df <- bind_rows(fixed_lineage.epas1_msa_df, tmp_species)
        }else{
            tmp_species <- tmp.lineage_msa %>% filter(species == species_oi)
            changes <- tmp_species %>% filter(mut_type != 'Match') %>% 
                    filter(mut2 != anc_mut) %>% dplyr::select(curpos, log_change = log2FC_MPRA) %>% 
                    mutate(match2 = 'Mutation')
            tmp_species <- tmp_species %>% left_join(changes) %>% 
                        mutate(match2 = case_when(
                            is.na(match2) ~ 'Match',
                            TRUE ~ match2),
                              log_change = case_when(
                             is.na(log_change) ~ ifelse(mut_type%in%c("DEL","INS"), NA, 0),
                              TRUE ~ log_change))
            if(!is.null(prev_mut_tbl) && nrow(prev_mut_tbl) > 0){
                tmp_species <- tmp_species %>%
                    left_join(prev_mut_tbl, by = "curpos") %>%
                    mutate(
                      match2 = case_when(
                        mut2 == prev_mut2 ~ 'Mutation_C',
                        TRUE ~ match2
                      )
                    ) %>%
                    dplyr::select(-prev_mut2)
            }
            fixed_lineage.epas1_msa_df <- bind_rows(fixed_lineage.epas1_msa_df, tmp_species)
            ### update mutations at each step
            mm <- match(tmp.lineage_msa$curpos, tmp_species$curpos)
            tmp.lineage_msa$anc_mut <-  tmp_species$mut2[mm] 
            
            ### update previous species' mutation set
            prev_mut_tbl <- tmp_species %>%
                  filter(match2%in%c('Mutation','Mutation_C')) %>%
                  dplyr::select(curpos, prev_mut2 = mut2)
        }
    }
    fixed_lineage.epas1_msa_df <- fixed_lineage.epas1_msa_df %>% mutate(log_change = case_when(
                                            mut_type%in%c('DEL','INS') ~ NA,
                                            TRUE ~ log_change))
    common_name <- common_name.df %>% filter(species == species_input) %>% pull(common_name)
    fixed_lineage.epas1_msa_df <- fixed_lineage.epas1_msa_df %>% 
                mutate(species2 = ifelse(species == species_input, common_name, str_replace(species, "fullTree","")))
    if(add_mouse){
        epas1_mouse_msa_df <- epas1_msa_df %>% filter(species == 'Mus_musculus')
        epas1_mouse_msa_df$log_change = 0
        epas1_mouse_msa_df$match2 = 'WT'
        epas1_mouse_msa_df$species2 <- 'Mouse'
        fixed_lineage.epas1_msa_df <- bind_rows(fixed_lineage.epas1_msa_df, epas1_mouse_msa_df)
        species.nodes2 <- c(species.nodes, 'Mus_musculus')
        fixed_lineage.epas1_msa_df$species2 <- factor(fixed_lineage.epas1_msa_df$species2, levels = c(
            ifelse(species.nodes2 == species_input, common_name, 
                   ifelse(species.nodes2 == 'Mus_musculus','Mouse' ,str_replace(species.nodes2, "fullTree","")))
          ))
        
    }else{
        fixed_lineage.epas1_msa_df$species2 <- factor(fixed_lineage.epas1_msa_df$species2, levels = c(
            ifelse(species.nodes == species_input, common_name, str_replace(species.nodes, "fullTree",""))
          ))
    }
    
    if(add_mouse){
        lineage_aff.mat <- epas1_msa_list$msa_df_TFBS %>% 
                filter(!is.na(TF_name2)) %>% filter(species%in%c(species.nodes,'Mus_musculus'))
    }else{
        lineage_aff.mat <- epas1_msa_list$msa_df_TFBS %>% 
                filter(!is.na(TF_name2)) %>% filter(species%in%species.nodes)
    }
    lineage_aff.mat <- lineage_aff.mat %>% left_join(aff.mat %>% dplyr::select(species,TF_name,TF_name2, hue_color), 
                                                     by = c('species','TF_name2','TF_name'))
    new_lineage_aff.mat <- data.frame()

    for(tf2 in unique(lineage_aff.mat$TF_name2)){
        tmp.mat <- lineage_aff.mat %>% filter(TF_name2 == tf2)
        fc_steps <- c(0)
        delta_steps <- c(0)
        initial_step <- tmp.mat %>% filter(species == anc_root) %>% pull(norm_affinity)
        for(species_oi in rev(species.nodes)){
            if(species_oi == anc_root){
            }else{
                current_step <- tmp.mat %>% filter(species == species_oi) %>% pull(norm_affinity)
                fc = current_step/initial_step
                delta <- current_step - initial_step
                fc_steps <- c(fc_steps, fc)
                delta_steps <- c(delta_steps, delta)
                initial_step <- current_step
            }
        }
        names(delta_steps) <- rev(species.nodes)
        names(fc_steps) <- rev(species.nodes)
        mm <- match(tmp.mat$species, names(delta_steps))
        tmp.mat$delta_steps <- delta_steps[mm]
        tmp.mat$fc_steps <- fc_steps[mm]
        start_aff <- tmp.mat %>% filter(species == anc_root) %>% pull(norm_affinity)
        wt_aff <- tmp.mat %>% filter(species == species_input) %>% pull(norm_affinity)
        tmp.mat <- tmp.mat %>% mutate(frac_delta = delta_steps/(wt_aff - start_aff))
        new_lineage_aff.mat <- bind_rows(new_lineage_aff.mat, tmp.mat)
    }
    new_lineage_aff.mat <- new_lineage_aff.mat %>% mutate(log2_fc_steps = ifelse(fc_steps == 0, 0, log2(fc_steps)))
    new_lineage_aff.mat <- new_lineage_aff.mat %>% 
                mutate(species2 = ifelse(species == species_input, common_name, str_replace(species, "fullTree","")))
    if(add_mouse){
        mouse_aff.mat <- lineage_aff.mat %>% filter(species == 'Mus_musculus')
        mouse_aff.mat$species2 <- 'Mouse'
        mouse_aff.mat$log2_fc_steps <- 0
        new_lineage_aff.mat <- bind_rows(new_lineage_aff.mat, mouse_aff.mat)
        new_lineage_aff.mat$species2 <- factor(new_lineage_aff.mat$species2, levels = c(
            ifelse(species.nodes2 == species_input, common_name, 
                   ifelse(species.nodes2 == 'Mus_musculus','Mouse' ,str_replace(species.nodes2, "fullTree","")))
            ))
        new_lineage_aff.mat <- new_lineage_aff.mat %>% filter(!is.na(species2))
    } else{
        new_lineage_aff.mat$species2 <- factor(new_lineage_aff.mat$species2, levels = c(
            ifelse(species.nodes == species_input, common_name, 
                   ifelse(species.nodes == 'Mus_musculus',"",str_replace(species.nodes, "fullTree","")))
          ))
    }
    
    return(list(msa_df = fixed_lineage.epas1_msa_df, msa_func_aff_TFBS = new_lineage_aff.mat))
}

### Giraffe lineage msa
giraffe.epas1_msa_info <- get_msa_aff_lineage("Giraffa_tippelskirchi",'fullTreeAnc238', add_mouse = T)
### primate lineage msa
species.pair <- c('Homo_sapiens', 'fullTreeAnc110')
primates.nodes <- make.nodepath(tree, species.pair)
primates.nodes2 <- c(primates.nodes, 'Mus_musculus')
primate.epas1_msa_info <- get_msa_aff_lineage("Homo_sapiens",'fullTreeAnc238', add_mouse = T)
primate.epas1_msa_info$msa_df <- primate.epas1_msa_info$msa_df %>% 
                        filter(species%in%primates.nodes2) %>% 
                        mutate(species2 = factor(species2, levels = 
                                                 ifelse(primates.nodes2 == 'Homo_sapiens', 'Human',
                                                        ifelse(primates.nodes2 == 'Mus_musculus', 'Mouse',
                                                              str_replace(primates.nodes2, "fullTree","")))))
primate.epas1_msa_info$msa_func_aff_TFBS <- primate.epas1_msa_info$msa_func_aff_TFBS %>% 
                        filter(species%in%primates.nodes2) %>% 
                        mutate(species2 = factor(species2, levels = 
                                                 ifelse(primates.nodes2 == 'Homo_sapiens', 'Human',
                                                        ifelse(primates.nodes2 == 'Mus_musculus', 'Mouse',
                                                              str_replace(primates.nodes2, "fullTree","")))))

### lemur lineage msa
lemur.epas1_msa_info <- get_msa_aff_lineage("Eulemur_fulvus",'fullTreeAnc110', add_mouse = T)
### dog lineage
dog.epas1_msa_info <- get_msa_aff_lineage("Lycaon_pictus",'fullTreeAnc239', add_mouse = T)
species.pair <- c('Lycaon_pictus', 'fullTreeAnc239')
dog.nodes <- make.nodepath(tree, species.pair)
dog.nodes2 <- c(dog.nodes, 'Mus_musculus')
dog.epas1_msa_info$msa_df <- dog.epas1_msa_info$msa_df %>% 
                        filter(species%in%dog.nodes2) %>% 
                        mutate(species2 = factor(species2, levels = 
                                                 ifelse(dog.nodes2 == 'Lycaon_pictus', 'African_dog',
                                                        ifelse(dog.nodes2 == 'Mus_musculus', 'Mouse',
                                                              str_replace(dog.nodes2, "fullTree","")))))
dog.epas1_msa_info$msa_func_aff_TFBS <- dog.epas1_msa_info$msa_func_aff_TFBS %>% 
                        filter(species%in%dog.nodes2) %>% 
                        mutate(species2 = factor(species2, levels = 
                                                 ifelse(dog.nodes2 == 'Lycaon_pictus', 'African_dog',
                                                        ifelse(dog.nodes2 == 'Mus_musculus', 'Mouse',
                                                              str_replace(dog.nodes2, "fullTree","")))))
###sand rat lineage
sand_rat.epas1_msa_info <- get_msa_aff_lineage("Psammomys_obesus",'fullTreeAnc239', add_mouse = T)
species.pair <- c('Psammomys_obesus', 'fullTreeAnc239')
sand_rat.nodes <- make.nodepath(tree, species.pair)
sand_rat.nodes2 <- c(sand_rat.nodes, 'Mus_musculus')
sand_rat.epas1_msa_info$msa_df <- sand_rat.epas1_msa_info$msa_df %>% 
                        filter(species%in%sand_rat.nodes2) %>% 
                        mutate(species2 = factor(species2, levels = 
                                                 ifelse(sand_rat.nodes2 == 'Psammomys_obesus', 'Sand_rat',
                                                        ifelse(sand_rat.nodes2 == 'Mus_musculus', 'Mouse',
                                                              str_replace(sand_rat.nodes2, "fullTree","")))))
sand_rat.epas1_msa_info$msa_func_aff_TFBS <- sand_rat.epas1_msa_info$msa_func_aff_TFBS %>% 
                        filter(species%in%sand_rat.nodes2) %>% 
                        mutate(species2 = factor(species2, levels = 
                                                 ifelse(sand_rat.nodes2 == 'Psammomys_obesus', 'Sand_rat',
                                                        ifelse(sand_rat.nodes2 == 'Mus_musculus', 'Mouse',
                                                              str_replace(sand_rat.nodes2, "fullTree","")))))
### make dog lineage 3' AP-1 affinity
species.pair <- c("Lycaon_pictus",'fullTreeAnc239')
species.nodes <- make.nodepath(tree, species.pair)
dog_ap1.mat <- ap1.mat  %>% filter(species%in%species.nodes)
new_dog_ap1.mat <- data.frame()

for(tf2 in unique(dog_ap1.mat$TF_name2)){
    tmp.mat <- dog_ap1.mat %>% filter(TF_name2 == tf2)
    fc_steps <- c(0)
    delta_steps <- c(0)
    initial_step <- tmp.mat %>% filter(species == 'fullTreeAnc239') %>% pull(norm_affinity)
    for(species_oi in rev(species.nodes)){
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
    names(delta_steps) <- rev(species.nodes)
    names(fc_steps) <- rev(species.nodes)
    mm <- match(tmp.mat$species, names(delta_steps))
    tmp.mat$delta_steps <- delta_steps[mm]
    tmp.mat$fc_steps <- fc_steps[mm]
    start_aff <- tmp.mat %>% filter(species == 'fullTreeAnc239') %>% pull(norm_affinity)
    wt_aff <- tmp.mat %>% filter(species == 'Lycaon_pictus') %>% pull(norm_affinity)
    tmp.mat <- tmp.mat %>% mutate(frac_delta = delta_steps/(wt_aff - start_aff))
    new_dog_ap1.mat <- bind_rows(new_dog_ap1.mat, tmp.mat)
}
new_dog_ap1.mat <- new_dog_ap1.mat %>% mutate(log2_fc_steps = ifelse(fc_steps == 0, 0, log2(fc_steps)))
new_dog_ap1.mat <- new_dog_ap1.mat %>% 
            mutate(species2 = ifelse(species == 'Lycaon_pictus', 'African_dog', str_replace(species, "fullTree","")))
species.nodes2 <- c(species.nodes, 'Mus_musculus')
mouse_aff.mat <- ap1.mat %>% filter(species == 'Mus_musculus')
mouse_aff.mat$species2 <- 'Mouse'
mouse_aff.mat$log2_fc_steps <- 0
new_dog_ap1.mat <- bind_rows(new_dog_ap1.mat, mouse_aff.mat)
new_dog_ap1.mat$species2 <- factor(new_dog_ap1.mat$species2, levels = c(
    ifelse(species.nodes2 == 'Lycaon_pictus', 'African_dog', 
           ifelse(species.nodes2 == 'Mus_musculus','Mouse' ,str_replace(species.nodes2, "fullTree","")))
    ))
new_dog_ap1.mat <- new_dog_ap1.mat %>% filter(!is.na(species2)) %>% filter(species%in%dog.nodes2)

# dog_klf4.mat <- klf.mat  %>% filter(species%in%species.nodes)
# new_dog_klf4.mat <- data.frame()

# for(tf2 in unique(dog_klf4.mat$TF_name2)){
#     tmp.mat <- dog_klf4.mat %>% filter(TF_name2 == tf2)
#     fc_steps <- c(0)
#     delta_steps <- c(0)
#     initial_step <- tmp.mat %>% filter(species == 'fullTreeAnc239') %>% pull(norm_affinity)
#     for(species_oi in rev(species.nodes)){
#         if(species_oi == 'fullTreeAnc239'){
#         }else{
#             current_step <- tmp.mat %>% filter(species == species_oi) %>% pull(norm_affinity)
#             fc = current_step/initial_step
#             delta <- current_step - initial_step
#             fc_steps <- c(fc_steps, fc)
#             delta_steps <- c(delta_steps, delta)
#             initial_step <- current_step
#         }
#     }
#     names(delta_steps) <- rev(species.nodes)
#     names(fc_steps) <- rev(species.nodes)
#     mm <- match(tmp.mat$species, names(delta_steps))
#     tmp.mat$delta_steps <- delta_steps[mm]
#     tmp.mat$fc_steps <- fc_steps[mm]
#     start_aff <- tmp.mat %>% filter(species == 'fullTreeAnc239') %>% pull(norm_affinity)
#     wt_aff <- tmp.mat %>% filter(species == 'Lycaon_pictus') %>% pull(norm_affinity)
#     tmp.mat <- tmp.mat %>% mutate(frac_delta = delta_steps/(wt_aff - start_aff))
#     new_dog_klf4.mat <- bind_rows(new_dog_klf4.mat, tmp.mat)
# }
# new_dog_klf4.mat <- new_dog_klf4.mat %>% mutate(log2_fc_steps = ifelse(fc_steps == 0, 0, log2(fc_steps)))
# new_dog_klf4.mat <- new_dog_klf4.mat %>% 
#             mutate(species2 = ifelse(species == 'Lycaon_pictus', 'African_dog', str_replace(species, "fullTree","")))
# species.nodes2 <- c(species.nodes, 'Mus_musculus')
# mouse_aff.mat <- klf.mat %>% filter(species == 'Mus_musculus')
# mouse_aff.mat$species2 <- 'Mouse'
# mouse_aff.mat$log2_fc_steps <- 0
# new_dog_klf4.mat <- bind_rows(new_dog_klf4.mat, mouse_aff.mat)
# new_dog_klf4.mat$species2 <- factor(new_dog_klf4.mat$species2, levels = c(
#     ifelse(species.nodes2 == 'Lycaon_pictus', 'African_dog', 
#            ifelse(species.nodes2 == 'Mus_musculus','Mouse' ,str_replace(species.nodes2, "fullTree","")))
#     ))
# new_dog_klf4.mat <- new_dog_klf4.mat %>% filter(!is.na(species2)) %>% filter(species%in%dog.nodes2)

species.pair <- c("Psammomys_obesus",'fullTreeAnc239')
species.nodes <- make.nodepath(tree, species.pair)
sand_rat_ap1.mat <- ap1.mat  %>% filter(species%in%species.nodes)
new_sand_rat_ap1.mat <- data.frame()
species.nodes <- species.nodes[species.nodes%in%unique(ap1.mat$species)]

for(tf2 in unique(sand_rat_ap1.mat$TF_name2)){
    tmp.mat <- sand_rat_ap1.mat %>% filter(TF_name2 == tf2)
    fc_steps <- c(0)
    delta_steps <- c(0)
    initial_step <- tmp.mat %>% filter(species == 'fullTreeAnc239') %>% pull(norm_affinity)
    for(species_oi in rev(species.nodes)){
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
    names(delta_steps) <- rev(species.nodes)
    names(fc_steps) <- rev(species.nodes)
    mm <- match(tmp.mat$species, names(delta_steps))
    tmp.mat$delta_steps <- delta_steps[mm]
    tmp.mat$fc_steps <- fc_steps[mm]
    start_aff <- tmp.mat %>% filter(species == 'fullTreeAnc239') %>% pull(norm_affinity)
    wt_aff <- tmp.mat %>% filter(species == 'Psammomys_obesus') %>% pull(norm_affinity)
    tmp.mat <- tmp.mat %>% mutate(frac_delta = delta_steps/(wt_aff - start_aff))
    new_sand_rat_ap1.mat <- bind_rows(new_sand_rat_ap1.mat, tmp.mat)
}
new_sand_rat_ap1.mat <- new_sand_rat_ap1.mat %>% mutate(log2_fc_steps = ifelse(fc_steps == 0, 0, log2(fc_steps)))
new_sand_rat_ap1.mat <- new_sand_rat_ap1.mat %>% 
            mutate(species2 = ifelse(species == 'Psammomys_obesus', 'Sand_rat', str_replace(species, "fullTree","")))
species.nodes2 <- c(species.nodes, 'Mus_musculus')
mouse_aff.mat <- ap1.mat %>% filter(species == 'Mus_musculus')
mouse_aff.mat$species2 <- 'Mouse'
mouse_aff.mat$log2_fc_steps <- 0
new_sand_rat_ap1.mat <- bind_rows(new_sand_rat_ap1.mat, mouse_aff.mat)
new_sand_rat_ap1.mat$species2 <- factor(new_sand_rat_ap1.mat$species2, levels = c(
    ifelse(species.nodes2 == 'Psammomys_obesus', 'Sand_rat', 
           ifelse(species.nodes2 == 'Mus_musculus','Mouse' ,str_replace(species.nodes2, "fullTree","")))
    ))
new_sand_rat_ap1.mat <- new_sand_rat_ap1.mat %>% filter(!is.na(species2)) %>% filter(species%in%sand_rat.nodes2)

### activity across species
CRE_oi <- "Epas1_chr17_10063"
big.tree <- as_tibble(tree)
big.tree <- groupClade(big.tree, c(293, 247,353,383,429,346,459))
anc.tips <- big.tree %>% filter(node%in%c(293, 247,353,383,429,346,459)) %>% pull(label)
anc.labels <- data.frame(species = anc.tips, anc_name = c('PRIMATES','RODENTIA','PERISSODACTYLA','CHIROPTERA',
                                                         'CETARTIODACTYLA','CARNIVORA','EULIPOTYPHLA'), group = c(2,1,3,4,5,6,7))

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
          species == anc_oi  ~ "Ancestor",
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
          species == anc_oi  ~ "Ancestor",
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
          species == anc_oi  ~ "Ancestor",
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

# epas1_avg_lineage_mpra$lineage <- factor(epas1_avg_lineage_mpra$lineage, levels = c('Mus_musculus','Lycaon_pictus','Murina_feae'))

# mouse_mpra.plot <- ggplot() + 
#                     geom_rect(data = epas1_avg_lineage_mpra %>% filter(species == 'Mus_musculus'), aes(ymin = mean_MPRA_act-sd_MPRA_act, ymax = mean_MPRA_act+sd_MPRA_act, color = NULL), xmin = -Inf, xmax = Inf,
#                    alpha = .2, fill = "#008000") + 
#                     geom_hline(data = epas1_avg_lineage_mpra %>% filter(species == 'Mus_musculus') , aes(yintercept = mean_MPRA_act), color = '#008000', linetype = 'dashed') +
#                     geom_rect(data = control_mean_act %>% filter(CRE_id == 'minP'), aes(ymin = mean_MPRA_act-sd_MPRA_act, ymax = mean_MPRA_act+sd_MPRA_act, color = NULL), xmin = -Inf, xmax = Inf,
#                    alpha = .2,fill = "#0057e7") +
#                     geom_hline(data = control_mean_act %>% filter(CRE_id == 'minP'), aes(yintercept = mean_MPRA_act), color = '#0057e7', linetype = 'dashed') +
#                     geom_line(data = epas1_avg_lineage_mpra %>% filter(lineage == 'Mus_musculus'), aes(x = species2, y = mean_MPRA_act, group = group)) + 
#                     geom_errorbar(data = epas1_avg_lineage_mpra %>% filter(lineage == 'Mus_musculus'), aes(x = species2, ymin=mean_MPRA_act-sd_MPRA_act, 
#                                                            ymax=mean_MPRA_act+sd_MPRA_act), 
#                                   width=.2, position=position_dodge(0.05)) + 
#                     geom_star(data = epas1_avg_lineage_mpra %>% filter(lineage == 'Mus_musculus'), aes(x = species2, y = mean_MPRA_act, starshape = marker, size = marker, fill = marker)) +  
#                     scale_starshape_manual(values = c('Other' = 15, "Ancestor" = 9, 'Extant' = 1), 
#                                            guide = guide_legend(override.aes = list(size = 2, fill = 'black')), name = "") +
#                     scale_size_manual(values = c('Other' = 1, "Ancestor" = 3, 'Extant' = 3), guide = 'none') +
#                     scale_fill_manual(values = c('Other' = 'black', "Ancestor" = 'red', 'Extant' = 'red'), guide = 'none') +
# #                     scale_y_log10(
# #                        breaks = c(10^-1, 10^0, 10^1),
# #                        labels = scales::trans_format("log10", scales::math_format(10^.x))
# #                      ) + annotation_logticks(sides = 'l') +
#                     labs(y = "MPRA activity", x = "") + coord_cartesian(clip="off") + ylim(c(0, 9)) + 
#                     theme_classic() + theme(legend.position = 'none', legend.box = "horizontal", 
#                                             axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# dog_mpra.plot <- ggplot() + 
#                     geom_rect(data = epas1_avg_lineage_mpra %>% filter(species == 'Mus_musculus'), aes(ymin = mean_MPRA_act-sd_MPRA_act, ymax = mean_MPRA_act+sd_MPRA_act, color = NULL), xmin = -Inf, xmax = Inf,
#                    alpha = .2, fill = "#008000") + 
#                     geom_hline(data = epas1_avg_lineage_mpra %>% filter(species == 'Mus_musculus') , aes(yintercept = mean_MPRA_act), color = '#008000', linetype = 'dashed') +
#                     geom_rect(data = control_mean_act %>% filter(CRE_id == 'minP'), aes(ymin = mean_MPRA_act-sd_MPRA_act, ymax = mean_MPRA_act+sd_MPRA_act, color = NULL), xmin = -Inf, xmax = Inf,
#                    alpha = .2,fill = "#0057e7") +
#                     geom_hline(data = control_mean_act %>% filter(CRE_id == 'minP'), aes(yintercept = mean_MPRA_act), color = '#0057e7', linetype = 'dashed') +
#                     geom_line(data = epas1_avg_lineage_mpra %>% filter(lineage == 'Lycaon_pictus'), aes(x = species2, y = mean_MPRA_act, group = group)) + 
#                     geom_errorbar(data = epas1_avg_lineage_mpra %>% filter(lineage == 'Lycaon_pictus'), aes(x = species2, ymin=mean_MPRA_act-sd_MPRA_act, 
#                                                            ymax=mean_MPRA_act+sd_MPRA_act), 
#                                   width=.2, position=position_dodge(0.05)) + 
#                     geom_star(data = epas1_avg_lineage_mpra %>% filter(lineage == 'Lycaon_pictus'), aes(x = species2, y = mean_MPRA_act, starshape = marker, size = marker, fill = marker)) +  
#                     scale_starshape_manual(values = c('Other' = 15, "Ancestor" = 9, 'Extant' = 1), 
#                                            guide = guide_legend(override.aes = list(size = 2, fill = 'black')), name = "") +
#                     scale_size_manual(values = c('Other' = 1, "Ancestor" = 3, 'Extant' = 3), guide = 'none') +
#                     scale_fill_manual(values = c('Other' = 'black', "Ancestor" = 'red', 'Extant' = 'red'), guide = 'none') +
# #                     scale_y_log10(
# #                        breaks = c(10^-1, 10^0, 10^1),
# #                        labels = scales::trans_format("log10", scales::math_format(10^.x))
# #                      ) + annotation_logticks(sides = 'l') +
#                     labs(y = "MPRA activity", x = "") + coord_cartesian(clip="off") + ylim(c(0, 9)) + 
#                     theme_classic() + theme(legend.position = 'none', legend.box = "horizontal", 
#                                             axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# bat_mpra.plot <- ggplot() + 
#                     geom_rect(data = epas1_avg_lineage_mpra %>% filter(species == 'Mus_musculus'), aes(ymin = mean_MPRA_act-sd_MPRA_act, ymax = mean_MPRA_act+sd_MPRA_act, color = NULL), xmin = -Inf, xmax = Inf,
#                    alpha = .2, fill = "#008000") + 
#                     geom_hline(data = epas1_avg_lineage_mpra %>% filter(species == 'Mus_musculus') , aes(yintercept = mean_MPRA_act), color = '#008000', linetype = 'dashed') +
#                     geom_rect(data = control_mean_act %>% filter(CRE_id == 'minP'), aes(ymin = mean_MPRA_act-sd_MPRA_act, ymax = mean_MPRA_act+sd_MPRA_act, color = NULL), xmin = -Inf, xmax = Inf,
#                    alpha = .2,fill = "#0057e7") +
#                     geom_hline(data = control_mean_act %>% filter(CRE_id == 'minP'), aes(yintercept = mean_MPRA_act), color = '#0057e7', linetype = 'dashed') +
#                     geom_line(data = epas1_avg_lineage_mpra %>% filter(lineage == 'Murina_feae'), aes(x = species2, y = mean_MPRA_act, group = group)) + 
#                     geom_errorbar(data = epas1_avg_lineage_mpra %>% filter(lineage == 'Murina_feae'), aes(x = species2, ymin=mean_MPRA_act-sd_MPRA_act, 
#                                                            ymax=mean_MPRA_act+sd_MPRA_act), 
#                                   width=.2, position=position_dodge(0.05)) + 
#                     geom_star(data = epas1_avg_lineage_mpra %>% filter(lineage == 'Murina_feae'), aes(x = species2, y = mean_MPRA_act, starshape = marker, size = marker, fill = marker)) +  
#                     scale_starshape_manual(values = c('Other' = 15, "Ancestor" = 9, 'Extant' = 1), 
#                                            guide = guide_legend(override.aes = list(size = 2, fill = 'black')), name = "") +
#                     scale_size_manual(values = c('Other' = 1, "Ancestor" = 3, 'Extant' = 3), guide = 'none') +
#                     scale_fill_manual(values = c('Other' = 'black', "Ancestor" = 'red', 'Extant' = 'red'), guide = 'none') +
# #                     scale_y_log10(
# #                        breaks = c(10^-1, 10^0, 10^1),
# #                        labels = scales::trans_format("log10", scales::math_format(10^.x))
# #                      ) + annotation_logticks(sides = 'l') +
#                     labs(y = "MPRA activity", x = "") + coord_cartesian(clip="off") + ylim(c(0, 9)) + 
#                     theme_classic() + theme(legend.position = 'none', legend.box = "horizontal", 
#                                             axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# layout <- 'AAAAABBBB
#            AAAAABBBB
#            AAAAABBBB
#            AAAAABBBB
#            CCCDDDEEE
#            CCCDDDEEE
#            FFFF#####
#            FFFF#####
#            FFFF#####'
# p <- epas.tree + aff.plot + mouse_mpra.plot  + bat_mpra.plot + dog_mpra.plot + 
#      no_aff.plot + plot_layout(design = layout)
# ggsave('test.pdf', plot = p, width = 210, height = 220, useDingbats = F, units = 'mm', device = 'pdf')

# labels <- levels(epas1_msa_list$msa_df_TFBS$species)
# epas1_gg_tfbs <- plot_tfbs_trajectory(epas1_msa_list$msa_df_TFBS, epas1_msa_list$msa_func_TFBS, CRE_oi, 'Mus_musculus',labels, 
#                                       tfbs_size = 0.1)


# # species.pair <- c("fullTreeAnc68",'fullTreeAnc239')
# # anc.cutoff <- make.nodepath(tree, species.pair)

# ### plot Gata4 MSA for five ancestors
# mouse_DMS_MPRA$mut <- toupper(mouse_DMS_MPRA$mut)
big.tree <- as_tibble(tree)
big.tree <- groupClade(big.tree, c(293, 247,353,383,429,346,459))
anc.tips <- big.tree %>% filter(node%in%c(293, 247,353,383,429,346,459)) %>% pull(label)
anc.labels <- data.frame(species = anc.tips, anc_name = c('PRIMATES','RODENTIA','PERISSODACTYLA','CHIROPTERA',
                                                         'CETARTIODACTYLA','CARNIVORA','EULIPOTYPHLA'), group = c(2,1,3,4,5,6,7)) 

# epas1_msa_df <- make_alignment_dataframe(fasta_file, c(anc.tips,'Mus_musculus'), 'Mus_musculus',CRE_oi)
# anc.tips <- anc.tips[anc.tips%in%unique(epas1_msa_df$species)]
# mouse_msa_df <- epas1_msa_df %>% filter(species == 'Mus_musculus')

# epas1_msa_list <- make_msa_tfbs_traj_df(fasta_file, c(anc.tips,'Mus_musculus'), 'Mus_musculus', CRE_oi, cactus.tfbs, mouse_func_tfbs)
# labels <- collect_labels(epas1_msa_list$msa_df_TFBS, side = T)
# epas1_func_tfbs <- mouse_func_tfbs %>% filter(CRE == CRE_oi)
# mm.start <- match(epas1_func_tfbs$TFBS_start, mouse_msa_df$ref_pos)
# mm.end <- match(epas1_func_tfbs$TFBS_end, mouse_msa_df$ref_pos)
# epas1_func_tfbs$msa_start <- mouse_msa_df$curpos[mm.start]
# epas1_func_tfbs$msa_end <- mouse_msa_df$curpos[mm.end]
# epas1_msa_df <- left_join(epas1_msa_df, mouse_DMS_MPRA %>% filter(CRE == CRE_oi) %>%
#                           mutate(log2FC_MPRA = ifelse(is.infinite(log2FC_MPRA) & log2FC_MPRA < 0, NA, log2FC_MPRA)) %>% ### replace -Inf with NA
#                           dplyr::select(ref_pos = pos, mut, log2FC_MPRA, sd_log2FC_MPRA, MPRA_act, sd_MPRA_act, min_act,max_act, wilcox_BH),
#                          by = c("ref_pos","mut"))

# epas1_msa_df <- epas1_msa_df %>% mutate(mut2 = case_when(
#                         mut == 'DEL' ~ '-',
#                         TRUE ~ mut))

# anc_tfbs_map <- data.frame(msa_start = c(162,243),
#                       msa_end = c(204, 283))
# epas1_msa_list$msa_df_TFBS <- epas1_msa_list$msa_df_TFBS %>% 
#                     mutate(species = factor(species, levels = c('Mus_musculus',"fullTreeAnc68","fullTreeAnc110",
#                                                                'fullTreeAnc121','fullTreeAnc150','fullTreeAnc194',
#                                                                'fullTreeAnc223','fullTreeAnc233')))
# ### now align every species
# #### generate affinity matrix
# epas1_oCRE_mean_act <- oCRE_mean_act %>%
#   filter(full_CRE_id == CRE_oi)
# all.species <- epas1_oCRE_mean_act %>% pull(species)
# all_epas1_msa_df <- make_alignment_dataframe(fasta_file, all.species, 'Mus_musculus',CRE_oi)
# mouse_msa_df <- all_epas1_msa_df %>% filter(species == 'Mus_musculus')

# all_epas1_msa_list <- make_msa_tfbs_traj_df(fasta_file, all.species, 'Mus_musculus', CRE_oi, cactus.tfbs, mouse_func_tfbs)
# labels <- collect_labels(all_epas1_msa_list$msa_df_TFBS, side = T)
# all_epas1_func_tfbs <- mouse_func_tfbs %>% filter(CRE == CRE_oi)
# mm.start <- match(all_epas1_func_tfbs$TFBS_start, mouse_msa_df$ref_pos)
# mm.end <- match(all_epas1_func_tfbs$TFBS_end, mouse_msa_df$ref_pos)
# all_epas1_func_tfbs$msa_start <- mouse_msa_df$curpos[mm.start]
# all_epas1_func_tfbs$msa_end <- mouse_msa_df$curpos[mm.end]
# all_epas1_msa_df <- left_join(all_epas1_msa_df, mouse_DMS_MPRA %>% filter(CRE == CRE_oi) %>%
#                           mutate(log2FC_MPRA = ifelse(is.infinite(log2FC_MPRA) & log2FC_MPRA < 0, NA, log2FC_MPRA)) %>% ### replace -Inf with NA
#                           dplyr::select(ref_pos = pos, mut, log2FC_MPRA, sd_log2FC_MPRA, MPRA_act, sd_MPRA_act, min_act,max_act, wilcox_BH),
#                          by = c("ref_pos","mut"))

# all_epas1_msa_df <- all_epas1_msa_df %>% mutate(mut2 = case_when(
#                         mut == 'DEL' ~ '-',
#                         TRUE ~ mut))

# all_epas1_msa_list$msa_df_TFBS <- all_epas1_msa_list$msa_df_TFBS %>% 
#                         left_join(big.tree %>% dplyr::select(species = label, group) %>% mutate(group = factor(group))) %>% 
#                         left_join(anc.labels %>% dplyr::select(group, Order = anc_name) %>% mutate(group = factor(group)))


# aff.mat <- all_epas1_msa_list$msa_df_TFBS %>% filter(!is.na(TF_name2)) %>% 
#             mutate(Order = factor(Order, levels = c('RODENTIA','PRIMATES','PERISSODACTYLA','CHIROPTERA','CETARTIODACTYLA','CARNIVORA','EULIPOTYPHLA'))) %>% 
#             mutate(genome = case_when(
#                         startsWith(species, 'full') ~ 'ancestor',
#                         TRUE ~ 'extant'))
# triple.mat <- aff.mat %>% filter(TF_name2%in%c('Foxa2_3','Sox17_1','Gata4/6_4')) %>% 
#                 mutate(TF_name = factor(TF_name, levels = c('Foxa2','Sox17','Gata4/6'))) %>% 
#                 mutate(log2_affFC = norm_affinity/norm_affinity_ref) %>% filter(!is.na(Order))
# library(ggbeeswarm)
# triplet.plot <- ggplot() + geom_quasirandom(data = triple.mat, aes(x = Order, y = log2_affFC, color = Order, 
#                                                                           shape = genome, alpha = genome), 
#             dodge.width=1, size = 1, stroke = 1.1) +  # Adjust the width for separation 
#             geom_hline(data = triple.mat %>% filter(species == 'Mus_musculus'), 
#                        aes(yintercept =  log2_affFC),color = "#008000", linetype = 'dashed') + 
# #             scale_y_log10(
# #                labels = scales::trans_format("log10", scales::math_format(10^.x))
# #              ) + #annotation_logticks(sides="l") + 
#             scale_color_manual(values = order.colors, guide = 'none') + 
# #                 scale_fill_manual(values = c("ancestor" = "white", "extant" = "black")) +
#             scale_shape_manual(values = c(23, 21), guide = guide_legend(override.aes = list(size = 1, color = 'black'))) +
#             scale_alpha_manual(values = c(0.3, 0.8)) + theme_classic() + scale_x_discrete(drop = FALSE)  +
#                             labs(y = 'log2FC predicted TFBS affinity (species/mouse)', x = '') + 
#                 facet_wrap(~TF_name, ncol = 1, scales = "free") + 
#                 theme(legend.position = 'bottom', 
#                       axis.text.x = element_text(angle = 45, hjust=1, color = "black"),
#                       axis.text.y = element_text(color = "black"))

# ap1_triplet.mat <- aff.mat %>% filter(TF_name2%in%c('Jun_Atf3_8','Jun_Atf3_2','Jun_Atf3_5')) %>% 
#                 mutate(TF_name2 = factor(TF_name2, levels = c('Jun_Atf3_8','Jun_Atf3_2','Jun_Atf3_5'))) %>% 
#                 mutate(log2_affFC = norm_affinity/norm_affinity_ref) %>% filter(!is.na(Order))

# ap1.plot <- ggplot() + geom_quasirandom(data = ap1_triplet.mat, aes(x = Order, y = log2_affFC, color = Order, 
#                                                                           shape = genome, alpha = genome), 
#             dodge.width=1, size = 1, stroke = 1.1) +  # Adjust the width for separation 
#             geom_hline(data = ap1_triplet.mat %>% filter(species == 'Mus_musculus'), 
#                        aes(yintercept =  log2_affFC),color = "#008000", linetype = 'dashed') + 
# #             scale_y_log10(
# #                labels = scales::trans_format("log10", scales::math_format(10^.x))
# #              ) + #annotation_logticks(sides="l") + 
#             scale_color_manual(values = order.colors, guide = 'none') + 
# #                 scale_fill_manual(values = c("ancestor" = "white", "extant" = "black")) +
#             scale_shape_manual(values = c(23, 21), guide = guide_legend(override.aes = list(size = 1, color = 'black'))) +
#             scale_alpha_manual(values = c(0.3, 0.8)) + theme_classic() + scale_x_discrete(drop = FALSE)  +
#                             labs(y = 'log2FC predicted TFBS affinity (species/mouse)', x = '') + 
#                 facet_wrap(~TF_name2, ncol = 1, scales = "free") + 
#                 theme(legend.position = 'bottom', 
#                       axis.text.x = element_text(angle = 45, hjust=1, color = "black"),
#                       axis.text.y = element_text(color = "black"))






# CRE_oi <- "Epas1_chr17_10063"
# ### plot Epas1 transitions
# species.pair <- c("Mus_musculus",'fullTreeAnc239')
# mouse.nodes <- make.nodepath(tree, species.pair)
# ### filter nodes if in activity
# mouse.nodes <- mouse.nodes[mouse.nodes%in%unique(oCRE_mean_act$species)]

# start_node <- "fullTreeAnc239"
# order_species <- rev(mouse.nodes)  # traverse in your chosen order

# ### get replicate
# oCRE_replicate_act <- df_mpra_act %>% filter(class == 'oCRE')
# oCRE_replicate_act <- oCRE_replicate_act %>% separate(CRE_id, c('tile_id','something'), sep = "_20240416") %>% dplyr::select(-something)
# oCRE_replicate_act <- oCRE_replicate_act %>% separate(tile_id, c('full_CRE_id','species'), sep = "__")
# mm <- match(oCRE_replicate_act$species, cactus.annotation$Species)
# oCRE_replicate_act$common_name <- cactus.annotation$`Common Name`[mm]
# oCRE_replicate_act <- left_join(oCRE_replicate_act, max.dist, by =c("full_CRE_id","species"))

# # base table (ordered once)
# epas1_avg_lineage_mpra <- oCRE_mean_act %>%
#   filter(full_CRE_id == CRE_oi, species %in% mouse.nodes) %>%
#   mutate(
#     species2  = ifelse(species == 'Mus_musculus', 'Mouse', str_replace(species, "fullTree","")),
#     marker = case_when(
#       species2 %in% c("Anc68", "Anc59") ~ "Transitions",
#       species2 == "Mouse"                 ~ "Mouse",
#       TRUE                                       ~ "Other"
#     ),
#     species2 = factor(species2, levels = c(
#     ifelse(order_species == "Mus_musculus", "Mouse", str_replace(order_species, "fullTree",""))
#   ))
#   ) %>%
#   arrange(species)

# epas1_lineage_mpra <- oCRE_replicate_act %>%
#   filter(full_CRE_id == CRE_oi, species %in% mouse.nodes) %>%
#   mutate(
#     species2  = ifelse(species == 'Mus_musculus', 'Mouse', str_replace(species, "fullTree","")),
#     marker = case_when(
#       species2 %in% c("Anc68", "Anc59") ~ "Transitions",
#       species2 == "Mouse"                 ~ "Mouse",
#       TRUE                                       ~ "Other"
#     ),
#     species = factor(species, levels = order_species)
#   ) %>%
#   arrange(species)

# # constants
# start_act_df <- epas1_lineage_mpra %>%
#   filter(species == start_node) %>%
#   dplyr::select(biol_rep, start_act = MPRA_act)

# # Get mouse wild-type activity per replicate
# wt_act_df <- epas1_lineage_mpra %>%
#   filter(species == 'Mus_musculus') %>%
#   dplyr::select(biol_rep, wt_act = MPRA_act)

# # stepwise deltas & fold-changes
# epas1_lineage_mpra_steps <- epas1_lineage_mpra %>%
#   arrange(biol_rep, match(species, order_species)) %>%
#   group_by(biol_rep) %>%
#   mutate(
#     prev_act     = lag(MPRA_act, default = NA),
#     delta_steps  = MPRA_act - prev_act,
#     fc_steps     = MPRA_act / prev_act,
#     log2FC_steps = log2(fc_steps)
#   ) %>%
#   ungroup() %>%
#   left_join(start_act_df, by = "biol_rep") %>%
#   left_join(wt_act_df, by = "biol_rep") %>%
#   mutate(
#     delta_steps  = if_else(species == start_node, 0, delta_steps),
#     fc_steps     = if_else(species == start_node, 1, fc_steps),
#     log2FC_steps = if_else(species == start_node, 0, log2FC_steps),
#     frac_delta   = delta_steps / (wt_act - start_act),
#     species2     = factor(species2, levels = str_replace(order_species, "^fullTree", ""))
#   )

# stepwise_summary <- epas1_lineage_mpra_steps %>%
#   group_by(species) %>%
#   summarise(
#     avg_delta_steps     = mean(delta_steps, na.rm = TRUE),
#     sd_delta_steps      = sd(delta_steps, na.rm = TRUE),
#     avg_frac_delta      = mean(frac_delta, na.rm = TRUE),
#     sd_frac_delta       = sd(frac_delta, na.rm = TRUE),
#     avg_fc_steps        = mean(fc_steps, na.rm = TRUE),
#     sd_fc_steps         = sd(fc_steps, na.rm = TRUE),
#     avg_log2FC_steps    = mean(log2FC_steps, na.rm = TRUE),
#     sd_log2FC_steps     = sd(log2FC_steps, na.rm = TRUE),
#     p_delta_steps       = ifelse(n() > 1, t.test(delta_steps)$p.value, NA_real_),
#     p_frac_delta        = ifelse(n() > 1, t.test(frac_delta)$p.value, NA_real_),
#     p_log2FC_steps      = ifelse(n() > 1, t.test(log2FC_steps)$p.value, NA_real_)
#   ) %>%
#   ungroup()
# stepwise_summary$p_adj_log2FC_steps <- p.adjust(stepwise_summary$p_log2FC_steps, method = 'fdr')
# stepwise_summary <- stepwise_summary %>% 
#     mutate(
#         species2  = ifelse(species == 'Mus_musculus', 'Mouse', str_replace(species, "fullTree","")),
#         marker = case_when(
#           species2 %in% c("Anc68", "Anc59") ~ "Transitions",
#           species2 == "Mouse"                 ~ "Mouse",
#           TRUE                                       ~ "Other"
#         ),
#         species2 = factor(species2, levels = c(
#     ifelse(order_species == "Mus_musculus", "Mouse", str_replace(order_species, "fullTree",""))
#   ))
#     )
# epas1_avg_lineage_mpra$group <- 1
# stepwise_summary$group <- 1


# CRE_oi <- "Epas1_chr17_10063"
# fasta_file <- '/net/shendure/vol10/projects/tli/nobackup/EB_dev_project/oCREs_v2/nobackup/reoriented_fasta/Epas1_chr17_10063_oCRE.fasta'
# # species.pair <- c("fullTreeAnc68",'fullTreeAnc239')
# # anc.cutoff <- make.nodepath(tree, species.pair)

# ### plot Gata4 MSA for five ancestors
# mouse_DMS_MPRA$mut <- toupper(mouse_DMS_MPRA$mut)
# big.tree <- as_tibble(tree)
# big.tree <- groupClade(big.tree, c(293, 247,353,383,429,346,459))
# anc.tips <- big.tree %>% filter(node%in%c(293, 247,353,383,429,346,459)) %>% pull(label)
# anc.labels <- data.frame(species = anc.tips, anc_name = c('PRIMATES','RODENTIA','PERISSODACTYLA','CHIROPTERA',
#                                                          'CETARTIODACTYLA','CARNIVORA','EULIPOTYPHLA'))

# epas1_msa_df <- make_alignment_dataframe(fasta_file, c(anc.tips,'Mus_musculus'), 'Mus_musculus',CRE_oi)
# anc.tips <- anc.tips[anc.tips%in%unique(epas1_msa_df$species)]
# mouse_msa_df <- epas1_msa_df %>% filter(species == 'Mus_musculus')

# epas1_msa_list <- make_msa_tfbs_traj_df(fasta_file, c(anc.tips,'Mus_musculus'), 'Mus_musculus', CRE_oi, cactus.tfbs, mouse_func_tfbs)
# labels <- collect_labels(epas1_msa_list$msa_df_TFBS, side = T)
# epas1_func_tfbs <- mouse_func_tfbs %>% filter(CRE == CRE_oi)
# mm.start <- match(epas1_func_tfbs$TFBS_start, mouse_msa_df$ref_pos)
# mm.end <- match(epas1_func_tfbs$TFBS_end, mouse_msa_df$ref_pos)
# epas1_func_tfbs$msa_start <- mouse_msa_df$curpos[mm.start]
# epas1_func_tfbs$msa_end <- mouse_msa_df$curpos[mm.end]
# epas1_msa_df <- left_join(epas1_msa_df, mouse_DMS_MPRA %>% filter(CRE == CRE_oi) %>%
#                           mutate(log2FC_MPRA = ifelse(is.infinite(log2FC_MPRA) & log2FC_MPRA < 0, NA, log2FC_MPRA)) %>% ### replace -Inf with NA
#                           dplyr::select(ref_pos = pos, mut, log2FC_MPRA, sd_log2FC_MPRA, MPRA_act, sd_MPRA_act, min_act,max_act, wilcox_BH),
#                          by = c("ref_pos","mut"))

# epas1_msa_df <- epas1_msa_df %>% mutate(mut2 = case_when(
#                         mut == 'DEL' ~ '-',
#                         TRUE ~ mut))
# epas1_msa_df <- epas1_msa_df %>% mutate(species2 = ifelse(species == 'Mus_musculus', 'Mouse', str_replace(species, "fullTree","")))
# epas1_msa_df$species2 <- factor(epas1_msa_df$species2, levels = c(
#     ifelse(keep.tips == "Mus_musculus", "Mouse", str_replace(keep.tips, "fullTree",""))
#   ))

# epas1_msa_df <- epas1_msa_df %>%
#   mutate(fill_type = case_when(
#     mut2 == "-" & is.na(log2FC_MPRA) ~ "black",
#     mut == ref_nuc ~ "grey",
#     TRUE ~ "gradient"
#   ))

# fixed_epas1_msa_df <- data.frame()

# for(pos in seq(1:max(epas1_msa_df$curpos))){
#     tmp_msa_df <- epas1_msa_df %>% filter(curpos == pos)
#     all_mut_types <- unique(tmp_msa_df$mut_type)
#     if(any(all_mut_types %in% c('Mismatch', 'DEL', 'INS'))){
#         current_log2fc <- 0
#         have_not_found_match <- TRUE
#         log_change <- c()
#         match2 <- c()
#         for(species_oi in rev(keep.tips)){
#             tmp_species <- tmp_msa_df %>% filter(species == species_oi)
#             if(have_not_found_match){
#                 check_match <- tmp_species %>% filter(mut2 == ref_nuc)
#                 current_muttype <- tmp_species %>% pull(mut_type)
#                 if(nrow(check_match) > 0){
#                     log_change <- c(log_change, -1*current_log2fc)
#                     have_not_found_match <- FALSE
#                     match2 <- c(match2, 'First Match')
#                 }else{
#                     current_log2fc <- tmp_species %>% pull(log2FC_MPRA)
#                     log_change <- c(log_change, 0)
#                     match2 <- c(match2, current_muttype)
#                 }
# #                 current_muttype <- tmp_species %>% pull(mut_type)
# #                 if(current_muttype != 'Match'){
# #                     current_log2fc <- tmp_species %>% pull(log2FC_MPRA)
# #                     log_change <- c(log_change, 0)
# #                     match2 <- c(match2, current_muttype)
# #                 }else{
# #                     log_change <- c(log_change, -1*current_log2fc)
# #                     have_not_found_match <- FALSE
# #                     match2 <- c(match2, 'First Match')
# #                 }
#             }else{
#                 log_change <- c(log_change, 0)
#                 match2 <- c(match2, 'Match')
#             }
#         }
#         names(log_change) <- rev(keep.tips)
#         names(match2) <- rev(keep.tips)
#         mm <- match(tmp_msa_df$species, names(log_change))
#         tmp_msa_df$log_change <- log_change[mm]
#         tmp_msa_df$match2 <- match2[mm]
#         fixed_epas1_msa_df <- bind_rows(fixed_epas1_msa_df, tmp_msa_df)
#     }else{
#         tmp_msa_df$log_change <- 0
#         tmp_msa_df$match2 <- 'Match'
#         fixed_epas1_msa_df <- bind_rows(fixed_epas1_msa_df, tmp_msa_df)
#     }
# }

# fixed_epas1_msa_df <- fixed_epas1_msa_df %>% mutate(log_change = case_when(
#                                         mut_type%in%c('DEL','INS') ~ NA,
#                                         TRUE ~ log_change))
# fixed_epas1_msa_df <- fixed_epas1_msa_df %>% 
#             mutate(species2 = ifelse(species == 'Mus_musculus', 'Mouse', str_replace(species, "fullTree","")))
# fixed_epas1_msa_df$species2 <- factor(fixed_epas1_msa_df$species2, levels = c(
#     ifelse(keep.tips == "Mus_musculus", "Mouse", str_replace(keep.tips, "fullTree",""))
#   ))


# #### calculate number of SNPs that become the original mouse sequence at each step and mark which original mouse SNP confers a TFBS
# ### calculate avg log2FC for each mutation gained
# all_mouse_tfbs_positions <- unlist(mapply(seq, epas1_func_tfbs$msa_start, epas1_func_tfbs$msa_end))
# epas1_current_snps <- epas1_msa_df %>% filter(species == 'fullTreeAnc239') %>% filter(mut2 != ref_nuc) %>% pull(curpos)
# gain_mouse_snps <- fixed_epas1_msa_df %>% filter(species == 'fullTreeAnc239') %>% filter(mut2 == ref_nuc) %>%
#                             arrange(curpos) %>% mutate(group_id = cumsum(c(TRUE, diff(curpos) > 1)),
#                                                       in_TFBS = case_when(
#                                                           curpos%in%all_mouse_tfbs_positions ~ TRUE,
#                                                           TRUE ~ FALSE))
# fc_calc <- gain_mouse_snps %>% group_by(in_TFBS) %>% 
#                         summarise(avg_log_change = mean(log_change, na.rm = TRUE), sd_log_change = sd(log_change, na.rm = TRUE))
# epas1_all_fc <- gain_mouse_snps

# epas1_snp_counts <- gain_mouse_snps %>%
#                   group_by(species, group_id) %>%
#                   summarise(
#                     num_snps = 1,                     # each contiguous run counts as 1
#                     tfbs_count = any(in_TFBS),        # TRUE if any in the run is in TFBS
#                     .groups = "drop_last"
#                   ) %>%
#                   summarise(
#                     n_snp_count = sum(num_snps),
#                     total_tfbs_snps = sum(tfbs_count),# same as above unless you want total TRUEs
#                     .groups = "drop"
#                   ) %>% mutate(total_nontfbs_snps = n_snp_count - total_tfbs_snps)
# epas1_avg_fc <- data.frame(
#                                 species         = 'fullTreeAnc239',
#                                 avg_tfbs_FC     = pull_or_zero(fc_calc %>% filter(in_TFBS), "avg_log_change"),
#                                 avg_nontfbs_FC  = pull_or_zero(fc_calc %>% filter(!in_TFBS), "avg_log_change"),
#                                 sd_tfbs_FC      = pull_or_zero(fc_calc %>% filter(in_TFBS), "sd_log_change"),
#                                 sd_nontfbs_FC   = pull_or_zero(fc_calc %>% filter(!in_TFBS), "sd_log_change")
#                               )
# for(species_oi in rev(keep.tips)){
#     if(species_oi == 'fullTreeAnc239'){
#     } else{
#         temp_msa_df <- fixed_epas1_msa_df %>% filter(species == species_oi)
#         current_snps <- temp_msa_df %>% filter(mut2 != ref_nuc) %>% pull(curpos)
#         if(length(epas1_current_snps) > 0){
#             gain_mouse_snps <- temp_msa_df %>% filter(curpos%in%epas1_current_snps) %>% filter(mut2 == ref_nuc) %>%
#                             arrange(curpos) %>% mutate(group_id = cumsum(c(TRUE, diff(curpos) > 1)),
#                                                       in_TFBS = case_when(
#                                                           curpos%in%all_mouse_tfbs_positions ~ TRUE,
#                                                           TRUE ~ FALSE))
#             if (nrow(gain_mouse_snps) > 0) {
#               fc_calc <- gain_mouse_snps %>% group_by(in_TFBS) %>% 
#                         summarise(avg_log_change = mean(log_change, na.rm = TRUE), sd_log_change = sd(log_change, na.rm = TRUE))
#               epas1_all_fc <- bind_rows(epas1_all_fc, gain_mouse_snps)
#               snp_counts <- gain_mouse_snps %>%
#                   group_by(species, group_id) %>%
#                   summarise(
#                     num_snps = 1,                     # each contiguous run counts as 1
#                     tfbs_count = any(in_TFBS),        # TRUE if any in the run is in TFBS
#                     .groups = "drop_last"
#                   ) %>%
#                   summarise(
#                     num_snps = sum(num_snps),
#                     total_tfbs_snps = sum(tfbs_count),# same as above unless you want total TRUEs
#                     .groups = "drop"
#                   ) %>% mutate(total_nontfbs_snps = num_snps - total_tfbs_snps)
#                 epas1_snp_counts <- bind_rows(epas1_snp_counts,
#                                      data.frame(species = c(species_oi), 
#                                     n_snp_count = c(snp_counts %>% pull(num_snps)),
#                                     total_tfbs_snps = c(snp_counts %>% pull(total_tfbs_snps)),
#                                      total_nontfbs_snps = c(snp_counts %>% pull(total_nontfbs_snps))))
#                 epas1_avg_fc <- bind_rows(epas1_avg_fc,
#                               data.frame(
#                                 species         = species_oi,
#                                 avg_tfbs_FC     = pull_or_zero(fc_calc %>% filter(in_TFBS), "avg_log_change"),
#                                 avg_nontfbs_FC  = pull_or_zero(fc_calc %>% filter(!in_TFBS), "avg_log_change"),
#                                 sd_tfbs_FC      = pull_or_zero(fc_calc %>% filter(in_TFBS), "sd_log_change"),
#                                 sd_nontfbs_FC   = pull_or_zero(fc_calc %>% filter(!in_TFBS), "sd_log_change")
#                               )
#                             )
#             }else{
#                 epas1_snp_counts <- bind_rows(epas1_snp_counts,
#                                      data.frame(species = c(species_oi), 
#                                     n_snp_count = 0,
#                                     total_tfbs_snps = 0,
#                                      total_nontfbs_snps = 0))
#                 epas1_avg_fc <- bind_rows(epas1_avg_fc, 
#                                     data.frame(species = species_oi, 
#                                      avg_tfbs_FC = 0, 
#                                      avg_nontfbs_FC = 0, 
#                                      sd_tfbs_FC = 0,
#                                      sd_nontfbs_FC = 0))
#                 epas1_all_fc <- bind_rows(epas1_all_fc, data.frame(species = species_oi))
#             }
#         }else{
#             epas1_snp_counts <- bind_rows(epas1_snp_counts,
#                                      data.frame(species = c(species_oi), 
#                                     n_snp_count = 0,
#                                     total_tfbs_snps = 0,
#                                      total_nontfbs_snps = 0))
#             epas1_avg_fc <- bind_rows(epas1_avg_fc, 
#                                     data.frame(species = species_oi, 
#                                      avg_tfbs_FC = 0, 
#                                      avg_nontfbs_FC = 0, 
#                                      sd_tfbs_FC = 0,
#                                      sd_nontfbs_FC = 0))
#             epas1_all_fc <- bind_rows(epas1_all_fc, data.frame(species = species_oi))
#         }
#         epas1_current_snps <- current_snps
#     }
# }

# epas1_avg_fc <- epas1_avg_fc %>%
#   mutate(across(
#     .cols = where(is.numeric),
#     .fns = ~ replace_na(.x, 0))
#   )


# epas1_snp_counts <- epas1_snp_counts %>% mutate(species2 = ifelse(species == 'Mus_musculus', 'Mouse', str_replace(species, "fullTree","")))
# epas1_snp_counts$species2 <- factor(epas1_snp_counts$species2, levels = c(
#     ifelse(keep.tips == "Mus_musculus", "Mouse", str_replace(keep.tips, "fullTree",""))
#   ))
# epas1_avg_fc <- epas1_avg_fc %>% mutate(species2 = ifelse(species == 'Mus_musculus', 'Mouse', str_replace(species, "fullTree","")))
# epas1_avg_fc$species2 <- factor(epas1_avg_fc$species2, levels = c(
#     ifelse(keep.tips == "Mus_musculus", "Mouse", str_replace(keep.tips, "fullTree",""))))
# epas1_all_fc <- epas1_all_fc %>% mutate(species2 = ifelse(species == 'Mus_musculus', 'Mouse', str_replace(species, "fullTree",""))) %>% 
#                 mutate(tf_call = case_when(
#                             in_TFBS ~ 'func. TFBS',
#                             TRUE ~ 'non func. TFBS'))
# epas1_all_fc$species2 <- factor(epas1_all_fc$species2, levels = c(
#     ifelse(keep.tips == "Mus_musculus", "Mouse", str_replace(keep.tips, "fullTree",""))))

# epas1_first_match_group_fc <- epas1_all_fc %>% filter(match2 == 'First Match') %>% group_by(species2, group_id, tf_call) %>% 
#                         summarise(avg_log_change = mean(log_change, na.rm = TRUE), sd_log_change = sd(log_change, na.rm = TRUE))

# #### generate affinity matrix
# aff.mat <- epas1_msa_list$msa_df_TFBS %>% filter(!is.na(TF_name2)) %>% 
#             mutate(species2 = ifelse(species == 'Mus_musculus', 'Mouse', str_replace(species, "fullTree",""))) %>% 
#                 mutate(marker = case_when(
#                         species2%in%c("Anc68",'Anc59') ~ 'Transitions',
#                         species2 == 'Mouse' ~ 'Mouse',
#                         TRUE ~ 'Other'))
# aff.mat$species2 <- factor(aff.mat$species2, levels = c(
#     ifelse(keep.tips == "Mus_musculus", "Mouse", str_replace(keep.tips, "fullTree",""))
#   ))
# new_aff.mat <- data.frame()

# for(tf2 in unique(aff.mat$TF_name2)){
#     tmp.mat <- aff.mat %>% filter(TF_name2 == tf2)
#     fc_steps <- c(0)
#     delta_steps <- c(0)
#     initial_step <- tmp.mat %>% filter(species == 'fullTreeAnc239') %>% pull(norm_affinity)
#     for(species_oi in rev(keep.tips)){
#         if(species_oi == 'fullTreeAnc239'){
#         }else{
#             current_step <- tmp.mat %>% filter(species == species_oi) %>% pull(norm_affinity)
#             fc = current_step/initial_step
#             delta <- current_step - initial_step
#             fc_steps <- c(fc_steps, fc)
#             delta_steps <- c(delta_steps, delta)
#             initial_step <- current_step
#         }
#     }
#     names(delta_steps) <- rev(keep.tips)
#     names(fc_steps) <- rev(keep.tips)
#     mm <- match(tmp.mat$species, names(delta_steps))
#     tmp.mat$delta_steps <- delta_steps[mm]
#     tmp.mat$fc_steps <- fc_steps[mm]
#     start_aff <- tmp.mat %>% filter(species == 'fullTreeAnc239') %>% pull(norm_affinity)
#     wt_aff <- tmp.mat %>% filter(species == 'Mus_musculus') %>% pull(norm_affinity)
#     tmp.mat <- tmp.mat %>% mutate(frac_delta = delta_steps/(wt_aff - start_aff))
#     new_aff.mat <- bind_rows(new_aff.mat, tmp.mat)
# }
# new_aff.mat <- new_aff.mat %>% mutate(log2_fc_steps = ifelse(fc_steps == 0, 0, log2(fc_steps)))
# new_aff.mat <- new_aff.mat %>% 
#             mutate(species2 = ifelse(species == 'Mus_musculus', 'Mouse', str_replace(species, "fullTree","")))
# new_aff.mat$species2 <- factor(new_aff.mat$species2, levels = c(
#     ifelse(keep.tips == "Mus_musculus", "Mouse", str_replace(keep.tips, "fullTree",""))
#   ))
# new_aff.mat$TF_name <- factor(new_aff.mat$TF_name, levels = c('Foxa2',"Sox17","Gata4/6","Jun_Atf3","Klf4"))
# # final_aff.mat <- new_aff.mat %>% filter(species%in%anc.cutoff)

# tfbs_map <- data.frame(msa_start = c(114,173,239, 296),
#                       msa_end = c(126,200, 279, 307))

# # new_aff.mat %>% filter(log2_fc_steps > 0) %>% dplyr::select(species, TF_name, rank_TF, log2_fc_steps)
# # epas1_msa_list$msa_func_TFBS %>% arrange(msa_start) %>% dplyr::select(TF_name, rank_TF, msa_start, msa_end)
# ### get activity trajectory from Anc41 to Sand Rat
# species.pair <- c("Psammomys_obesus",'fullTreeAnc41')
# sand_rat.nodes <- make.nodepath(tree, species.pair)
# ### filter nodes if in activity
# sand_rat.nodes <- sand_rat.nodes[sand_rat.nodes%in%unique(oCRE_mean_act$species)]

# start_node <- "fullTreeAnc41"
# order_species <- rev(sand_rat.nodes)  # traverse in your chosen order

# # base table (ordered once)
# sandrat_avg_lineage_mpra <- oCRE_mean_act %>%
#   filter(full_CRE_id == CRE_oi, species %in% sand_rat.nodes) %>%
#   mutate(
#     species2  = ifelse(species == 'Psammomys_obesus', 'Fat sand rat', str_replace(species, "fullTree","")),
#     marker = case_when(
#       species2 %in% c("Anc68", "Anc59") ~ "Transitions",
#       species2 == "Fat sand rat"                 ~ "Fat sand rat",
#       TRUE                                       ~ "Other"
#     ),
#     species2 = factor(species2, levels = c(
#     ifelse(order_species == "Psammomys_obesus", "Fat sand rat", str_replace(order_species, "fullTree",""))
#   ))
#   ) %>%
#   arrange(species)

# sandrat_lineage_mpra <- oCRE_replicate_act %>%
#   filter(full_CRE_id == CRE_oi, species %in% sand_rat.nodes) %>%
#   mutate(
#     species2  = ifelse(species == 'Psammomys_obesus', 'Fat sand rat', str_replace(species, "fullTree","")),
#     marker = case_when(
#       species2 %in% c("Anc68", "Anc59") ~ "Transitions",
#       species2 == "Psammomys_obesus"                 ~ "Sand Rat",
#       TRUE                                       ~ "Other"
#     ),
#     species = factor(species, levels = order_species)
#   ) %>%
#   arrange(species)

# # constants
# start_act_df <- sandrat_lineage_mpra %>%
#   filter(species == start_node) %>%
#   dplyr::select(biol_rep, start_act = MPRA_act)

# # Get sand rat wild-type activity per replicate
# wt_act_df <- sandrat_lineage_mpra %>%
#   filter(species == 'Psammomys_obesus') %>%
#   dplyr::select(biol_rep, wt_act = MPRA_act)

# # stepwise deltas & fold-changes
# sandrat_lineage_mpra_steps <- sandrat_lineage_mpra %>%
#   arrange(biol_rep, match(species, order_species)) %>%
#   group_by(biol_rep) %>%
#   mutate(
#     prev_act     = lag(MPRA_act, default = NA),
#     delta_steps  = MPRA_act - prev_act,
#     fc_steps     = MPRA_act / prev_act,
#     log2FC_steps = log2(fc_steps)
#   ) %>%
#   ungroup() %>%
#   left_join(start_act_df, by = "biol_rep") %>%
#   left_join(wt_act_df, by = "biol_rep") %>%
#   mutate(
#     delta_steps  = if_else(species == start_node, 0, delta_steps),
#     fc_steps     = if_else(species == start_node, 1, fc_steps),
#     log2FC_steps = if_else(species == start_node, 0, log2FC_steps),
#     frac_delta   = delta_steps / (wt_act - start_act),
#     species2     = factor(species2, levels = str_replace(order_species, "^fullTree", ""))
#   )

# sandrat_stepwise_summary <- sandrat_lineage_mpra_steps %>%
#   group_by(species) %>%
#   summarise(
#     avg_delta_steps     = mean(delta_steps, na.rm = TRUE),
#     sd_delta_steps      = sd(delta_steps, na.rm = TRUE),
#     avg_frac_delta      = mean(frac_delta, na.rm = TRUE),
#     sd_frac_delta       = sd(frac_delta, na.rm = TRUE),
#     avg_fc_steps        = mean(fc_steps, na.rm = TRUE),
#     sd_fc_steps         = sd(fc_steps, na.rm = TRUE),
#     avg_log2FC_steps    = mean(log2FC_steps, na.rm = TRUE),
#     sd_log2FC_steps     = sd(log2FC_steps, na.rm = TRUE),
#     p_delta_steps       = ifelse(n() > 1, t.test(delta_steps)$p.value, NA_real_),
#     p_frac_delta        = ifelse(n() > 1, t.test(frac_delta)$p.value, NA_real_),
#     p_log2FC_steps      = ifelse(n() > 1, t.test(log2FC_steps)$p.value, NA_real_)
#   ) %>%
#   ungroup()
# sandrat_stepwise_summary$p_adj_log2FC_steps <- p.adjust(sandrat_stepwise_summary$p_log2FC_steps, method = 'fdr')
# sandrat_stepwise_summary <- sandrat_stepwise_summary %>% 
#     mutate(
#         species2  = ifelse(species == 'Psammomys_obesus', 'Fat sand rat', str_replace(species, "fullTree","")),
#         marker = case_when(
#           species2 %in% c("Anc68", "Anc59") ~ "Transitions",
#           species2 == "Psammomys_obesus"                 ~ "Sand Rat",
#           TRUE                                       ~ "Other"
#         ),
#         species2 = factor(species2, levels = c(
#     ifelse(order_species == "Psammomys_obesus", "Fat sand rat", str_replace(order_species, "fullTree",""))
#   ))
#     )
# sandrat_avg_lineage_mpra$group <- 2
# sandrat_stepwise_summary$group <- 2

# ### plot Epas1 transitions of sand-rat expression
# rodent.tree <- tree_subset(tree, node = 'fullTreeAnc41', levels_back = 0)
# rodent.tree <- as_tibble(rodent.tree)
# cre_mean_act <- oCRE_mean_act %>% filter(full_CRE_id == CRE_oi)
# shape.values <- c(1,24,13,15)
# names(shape.values) <- c('mouse','no ortholog','not sequenced', 'sequenced')
# sub.seq <- cre_mean_act[,c('species','mean_MPRA_act','common_name')]
# sub.seq$label <- sub.seq$species
# sub.seq <- sub.seq %>% mutate(capture = case_when(
#                         species == 'Mus_musculus' ~ 'mouse',
#                         species == 'Psammomys_obesus' ~ 'sand rat',
#                         is.na(mean_MPRA_act) ~ 'not sequenced',
#                         TRUE ~ 'sequenced'))

# sub.tree <- left_join(rodent.tree, sub.seq, by = 'label')
# sub.tree <- sub.tree %>% mutate(capture = case_when(
#                         is.na(capture) ~ 'no ortholog',
#                         TRUE ~ capture)) %>% 
#             mutate(label2 = case_when(
#                 label == 'Mus_musculus' ~ 'Mouse',
#                 label == 'Psammomys_obesus' ~ 'Fat sand rat',
#                 label%in%c(mouse.nodes,sand_rat.nodes) ~ str_replace(label, "fullTree",""),
#                 TRUE ~ NA))
# sub.tips <- sub.tree %>% filter(label%in%c(mouse.nodes,sand_rat.nodes)) %>%
#             mutate(Shape = case_when(
#                 label == 'Mus_musculus' ~ 1,
#                 label == 'Psammomys_obesus' ~ 1,
#                 label%in%c(mouse.nodes, sand_rat.nodes) ~ 9,
#                 TRUE ~ 15)) %>% 
#             mutate(Size = case_when(
#                 label == 'Mus_musculus' ~ 5,
#                 label == 'Psammomys_obesus' ~ 5,
#                 label%in%c(mouse.nodes, sand_rat.nodes) ~ 5,
#                 TRUE ~ 0))
# sub.tips$Shape <- factor(sub.tips$Shape)
# sub.tips <- as.tibble(sub.tips)
# sub.tree <- as.treedata(sub.tree)



# #### get mutation series for Epas1 mouse vs sand-rat
# ### plot Epas1 MSA for five rodents
# CRE_oi <- "Epas1_chr17_10063"
# fasta_file <- '/net/shendure/vol10/projects/tli/nobackup/EB_dev_project/oCREs_v2/nobackup/reoriented_fasta/Epas1_chr17_10063_oCRE.fasta'
# mouse_DMS_MPRA$mut <- toupper(mouse_DMS_MPRA$mut)
# sandrat_epas1_msa_df <- make_alignment_dataframe(fasta_file, c(sand_rat.nodes,'Mus_musculus'), 'Mus_musculus',CRE_oi)
# mouse_msa_df <- sandrat_epas1_msa_df %>% filter(species == 'Mus_musculus')

# sandrat_epas1_msa_list <- make_msa_tfbs_traj_df(fasta_file, c(sand_rat.nodes,'Mus_musculus'), 'Mus_musculus', CRE_oi, cactus.tfbs, mouse_func_tfbs)
# labels <- collect_labels(sandrat_epas1_msa_df$msa_df_TFBS, side = T)
# sandrat_epas1_func_tfbs <- mouse_func_tfbs %>% filter(CRE == CRE_oi)
# mm.start <- match(sandrat_epas1_func_tfbs$TFBS_start, mouse_msa_df$ref_pos)
# mm.end <- match(sandrat_epas1_func_tfbs$TFBS_end, mouse_msa_df$ref_pos)
# sandrat_epas1_func_tfbs$msa_start <- mouse_msa_df$curpos[mm.start]
# sandrat_epas1_func_tfbs$msa_end <- mouse_msa_df$curpos[mm.end]
# sandrat_epas1_msa_df <- left_join(sandrat_epas1_msa_df, mouse_DMS_MPRA %>% filter(CRE == CRE_oi) %>%
#                           dplyr::select(ref_pos = pos, mut, log2FC_MPRA, sd_log2FC_MPRA, MPRA_act, sd_MPRA_act, min_act,max_act, wilcox_BH),
#                          by = c("ref_pos","mut"))

# sandrat_epas1_msa_df <- sandrat_epas1_msa_df %>% mutate(mut2 = case_when(
#                         mut == 'DEL' ~ '-',
#                         TRUE ~ mut))
# sandrat_epas1_msa_df <- sandrat_epas1_msa_df %>% 
#                 mutate(species2 = ifelse(species == 'Mus_musculus', 'Mouse', 
#                         ifelse(species == 'Psammomys_obesus','Fat sand rat', str_replace(species, "fullTree",""))))
# sandrat_epas1_msa_df$species2 <- factor(sandrat_epas1_msa_df$species2, levels = rev(c('Anc41','Anc40','Anc39','Fat sand rat', 'Mouse')))

# sandrat_epas1_msa_df <- sandrat_epas1_msa_df %>%
#   mutate(fill_type = case_when(
#     mut2 == "-" & is.na(log2FC_MPRA) ~ "black",
#     mut == ref_nuc ~ "grey",
#     TRUE ~ "gradient"
#   ))

# ### mark stable mutations from Anc41 to SandRat
# # Step 1: filter for mut_type not equal to "Match" and not Mus_musculus
# filtered_df <- sandrat_epas1_msa_df %>%
#   filter(mut_type != "Match", species != "Mus_musculus")

# # Step 2: For each curpos, check if mut2 is the same across all species
# uniform_mut_df <- filtered_df %>%
#   group_by(curpos) %>%
#   filter(n_distinct(mut2) == 1 & n_distinct(species2) == 4) %>%
#   ungroup()

# sandrat_epas1_msa_df <- sandrat_epas1_msa_df %>%
#   mutate(
#     is_uniform_nonmatch = curpos %in% uniform_mut_df$curpos &
#                           mut_type != "Match" &
#                           species != "Mus_musculus"
#   )

# #### calculate number of mutations and mark which mutations links with a TFBS
# ### calculate avg log2FC for each mutation gained
# all_mouse_tfbs_positions <- unlist(mapply(seq, sandrat_epas1_func_tfbs$msa_start, sandrat_epas1_func_tfbs$msa_end))

# sandrat_epas1_all_fc <- sandrat_epas1_msa_df %>% filter(species2 == 'Fat sand rat') %>% 
#                     filter(is_uniform_nonmatch) %>% arrange(curpos) %>% 
#                     mutate(group_id = cumsum(c(TRUE, diff(curpos) > 1)),
#                       in_TFBS = case_when(
#                           curpos%in%all_mouse_tfbs_positions ~ TRUE,
#                           TRUE ~ FALSE)) %>% 
#                     mutate(tf_call = case_when(
#                         in_TFBS ~ 'func. TFBS',
#                         TRUE ~ 'non func. TFBS'))
# sandrat_epas1_avg_fc <- sandrat_epas1_all_fc %>% group_by(species2, in_TFBS) %>% 
#                     summarise(avg_log2FC = mean(log2FC_MPRA, na.rm = TRUE), sd_log2FC = sd(log2FC_MPRA, na.rm = TRUE))
# sandrat_epas1_snp_counts <- sandrat_epas1_all_fc %>% filter(species2 == 'Fat sand rat') %>%
#               group_by(species2, group_id) %>%
#               summarise(
#                 num_snps = 1,                     # each contiguous run counts as 1
#                 tfbs_count = any(in_TFBS),        # TRUE if any in the run is in TFBS
#                 .groups = "drop_last"
#               ) %>%
#               summarise(
#                 n_snp_count = sum(num_snps),
#                 total_tfbs_snps = sum(tfbs_count),# same as above unless you want total TRUEs
#                 .groups = "drop"
#               ) %>% mutate(total_nontfbs_snps = n_snp_count - total_tfbs_snps) 
# sandrat_epas1_group_fc <- sandrat_epas1_all_fc %>% group_by(tf_call, group_id) %>% 
#                         summarise(avg_log2FC = mean(log2FC_MPRA, na.rm = TRUE), sd_log2FC = sd(log2FC_MPRA, na.rm = TRUE))


# fixed_sandrat_msa_df <- sandrat_epas1_msa_df %>% filter(species != 'Mus_musculus')
# fixed_sandrat_all_fc <- sandrat_epas1_all_fc %>% filter(species2 != 'Mouse')
# fixed_sandrat_snp_counts <- sandrat_epas1_snp_counts %>% filter(species2 != 'Mouse')


# # # Get the set of curpos values you care about
# # valid_curpos <- unique(epas1_all_fc$curpos)

# # epas1_sandrat_msa_snp_tfbs <- epas1_msa_list$msa_df_TFBS %>% filter(species == 'Psammomys_obesus') %>% 
# #   filter(norm_affinity > 0.05 | norm_affinity_ref > 0.05) %>% 
# #   rowwise() %>%
# #   filter(any(seq(msa_start, msa_end) %in% valid_curpos)) %>%
# #   ungroup() %>% mutate(fc = norm_affinity/norm_affinity_ref, tf_call = case_when(
# #                           is.na(TF_name2) ~ 'non func. TFBS',
# #                           TRUE ~ 'func. TFBS')) %>% mutate(log2_aff = log2(fc))


# # ### make dog lineage
# # species.pair <- c("Canis_lupus_familiaris",'fullTreeAnc239')
# # dog.nodes <- make.nodepath(tree, species.pair)
# # ### filter nodes if in activity
# # dog.nodes <- dog.nodes[dog.nodes%in%unique(oCRE_mean_act$species)]

# # start_node <- "fullTreeAnc239"
# # order_species <- rev(dog.nodes)  # traverse in your chosen order

# # # base table (ordered once)
# # dog_avg_lineage_mpra <- oCRE_mean_act %>%
# #   filter(full_CRE_id == CRE_oi, species %in% dog.nodes) %>%
# #   mutate(
# #     species2  = ifelse(species == 'Canis_lupus_familiaris', 'Dog', str_replace(species, "fullTree","")),
# #     marker = case_when(
# #       species2 %in% c("Anc222") ~ "Transitions",
# #       species2 == "Dog"                 ~ "Dog",
# #       TRUE                                       ~ "Other"
# #     ),
# #     species2 = factor(species2, levels = c(
# #     ifelse(order_species == "Canis_lupus_familiaris", "Dog", str_replace(order_species, "fullTree",""))
# #   ))
# #   ) %>%
# #   arrange(species)

# # dog_lineage_mpra <- oCRE_replicate_act %>%
# #   filter(full_CRE_id == CRE_oi, species %in% dog.nodes) %>%
# #   mutate(
# #     species2  = ifelse(species == 'Canis_lupus_familiaris', 'Dog', str_replace(species, "fullTree","")),
# #     marker = case_when(
# #       species2 %in% c("Anc222") ~ "Transitions",
# #       species2 == "Dog"                 ~ "Dog",
# #       TRUE                                       ~ "Other"
# #     ),
# #     species = factor(species, levels = order_species)
# #   ) %>%
# #   arrange(species)

# # # constants
# # start_act_df <- dog_lineage_mpra %>%
# #   filter(species == start_node) %>%
# #   dplyr::select(biol_rep, start_act = MPRA_act)

# # # Get sand rat wild-type activity per replicate
# # wt_act_df <- dog_lineage_mpra %>%
# #   filter(species == 'Canis_lupus_familiaris') %>%
# #   dplyr::select(biol_rep, wt_act = MPRA_act)

# # # stepwise deltas & fold-changes
# # dog_lineage_mpra_steps <- dog_lineage_mpra %>%
# #   arrange(biol_rep, match(species, order_species)) %>%
# #   group_by(biol_rep) %>%
# #   mutate(
# #     prev_act     = lag(MPRA_act, default = NA),
# #     delta_steps  = MPRA_act - prev_act,
# #     fc_steps     = MPRA_act / prev_act,
# #     log2FC_steps = log2(fc_steps)
# #   ) %>%
# #   ungroup() %>%
# #   left_join(start_act_df, by = "biol_rep") %>%
# #   left_join(wt_act_df, by = "biol_rep") %>%
# #   mutate(
# #     delta_steps  = if_else(species == start_node, 0, delta_steps),
# #     fc_steps     = if_else(species == start_node, 1, fc_steps),
# #     log2FC_steps = if_else(species == start_node, 0, log2FC_steps),
# #     frac_delta   = delta_steps / (wt_act - start_act),
# #     species2     = factor(species2, levels = str_replace(order_species, "^fullTree", ""))
# #   )

# # dog_stepwise_summary <- dog_lineage_mpra_steps %>%
# #   group_by(species) %>%
# #   summarise(
# #     avg_delta_steps     = mean(delta_steps, na.rm = TRUE),
# #     sd_delta_steps      = sd(delta_steps, na.rm = TRUE),
# #     avg_frac_delta      = mean(frac_delta, na.rm = TRUE),
# #     sd_frac_delta       = sd(frac_delta, na.rm = TRUE),
# #     avg_fc_steps        = mean(fc_steps, na.rm = TRUE),
# #     sd_fc_steps         = sd(fc_steps, na.rm = TRUE),
# #     avg_log2FC_steps    = mean(log2FC_steps, na.rm = TRUE),
# #     sd_log2FC_steps     = sd(log2FC_steps, na.rm = TRUE),
# #     p_delta_steps       = ifelse(n() > 1, t.test(delta_steps)$p.value, NA_real_),
# #     p_frac_delta        = ifelse(n() > 1, t.test(frac_delta)$p.value, NA_real_),
# #     p_log2FC_steps      = ifelse(n() > 1, t.test(log2FC_steps)$p.value, NA_real_)
# #   ) %>%
# #   ungroup()
# # dog_stepwise_summary$p_adj_log2FC_steps <- p.adjust(dog_stepwise_summary$p_log2FC_steps, method = 'fdr')
# # dog_stepwise_summary <- dog_stepwise_summary %>% 
# #     mutate(
# #         species2  = ifelse(species == 'Canis_lupus_familiaris', 'Dog', str_replace(species, "fullTree","")),
# #         marker = case_when(
# #           species2 %in% c("Anc222") ~ "Transitions",
# #           species2 == "Dog"                 ~ "Dog",
# #           TRUE                                       ~ "Other"
# #         ),
# #         species2 = factor(species2, levels = c(
# #     ifelse(order_species == "Canis_lupus_familiaris", "Dog", str_replace(order_species, "fullTree",""))
# #   ))
# #     )
# # dog_avg_lineage_mpra$group <- 3
# # dog_stepwise_summary$group <- 3


# # dog_mpra.plot <- ggplot() + 
# #                     geom_rect(data = epas1_avg_lineage_mpra %>% filter(species == 'Mus_musculus'), aes(ymin = mean_MPRA_act-sd_MPRA_act, ymax = mean_MPRA_act+sd_MPRA_act, color = NULL), xmin = -Inf, xmax = Inf,
# #                    alpha = .2, fill = "#008000") + 
# #                     geom_hline(data = epas1_avg_lineage_mpra %>% filter(species == 'Mus_musculus') , aes(yintercept = mean_MPRA_act), color = '#008000', linetype = 'dashed') +
# #                     geom_rect(data = control_mean_act %>% filter(CRE_id == 'minP'), aes(ymin = mean_MPRA_act-sd_MPRA_act, ymax = mean_MPRA_act+sd_MPRA_act, color = NULL), xmin = -Inf, xmax = Inf,
# #                    alpha = .2,fill = "#0057e7") +
# #                     geom_hline(data = control_mean_act %>% filter(CRE_id == 'minP'), aes(yintercept = mean_MPRA_act), color = '#0057e7', linetype = 'dashed') +
# #                     geom_line(data = dog_avg_lineage_mpra, aes(x = species2, y = mean_MPRA_act, group = group)) + 
# #                     geom_errorbar(data = dog_avg_lineage_mpra, aes(x = species2, ymin=mean_MPRA_act-sd_MPRA_act, 
# #                                                            ymax=mean_MPRA_act+sd_MPRA_act), 
# #                                   width=.2, position=position_dodge(0.05)) + 
# #                     geom_star(data = dog_avg_lineage_mpra, aes(x = species2, y = mean_MPRA_act, starshape = marker, size = marker, fill = marker)) +  
# #                     scale_starshape_manual(values = c('Other' = 15, "Transitions" = 9, "Dog" = 1, 'Rat' = 1), guide = guide_legend(override.aes = list(size = 2, fill = 'black')), name = "") +
# #                     scale_size_manual(values = c('Other' = 1, "Transitions" = 3, "Dog" = 3, 'Rat' = 3), guide = 'none') +
# #                     scale_fill_manual(values = c('Other' = 'black', "Transitions" = 'red', "Dog" = 'gold'), guide = 'none') +
# # #                     scale_y_log10(
# # #                        breaks = c(10^-1, 10^0, 10^1),
# # #                        labels = scales::trans_format("log10", scales::math_format(10^.x))
# # #                      ) + annotation_logticks(sides = 'l') +
# #                     labs(y = "MPRA activity", x = "") + coord_cartesian(clip="off") + 
# #                     theme_classic() + theme(legend.position = 'bottom', legend.box = "horizontal", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


# # ### make seal lineage
# # species.pair <- c("Neomonachus_schauinslandi",'fullTreeAnc239')
# # seal.nodes <- make.nodepath(tree, species.pair)
# # ### filter nodes if in activity
# # seal.nodes <- seal.nodes[seal.nodes%in%unique(oCRE_mean_act$species)]

# # start_node <- "fullTreeAnc239"
# # order_species <- rev(seal.nodes)  # traverse in your chosen order

# # # base table (ordered once)
# # seal_avg_lineage_mpra <- oCRE_mean_act %>%
# #   filter(full_CRE_id == CRE_oi, species %in% seal.nodes) %>%
# #   mutate(
# #     species2  = ifelse(species == 'Neomonachus_schauinslandi', 'Hawaiian monk seal', str_replace(species, "fullTree","")),
# #     marker = case_when(
# #       species2 %in% c("Anc222",'Anc217') ~ "Transitions",
# #       species2 == "Hawaiian monk seal"                 ~ "Hawaiian monk seal",
# #       TRUE                                       ~ "Other"
# #     ),
# #     species2 = factor(species2, levels = c(
# #     ifelse(order_species == "Neomonachus_schauinslandi", "Hawaiian monk seal", str_replace(order_species, "fullTree",""))
# #   ))
# #   ) %>%
# #   arrange(species)

# # seal_lineage_mpra <- oCRE_replicate_act %>%
# #   filter(full_CRE_id == CRE_oi, species %in% seal.nodes) %>%
# #   mutate(
# #     species2  = ifelse(species == 'Neomonachus_schauinslandi', 'Hawaiian monk seal', str_replace(species, "fullTree","")),
# #     marker = case_when(
# #       species2 %in% c("Anc222",'Anc217') ~ "Transitions",
# #       species2 == "Hawaiian monk seal"                 ~ "Hawaiian monk seal",
# #       TRUE                                       ~ "Other"
# #     ),
# #     species2 = factor(species2, levels = c(
# #     ifelse(order_species == "Neomonachus_schauinslandi", "Hawaiian monk seal", str_replace(order_species, "fullTree",""))
# #   ))
# #   ) %>%
# #   arrange(species)

# # # constants
# # start_act_df <- seal_lineage_mpra %>%
# #   filter(species == start_node) %>%
# #   dplyr::select(biol_rep, start_act = MPRA_act)

# # # Get sand rat wild-type activity per replicate
# # wt_act_df <- seal_lineage_mpra %>%
# #   filter(species == 'Neomonachus_schauinslandi') %>%
# #   dplyr::select(biol_rep, wt_act = MPRA_act)

# # # stepwise deltas & fold-changes
# # seal_lineage_mpra_steps <- seal_lineage_mpra %>%
# #   arrange(biol_rep, match(species, order_species)) %>%
# #   group_by(biol_rep) %>%
# #   mutate(
# #     prev_act     = lag(MPRA_act, default = NA),
# #     delta_steps  = MPRA_act - prev_act,
# #     fc_steps     = MPRA_act / prev_act,
# #     log2FC_steps = log2(fc_steps)
# #   ) %>%
# #   ungroup() %>%
# #   left_join(start_act_df, by = "biol_rep") %>%
# #   left_join(wt_act_df, by = "biol_rep") %>%
# #   mutate(
# #     delta_steps  = if_else(species == start_node, 0, delta_steps),
# #     fc_steps     = if_else(species == start_node, 1, fc_steps),
# #     log2FC_steps = if_else(species == start_node, 0, log2FC_steps),
# #     frac_delta   = delta_steps / (wt_act - start_act),
# #     species2     = factor(species2, levels = str_replace(order_species, "^fullTree", ""))
# #   )

# # seal_stepwise_summary <- seal_lineage_mpra_steps %>%
# #   group_by(species) %>%
# #   summarise(
# #     avg_delta_steps     = mean(delta_steps, na.rm = TRUE),
# #     sd_delta_steps      = sd(delta_steps, na.rm = TRUE),
# #     avg_frac_delta      = mean(frac_delta, na.rm = TRUE),
# #     sd_frac_delta       = sd(frac_delta, na.rm = TRUE),
# #     avg_fc_steps        = mean(fc_steps, na.rm = TRUE),
# #     sd_fc_steps         = sd(fc_steps, na.rm = TRUE),
# #     avg_log2FC_steps    = mean(log2FC_steps, na.rm = TRUE),
# #     sd_log2FC_steps     = sd(log2FC_steps, na.rm = TRUE),
# #     p_delta_steps       = ifelse(n() > 1, t.test(delta_steps)$p.value, NA_real_),
# #     p_frac_delta        = ifelse(n() > 1, t.test(frac_delta)$p.value, NA_real_),
# #     p_log2FC_steps      = ifelse(n() > 1, t.test(log2FC_steps)$p.value, NA_real_)
# #   ) %>%
# #   ungroup()
# # seal_stepwise_summary$p_adj_log2FC_steps <- p.adjust(seal_stepwise_summary$p_log2FC_steps, method = 'fdr')
# # seal_stepwise_summary <- seal_stepwise_summary %>% 
# #     mutate(
# #         species2  = ifelse(species == 'Neomonachus_schauinslandi', 'Hawaiian monk seal', str_replace(species, "fullTree","")),
# #     marker = case_when(
# #       species2 %in% c("Anc222", 'Anc217') ~ "Transitions",
# #       species2 == "Hawaiian monk seal"                 ~ "Hawaiian monk seal",
# #       TRUE                                       ~ "Other"
# #     ),
# #     species2 = factor(species2, levels = c(
# #     ifelse(order_species == "Neomonachus_schauinslandi", "Hawaiian monk seal", str_replace(order_species, "fullTree",""))
# #   ))
# #     )
# # seal_avg_lineage_mpra$group <- 4
# # seal_stepwise_summary$group <- 4

# # #### NOTE maybe focus on seals as there is representation of consistent activity between extant species 
# # ## Neomonachus_schauinslandi, Mirounga_angustirostris, and Leptonychotes_weddellii all within Anc217
# # seal_epas1_msa_df <- make_alignment_dataframe(fasta_file, c('Neomonachus_schauinslandi','Mirounga_angustirostris',
# #                                                             'Leptonychotes_weddellii','Mus_musculus'), 'Mus_musculus',CRE_oi)
# # mouse_msa_df <- seal_epas1_msa_df %>% filter(species == 'Mus_musculus')

# # seal_epas1_msa_list <- make_msa_tfbs_traj_df(fasta_file, c('Neomonachus_schauinslandi','Mirounga_angustirostris',
# #                                                             'Leptonychotes_weddellii','Mus_musculus'), 'Mus_musculus', CRE_oi, cactus.tfbs, mouse_func_tfbs)
# # labels <- collect_labels(seal_epas1_msa_list$msa_df_TFBS, side = T)
# # seal_epas1_func_tfbs <- mouse_func_tfbs %>% filter(CRE == CRE_oi)
# # mm.start <- match(seal_epas1_func_tfbs$TFBS_start, mouse_msa_df$ref_pos)
# # mm.end <- match(seal_epas1_func_tfbs$TFBS_end, mouse_msa_df$ref_pos)
# # seal_epas1_func_tfbs$msa_start <- mouse_msa_df$curpos[mm.start]
# # seal_epas1_func_tfbs$msa_end <- mouse_msa_df$curpos[mm.end]
# # seal_epas1_msa_df <- left_join(seal_epas1_msa_df, mouse_DMS_MPRA %>% filter(CRE == CRE_oi) %>%
# #                           dplyr::select(ref_pos = pos, mut, log2FC_MPRA, sd_log2FC_MPRA, MPRA_act, sd_MPRA_act, min_act,max_act, wilcox_BH),
# #                          by = c("ref_pos","mut"))

# # seal_epas1_msa_df <- seal_epas1_msa_df %>% mutate(mut2 = case_when(
# #                         mut == 'DEL' ~ '-',
# #                         TRUE ~ mut))
# # seal_epas1_msa_df <- seal_epas1_msa_df %>% mutate(species2 = case_when(
# #                                 species == 'Mus_musculus' ~ 'Mouse',
# #                                 species == 'Neomonachus_schauinslandi' ~ 'Hawaiian monk seal',
# #                                 species == 'Mirounga_angustirostris' ~ 'Northen elephant seal',
# #                                 species == 'Leptonychotes_weddellii' ~ 'Weddell seal'))
# # seal_epas1_msa_df$species2 <- factor(seal_epas1_msa_df$species2, levels = c('Hawaiian monk seal','Northen elephant seal','Weddell seal', 'Mouse'))

# # seal_epas1_msa_df <- seal_epas1_msa_df %>%
# #   mutate(fill_type = case_when(
# #     mut2 == "-" & is.na(log2FC_MPRA) ~ "black",
# #     mut == ref_nuc ~ "grey",
# #     TRUE ~ "gradient"
# #   ))

# # #### calculate number of mutations and mark which mutations links with a TFBS
# # ### calculate avg log2FC for each mutation gained
# # all_mouse_tfbs_positions <- unlist(mapply(seq, seal_epas1_func_tfbs$msa_start, seal_epas1_func_tfbs$msa_end))

# # seal_epas1_all_fc <- seal_epas1_msa_df %>% filter(species != 'Mus_musculus') %>% 
# #                 filter(mut_type != 'Match') %>% arrange(curpos) %>% 
# #                     mutate(group_id = cumsum(c(TRUE, diff(curpos) > 1)),
# #                       in_TFBS = case_when(
# #                           curpos%in%all_mouse_tfbs_positions ~ TRUE,
# #                           TRUE ~ FALSE)) %>% 
# #                     mutate(tf_call = case_when(
# #                         in_TFBS ~ 'func. TFBS',
# #                         TRUE ~ 'non func. TFBS'))
# # seal_epas1_avg_fc <- seal_epas1_all_fc %>%
# #   filter(is.finite(log2FC_MPRA)) %>%              # remove -Inf and Inf
# #   group_by(species, in_TFBS) %>%
# #   summarise(
# #     avg_log2FC = mean(log2FC_MPRA, na.rm = TRUE),
# #     sd_log2FC  = sd(log2FC_MPRA, na.rm = TRUE),
# #     .groups = "drop"
# #   )
# # seal_epas1_snp_counts <- seal_epas1_all_fc %>%
# #               group_by(species, group_id) %>%
# #               summarise(
# #                 num_snps = 1,                     # each contiguous run counts as 1
# #                 tfbs_count = any(in_TFBS),        # TRUE if any in the run is in TFBS
# #                 .groups = "drop_last"
# #               ) %>%
# #               summarise(
# #                 n_snp_count = sum(num_snps),
# #                 total_tfbs_snps = sum(tfbs_count),# same as above unless you want total TRUEs
# #                 .groups = "drop"
# #               ) %>% mutate(total_nontfbs_snps = n_snp_count - total_tfbs_snps)

# # # Get the set of curpos values you care about
# # valid_curpos <- unique(seal_epas1_all_fc$curpos)

# # epas1_seal_msa_snp_tfbs <- seal_epas1_msa_list$msa_df_TFBS %>% filter(species != 'Mus_musculus') %>% 
# #   filter(norm_affinity > 0.05 | !is.na(TF_name2)) %>% 
# #   rowwise() %>%
# #   filter(any(seq(msa_start, msa_end) %in% valid_curpos)) %>%
# #   ungroup() %>% mutate(fc = norm_affinity/norm_affinity_ref, tf_call = case_when(
# #                           is.na(TF_name2) ~ 'non func. TFBS',
# #                           TRUE ~ 'func. TFBS')) %>% mutate(log2_aff = log2(fc))


# # seal_epas1_fc_plot <- ggplot(seal_epas1_all_fc, aes(y=tf_call, x = log2FC_MPRA, color = tf_call)) +
# #   geom_vline(xintercept = 0, color = "black", linetype = 'dashed') +
# #   geom_quasirandom_rast(size = 1, alpha = 0.6) +
# #   scale_color_manual(values = c('func. TFBS' = 'red', 'non func. TFBS' = 'black')) +
# #   labs(x = expression(atop(log[2] ~ "FC (MPRA)")), y = "") + facet_wrap(~species2) + 
# #   theme_classic() +  coord_cartesian(clip = 'off') +
# #   theme(
# #     legend.position = 'none'
# #   )

# # seal_epas1_aff_plot <- ggplot(epas1_seal_msa_snp_tfbs, aes(x=log2_aff, y = TF_name, color = TF_name)) +
# #   geom_vline(xintercept = 0, color = "black", linetype = 'dashed') +
# #   geom_quasirandom_rast(size = 1, alpha = 0.6) +
# #   scale_color_manual(values = col_TFs) +
# #   labs(x = expression(log[2]~"FC predicted TFBS affinity"), y = "") + facet_grid(species~tf_call) +
# #   theme_classic() +  coord_cartesian(clip = 'off') +
# #   theme(
# #     legend.position = 'none'
# #   )


# # #
# # species.pair <- c("Canis_lupus_familiaris",'fullTreeAnc239')
# # dog.nodes <- make.nodepath(tree, species.pair)
# # ### filter nodes if in activity
# # dog.nodes <- dog.nodes[dog.nodes%in%unique(oCRE_mean_act$species)]
# # dog.nodes <- c(dog.nodes, 'Mus_musculus')

# # rm.nodes <- intersect(dog.nodes, mouse.nodes)

# # ### dog figure plotting
# # epas_msa_list <- make_msa_tfbs_traj_df(fasta_file, dog.nodes, 'Mus_musculus', CRE_oi, cactus.tfbs, mouse_func_tfbs)
# # labels <- collect_labels(epas_msa_list$msa_df_TFBS %>% filter(!species%in%rm.nodes))[1:11]
# # epas1_dog_gg_tfbs <- plot_tfbs_trajectory(epas_msa_list$msa_df_TFBS, epas_msa_list$msa_func_TFBS, CRE_oi, 'Mus_musculus',labels, 
# #                                      tfbs_size = 0.1, rm_species = rm.nodes)
# # epas1_dog_mpra_bar <- plot_mpra_act_trajectory(oCRE_mean_act, df_mouse_tiles_mpra, control_mean_act, CRE_oi, setdiff(dog.nodes, mouse.nodes)) 
# # epas1_dog_seq_bar <- plot_seqID_trajectory(oCRE_mean_act, df_mouse_tiles_mpra, control_mean_act, CRE_oi, setdiff(dog.nodes, mouse.nodes)) 

# # epas1_dog.df <- make_aff_corr_df(epas_msa_list$msa_df_TFBS, epas_msa_list$msa_func_TFBS, 'Canis_lupus_familiaris', 'Mus_musculus')

# # epas1_dog.aff <- plot_aff_corr_plot(epas1_dog.df) + 
# #                 labs(x = 'Mouse binding affinity', y = 'Dog Anc.\nbinding affinity') + theme(legend.position = 'bottom')
# # written.labels <- collect_labels(epas_msa_list$msa_df_TFBS, written = T)

# # species.pair <- c("Psammomys_obesus",'fullTreeAnc239')
# # sand.nodes <- make.nodepath(tree, species.pair)
# # ### filter nodes if in activity
# # sand.nodes <- sand.nodes[dog.nodes%in%unique(oCRE_mean_act$species)]
# # sand.nodes <- c(sand.nodes, 'Mus_musculus')

# # rm.nodes <- intersect(sand.nodes, mouse.nodes)

# # epas_msa_list <- make_msa_tfbs_traj_df(fasta_file, sand.nodes, 'Mus_musculus', CRE_oi, cactus.tfbs, mouse_func_tfbs)
# # labels <- collect_labels(epas_msa_list$msa_df_TFBS %>% filter(!species%in%rm.nodes))[1:3]
# # epas1_sand_gg_tfbs <- plot_tfbs_trajectory(epas_msa_list$msa_df_TFBS, epas_msa_list$msa_func_TFBS, CRE_oi, 'Mus_musculus',labels, 
# #                                      tfbs_size = 0.1, rm_species = rm.nodes)
# # epas1_sand_mpra_bar <- plot_mpra_act_trajectory(oCRE_mean_act, df_mouse_tiles_mpra, control_mean_act, CRE_oi, setdiff(sand.nodes, mouse.nodes))
# # epas1_sand_seq_bar <- plot_seqID_trajectory(oCRE_mean_act, df_mouse_tiles_mpra, control_mean_act, CRE_oi, setdiff(sand.nodes, mouse.nodes)) 

# # order.colors <- c("steelblue","#cc0000","#ffc425",'#905aa3','#ae5a41','#ff7f50','#1c9e48')
# # names(order.colors) <- c('PRIMATES','RODENTIA','EULIPOTYPHLA','CHIROPTERA',
# #                          'PERISSODACTYLA','CARNIVORA','CETARTIODACTYLA')

# # layout <- '
# #     ###AAAABC
# #     ###AAAABC
# #     ###AAAABC
# #     ###AAAABC
# #     MM#DDDDEF
# #     MM#DDDDEF
# #     MM#DDDDEF
# #     MM#DDDDEF
# #     MM#GGGGHI
# #     MM#JJJJKL
# #     MM#JJJJKL
# #     MM#JJJJKL
# #     '
# # #     GAAAAAB
# # #     GAAAAAB
# # #     GCCCCCD
# # #     GEEEEEF
# # #     GEEEEEF'
# # p <- (gata4_gg_tfbs + ggtitle("Gata4:chr14_5729") + theme(plot.title = element_text(size = 12, face = "bold"))) + gata4_mpra_bar + gata4_seq_bar + (epas1_gg_tfbs + labs(x = "") + ggtitle("Epas1:chr17_10063") + theme(plot.title = element_text(size = 12, face = "bold"))) + epas1_mpra_bar + epas1_seq_bar + (epas1_sand_gg_tfbs + labs(x = "") ) + epas1_sand_mpra_bar + epas1_sand_seq_bar + epas1_dog_gg_tfbs + epas1_dog_mpra_bar + epas1_dog_seq_bar + sub.tree + plot_layout(design = layout)
# # ggsave('fig2_v2.pdf', plot = p, width = 230, height = 240, useDingbats = F, units = 'mm', device = 'pdf')
