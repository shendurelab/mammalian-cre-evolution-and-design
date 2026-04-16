library(tidyverse)
library(castor)
suppressMessages(suppressWarnings(library(GenomicFeatures)))
suppressMessages(suppressWarnings(library(GenomicRanges)))
suppressMessages(suppressWarnings(library(dplyr)))
library(reshape2)
require(grid)
library(stringr)
library(IRanges)
library(seqinr)
library(Biostrings)

make.nodepath <- function(tree, pair){
    ### get all ancestral node trajectory from x species to y species (here x = mouse and y = human)
    all.labels <- c(tree$tip.label, tree$node.label)
    tip.numbers <- match(pair, all.labels)
    pathway <- nodepath(tree, tip.numbers[1], tip.numbers[2])
    return(all.labels[pathway])
}

make_alignment_dataframe_from_seq <- function(seqs, seq_names, ref_name, CRE_oi){
    sequences <- DNAStringSet(seqs)
    
    # set headers
    names(sequences) <- seq_names
    # Create a temporary file to write the FASTA content
    tmp_fasta <- tempfile(fileext = ".fasta")
    # Write the sequences to a new FASTA file
    writeXStringSet(sequences, tmp_fasta)
    
    tmp_align <- tempfile(fileext = ".afa")
    command <- paste("mafft --adjustdirection --thread 4", tmp_fasta,">",tmp_align)
    system(command)
    
    # Function to assign the value based on character
    assign_value <- function(char) {
      if (char == '-') {
        return(NA)
      } else {
        non_dash_count <<- non_dash_count + 1
        return(non_dash_count)
      }
    }
    
    cactus.seq <- data.frame()
    msa <- read.fasta(tmp_align,  as.string = TRUE)
    ref.seq <- toupper(msa[[ref_name]][1])
    
    for(tile in names(msa)){
        seq <- toupper(msa[[tile]][1])
        species <- tile
        tmp <- data.frame(ori_CRE_id = CRE_oi, species = species, seq = seq)
        cactus.seq <- bind_rows(cactus.seq, tmp)
    }
    
    msa.df <- data.frame()
    ref.vec <- unlist(strsplit(ref.seq,split = ""))
    ref.df <- data.frame(ref_nuc = ref.vec, ori_CRE_id = CRE_oi, curpos = seq(1:length(ref.vec)))
    non_dash_count <- 0
    # Apply the function to each character
    values <- sapply(ref.vec, assign_value)
    ref.df$ref_pos <- values
    
    
    for(species_oi in unique(cactus.seq$species)){
        species.seq <- cactus.seq %>% filter(species == species_oi) %>% pull(seq)
        target.vec <- unlist(strsplit(species.seq,split = ""))

        target.df <- data.frame(mut = target.vec, species = species_oi, ori_CRE_id = CRE_oi, curpos = seq(1:length(target.vec)))
        target.df <- left_join(target.df, ref.df, by = c("curpos",'ori_CRE_id'))
        target.df <- target.df %>% mutate(mut_type = case_when(
                              mut == "-" ~ "DEL",
                              ref_nuc == "-" ~ "INS",  
                              mut == ref_nuc ~ "Match",
                              TRUE ~ "Mismatch")) %>%
                    mutate(mut = case_when(
                        mut_type == 'DEL' ~ 'DEL',
                        TRUE ~ mut))
        non_dash_count <- 0
        # Apply the function to each character
        values <- sapply(target.vec, assign_value)
        target.df$query_pos <- values 
        msa.df <- bind_rows(msa.df, target.df)
    }
    
    unlink(tmp_align)
    unlink(tmp_fasta)
    
    return(msa.df)
}


make_alignment_from_fasta <- function(fasta, seq_names, ref_name, in_reference = FALSE){
    sequences <- readDNAStringSet(fasta)
    if(in_reference){ ### align to reference sequence
        non_ref_seqnames <- seq_names[!seq_names%in%c(ref_name)]
        
        non_ref_sequences <- sequences[names(sequences) %in% non_ref_seqnames]
        ref_sequence <- sequences[names(sequences) %in% c(ref_name)]
        
        # Create a temporary file to write the FASTA content
        tmp_non_ref_fasta <- tempfile(fileext = ".fasta")
        # Write the filtered sequences to a new FASTA file
        writeXStringSet(non_ref_sequences, tmp_non_ref_fasta)
        
        # Create a temporary file to write the FASTA content
        ref_fasta <- tempfile(fileext = ".fasta")
        # Write the filtered sequences to a new FASTA file
        writeXStringSet(ref_sequence, ref_fasta)
        
        tmp_align <- tempfile(fileext = ".afa")
        command <- paste("mafft --6merpair --keeplength --addfragments", tmp_non_ref_fasta, ref_fasta,">",tmp_align)
        system(command)
        unlink(tmp_non_ref_fasta)
        unlink(ref_fasta)
    } else{  
        # Filter by selected headers
        filtered_sequences <- sequences[names(sequences)%in%seq_names]
        # Create a temporary file to write the FASTA content
        tmp_fasta <- tempfile(fileext = ".fasta")
        # Write the filtered sequences to a new FASTA file
        writeXStringSet(filtered_sequences, tmp_fasta)
        
        tmp_align <- tempfile(fileext = ".afa")
        command <- paste("mafft --adjustdirection --thread 4", tmp_fasta,">",tmp_align)
        system(command)
        unlink(tmp_fasta)
    }
    
    msa <- read.fasta(tmp_align,  as.string = TRUE)
    names(msa) <- gsub("_R_", "", names(msa))
    unlink(tmp_align)
    return(msa)
}



make_alignment_dataframe <- function(fasta, seq_names, ref_name, CRE_oi, in_reference = FALSE){    
    msa <- make_alignment_from_fasta(fasta, seq_names, ref_name, in_reference = in_reference)

    # Function to assign the value based on character
    assign_value <- function(char) {
      if (char == '-') {
        return(NA)
      } else {
        non_dash_count <<- non_dash_count + 1
        return(non_dash_count)
      }
    }
    
    cactus.seq <- data.frame()
    ref.seq <- toupper(msa[[ref_name]][1])
    
    for(tile in names(msa)){
        seq <- toupper(msa[[tile]][1])
        species <- tile
        tmp <- data.frame(ori_CRE_id = CRE_oi, species = species, seq = seq)
        cactus.seq <- bind_rows(cactus.seq, tmp)
    }
    
    msa.df <- data.frame()
    ref.vec <- unlist(strsplit(ref.seq,split = ""))
    ref.df <- data.frame(ref_nuc = ref.vec, ori_CRE_id = CRE_oi, curpos = seq(1:length(ref.vec)))
    non_dash_count <- 0
    # Apply the function to each character
    values <- sapply(ref.vec, assign_value)
    ref.df$ref_pos <- values
    
    
    for(species_oi in unique(cactus.seq$species)){
        species.seq <- cactus.seq %>% filter(species == species_oi) %>% pull(seq)
        target.vec <- unlist(strsplit(species.seq,split = ""))

        target.df <- data.frame(mut = target.vec, species = species_oi, ori_CRE_id = CRE_oi, curpos = seq(1:length(target.vec)))
        target.df <- left_join(target.df, ref.df, by = c("curpos",'ori_CRE_id'))
        target.df <- target.df %>% mutate(mut_type = case_when(
                              mut == "-" ~ "DEL",
                              ref_nuc == "-" ~ "INS",  
                              mut == ref_nuc ~ "Match",
                              TRUE ~ "Mismatch")) %>%
                    mutate(mut = case_when(
                        mut_type == 'DEL' ~ 'DEL',
                        TRUE ~ mut))
        non_dash_count <- 0
        # Apply the function to each character
        values <- sapply(target.vec, assign_value)
        target.df$query_pos <- values 
        msa.df <- bind_rows(msa.df, target.df)
    }
    
    return(msa.df)
}

make_cre_tfbs_alignment_dataframe <- function(msa.df, tfbs_df, CRE_oi){
    msa_df_TFBS <- data.frame()

    cre.tfbs <- tfbs_df %>% filter(ori_CRE_id == CRE_oi)

    for(species_oi in unique(msa.df$species)){
        ### then filter by species
        print(paste0("processing on ", species_oi))
        species.tfbs <- cre.tfbs %>% filter(species == species_oi)
        species.map <- msa.df %>% filter(species == species_oi)
        
        mm.start <- match(species.tfbs$start, species.map$query_pos)
        mm.end <- match(species.tfbs$end, species.map$query_pos)
        
        species.tfbs$msa_start <- species.map$curpos[mm.start]
        species.tfbs$msa_end <- species.map$curpos[mm.end]
        msa_df_TFBS <- bind_rows(msa_df_TFBS, species.tfbs)
    }
    msa_df_TFBS <- msa_df_TFBS %>% mutate(center_pos = (msa_start + msa_end)/2)
    
    return(msa_df_TFBS)
    
}

make_delta_tfbs_alignment_dataframe <- function(msa_df_TFBS, ref_name, comp_name, 
                                                norm_affinity_cutoff = 0.1){
    ref_df_TFBS <- msa_df_TFBS %>% filter(species == ref_name) %>% filter(norm_affinity > norm_affinity_cutoff)
    comp_df_TFBS <- msa_df_TFBS %>% filter(species == comp_name) %>%
                select(msa_start, msa_end, TF_name, norm_affinity)
    
    delta_df_TFBS <- data.frame()
    
    for(tf in unique(ref_df_TFBS$TF_name)){
        ref_tfs <- ref_df_TFBS %>% filter(TF_name == tf)
        comp_tfs <- comp_df_TFBS %>% filter(TF_name == tf)
        
        len_substring <- unique(ref_tfs %>% pull(substring_l))
        
        
        ref.range <- IRanges(start = ref_tfs$msa_start, end = ref_tfs$msa_end)
        comp.range <- IRanges(start = comp_tfs$msa_start, end = comp_tfs$msa_end)
        
        overlaps <- findOverlaps(ref.range, comp.range, minoverlap = len_substring)
        overlapping_tfbs <- comp_tfs[unique(subjectHits(overlaps)), ]
        
        ref_tfs <- left_join(ref_tfs, overlapping_tfbs, by=c("msa_start","msa_end",'TF_name')) %>% 
                mutate(delta_affinity = norm_affinity.x - norm_affinity.y) %>% 
                select(start_pos, end_pos, affinity, TF_name, ori_CRE_id, species, norm_affinity = norm_affinity.x, 
                      msa_start, msa_end, delta_affinity)
        delta_df_TFBS <- bind_rows(delta_df_TFBS, ref_tfs)
    }
    
    return(delta_df_TFBS)
    
}

rolling_logFC_dms <- function(msa_dms_df, width_av = 8){
  df_median <- msa_dms_df %>% group_by(species, curpos) %>%
    summarize(av_FC=mean(log2FC_MPRA,na.rm=TRUE),
              median_FC=median(log2FC_MPRA,na.rm=TRUE),
              min_FC=min(log2FC_MPRA,na.rm=TRUE),
              max_FC=max(log2FC_MPRA,na.rm=TRUE))

  ma <- function(x, width_av = 6){stats::filter(x, rep(1 / width_av, width_av), sides = 2)}
  width_filter <- c(1,2,3,4,5,6,8,10,15,20)
  df_median_FC_ma <- data.frame()
  for (wf in width_filter){
    ma_oi <- ma(df_median$median_FC,width_av=wf) %>% as.vector()
    df_median_FC_ma <- rbind(df_median_FC_ma,
                             data.frame(curpos=df_median$curpos,
                                        median_FC_ma=ma_oi,
                                        filter_width=wf))
  }

  df_median$rollmean_median <- ma(df_median$median_FC,width_av=width_av)
  return(df_median)  
}