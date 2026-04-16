
library(dplyr)
library(tidyverse)


args <- commandArgs(trailingOnly=TRUE)
input_file <- args[1]
out_file_full_pileup <- args[2]


# read in BC-CRE merged table
print("reading merged BC+BAM file")
df_BC_CRE <- read.table(input_file, sep="\t",header=TRUE)
print(paste0("Number of paired BC + CRE_id: ", nrow(df_BC_CRE)))
print(paste0("Number of unique headers: ", length(unique(df_BC_CRE$read_id))))
print(table(df_BC_CRE$is_read1))

# pile up, keeping unique tagmentation position information, thresholding on mapping quality
mapq_thresh <- 40

# pile up also over position information. calculate summary stats of position dispersion as a downstream flag 

print("filtering and grouping as by CRE and orientation")
df_BC_CRE_filtered <- df_BC_CRE %>% filter(mapq>mapq_thresh )  
print(paste0("Number of filtered paired BC + CRE_id: ", nrow(df_BC_CRE_filtered)))
df_BC_CRE_full_pileup <- df_BC_CRE_filtered %>% 
        group_by(BC1, CRE_id,read_forward_in_ref, read_start, read_end) %>% summarise(n_count = n()) %>% 
        group_by(BC1, read_forward_in_ref) %>% mutate(n_prct = n_count/sum(n_count))

print("writing to file pile up file")
write.table(df_BC_CRE_full_pileup,out_file_full_pileup,
            sep="\t", row.names=FALSE, quote=FALSE)






