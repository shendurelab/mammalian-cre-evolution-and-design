### module load pcre2/10.39; module load R/4.3.1

library(tidyverse)
library(ggplot2)
library(cowplot)
library(ggrastr)
library(castor)
library(reshape2)
library(readxl)
library(ggridges)
library(tidytree)
library(DescTools)
library(patchwork)
library(ape)
library(phytools)
library(gtools)
library(GenomicRanges)
library(Gviz)
library(ggplotify)
library(rtracklayer)
library(BSgenome.Mmusculus.UCSC.mm10)
library(reshape2)
library(edgeR)
source("utils/mpra_tablemaker.R")

# https://joeystanley.com/blog/custom-themes-in-ggplot2
pub_theme <- function() {
  theme_cowplot(font_size = 6) %+replace%
    theme(
      strip.background = element_rect(color='black', 
                                      fill='grey90',
                                      size=0.5, 
                                      linetype="solid",
                                      linewidth = 0.25),
      text = element_text(size=6, family='Helvetica'),
      strip.text = element_text(family = "Helvetica", size=6, face='bold', margin=margin(b=2.5, t=2.5)),
      axis.text.x = element_text(size=6, color = "black"),
      axis.text.y = element_text(size=6, color = "black"),
      legend.text = element_text(size=6),
      legend.title = element_text(size=6),
      axis.line = element_line(linewidth = 0.25),
      axis.title = element_text(size=7),
      legend.position = 'bottom',
      legend.title.align = 0.5,
      plot.tag = element_text(size = 10, face='bold'),
      plot.title = element_text(hjust = 0.5, face='plain', size=9,
                                margin=margin(b=10))
    )
}
theme_set(pub_theme())

######## plot ATAC seq profile correlations
# Define cell types and file paths
shared_fill_scale <- scale_fill_viridis_c(
  trans = "log10", 
  limits = c(1, 1e4),  # ensure both plots use same range
  breaks = c(10^0,10^2, 10^4), 
  labels = scales::trans_format("log10", scales::math_format(10^.x)),
  name = "peak count"
)


cell_types <- c("PYS2", "Parietal Endoderm", "Pluripotent", "Ectoderm", "Mesoderm")
file_paths <- c("PYS2_signal.bed", "endo_mEB_signal.bed", "pluri_mEB_signal.bed", 
                "ecto_mEB_signal.bed", "meso_mEB_signal.bed")
base_path <- "../data/ATAC_data/"

# Read and merge all data into a single dataframe
signal.df <- map2_dfr(file_paths, cell_types, ~{
  read.delim(file.path(base_path, .x), sep = "\t", header = FALSE) %>%
    setNames(c("chr", "start", "end", "peak_id", "signal")) %>%
    mutate(cell = .y)
})

# Reshape data to long format for ggplot facets
signal.long <- signal.df %>%
  pivot_wider(names_from = cell, values_from = signal) %>%
  pivot_longer(cols = -c(peak_id, chr, start, end, PYS2), names_to = "cell_type", values_to = "signal") %>%
  drop_na()  # Remove missing values

signal.filter <- signal.long %>% filter(PYS2 > 0) %>% filter(signal > 0)
dim(signal.filter)

# Compute correlations
cell.corr <- signal.filter %>%
  group_by(cell_type) %>%
  summarise(correlation = cor(signal, PYS2, method = "spearman"), r2 = correlation^2)

signal.filter$cell_type <- factor(signal.filter$cell_type, levels = c("Parietal Endoderm", "Pluripotent", "Ectoderm", "Mesoderm"))
cell.corr$cell_type <- factor(cell.corr$cell_type, levels = c("Parietal Endoderm", "Pluripotent", "Ectoderm", "Mesoderm"))
signal.filter <- signal.filter %>%
  left_join(cell.corr, by = "cell_type") %>%
  mutate(cell_label = paste0(cell_type, "\nrho = ", round(correlation, 3)))  # use for labeling only


acc.plot <- ggplot() +
  stat_bin2d(data = signal.filter, aes(x = signal, y = PYS2, fill = after_stat(count)), bins = 30) +
  scale_x_log10(breaks = c(10^-2, 10^0, 10^2), labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_log10(breaks = c(10^-2, 10^0, 10^2), labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  shared_fill_scale + 
  annotation_logticks() +
  theme_classic() +
  facet_wrap(~cell_type, ncol = 1, labeller = as_labeller(setNames(signal.filter$cell_label, signal.filter$cell_type))) +
  labs(x = "RPM accessibility in cell types (mEB)", y = "RPM accessibility in PYS-2") +
  theme(legend.position = 'bottom',
        panel.spacing = unit(1, "lines"))

######## plot DA profile correlations
# Define cell types and file paths
file_paths <- c("PYS2_DA_signal.bed", "endo_mEB_DA_signal.bed", "pluri_mEB_DA_signal.bed", 
                "ecto_mEB_DA_signal.bed", "meso_mEB_DA_signal.bed")

# Read and merge all data into a single dataframe
da_signal.df <- map2_dfr(file_paths, cell_types, ~{
  read.delim(file.path(base_path, .x), sep = "\t", header = FALSE) %>%
    setNames(c("chr", "start", "end", "peak_id", "signal")) %>%
    mutate(cell = .y)
})

# Reshape data to long format for ggplot facets
da_signal.long <- da_signal.df %>%
  pivot_wider(names_from = cell, values_from = signal) %>%
  pivot_longer(cols = -c(peak_id, chr, start, end, PYS2), names_to = "cell_type", values_to = "signal") %>%
  drop_na()  # Remove missing values

da_signal.filter <- da_signal.long %>% filter(PYS2 > 0 ) %>% filter(signal > 0)

# Compute correlations
da_cell.corr <- da_signal.long %>%
  group_by(cell_type) %>%
  summarise(correlation = cor(signal, PYS2, method = "spearman"), r2 = correlation^2)

da_signal.filter$cell_type <- factor(da_signal.filter$cell_type, levels = c("Parietal Endoderm", "Pluripotent", "Ectoderm", "Mesoderm"))
da_cell.corr$cell_type <- factor(da_cell.corr$cell_type, levels = c("Parietal Endoderm", "Pluripotent", "Ectoderm", "Mesoderm"))
da_signal.filter <- da_signal.filter %>%
  left_join(da_cell.corr, by = "cell_type") %>%
  mutate(cell_label = paste0(cell_type, "\nrho = ", round(correlation, 3)))  # use for labeling only

da_acc.plot <- ggplot() +
  stat_bin2d(data = da_signal.filter, aes(x = signal, y = PYS2, fill = after_stat(count)), bins = 30) +
  scale_x_log10(breaks = c(10^-2, 10^0, 10^2), labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_log10(breaks = c(10^-2, 10^0, 10^2), labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  shared_fill_scale + 
  annotation_logticks() +
  theme_classic() +
  facet_wrap(~cell_type, ncol = 1, labeller = as_labeller(setNames(da_signal.filter$cell_label, da_signal.filter$cell_type))) +
  labs(x = "RPM accessibility in cell types (mEB)", y = "RPM accessibility in PYS-2") +
  theme(legend.position = 'bottom',
        panel.spacing = unit(1, "lines"))   


#### in vivo correlations
file_paths <- list.files("../data/ATAC_data/Argelaguet_gastrulation_2022/", full.names = T)

# Read and merge all data into a single dataframe
in_vivo_signal.df <- map_dfr(file_paths, ~{
  # Extract cell type from the filename
  cell_type <- basename(.x) %>%
    str_remove("_in_vivo_signal.bed")

  read.delim(.x, sep = "\t", header = FALSE) %>%
    setNames(c("chr", "start", "end", "peak_id", "signal")) %>%
    mutate(cell = cell_type)
})

# Reshape data to long format for ggplot facets
in_vivo_signal.long <- in_vivo_signal.df %>%
  pivot_wider(names_from = cell, values_from = signal) %>%
  pivot_longer(cols = -c(peak_id, chr, start, end, PYS2), names_to = "cell_type", values_to = "signal") %>%
  drop_na()  # Remove missing values

in_vivo_signal.filter <- in_vivo_signal.long %>% filter(PYS2 > 0) %>% filter(signal > 0)
dim(in_vivo_signal.filter)

# Compute correlations
cell.corr <- in_vivo_signal.filter %>%
  group_by(cell_type) %>%
  summarise(correlation = cor(signal, PYS2, method = "spearman"), r2 = correlation^2)
cell.corr <- cell.corr %>% arrange(-correlation)
cell.corr$cell_type <- factor(cell.corr$cell_type, levels = rev(cell.corr$cell_type))

in_vivo.plot <- ggplot(cell.corr, aes(x = cell_type, y = correlation)) +
  geom_point_rast(aes(color = ifelse(cell_type == "Parietal_endoderm", "highlight", "other")), size = 2) +
  scale_color_manual(values = c("highlight" = "#FC8D62", "other" = "darkgrey")) +
  theme_classic() + scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0,1)) +
  labs(x = expression(italic("in vivo")~"cell types"), y = "PYS-2 Spearman correlation") + coord_flip() + 
  theme(legend.position = "none",
        panel.spacing = unit(1, "lines"))

### RNAseq correlations
gene.name <- read.csv("../data/references/mm10_biomart_GeneIDexport.csv")

### generate coutn matrices
eb.files <- data.frame(files = Sys.glob('../data/RNA_data/*_pseudobulk_count.txt'),
                       cell_name = c("Ectoderm", "Mesoderm", "Parietal Endoderm", "Pluripotent"))

eb.mean <- data.frame()
eb.sd <- data.frame()
for(mat in unique(eb.files$files)){
    cell_name <- eb.files %>% filter(files == mat) %>% pull(cell_name)
    tmp <- read.delim(mat, sep ='\t')
    tmp <- cpm(tmp) ### cpm normalize
    tmp_mean <- rowMeans(tmp, na.rm = T)
    tmp_sd <- apply(tmp, 1, sd, na.rm = TRUE)
    tmp_mean <- data.frame(exp=tmp_mean)
    tmp_sd <- data.frame(sd=tmp_sd)
    colnames(tmp_mean) <- cell_name
    colnames(tmp_sd) <- cell_name
    tmp_mean$geneID <- rownames(tmp_mean)
    tmp_sd$geneID <- rownames(tmp_sd)
    if(nrow(eb.mean) == 0){
        eb.mean <- tmp_mean
        eb.sd <- tmp_sd
    } else{
        eb.mean <- merge(eb.mean, tmp_mean, by = c('geneID'))
        eb.sd <- merge(eb.sd, tmp_sd, by = c('geneID'))
    }
}
keep.genes <- unique(eb.mean$geneID)

pys2.files <- Sys.glob('../data/RNA_data/*_counts.txt')
pys2.mat <- data.frame()

count <- 1

for(mat in pys2.files){
    rep.name <- paste0('PYS2_rep',count)
    tmp <- read.delim(mat, sep ='\t', skip = 1)
    tmp <- tmp[ , grepl( "Geneid|BYVG9T" , names( tmp ) ) ]
    colnames(tmp) <- c('geneID',rep.name)
    tmp$geneID <- sapply(str_split(tmp$geneID, "\\."), function(x) head(x, 1))
    mm <- match(tmp$geneID, gene.name$Gene.stable.ID)
    tmp$geneID <- gene.name$Gene.name[mm]
    tmp <- tmp %>% filter(geneID%in%keep.genes) %>%
          group_by(geneID) %>% 
          summarise(across(everything(), sum))
    if(nrow(pys2.mat) == 0){
        pys2.mat <- tmp
    } else{
        pys2.mat <- merge(pys2.mat, tmp, by = c('geneID'))
    }
    count <- count + 1
}
rownames(pys2.mat) <- pys2.mat$geneID
pys2.mat <- pys2.mat[,-1]
pys2.mat <- cpm(pys2.mat)
pys2.mean <- data.frame(PYS2 = rowMeans(pys2.mat))
pys2.sd <- data.frame(PYS2 = apply(pys2.mat, 1, sd, na.rm = TRUE))
pys2.mean$geneID <- rownames(pys2.mean)
pys2.sd$geneID <- rownames(pys2.sd)
                         
                         
cpm.mean <- merge(eb.mean, pys2.mean, by = c('geneID'))
cpm.melt <- melt(cpm.mean)
colnames(cpm.melt) <- c('geneID','cell', "cpm_mean")
cpm.long <- cpm.melt %>%
  pivot_wider(names_from = cell, values_from = cpm_mean) %>%
  pivot_longer(cols = -c(geneID, PYS2), names_to = "cell_type", values_to = "cpm_mean") %>%
  drop_na()  # Remove missing values

# Compute correlations
cpm.corr <- cpm.long %>%
  group_by(cell_type) %>%
  summarise(spearman = cor(cpm_mean, PYS2, method = "spearman"), 
            pearson = cor(log(cpm_mean+0.01), log(PYS2+0.01), method = "pearson"),
            pearson_r2 = pearson^2)

final_TFs <- c("Jun_Atf3","Foxa2","Gata4/6","Sox17","Klf4","Hnf1b")
col_TFs <- c("firebrick2","forestgreen","darkorange2","dodgerblue1","mediumorchid3",'goldenrod1')
names(col_TFs) <- final_TFs 
check.TFs <- data.frame(tf_name = c('Sox17','Gata4','Gata6','Foxa2', 'Klf4',
              'Hnf1b','Jun','Jund','Atf3'), group = c("Sox17","Gata4/6",'Gata4/6','Foxa2','Klf4','Hnf1b',
                                                     'Jun_Atf3','Jun_Atf3','Jun_Atf3'))   
                         
cpm.long$cell_type <- factor(cpm.long$cell_type, levels = c("Parietal Endoderm", "Pluripotent", "Ectoderm", "Mesoderm"))
cpm.corr$cell_type <- factor(cpm.corr$cell_type, levels = c("Parietal Endoderm", "Pluripotent", "Ectoderm", "Mesoderm"))
cpm.long <- cpm.long %>%
  left_join(cpm.corr, by = "cell_type") %>%
  mutate(cell_label = paste0(cell_type, "\nrho = ", round(spearman, 3)))  # use for labeling only
tf.means <- cpm.long %>% filter(geneID%in%unique(check.TFs$tf_name))
mm <- match(tf.means$geneID, check.TFs$tf_name)
tf.means$tf_group <- check.TFs$group[mm]

            
cpm.plot <- ggplot() +
  stat_bin2d(data = cpm.long, aes(x = cpm_mean, y = PYS2, fill = after_stat(count)), bins = 30) +
  geom_point(data = tf.means, aes(x = cpm_mean, y = PYS2, color = tf_group), size = 4, shape= 18) + 
  scale_color_manual(values = col_TFs, "endoderm TF group") + 
  scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_fill_viridis_c(
          trans = "log10", 
          limits = c(1, 1e4),  # ensure both plots use same range
          breaks = c(10^0,10^2, 10^4), 
          labels = scales::trans_format("log10", scales::math_format(10^.x)),
          name = "gene count"
        ) + 
  annotation_logticks() +
  theme_classic() +
  facet_wrap(~cell_type, ncol = 4, labeller = as_labeller(setNames(cpm.long$cell_label, cpm.long$cell_type))) +
  labs(x = "mean CPM gene expression in cell types (mEB)", y = "mean CPM gene expression\nin PYS-2") +
  theme(legend.position = 'bottom',
        panel.spacing = unit(1, "lines"))                         
                         
layout <- "AABB#CCC
           AABB#CCC
           AABB#CCC
           AABB#CCC
           DDDDDDDD"

p <- acc.plot + da_acc.plot + in_vivo.plot + cpm.plot + plot_layout(design = layout)
ggsave('supp_fig1_acc.pdf', plot  = p, width = 210, height = 210, useDingbats = F, units = 'mm')

###### plot tracks
max.tiles <- read.delim("../data/DMS_300bp_tiles_phyloP.bed",
                       sep = '\t', header = F)
colnames(max.tiles) <- c("chr",'start','end','CRE_ID','phyloP')

mm10 = BSgenome.Mmusculus.UCSC.mm10
col.cells <- c("#f1b873","#FC8D62","#95A7CF","#64C0A3","#E994C8")
##### download bigwigs from GEO 
bws <- c("/net/shendure/vol10/projects/tli/nobackup/EB_dev_project/EB_chrombpnet/compare_PYS2_ATAC/PYS2.bw",
        "/net/shendure/vol10/projects/tli/nobackup/EB_dev_project/EB_chrombpnet/merged_bigwigs/endo_normalized.bw",
        "/net/shendure/vol10/projects/tli/nobackup/EB_dev_project/EB_chrombpnet/merged_bigwigs/meso_normalized.bw",
        "/net/shendure/vol10/projects/tli/nobackup/EB_dev_project/EB_chrombpnet/merged_bigwigs/ecto_normalized.bw",
        "/net/shendure/vol10/projects/tli/nobackup/EB_dev_project/EB_chrombpnet/merged_bigwigs/pluri_normalized.bw")
clusters <- c('PYS-2','endo.','meso.','ecto.','pluri.')
atac_config <- data.frame(cluster = clusters, 
                            colour = col.cells,
                            bigwig = bws)
atac_config$colour = as.character(atac_config$colour)

rownames(atac_config) = as.factor(atac_config$cluster)
CLUSTERS <- rownames(atac_config)

gene.bed <- read.delim('../data/references/mm10_biomart_GeneExport.txt', sep = '\t')
gene.bed <- gene.bed[gene.bed$Gene.type == 'protein_coding',]

refGene.mm10 <- import.gff(gzfile("../data/references/mm10.refGene.gtf.gz"))
refGene.mm10 = keepStandardChromosomes(refGene.mm10, pruning.mode = "coarse")
gene_anno <- data.frame(refGene.mm10)
gene_anno <- gene_anno[!is.na(gene_anno$exon_number),]
gene_anno <- gene_anno[gene_anno$gene_name%in%gene.bed$Gene.name, ]

# rename some columns to match requirements
gene_anno <- gene_anno[,c('seqnames','start','end','strand','gene_id','gene_name','transcript_id','type','exon_id')]
colnames(gene_anno) <- c("chromosome", "start", "end", "strand", "gene", "symbol", "transcript", "feature", "exon")

# one bigwig at a time, any number of peaks
get_matrix_from_bigwig <- function(bigwig_path, peak_set) {
    # ensure peak set has fixed width
    stopifnot(length(unique(width(peak_set)))==1)
    
    as.matrix(import(bigwig_path, 
      which=peak_set, as="NumericList"))
}

# calculate top and bottom percentiles for obs/pred for each state
NUM_SAMP = 1000

atac_peaks = read.table("../data/ATAC_data/consensus_peaks.bed")[,1:3]
colnames(atac_peaks) = c("chr", "start", "end")
atac_peaks = GRanges(atac_peaks)
atac.upper_lims = c()
atac.lower_lims = c()

sampled_peaks = sample(atac_peaks, NUM_SAMP)
ATAC.BIGWIGS <- as.vector(atac_config$bigwig)

for (i in seq(1:length(ATAC.BIGWIGS)))  {
    
    vals = as.vector(get_matrix_from_bigwig(ATAC.BIGWIGS[[i]], resize(sampled_peaks, width=500, fix='center')))
    if(i == 1){ ### given PYS-2 is bulk 
        atac.upper_lims = c(atac.upper_lims, quantile(vals, 0.9999))
        atac.lower_lims = c(atac.lower_lims, quantile(vals, 0.01))
    }else{
        atac.upper_lims = c(atac.upper_lims, quantile(vals, 0.999))
        atac.lower_lims = c(atac.lower_lims, quantile(vals, 0.01))
    }
}

names(atac.upper_lims) <- CLUSTERS
names(atac.lower_lims) <- CLUSTERS

get_region_tracks <- function(chr, show_axis=T) {
    # NOTE: this does not perform any max aggregation
    
    bw_tracks = c()
    for (i in CLUSTERS) {
        atac_bw_path = atac_config[i, "bigwig"]
        atac_bw_path <- as.character(atac_bw_path)
        
        atac_track = DataTrack(atac_bw_path, 
                             genome="mm10", 
                             chromosome = chr, 
                             name=sprintf("%s", atac_config[i, "cluster"]),
                             ylim=c(0,atac.upper_lims[[i]]),
#                              transformation = function(x) log2(x + 1),
#                              window=-1, 
                             type="hist", col.title="black",
                             col.histogram=atac_config[i, "colour"],
                             fill.histogram=atac_config[i, "colour"],
                             background.title = atac_config[i, "colour"], 
                               background.panel = "transparent",
                               littleTicks=F)
        
        # Reduce title size
        displayPars(atac_track)$cex.title <- 0.8  # Change to desired size (e.g., 0.5, 0.8)
        
        # don't show axis ticks for each plot
        displayPars(atac_track)$showAxis = FALSE
        if (!show_axis) {
            displayPars(atac_track)$showTitle = FALSE
        } 
        bw_tracks = c(bw_tracks, atac_track) 
    } 
    return(bw_tracks)
}

get_gene_track <- function(chr, gene_anno, show_axis=T) {
    gene_track <- GeneRegionTrack(gene_anno, genome = "mm10", 
                                 chromosome = chr, 
                                 name = "", 
                                 collapseTranscripts="longest",
                                 transcriptAnnotation="symbol",
                                background.title = 'transparent',
                                  background.panel = "transparent",
                                fill='#000000',
                                stackHeight=0.5)
    if (!show_axis) {
            displayPars(gene_track)$showTitle = F
        }
    
    return(gene_track)
}

#### for Gata4 CRE
CHR = "chr14"
FROM = 63230939 - 500
TO = 63231783 + 500
HIGHLIGHT_COL = "black"
MAX_COL = "red"
MAX_COORDS = max.tiles %>% filter(CRE_ID == 'Gata4_chr14_5729')


bw_tracks = get_region_tracks(CHR, show_axis=T)
gene_track = get_gene_track(CHR, gene_anno, show_axis=T)
#### Gata4 CRE max tile
bw_highlighted = HighlightTrack(trackList = c(bw_tracks),
                               start=c(MAX_COORDS$start, 63230939), 
                               end=c(MAX_COORDS$end, 63231783), 
                               chromosome=CHR,
                               col=c(MAX_COL, HIGHLIGHT_COL), fill=c(MAX_COL, HIGHLIGHT_COL), alpha = c(0.2, 0.1))

GATA4 = as.ggplot(~plotTracks(c(bw_highlighted), from=FROM, to=TO,
                             sizes=c(rep(0.2,nrow(atac_config)))))+
       ggtitle("Gata4:chr14_5729") + theme(plot.title = element_text(size = 12, face = "bold"))

ggsave('test.pdf', height = 4, width = 8, plot = GATA4, useDingbats = F)
# Rasterize only the plot elements from plotTracks()
# GATA4_rasterized <- rasterize(GATA4, dpi = 500) # Adjust DPI as needed

#### for Epas1 CRE
CHR = "chr17"
FROM = 86740868 - 500
TO = 86741862 + 500
HIGHLIGHT_COL = "black"
MAX_COL = "red"
MAX_COORDS = max.tiles %>% filter(CRE_ID == 'Epas1_chr17_10063')

bw_tracks = get_region_tracks(CHR, show_axis=T)
gene_track = get_gene_track(CHR, gene_anno, show_axis=T)
#### Gata4 CRE max tile
bw_highlighted = HighlightTrack(trackList = c(bw_tracks),
                               start=c(MAX_COORDS$start,86740868), 
                               end=c(MAX_COORDS$end, 86741862), 
                               chromosome=CHR,
                               col=c(MAX_COL, HIGHLIGHT_COL), fill=c(MAX_COL, HIGHLIGHT_COL), alpha = c(0.2, 0.1))

EPAS1 = as.ggplot(~plotTracks(c(bw_highlighted), from=FROM, to=TO,
                             sizes=c(rep(0.2,nrow(atac_config)))))+
       ggtitle("Epas1:chr17_10063") + theme(plot.title = element_text(size = 12, face = "bold"))

# Rasterize only the plot elements from plotTracks()
# EPAS1_rasterized <- rasterize(EPAS1, dpi = 500) # Adjust DPI as needed
#### for Lama1 CRE
CHR = "chr17"
FROM = 67653831 - 500
TO = 67654710 + 500
HIGHLIGHT_COL = "black"
MAX_COL = "red"
MAX_COORDS = max.tiles %>% filter(CRE_ID == 'Lama1_chr17_7784')

bw_tracks = get_region_tracks(CHR, show_axis=T)
gene_track = get_gene_track(CHR, gene_anno, show_axis=T)
#### Lama1 CRE max tile
bw_highlighted = HighlightTrack(trackList = c(bw_tracks),
                               start=c(MAX_COORDS$start, 67653831), 
                               end=c(MAX_COORDS$end, 67654710), 
                               chromosome=CHR,
                               col=c(MAX_COL, HIGHLIGHT_COL), fill=c(MAX_COL, HIGHLIGHT_COL), alpha = c(0.2, 0.1))

LAMA1 = as.ggplot(~plotTracks(c(bw_highlighted), from=FROM, to=TO,
                             sizes=c(rep(0.2,nrow(atac_config)))))+
       ggtitle("Lama1:chr17_7784") + theme(plot.title = element_text(size = 12, face = "bold"))

#### for Sparc CRE
CHR = "chr11"
FROM = 55428316 - 500
TO = 55429205 + 500
HIGHLIGHT_COL = "black"
MAX_COL = "red"
MAX_COORDS = max.tiles %>% filter(CRE_ID == 'Sparc_chr11_7211')

bw_tracks = get_region_tracks(CHR, show_axis=T)
gene_track = get_gene_track(CHR, gene_anno, show_axis=T)
#### Sparc CRE max tile
bw_highlighted = HighlightTrack(trackList = c(bw_tracks),
                               start=c(MAX_COORDS$start, 55428316), 
                               end=c(MAX_COORDS$end, 55429205), 
                               chromosome=CHR,
                               col=c(MAX_COL, HIGHLIGHT_COL), fill=c(MAX_COL, HIGHLIGHT_COL), alpha = c(0.2, 0.1))

SPARC = as.ggplot(~plotTracks(c(bw_highlighted), from=FROM, to=TO,
                             sizes=c(rep(0.2,nrow(atac_config)))))+
       ggtitle("Sparc:chr11_7211") + theme(plot.title = element_text(size = 12, face = "bold"))

#### for Sox2_2007 CRE
CHR = "chr3"
FROM = 34757573 - 500
TO = 34758631 + 500
HIGHLIGHT_COL = "black"

bw_tracks = get_region_tracks(CHR, show_axis=T)
gene_track = get_gene_track(CHR, gene_anno, show_axis=T)
#### Sox2_2007 CRE max tile
bw_highlighted = HighlightTrack(trackList = c(bw_tracks),
                               start=c(34757573), 
                               end=c(34758631), 
                               chromosome=CHR,
                               col=HIGHLIGHT_COL, fill=HIGHLIGHT_COL, alpha = 0.1)

SOX2_2007 = as.ggplot(~plotTracks(c(bw_highlighted), from=FROM, to=TO,
                             sizes=c(rep(0.2,nrow(atac_config)))))+
       ggtitle("Sox2:chr3_2007") + theme(plot.title = element_text(size = 12, face = "bold"))

#### for Bend5 CRE
CHR = "chr4"
FROM = 111483057 - 500
TO = 111483881 + 500
HIGHLIGHT_COL = "black"
MAX_COL = "red"
MAX_COORDS = max.tiles %>% filter(CRE_ID == 'Bend5_chr4_8201')


bw_tracks = get_region_tracks(CHR, show_axis=T)
gene_track = get_gene_track(CHR, gene_anno, show_axis=T)
#### Bend5 CRE max tile
bw_highlighted = HighlightTrack(trackList = c(bw_tracks),
                               start=c(MAX_COORDS$start,111483057), 
                               end=c(MAX_COORDS$end, 111483881), 
                               chromosome=CHR,
                               col=c(MAX_COL, HIGHLIGHT_COL), fill=c(MAX_COL, HIGHLIGHT_COL), alpha = c(0.2, 0.1))

BEND = as.ggplot(~plotTracks(c(bw_highlighted), from=FROM, to=TO,
                             sizes=c(rep(0.2,nrow(atac_config)))))+
       ggtitle("Bend5:chr4_8201") + theme(plot.title = element_text(size = 12, face = "bold"))
# Rasterize only the plot elements from plotTracks()
# GATA4_rasterized <- rasterize(GATA4, dpi = 500) # Adjust DPI as needed

#### for Lamb CRE
CHR = "chr12"
FROM = 31221343 - 500
TO = 31221863 + 500
HIGHLIGHT_COL = "black"

bw_tracks = get_region_tracks(CHR, show_axis=T)
gene_track = get_gene_track(CHR, gene_anno, show_axis=T)
#### Lamb CRE max tile
bw_highlighted = HighlightTrack(trackList = c(bw_tracks),
                               start=c(31221343), 
                               end=c(31221863), 
                               chromosome=CHR,
                               col=HIGHLIGHT_COL, fill=HIGHLIGHT_COL, alpha = 0.1)

LAMB = as.ggplot(~plotTracks(c(bw_highlighted), from=FROM, to=TO,
                             sizes=c(rep(0.2,nrow(atac_config)))))+
       ggtitle("Lamb1:chr12_2183") + theme(plot.title = element_text(size = 12, face = "bold"))

# Rasterize only the plot elements from plotTracks()
# EPAS1_rasterized <- rasterize(EPAS1, dpi = 500) # Adjust DPI as needed
#### for Lamc1 CRE
CHR = "chr1"
FROM = 153376348 - 500
TO = 153377201 + 500
HIGHLIGHT_COL = "black"

bw_tracks = get_region_tracks(CHR, show_axis=T)
gene_track = get_gene_track(CHR, gene_anno, show_axis=T)
#### Lamc1 CRE max tile
bw_highlighted = HighlightTrack(trackList = c(bw_tracks),
                               start=c(153376348), 
                               end=c(153377201), 
                               chromosome=CHR,
                               col=HIGHLIGHT_COL, fill=HIGHLIGHT_COL, alpha = 0.1)

LAMC1 = as.ggplot(~plotTracks(c(bw_highlighted), from=FROM, to=TO,
                             sizes=c(rep(0.2,nrow(atac_config)))))+
       ggtitle("Lamc1:chr1_12189") + theme(plot.title = element_text(size = 12, face = "bold"))

#### for Foxa2 CRE
CHR = "chr2"
FROM = 148117185 - 500
TO = 148118133 + 500
HIGHLIGHT_COL = "black"

bw_tracks = get_region_tracks(CHR, show_axis=T)
gene_track = get_gene_track(CHR, gene_anno, show_axis=T)
#### Foxa2 CRE max tile
bw_highlighted = HighlightTrack(trackList = c(bw_tracks),
                               start=c(148117185), 
                               end=c(148118133), 
                               chromosome=CHR,
                               col=HIGHLIGHT_COL, fill=HIGHLIGHT_COL, alpha = 0.1)

FOXA = as.ggplot(~plotTracks(c(bw_highlighted), from=FROM, to=TO,
                             sizes=c(rep(0.2,nrow(atac_config)))))+
       ggtitle("Foxa2:chr2_13858") + theme(plot.title = element_text(size = 12, face = "bold"))

#### for Sox2_2009 CRE
CHR = "chr3"
FROM = 34760351 - 500
TO = 34762064 + 500
HIGHLIGHT_COL = "black"

bw_tracks = get_region_tracks(CHR, show_axis=T)
gene_track = get_gene_track(CHR, gene_anno, show_axis=T)
#### Sox2_2009 CRE max tile
bw_highlighted = HighlightTrack(trackList = c(bw_tracks),
                               start=c(34760351), 
                               end=c(34762064), 
                               chromosome=CHR,
                               col=HIGHLIGHT_COL, fill=HIGHLIGHT_COL, alpha = 0.1)

SOX2_2009 = as.ggplot(~plotTracks(c(bw_highlighted), from=FROM, to=TO,
                             sizes=c(rep(0.2,nrow(atac_config)))))+
       ggtitle("Sox2:chr3_2009") + theme(plot.title = element_text(size = 12, face = "bold"))


layout <- 'ABCD#E
           FGHI#J'

p <- GATA4 + EPAS1 + LAMA1 + SPARC + SOX2_2007 + BEND + LAMB + LAMC1 + FOXA + SOX2_2009 + plot_layout(design = layout)
ggsave('supp_fig1_DAsites.pdf', plot = p, width = 12, height = 5)

########## plot MPRA quality controls
df_counts_w_CREs2 <- read.table("../data/oCRE_and_CRE_optimization_MPRA_count_table.txt.gz",
                                header=TRUE)
dna_BC_recovery_per_lib <- df_counts_w_CREs2 %>% group_by(class,lib_name) %>% 
                    summarize(frac_BC_recovered=mean(!is.na(UMI)), frac_BC_recovered_2plus=sum(UMI>=2,na.rm=TRUE)/length(UMI)) %>% 
                    filter(startsWith(lib_name, 'DNA')) %>% filter(class %in% c("pos_control","neg_control","Engreitz_control","full_CRE"))
mean(dna_BC_recovery_per_lib$frac_BC_recovered_2plus)
# [1] 0.9935596

coverage_lib <- calculate_MPRA_coverage(df_counts_w_CREs2)
df_wide <- generate_MPRA_wide_table(df_counts_w_CREs2, coverage_lib)

shared_fill_scale <- scale_fill_viridis_c(
  trans = "log10", 
  limits = c(1, 1e4),  # ensure both plots use same range
  breaks = c(10^0,10^2, 10^4), 
  labels = scales::trans_format("log10", scales::math_format(10^.x)),
  name = "BC count"
)

generate_DNA_BC_R2_plots_facet <- function(data) {
  library(data.table)

  # Step 1: Cast to wide format
  bc_table_cast <- dcast(data, BC ~ biol_rep, value.var = 'UMI_DNA', fun.aggregate = mean)
  bc_table_cast <- bc_table_cast[complete.cases(bc_table_cast), ]
  bc_table_cast <- log2(bc_table_cast[ , -1] + 1)  # Remove BC column, log transform

  # Rename columns
  colnames(bc_table_cast) <- paste0("rep", colnames(bc_table_cast))

  # Step 2: Create long-form data for all replicate pairs
  pair_df <- bind_rows(
    tibble(x = bc_table_cast$rep1, y = bc_table_cast$rep2, pair = "rep1 vs rep2"),
    tibble(x = bc_table_cast$rep2, y = bc_table_cast$rep3, pair = "rep2 vs rep3"),
    tibble(x = bc_table_cast$rep3, y = bc_table_cast$rep1, pair = "rep3 vs rep1")
  )

  # Step 3: Calculate R² values
  rep.corrs <- pair_df %>%
    group_by(pair) %>%
    summarise(R = cor(x, y, method = "spearman"), .groups = "drop")

  # Step 4: Join R² back into data for plotting
  pair_df <- left_join(pair_df, rep.corrs, by = "pair")

  # Step 5: Plot with facet_wrap
  p <- ggplot(pair_df, aes(x = x, y = y)) +
    stat_bin2d(data = pair_df, aes(x = x, y = y, fill = after_stat(count)), bins = 30) +
    scale_fill_viridis_c(name = "BC count", 
                     trans = "log10", 
                      limits = c(1, 1e4),  # ensure both plots use same range
                      breaks = c(10^0,10^1,10^2, 10^3, 10^4), 
                      labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
    facet_wrap(~ pair, ncol = 3) +
    labs(x = "log2(DNA UMI reads)",
         y = "log2(DNA UMI reads)") +
    geom_text(data = rep.corrs,
              aes(x = Inf, y = -Inf, label = paste0("rho = ", round(R, 3))),
              inherit.aes = FALSE,
              hjust = 1.1, vjust = -0.6) +
    theme_classic() +
    theme(panel.spacing = unit(1, "lines"))

  # Optional: print mean R²
  print(paste0("Mean rho = ", round(mean(rep.corrs$R), 3)))

  return(p)
}


DNA_bc_corr.plot <- generate_DNA_BC_R2_plots_facet(df_wide %>% filter(class %in% c("pos_control","neg_control","Engreitz_control","full_CRE")))
control_bc_count <- df_wide %>% filter(class %in% c("pos_control","neg_control","Engreitz_control","full_CRE")) %>% 
                    group_by(class, CRE_id, biol_rep) %>% summarise(n_BC = n_distinct(BC)) %>% group_by(class) %>% 
                    summarise(mean_BC = mean(n_BC), sd_BC = sd(n_BC))
control_oi_wide <- df_wide %>% filter(class %in% c("pos_control","neg_control","Engreitz_control"))

df_mpra_act <- generate_winsorize_MPRA_table(df_wide)
control_rep_act <- df_mpra_act %>% filter(class%in%c('neg_control','pos_control','Engreitz_control')) 
control_rep_act <- dcast(control_rep_act, CRE_id + class ~ biol_rep, value.var = 'MPRA_act')
colnames(control_rep_act) <- c('CRE_id','class','rep1','rep2','rep3')
control_rep_act <- control_rep_act %>% mutate(class = case_when(
                    class == 'neg_control' ~ 'Negative control (minP)',
                    class == 'pos_control' ~ 'Positive control (EEF1A1p)',
                    class == 'Engreitz_control' ~ 'IGVF promoter series'))
control.plot <- ggplot(control_rep_act, aes(x = rep1, y = rep2, fill = class)) + 
        geom_point_rast(size = 2, pch = 21, color = 'black', alpha = 0.8) + scale_fill_manual(values = c("#008744","#0057e7","#d62d20"), name = '') +
            scale_x_log10(
               breaks = c(10^-1, 10^0, 10^1, 10^2),limits = c(0.05, 100),
               labels = scales::trans_format("log10", scales::math_format(10^.x))
             ) +
             scale_y_log10(
               breaks = c(10^-1, 10^0, 10^1, 10^2), limits = c(0.05, 100),
               labels = scales::trans_format("log10", scales::math_format(10^.x))
             )+ geom_abline(intercept=0, slope = 1,linetype="dashed", alpha = 0.5)+ theme_classic() +
              annotation_logticks() + labs(x = 'MPRA activity replicate 1', y = 'MPRA activity replicate 2') + 
            theme(legend.position = c(0.35, 0.9),
                 legend.background = element_blank(),   # Removes legend background
            legend.key = element_blank())           # Removes individual legend key background

control_mean_act <- df_mpra_act %>% filter(class%in%c('neg_control','pos_control')) %>% group_by(CRE_id) %>% 
        summarise(q95 = quantile(MPRA_act, probs = 0.95), mean_MPRA_act = mean(MPRA_act), sd_MPRA_act = sd(MPRA_act))
control_fc <- foldchange(control_mean_act %>% filter(CRE_id == 'EEF1A1p') %>% pull(mean_MPRA_act),
                        control_mean_act %>% filter(CRE_id == 'minP') %>% pull(mean_MPRA_act))
control_fc
# [1] 536.9103
CRE.list <- c("EEF1A1p","Gata4:chr14_5729","Epas1:chr17_10063","Lama1:chr17_7784",
              "Sparc:chr11_7211","Bend5:chr4_8201","Lamb1:chr12_2183","Lamc1:chr1_12189",
              "Foxa2:chr2_13858","Sox2:chr3_2007","Sox2:chr3_2009","minP")  
full_CRE_act <- df_mpra_act %>% filter(class == 'full_CRE') %>% #filter(CRE_id != "Sox2_chr3_2007") %>% filter(CRE_id != "Sox2_chr3_2009") %>% 
            separate(CRE_id, into = c("gene",'chr','id'), sep = '_') %>% mutate(CRE_id = paste0(gene,":",chr,"_",id)) %>% 
            select(-gene, -chr,-id) %>% mutate(cell_type = case_when(
                        CRE_id%in%c("Sox2:chr3_2007","Sox2:chr3_2009") ~ 'pluri.',
                        TRUE ~ 'endo.'))
control_act <- df_mpra_act %>% filter(class%in%c('neg_control','pos_control')) %>% 
                mutate(cell_type = CRE_id)
full_CRE_act <- bind_rows(full_CRE_act, control_act)
full_CRE_mean_act <- full_CRE_act %>% group_by(CRE_id) %>% 
        summarise(q95 = quantile(MPRA_act, probs = 0.95), mean_MPRA_act = mean(MPRA_act), sd_MPRA_act = sd(MPRA_act))
full_CRE_act$CRE_id <- factor(full_CRE_act$CRE_id, levels = CRE.list)
pos_mean_act <- control_mean_act %>% filter(CRE_id == 'EEF1A1p') %>% mutate(xmin = mean_MPRA_act - sd_MPRA_act, xmax = mean_MPRA_act + sd_MPRA_act)
neg_mean_act <- control_mean_act %>% filter(CRE_id == 'minP') %>% mutate(xmin = mean_MPRA_act - sd_MPRA_act, xmax = mean_MPRA_act + sd_MPRA_act)
library(ggbeeswarm)
x_breaks = c(10^-1, 10^0, 10^1, 10^2)

fullCRE.plot <- ggplot() + 
        geom_rect(data = neg_mean_act, aes(xmin = xmin, xmax = xmax), ymin = -Inf, ymax = Inf,
           alpha = .4,fill = "#0057e7")  +
        geom_vline(xintercept = control_mean_act %>% filter(CRE_id == 'minP') %>% pull(mean_MPRA_act) ,color = "#0057e7", linetype = 'dashed') + 
        annotate("text", 
           x = control_mean_act %>% filter(CRE_id == 'minP') %>% pull(mean_MPRA_act), 
           y = Inf, 
           label = "minP", 
           angle = 270, 
           vjust = -0.5,
           hjust = -0.5,
           size = 4, 
           color = "#0057e7") +
        geom_rect(data = pos_mean_act, aes(xmin = xmin, xmax = xmax), ymin = -Inf, ymax = Inf,
           alpha = .4,fill = "#d62d20")  +
        geom_vline(xintercept = control_mean_act %>% filter(CRE_id == 'EEF1A1p') %>% pull(mean_MPRA_act) ,color = "#d62d20", linetype = 'dashed') + 
        annotate("text", 
           x = control_mean_act %>% filter(CRE_id == 'EEF1A1p') %>% pull(mean_MPRA_act), 
           y = Inf, 
           label = "EEF1A1p", 
           angle = 270, 
           vjust = -0.5, 
           hjust = -0.5,
           size = 4, 
           color = "#d62d20") +
        geom_vline(xintercept = 0.3 ,color = "grey", linetype = 'dashed') + 
        geom_quasirandom(data = full_CRE_act, aes(y = CRE_id, x = MPRA_act, fill = cell_type), 
        dodge.width=1, size = 2, pch = 21, color = 'black', alpha = 0.8) +
        scale_fill_manual(values = c("#d62d20", "#FC8D62","#0057e7","#E994C8"), name = "") +
            scale_x_log10(
               breaks = x_breaks, limits = c(0.05, 100),
               labels = scales::trans_format("log10", scales::math_format(10^.x))
             )+ annotation_logticks(sides = 'b') + theme_classic() + 
            labs(x = 'MPRA activity', y = '')

layout <- '
        AAAAAAAAAAA
        BBBBBBCCCCC
        '
p <- DNA_bc_corr.plot + control.plot + fullCRE.plot + plot_layout(design = layout)
ggsave('supp_fig2_QC.pdf', plot  = p, width = 210, height = 150, useDingbats = F, units = 'mm')

