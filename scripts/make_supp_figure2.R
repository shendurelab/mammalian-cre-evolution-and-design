### module load pcre2/10.39; module load R/4.3.1
suppressMessages(library(tidyverse))
library(ggplot2)
library(patchwork)

to_cpm <- function(mat) sweep(mat, 2, colSums(mat), "/") * 1e6
read_pb <- function(path) as.matrix(read.delim(path, row.names = 1, check.names = FALSE))
mean_cpm <- function(path) rowMeans(to_cpm(read_pb(path)))

parietal.cpm <- mean_cpm("../data/RNA_data/parietal_endo_mEB_pseudobulk_count.txt")

### human cell types of interest, mapped to display labels
human.labels <- c(pluripotent = "Pluripotent", mesoderm = "Mesoderm",
                   early_ectoderm = "Ectoderm", endoderm = "Endoderm")
human.cpm <- setNames(lapply(names(human.labels), function(ct)
    mean_cpm(paste0("../data/RNA_data/", ct, "_hEB_pseudobulk_count.txt"))), human.labels)

### 1:1 mouse-human orthologs, restricted to TF genes (same filtering as hEB_mEB_ortholog_correlation_TF.R)
ortho <- read.csv("../data/references/mm10_hg38_geneID_homology.csv", check.names = FALSE)
colnames(ortho) <- c("mouse_ensembl", "mouse_symbol", "human_ensembl", "human_symbol")
ortho <- ortho %>% distinct(mouse_symbol, human_symbol)

mouse.ok <- ortho %>% count(mouse_symbol) %>% filter(n == 1) %>% pull(mouse_symbol)
human.ok <- ortho %>% count(human_symbol) %>% filter(n == 1) %>% pull(human_symbol)
ortho.1to1 <- ortho %>% filter(mouse_symbol %in% mouse.ok, human_symbol %in% human.ok)

tf.genes <- read.csv("../data/references/mm_ebi_go_dbTF_overlaps.csv")
ortho.1to1 <- ortho.1to1 %>%
    filter(mouse_symbol %in% tf.genes$Gene.name, mouse_symbol %in% names(parietal.cpm))
message("1:1 TF orthologs present in parietal_endo pseudobulk: ", nrow(ortho.1to1))

### build long-format table: one row per (TF ortholog pair, human cell type)
cpm.long <- map_dfr(names(human.cpm), function(lab) {
    hc <- human.cpm[[lab]]
    pairs <- ortho.1to1 %>% filter(human_symbol %in% names(hc))
    tibble(geneID = pairs$mouse_symbol,
           parietal_endo = parietal.cpm[pairs$mouse_symbol],
           cell_type = lab,
           cpm_mean = hc[pairs$human_symbol])
})

### per-facet correlations
cpm.corr <- cpm.long %>%
    group_by(cell_type) %>%
    summarise(spearman = cor(cpm_mean, parietal_endo, method = "spearman"),
              pearson = cor(log(cpm_mean + 0.01), log(parietal_endo + 0.01), method = "pearson"),
              pearson_r2 = pearson^2)
message("Correlations (mouse parietal_endo vs. human cell type, TF genes only):")
print(cpm.corr)

cpm.long$cell_type <- factor(cpm.long$cell_type, levels = c("Endoderm", "Pluripotent", "Ectoderm", "Mesoderm"))
cpm.corr$cell_type <- factor(cpm.corr$cell_type, levels = c("Endoderm", "Pluripotent", "Ectoderm", "Mesoderm"))
cpm.long <- cpm.long %>%
    left_join(cpm.corr, by = "cell_type") %>%
    mutate(cell_label = paste0("Human ", cell_type, "\nrho = ", round(spearman, 3)))

### highlight known endoderm-specifying TFs, same grouping/palette convention as
### make_supp_figure1_and_3.R::cpm.plot
final_TFs <- c("Jun_Atf3", "Foxa2", "Gata4/6", "Sox17", "Klf4", "Hnf1b")
col_TFs <- c("firebrick2", "forestgreen", "darkorange2", "dodgerblue1", "mediumorchid3", "goldenrod1")
names(col_TFs) <- final_TFs
check.TFs <- data.frame(tf_name = c('Sox17', 'Gata4', 'Gata6', 'Foxa2', 'Klf4',
                                     'Hnf1b', 'Jun', 'Jund', 'Atf3'),
                         group = c("Sox17", "Gata4/6", "Gata4/6", "Foxa2", "Klf4", "Hnf1b",
                                   "Jun_Atf3", "Jun_Atf3", "Jun_Atf3"))
tf.means <- cpm.long %>% filter(geneID %in% unique(check.TFs$tf_name))
mm <- match(tf.means$geneID, check.TFs$tf_name)
tf.means$tf_group <- check.TFs$group[mm]

cpm.plot <- ggplot() +
    stat_bin2d(data = cpm.long, aes(x = cpm_mean, y = parietal_endo, fill = after_stat(count)), bins = 30) +
    geom_point(data = tf.means, aes(x = cpm_mean, y = parietal_endo, color = tf_group), size = 4, shape = 18) +
    scale_color_manual(values = col_TFs, "endoderm TF group") +
    scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    scale_fill_viridis_c(trans = "log10", name = "gene count") +
    annotation_logticks() +
    theme_classic() +
    facet_wrap(~cell_type, ncol = 4, labeller = as_labeller(setNames(cpm.long$cell_label, cpm.long$cell_type))) +
    labs(x = "mean CPM in human EB cell type\n(TF genes)",
         y = "mean CPM in mouse parietal endoderm\n(mEB, TF genes)") +
    theme(legend.position = 'bottom', panel.spacing = unit(1, "lines"))

### ---------------------------------------------------------------------
### PYS2 vs. mouse mEB cell types (TF genes only, single species - no ortholog
### mapping needed). Adapted from
### make_supp_figure1_and_3.R::cpm.plot
### ---------------------------------------------------------------------

mouse.labels <- c(ecto = "Ectoderm", meso = "Mesoderm", parietal_endo = "Parietal Endoderm", pluri = "Pluripotent")
mouse.avg <- sapply(names(mouse.labels), function(ct) mean_cpm(paste0("../data/RNA_data/", ct, "_mEB_pseudobulk_count.txt")))
colnames(mouse.avg) <- mouse.labels[colnames(mouse.avg)]

gene.name <- read.csv("../data/references/mm10_biomart_GeneIDexport.csv")

pys2.files <- Sys.glob("../data/RNA_data/*_counts.txt")
pys2.mat <- data.frame()
for (mat in pys2.files) {
    rep.name <- str_split(basename(mat), "_counts.txt")[[1]][1]
    tmp <- read.delim(mat, sep = '\t', skip = 1)
    tmp <- tmp[, grepl("Geneid|BYVG9T", names(tmp))]
    colnames(tmp) <- c('geneID', rep.name)
    tmp$geneID <- sapply(str_split(tmp$geneID, "\\."), function(x) head(x, 1))
    mm <- match(tmp$geneID, gene.name$Gene.stable.ID)
    tmp$geneID <- gene.name$Gene.name[mm]
    tmp <- tmp %>%
        filter(!is.na(geneID)) %>%
        group_by(geneID) %>%
        summarise(across(everything(), sum))
    if (nrow(pys2.mat) == 0) {
        pys2.mat <- tmp
    } else {
        pys2.mat <- merge(pys2.mat, tmp, by = c('geneID'))
    }
}
rownames(pys2.mat) <- pys2.mat$geneID
pys2.mat <- pys2.mat[, -1]
pys2.mean <- rowMeans(to_cpm(as.matrix(pys2.mat)))

### restrict to TF genes shared between mEB pseudobulk and PYS2 -- no ortholog mapping needed, both mouse
mouse.tf <- mouse.avg[rownames(mouse.avg) %in% tf.genes$Gene.name & rownames(mouse.avg) %in% names(pys2.mean), , drop = FALSE]
message("TF genes shared between mEB pseudobulk and PYS2: ", nrow(mouse.tf))

pys2.long <- map_dfr(colnames(mouse.tf), function(ct) {
    tibble(geneID = rownames(mouse.tf),
           cpm_mean = mouse.tf[, ct],
           cell_type = ct,
           PYS2 = pys2.mean[rownames(mouse.tf)])
})

pys2.corr <- pys2.long %>%
    group_by(cell_type) %>%
    summarise(spearman = cor(cpm_mean, PYS2, method = "spearman"),
              pearson = cor(log(cpm_mean + 0.01), log(PYS2 + 0.01), method = "pearson"),
              pearson_r2 = pearson^2)
message("PYS2 vs. mEB cell type correlations (TF genes only):")
print(pys2.corr)

pys2.long$cell_type <- factor(pys2.long$cell_type, levels = c("Parietal Endoderm", "Pluripotent", "Ectoderm", "Mesoderm"))
pys2.corr$cell_type <- factor(pys2.corr$cell_type, levels = c("Parietal Endoderm", "Pluripotent", "Ectoderm", "Mesoderm"))
pys2.long <- pys2.long %>%
    left_join(pys2.corr, by = "cell_type") %>%
    mutate(cell_label = paste0("Mouse ", cell_type, "\nrho = ", round(spearman, 3)))

tf.means.pys2 <- pys2.long %>% filter(geneID %in% unique(check.TFs$tf_name))
mm <- match(tf.means.pys2$geneID, check.TFs$tf_name)
tf.means.pys2$tf_group <- check.TFs$group[mm]

pys2.cpm.plot <- ggplot() +
    stat_bin2d(data = pys2.long, aes(x = cpm_mean, y = PYS2, fill = after_stat(count)), bins = 30) +
    geom_point(data = tf.means.pys2, aes(x = cpm_mean, y = PYS2, color = tf_group), size = 4, shape = 18) +
    scale_color_manual(values = col_TFs, "endoderm TF group") +
    scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    scale_fill_viridis_c(trans = "log10", name = "gene count") +
    annotation_logticks() +
    theme_classic() +
    facet_wrap(~cell_type, ncol = 4, labeller = as_labeller(setNames(pys2.long$cell_label, pys2.long$cell_type))) +
    labs(x = "mean CPM gene expression in cell types (mEB, TF genes)",
         y = "mean CPM gene expression\nin PYS-2 (TF genes)") +
    theme(legend.position = 'bottom', panel.spacing = unit(1, "lines"))


layout <- "AAAAAAAA
           BBBBBBBB"

p <- pys2.cpm.plot + cpm.plot + plot_layout(design = layout) +
    plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(size = 12, face='bold'))
ggsave('supp_fig2_RNA.pdf', plot  = p, width = 210, height = 160, useDingbats = F, units = 'mm')
