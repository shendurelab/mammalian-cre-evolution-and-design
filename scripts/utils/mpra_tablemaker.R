### module load pcre2/10.39; module load R/4.3.1

library(tidyverse)
library(ggplot2)
library(cowplot)
library(ggrastr)
library(castor)
library(ape)
library(phytools)
library(reshape2)
library(DescTools)

calculate_MPRA_coverage <- function(data){
    df_metad <- data %>% filter(is_duplicate == 'False') %>% distinct(lib_name, material, biol_rep)
    
    df_BC_recovery_per_lib_simple <- data %>% filter(is_duplicate == 'False') %>% group_by(lib_name) %>% 
    summarize(frac_BC_recovered=mean(!is.na(UMI)), 
              on_target_reads=sum(reads,na.rm=TRUE),
              on_target_UMI=sum(UMI,na.rm=TRUE),
              frac_BC_recovered_2plus=sum(UMI>=2,na.rm=TRUE)/length(UMI))
  
    df_BC_recovery_per_lib_simple <- df_BC_recovery_per_lib_simple %>% transform(on_target_reads=on_target_reads/1e6,
                                                                                on_target_UMI=on_target_UMI/1e6,
                                                                                saturation=1- (on_target_UMI/on_target_reads))
    df_BC_recovery_per_lib_simple <- df_BC_recovery_per_lib_simple %>% left_join(df_metad)
    return(df_BC_recovery_per_lib_simple)
    
}

generate_MPRA_wide_table <- function(data, coverage_lib){
    df_read_coverage <- coverage_lib %>% filter(material != 'plasmid') %>% dplyr::select(biol_rep,material,on_target_UMI) %>% 
        pivot_wider(id_cols=biol_rep,names_from=material, values_from=on_target_UMI,names_prefix="on_target_UMI_")
    
    df_oi_wide <- data %>% filter(material != 'plasmid') %>% 
    pivot_wider(id_cols=c(BC,CRE_id,biol_rep,class),
                values_from=c(UMI,reads),
                names_from=material)
  df_oi_wide2 <- df_oi_wide %>% filter(!is.na(CRE_id)) %>% left_join(df_read_coverage)
  df_oi_wide2[is.na(df_oi_wide2)] <- 0
    return(df_oi_wide2)
}

generate_winsorize_MPRA_table <- function(data_wide, DNA_UMI_thresh_oi = 2, winsor_cut = 0.01){
    df_mpra_act_oi <- data_wide %>% filter(UMI_DNA>=DNA_UMI_thresh_oi) %>% 
      group_by(biol_rep,CRE_id,class) %>% 
      summarize(n_BC=length(BC),
                MPRA_act=sum(Winsorize(UMI_RNA/on_target_UMI_RNA,probs=c(0,1-winsor_cut)))/sum(Winsorize(UMI_DNA/on_target_UMI_DNA,probs=c(0,1-winsor_cut))),
                sum_DNA_umi=sum(UMI_DNA),
                sum_RNA_umi=sum(UMI_RNA),
                norm_DNA_umi = sum(UMI_DNA/on_target_UMI_DNA),
                norm_RNA_umi = sum(UMI_RNA/on_target_UMI_RNA),
                std_DNA_umi=sd(UMI_DNA),
                std_RNA_umi=sd(UMI_RNA))
    return(df_mpra_act_oi)
}

generate_DNA_BC_R2_plots <- function(data){
    bc_table_cast <- dcast(data, BC ~ biol_rep, value.var = 'UMI_DNA', fun.aggregate = mean) ### have to take mean given duplicate CREs
    bc_table_cast <- bc_table_cast[,-1]
    bc_table_cast <- bc_table_cast[complete.cases(bc_table_cast),]
    bc_table_cast <- log(bc_table_cast + 1)
    colnames(bc_table_cast) <- paste0('rep',colnames(bc_table_cast))
    rep.corrs <- data.frame(replicates = c("rep1-rep2",'rep2-rep3','rep3-rep1'),
                           corr_sq = c(cor(bc_table_cast$rep1, bc_table_cast$rep2, method = 'spearman')^2,
                            cor(bc_table_cast$rep2, bc_table_cast$rep3, method = 'spearman')^2,
                            cor(bc_table_cast$rep1, bc_table_cast$rep3, method = 'spearman')^2))
    
    plot_pair <- function(replicate1, replicate2, corr) {
      rep.comp <- paste0(replicate1, "-",replicate2)
      corr <- corr %>% filter(replicates == rep.comp)  
      ggplot(bc_table_cast, aes_string(x = replicate1, y = replicate2)) +
        geom_point_rast(alpha = 0.5) +
        geom_smooth(method = "lm", col = "blue", se = FALSE) +
        labs(x = paste(replicate1, "\nlog(DNA UMI reads)"),
             y = paste(replicate2, "\nlog(DNA UMI reads)"),
             subtitle = paste("R2 =", round(corr$corr_sq, 3))) +
#         scale_x_log10(
#            breaks = scales::trans_breaks("log10", function(x) 10^x),
#            labels = scales::trans_format("log10", scales::math_format(10^.x))
#          ) +
#          scale_y_log10(
#            breaks = scales::trans_breaks("log10", function(x) 10^x),
#            labels = scales::trans_format("log10", scales::math_format(10^.x))
#          ) + annotation_logticks() + 
        theme_classic() + theme(text = element_text(size = 14))
    }

    p1 <- plot_pair('rep1','rep2',rep.corrs)
    p2 <- plot_pair('rep2','rep3',rep.corrs)
    p3 <- plot_pair('rep3','rep1',rep.corrs)
    
    print(paste0("Mean R2 = ",mean(rep.corrs$corr_sq)))
    corr.plots <- plot_grid(p1, p2, p3, ncol = 3, align = 'vh')
    return(corr.plots)
}
