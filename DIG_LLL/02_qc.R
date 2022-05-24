library(dplyr)
library(data.table)
library(stringr)
library(BuenColors)

setwd('~/RWorkSpace_local/CITE-seq/Duerr/DOGMA-seq/DIG_LLL/code/')

source('../../global_functions/estimateLibraryComplexity.R')

process_qc <- function(dt, what ){
  dt %>% filter(is_cell == 1) %>% arrange(desc(atac_fragments)) %>% 
    mutate(pct_mito = atac_mitochondrial_reads/ atac_raw_reads) %>%
    mutate(assay = what)  %>%
    mutate(atac_total = atac_dup_reads + atac_fragments) %>%
    mutate(gex_genes_count = nFeature_RNA) -> odf
  odf$atac_complexity <- sapply(1:dim(odf)[1], function(i){
   estimateLibrarySize(odf$atac_total[i] ,odf$atac_fragments[i])
  })
  odf[,c("pct_mito", "atac_complexity", "gex_genes_count","assay", "condition")]
}

load('../data/DIG_data.RData')
DIG <- process_qc(data@meta.data, "DIG")
rm(data)
load('../data/LLL_data.RData')
LLL <-  process_qc(data@meta.data, "LLL")
rm(data)

all_df_qc <- rbind(DIG, LLL)

pMito <- ggplot(all_df_qc, aes(x = assay, y =pct_mito*100, color = assay)) +
  geom_boxplot(outlier.shape = NA) + 
  pretty_plot(fontsize = 7) + L_border() + scale_y_continuous(limits = c(0,60)) +
  scale_color_manual(values = c("dodgerblue3", "firebrick", "darkgrey")) +
  labs(x = "Condition", y = "%mtDNA (ATAC)", color = "") +
  theme(legend.position = "none")

pATAC <- ggplot(all_df_qc, aes(x = assay, y = log10(atac_complexity), color = assay)) +
  geom_boxplot(outlier.shape = NA) + 
  pretty_plot(fontsize = 7) + L_border() + scale_y_continuous(limits = c(2, 5.5)) +
  scale_color_manual(values = c("dodgerblue3", "firebrick", "darkgrey")) +
  labs(x = "Condition", y = "log10 ATAC Complexity", color = "") +
  theme(legend.position = "none")

pRNA <- ggplot(all_df_qc, aes(x = assay, y = (gex_genes_count), color = assay)) +
  geom_boxplot(outlier.shape = NA) + 
  pretty_plot(fontsize = 7) + L_border() + scale_y_log10(limits = c(10^2, 10^4), breaks = c(100, 500, 1000, 5000, 10000) ) +
  scale_color_manual(values = c("dodgerblue3", "firebrick", "darkgrey")) +
  labs(x = "Condition", y = "# of genes detected (RNA)", color = "") +
  theme(legend.position = "none")

cowplot::ggsave2(cowplot::plot_grid(
  pATAC, pMito, pRNA, nrow = 1
), file = "../plots/all_QC.pdf", width = 3, height = 1.5)

all_df_qc %>% group_by(assay) %>%
  summarize(mito = median(pct_mito*100), median(atac_complexity), 
            median(gex_genes_count), count = n())

pMito <- ggplot(all_df_qc, aes(x = assay, y =pct_mito*100, color = condition)) +
  geom_boxplot(outlier.shape = NA) + 
  pretty_plot(fontsize = 7) + L_border() + scale_y_continuous(limits = c(0,60)) +
  labs(x = "Condition", y = "%mtDNA (ATAC)", color = "")

pATAC <- ggplot(all_df_qc, aes(x = assay, y = log10(atac_complexity), color = condition)) +
  geom_boxplot(outlier.shape = NA) + 
  pretty_plot(fontsize = 7) + L_border() + scale_y_continuous(limits = c(2, 5.5)) +
  labs(x = "Condition", y = "log10 ATAC Complexity", color = "")

pRNA <- ggplot(all_df_qc, aes(x = assay, y = (gex_genes_count), color = condition)) +
  geom_boxplot(outlier.shape = NA) + 
  pretty_plot(fontsize = 7) + L_border() + scale_y_log10(limits = c(10^2, 10^4), breaks = c(100, 500, 1000, 5000, 10000) ) +
  labs(x = "Condition", y = "log10 genes detected (RNA)", color = "")

cowplot::ggsave2(cowplot::plot_grid(
  pATAC, pMito, pRNA, nrow = 1
), file = "../plots/all_QC_split.pdf", width = 9, height = 1.5)

all_df_qc %>% group_by(assay, condition) %>%
  summarize(mito = median(pct_mito*100), median(atac_complexity), 
            median(gex_genes_count), count = n())
