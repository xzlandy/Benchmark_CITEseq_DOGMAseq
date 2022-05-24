library(dplyr)
library(data.table)
library(stringr)
library(BuenColors)

setwd('~/CITE-seq/Duerr/DOGMA-seq/DIG_CITE_rerun_1/code/')

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
load('../data/CITE_data_wt_intron.RData')
CITE <- data@meta.data
CITE <- CITE %>% mutate(assay = 'CITE') %>%
  mutate(gex_genes_count = nFeature_RNA) %>%
  mutate(condition = condition)
CITE <- CITE[,c("gex_genes_count","assay", "condition")]
rm(data)
DIG <- DIG[,c("gex_genes_count","assay", "condition")]

all_df_qc <- rbind(DIG, CITE)
all_df_qc$sample <- str_split_fixed(all_df_qc$condition, '_Act', 2)[,1]

# pMito <- ggplot(all_df_qc, aes(x = assay, y =pct_mito*100, color = assay)) +
#   geom_boxplot(outlier.shape = NA) + 
#   pretty_plot(fontsize = 7) + L_border() + scale_y_continuous(limits = c(0,60)) +
#   scale_color_manual(values = c("dodgerblue3", "firebrick", "darkgrey")) +
#   labs(x = "Condition", y = "%mtDNA (ATAC)", color = "") +
#   theme(legend.position = "none")

# pATAC <- ggplot(all_df_qc, aes(x = assay, y = log10(atac_complexity), color = assay)) +
#   geom_boxplot(outlier.shape = NA) + 
#   pretty_plot(fontsize = 7) + L_border() + scale_y_continuous(limits = c(2, 5.5)) +
#   scale_color_manual(values = c("dodgerblue3", "firebrick", "darkgrey")) +
#   labs(x = "Condition", y = "log10 ATAC Complexity", color = "") +
#   theme(legend.position = "none")

pRNA <- ggplot(all_df_qc, aes(x = assay, y = (gex_genes_count), color = assay)) +
  geom_boxplot(outlier.shape = NA) + 
  pretty_plot(fontsize = 7) + L_border() + scale_y_log10(limits = c(10^2, 10^4), breaks = c(100, 500, 1000, 5000, 10000) ) +
  scale_color_manual(values = c("dodgerblue3", "firebrick", "darkgrey")) +
  labs(x = "Condition", y = "# of genes detected (RNA)", color = "") +
  theme(legend.position = "none")

cowplot::ggsave2(cowplot::plot_grid(
  pRNA, nrow = 1
), file = "../plots/all_QC_wt_intron.pdf", width = 1, height = 1.5)

all_df_qc %>% group_by(assay) %>%
  summarize(median(gex_genes_count), count = n())

# # A tibble: 2 × 3
# assay `median(gex_genes_count)` count
# <chr>                     <dbl> <int>
# 1 CITE                       3401 13206
# 2 DIG                        2106 14172

# pMito <- ggplot(all_df_qc, aes(x = assay, y =pct_mito*100, color = condition)) +
#   geom_boxplot(outlier.shape = NA) + 
#   pretty_plot(fontsize = 7) + L_border() + scale_y_continuous(limits = c(0,60)) +
#   labs(x = "Condition", y = "%mtDNA (ATAC)", color = "")

# pATAC <- ggplot(all_df_qc, aes(x = assay, y = log10(atac_complexity), color = condition)) +
#   geom_boxplot(outlier.shape = NA) + 
#   pretty_plot(fontsize = 7) + L_border() + scale_y_continuous(limits = c(2, 5.5)) +
#   labs(x = "Condition", y = "log10 ATAC Complexity", color = "")

pRNA <- ggplot(all_df_qc, aes(x = assay, y = (gex_genes_count), color = condition)) +
  geom_boxplot(outlier.shape = NA) + 
  pretty_plot(fontsize = 7) + L_border() + scale_y_log10(limits = c(10^2, 10^4), breaks = c(100, 500, 1000, 5000, 10000) ) +
  labs(x = "Condition", y = "log10 genes detected (RNA)", color = "")

cowplot::ggsave2(cowplot::plot_grid(
  pRNA, nrow = 1
), file = "../plots/all_QC_split_wt_intron.pdf", width = 3, height = 1.5)

# all_df_qc %>% group_by(assay, condition) %>%
#   summarize(mito = median(pct_mito*100), median(atac_complexity), 
#             median(gex_genes_count), count = n())

pRNA <- ggplot(all_df_qc, aes(x = assay, y = (gex_genes_count), color = sample)) +
  geom_boxplot(outlier.shape = NA) + 
  pretty_plot(fontsize = 7) + L_border() + scale_y_log10(limits = c(10^2, 10^4), breaks = c(100, 500, 1000, 5000, 10000) ) +
  scale_color_manual(values = c("dodgerblue3", "firebrick", "darkgrey")) +
  labs(x = "Condition", y = "log10 genes detected (RNA)", color = "") +
  theme(legend.position = "bottom", legend.direction = 'horizontal')

cowplot::ggsave2(cowplot::plot_grid(
  pRNA, nrow = 1
), file = "../plots/all_QC_split_sample_wt_intron.pdf", width = 2, height = 2)

all_df_qc %>% group_by(assay, sample) %>%
  summarize(median(gex_genes_count), count = n())

# # A tibble: 4 × 4
# # Groups:   assay [2]
# assay sample   `median(gex_genes_count)` count
# <chr> <chr>                        <int> <int>
# 1 CITE  SB775372                      3448  7171
# 2 CITE  SB775393                      3361  6035
# 3 DIG   SB775372                      2152  7943
# 4 DIG   SB775393                      2062  6229