library(Seurat)
library(BuenColors)
library(data.table)
library(dplyr)
library(stringr)

setwd('~/CITE-seq/Duerr/DOGMA-seq/DIG_CITE_rerun_1/code/')

# the 10x hdf5 file contains both data types. 
import_pctMT_exon_RNAseq <- function(file,qc_file,  condition, what){
  inputdata.10x <- Read10X_h5(file)
  pbmc <- CreateSeuratObject(counts = inputdata.10x$`Gene Expression`)
  initial_df <- merge(data.frame(
    barcode = colnames(pbmc),
    what, 
    condition, 
    count_mt = PercentageFeatureSet(pbmc, pattern = "^MT-")[[1]] * pbmc@meta.data$nCount_RNA /100),
    fread(qc_file) %>% filter(is_cell == 1), by = "barcode")
  
  initial_df[,c("count_mt", "gex_exonic_umis", "gex_intronic_umis", "what", "condition")] %>%
    mutate(exon_no_mt = gex_exonic_umis-count_mt) %>% 
    mutate(pct_exon = gex_exonic_umis/(gex_exonic_umis + gex_intronic_umis)*100) %>%
    mutate(pct_exon_no_mt = exon_no_mt/(exon_no_mt + gex_intronic_umis)*100) 
  
}

load("../data/DIG_data.RData")
dig <- data
rm(data)

all_mtRNA_df <- dig@meta.data
all_mtRNA_df <- all_mtRNA_df %>% mutate(count_mt = percent.mt * nCount_RNA /100) %>%
  mutate(exon_no_mt = gex_exonic_umis-count_mt) %>% 
  mutate(pct_exon = gex_exonic_umis/(gex_exonic_umis + gex_intronic_umis)*100) %>%
  mutate(pct_exon_no_mt = exon_no_mt/(exon_no_mt + gex_intronic_umis)*100)
all_mtRNA_df$what <- 'DIG'
all_mtRNA_df$sample <- str_split_fixed(all_mtRNA_df$condition, '_Act', 2)[,1]

propRNA1 <- ggplot(all_mtRNA_df, aes(x = what, y = (pct_exon), color = what)) +
  geom_boxplot(outlier.shape = NA,  position = position_dodge(preserve = "single")) + 
  pretty_plot(fontsize = 7) + L_border() + scale_y_continuous(limits = c(0,100)) +
  scale_color_manual(values = c("dodgerblue3", "firebrick", "darkgrey")) +
  labs(x = "Condition", y = "%UMIs mapped to exons", color = "") +
  theme(legend.position = "none")

propRNA2 <- ggplot(all_mtRNA_df, aes(x = what, y = (pct_exon_no_mt), color = what)) +
  geom_boxplot(outlier.shape = NA,  position = position_dodge(preserve = "single")) + 
  pretty_plot(fontsize = 7) + L_border() + scale_y_continuous(limits = c(0,100)) +
  scale_color_manual(values = c("dodgerblue3", "firebrick", "darkgrey")) +
  labs(x = "Condition", y = "%UMIs mappted exons - no mito", color = "") +
  theme(legend.position = "none")


cowplot::ggsave2(cowplot::plot_grid(
  propRNA1, propRNA2, nrow =1
), file = "../plots/RNA_features_pct_QC_DIG_only.pdf", width = 2, height = 1.5)

all_mtRNA_df %>% group_by(what) %>% summarize(median(pct_exon), median(pct_exon_no_mt))

# # A tibble: 1 × 3
# what  `median(pct_exon)` `median(pct_exon_no_mt)`
# <chr>              <dbl>                    <dbl>
# 1 DIG                 42.2                     37.3

propRNA1 <- ggplot(all_mtRNA_df, aes(x = what, y = (pct_exon), color = sample)) +
  geom_boxplot(outlier.shape = NA,  position = position_dodge(preserve = "single")) + 
  pretty_plot(fontsize = 7) + L_border() + scale_y_continuous(limits = c(0,100)) +
  scale_color_manual(values = c("dodgerblue3", "firebrick", "darkgrey")) +
  labs(x = "Condition", y = "%UMIs mapped to exons", color = "") +
  theme(legend.position = "none")

propRNA2 <- ggplot(all_mtRNA_df, aes(x = what, y = (pct_exon_no_mt), color = sample)) +
  geom_boxplot(outlier.shape = NA,  position = position_dodge(preserve = "single")) + 
  pretty_plot(fontsize = 7) + L_border() + scale_y_continuous(limits = c(0,100)) +
  scale_color_manual(values = c("dodgerblue3", "firebrick", "darkgrey")) +
  labs(x = "Condition", y = "%UMIs mappted exons - no mito", color = "") +
  theme(legend.position = "none")


cowplot::ggsave2(cowplot::plot_grid(
  propRNA1, propRNA2, nrow =1
), file = "../plots/RNA_features_pct_QC_split_sample_DIG_only.pdf", width = 2, height = 1.5)

all_mtRNA_df %>% group_by(sample, what) %>% summarize(median(pct_exon), median(pct_exon_no_mt))

# # A tibble: 2 × 4
# # Groups:   sample [2]
# sample   what  `median(pct_exon)` `median(pct_exon_no_mt)`
# <chr>    <chr>              <dbl>                    <dbl>
# 1 SB775372 DIG                 41.8                     37.3
# 2 SB775393 DIG                 42.7                     37.4