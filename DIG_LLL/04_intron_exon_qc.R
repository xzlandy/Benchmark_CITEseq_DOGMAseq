library(Seurat)
library(BuenColors)
library(data.table)
library(dplyr)

setwd('~/RWorkSpace_local/CITE-seq/Duerr/DOGMA-seq/DIG_LLL/code/')

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
load("../data/LLL_data.RData")
lll <- data
rm(data)

all_mtRNA_df <- rbind(dig@meta.data, lll@meta.data)
all_mtRNA_df <- all_mtRNA_df %>% mutate(count_mt = percent.mt * nCount_RNA /100) %>%
  mutate(exon_no_mt = gex_exonic_umis-count_mt) %>% 
  mutate(pct_exon = gex_exonic_umis/(gex_exonic_umis + gex_intronic_umis)*100) %>%
  mutate(pct_exon_no_mt = exon_no_mt/(exon_no_mt + gex_intronic_umis)*100)
all_mtRNA_df$what <- all_mtRNA_df$orig.ident

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
), file = "../plots/RNA_features_pct_QC.pdf", width = 2, height = 1.5)

all_mtRNA_df %>% group_by(what) %>% summarize(median(pct_exon), median(pct_exon_no_mt))

propRNA1 <- ggplot(all_mtRNA_df, aes(x = what, y = (pct_exon), color = condition)) +
  geom_boxplot(outlier.shape = NA,  position = position_dodge(preserve = "single")) + 
  pretty_plot(fontsize = 7) + L_border() + scale_y_continuous(limits = c(0,100)) +
  labs(x = "Condition", y = "%UMIs mapped to exons", color = "")

propRNA2 <- ggplot(all_mtRNA_df, aes(x = what, y = (pct_exon_no_mt), color = condition)) +
  geom_boxplot(outlier.shape = NA,  position = position_dodge(preserve = "single")) + 
  pretty_plot(fontsize = 7) + L_border() + scale_y_continuous(limits = c(0,100)) +
  labs(x = "Condition", y = "%UMIs mappted exons - no mito", color = "")


cowplot::ggsave2(cowplot::plot_grid(
  propRNA1, propRNA2, nrow =1
), file = "../plots/RNA_features_pct_QC_split.pdf", width = 6, height = 1.5)

all_mtRNA_df %>% group_by(condition, what) %>% summarize(median(pct_exon), median(pct_exon_no_mt))
