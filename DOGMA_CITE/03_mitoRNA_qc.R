library(Seurat)
library(BuenColors)
library(dplyr)
library(stringr)

setwd('~/RWorkSpace/CITE-seq/Duerr/DOGMA-seq/DIG_CITE_rerun_1/code/')

# the 10x hdf5 file contains both data types. 
import_pctMT_RNAseq <- function(file, condition, what){
  inputdata.10x <- Read10X_h5(file)
  pbmc <- CreateSeuratObject(counts = inputdata.10x$`Gene Expression`)
  data.frame(
    what, 
    condition, 
    pct_mt = PercentageFeatureSet(pbmc, pattern = "^MT-")[[1]])
}

load("../data/DIG_data.RData")
dig <- data
load("../data/CITE_data.RData")
cite <- data
rm(data)

dig$what <- 'DIG'
cite$what <- 'CITE'
# cite$condition <- cite$orig.ident
all_mtRNA_df <- rbind(dig@meta.data[,c('percent.mt', 'what', 'condition')], cite@meta.data[,c('percent.mt', 'what', 'condition')])
all_mtRNA_df$pct_mt <- all_mtRNA_df$percent.mt
all_mtRNA_df$sample <- str_split_fixed(all_mtRNA_df$condition, '_Act', 2)[,1]

pMitoRNA <- ggplot(all_mtRNA_df, aes(x = what, y = (pct_mt), color = what)) +
  geom_boxplot(outlier.shape = NA,  position = position_dodge(preserve = "single")) + 
  pretty_plot(fontsize = 7) + L_border() + scale_y_continuous(limits = c(0,50)) +
  scale_color_manual(values = c("dodgerblue3", "firebrick", "darkgrey")) +
  labs(x = "Condition", y = "%UMIs from mtRNA", color = "") +
  theme(legend.position = "none")

pMitoRNA


cowplot::ggsave2(cowplot::plot_grid(
  pMitoRNA
), file = "../plots/mitoRNApct_QC.pdf", width = 1, height = 1.5)

all_mtRNA_df %>% group_by(what) %>% summarize(median(pct_mt))

# # A tibble: 2 × 2
# what  `median(pct_mt)`
# <chr>            <dbl>
# 1 CITE              6.32
# 2 DIG               7.56

pMitoRNA <- ggplot(all_mtRNA_df, aes(x = what, y = (pct_mt), color = condition)) +
  geom_boxplot(outlier.shape = NA,  position = position_dodge(preserve = "single")) + 
  pretty_plot(fontsize = 7) + L_border() + scale_y_continuous(limits = c(0,50)) +
  labs(x = "Condition", y = "%UMIs from mtRNA", color = "")

pMitoRNA


cowplot::ggsave2(cowplot::plot_grid(
  pMitoRNA
), file = "../plots/mitoRNApct_QC_split.pdf", width = 3, height = 1.5)

all_mtRNA_df %>% group_by(condition, what) %>% summarize(median(pct_mt))

pMitoRNA <- ggplot(all_mtRNA_df, aes(x = what, y = (pct_mt), color = sample)) +
  geom_boxplot(outlier.shape = NA,  position = position_dodge(preserve = "single")) + 
  pretty_plot(fontsize = 7) + L_border() + scale_y_continuous(limits = c(0,50)) +
  scale_color_manual(values = c("dodgerblue3", "firebrick", "darkgrey")) +
  labs(x = "Condition", y = "%UMIs from mtRNA", color = "") +
  theme(legend.position = "bottom", legend.direction = 'horizontal')

pMitoRNA


cowplot::ggsave2(cowplot::plot_grid(
  pMitoRNA
), file = "../plots/mitoRNApct_QC_split_sample.pdf", width = 2, height = 2)

all_mtRNA_df %>% group_by(sample, what) %>% summarize(median(pct_mt))

# # A tibble: 4 × 3
# # Groups:   sample [2]
# sample   what  `median(pct_mt)`
# <chr>    <chr>            <dbl>
# 1 SB775372 CITE              5.77
# 2 SB775372 DIG               7.04
# 3 SB775393 CITE              6.90
# 4 SB775393 DIG               8.30