library(Seurat)
library(BuenColors)
library(dplyr)

setwd('~/RWorkSpace/CITE-seq/Duerr/DOGMA-seq/DIG_LLL/code/')

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
load("../data/LLL_data.RData")
lll <- data
rm(data)

all_mtRNA_df <- rbind(dig@meta.data, lll@meta.data)
all_mtRNA_df$pct_mt <- all_mtRNA_df$percent.mt
all_mtRNA_df$what <- all_mtRNA_df$orig.ident

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

pMitoRNA <- ggplot(all_mtRNA_df, aes(x = what, y = (pct_mt), color = condition)) +
  geom_boxplot(outlier.shape = NA,  position = position_dodge(preserve = "single")) + 
  pretty_plot(fontsize = 7) + L_border() + scale_y_continuous(limits = c(0,50)) +
  labs(x = "Condition", y = "%UMIs from mtRNA", color = "")

pMitoRNA


cowplot::ggsave2(cowplot::plot_grid(
  pMitoRNA
), file = "../plots/mitoRNApct_QC_split.pdf", width = 3, height = 1.5)

all_mtRNA_df %>% group_by(condition, what) %>% summarize(median(pct_mt))
