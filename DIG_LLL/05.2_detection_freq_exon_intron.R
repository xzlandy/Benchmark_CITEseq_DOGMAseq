library(Seurat)
library(BuenColors)
library(data.table)
library(dplyr)
library(ggsci)

setwd('~/RWorkSpace/CITE-seq/Duerr/DOGMA-seq/DIG_LLL/code/')

load("../data/DIG_data.RData")
dig <- data

counts <- Read10X_h5('../../../Duerr_20210419_DOGMAseq_DIG_gex_exclude_introns/outs/raw_feature_bc_matrix.h5')
rna_counts <- counts$`Gene Expression`
atac_counts <- counts$Peaks
# create a Seurat object containing the RNA data
metadata <- read.csv(
  file = "../../../Duerr_20210419_DOGMAseq_DIG_gex_exclude_introns/outs/per_barcode_metrics.csv",
  header = TRUE,
  row.names = 1
)

dig_exon <- CreateSeuratObject(
  counts = rna_counts,
  assay = "RNA",
  meta.data = metadata
)

dig_exon <- dig_exon[,colnames(dig)]

DIG <- dig@assays$RNA@counts
DIG_exon <- dig_exon@assays$RNA@counts
DIG_exon <- DIG_exon[rownames(DIG),colnames(DIG)]
DIG_intron <- DIG - DIG_exon

save(DIG_exon, DIG_intron, file = '../output/DIG_exon_intron.RData')

load("../data/LLL_data.RData")
lll <- data

counts <- Read10X_h5('../../../Duerr_20210419_DOGMAseq_PFA_LLL_gex_exclude_introns/outs/raw_feature_bc_matrix.h5')
rna_counts <- counts$`Gene Expression`
atac_counts <- counts$Peaks
# create a Seurat object containing the RNA data
metadata <- read.csv(
  file = "../../../Duerr_20210419_DOGMAseq_PFA_LLL_gex_exclude_introns/outs/per_barcode_metrics.csv",
  header = TRUE,
  row.names = 1
)

lll_exon <- CreateSeuratObject(
  counts = rna_counts,
  assay = "RNA",
  meta.data = metadata
)

lll_exon <- lll_exon[,colnames(lll)]

LLL <- lll@assays$RNA@counts
LLL_exon <- lll_exon@assays$RNA@counts
LLL_exon <- LLL_exon[rownames(LLL),colnames(LLL)]
LLL_intron <- LLL - LLL_exon

save(LLL_exon, LLL_intron, file = '../output/LLL_exon_intron.RData')

load('../output/DIG_exon_intron.RData')
load('../output/LLL_exon_intron.RData')

DIG <- apply(DIG_exon, 1, function(x){sum(ifelse(x > 0, 1, 0))/length(x)})
LLL <- apply(LLL_exon, 1, function(x){sum(ifelse(x > 0, 1, 0))/length(x)})
table(names(DIG) == names(LLL))

df <- as.data.frame(cbind(DIG, LLL))
df$Type <- 'Other Genes'
df$Type[grep('^MT-', rownames(df))] <- 'MT- Genes'
df$Type[grep('^RPL|^RPS', rownames(df))] <- 'RPL/RPS Genes'
df$Type <- factor(df$Type, levels = c('Other Genes', 'RPL/RPS Genes', 'MT- Genes'))

p1 <- ggplot(df, aes(x = LLL, y = DIG, color = Type)) +
  geom_point() +
  pretty_plot(fontsize = 10) + L_border() + 
  scale_color_manual(values = c("darkgrey", "dodgerblue3", "firebrick")) +
  labs(x = "LLL Fraction > 0", y = "DIG Fraction > 0", color = "") +
  theme(legend.position = "bottom", legend.direction = 'horizontal')
cowplot::ggsave2(p1, file = "../plots/RNA_detection_freq_exon_only.pdf", width = 4, height = 4)

#################
DIG <- apply(DIG_intron, 1, function(x){sum(ifelse(x > 0, 1, 0))/length(x)})
LLL <- apply(LLL_intron, 1, function(x){sum(ifelse(x > 0, 1, 0))/length(x)})
table(names(DIG) == names(LLL))

df <- as.data.frame(cbind(DIG, LLL))
df$Type <- 'Other Genes'
df$Type[grep('^MT-', rownames(df))] <- 'MT- Genes'
df$Type[grep('^RPL|^RPS', rownames(df))] <- 'RPL/RPS Genes'
df$Type <- factor(df$Type, levels = c('Other Genes', 'RPL/RPS Genes', 'MT- Genes'))

p1 <- ggplot(df, aes(x = LLL, y = DIG, color = Type)) +
  geom_point() +
  pretty_plot(fontsize = 10) + L_border() + 
  scale_color_manual(values = c("darkgrey", "dodgerblue3", "firebrick")) +
  labs(x = "LLL Fraction > 0", y = "DIG Fraction > 0", color = "") +
  theme(legend.position = "bottom", legend.direction = 'horizontal')
cowplot::ggsave2(p1, file = "../plots/RNA_detection_freq_intron_only.pdf", width = 4, height = 4)


