library(Seurat)
library(BuenColors)
library(data.table)
library(dplyr)

setwd('~/CITE-seq/Duerr/DOGMA-seq/DIG_CITE_rerun_1/code/')

load("../data/DIG_data.RData")
dig <- data

counts <- Read10X_h5('../../../Duerr_20210610_DOGMAseq_1_gex_exclude_introns/outs/raw_feature_bc_matrix.h5')
rna_counts <- counts$`Gene Expression`
atac_counts <- counts$Peaks
# create a Seurat object containing the RNA data
metadata <- read.csv(
  file = "../../../Duerr_20210610_DOGMAseq_1_gex_exclude_introns/outs/per_barcode_metrics.csv",
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

exon <- rowSums(DIG_exon)
exon_intron <- rowSums(DIG)
exon[exon > exon_intron] <- exon_intron[exon > exon_intron]
exon_rate <- exon/exon_intron

save(exon_rate, file = '../output/exon_rate.RData')

load('../output/exon_rate.RData')
load("../data/DIG_data.RData")
dig <- data
load("../data/CITE_data.RData")
cite <- data
rm(data)

DIG <- apply(dig@assays$RNA@counts, 1, function(x){sum(ifelse(x > 0, 1, 0))/length(x)})
CITE <- apply(cite@assays$RNA@counts, 1, function(x){sum(ifelse(x > 0, 1, 0))/length(x)})
table(names(DIG) == names(CITE))

df <- as.data.frame(cbind(DIG, CITE))
df$Type <- 'Other Genes'
df$Type[grep('^MT-', rownames(df))] <- 'MT- Genes'
df$Type[grep('^RPL|^RPS', rownames(df))] <- 'RPL/RPS Genes'
df$Type <- factor(df$Type, levels = c('Other Genes', 'RPL/RPS Genes', 'MT- Genes'))
table(rownames(df) == names(exon_rate))
df$Exon <- exon_rate
df$Exon_bi <- ifelse(df$Exon > 0.5, 'Exon-dominated Genes', 'Intron-dominated Genes')

p1 <- ggplot(df, aes(x = CITE, y = DIG, color = Exon)) +
  geom_point() +
  pretty_plot(fontsize = 10) + L_border() + 
  scale_color_gradient(low = "dodgerblue3", high = "firebrick") +
  labs(x = "CITE Fraction > 0", y = "DIG Fraction > 0", color = "Proportion of Exonic UMIs") +
  theme(legend.position = "bottom", legend.direction = 'horizontal')
cowplot::ggsave2(p1, file = "../plots/RNA_detection_freq_exon_rate.pdf", width = 4, height = 4)

p1 <- ggplot(df, aes(x = CITE, y = DIG, color = Exon_bi)) +
  geom_point() +
  pretty_plot(fontsize = 10) + L_border() + 
  scale_color_manual(values = c("dodgerblue3", "firebrick")) +
  labs(x = "CITE Fraction > 0", y = "DIG Fraction > 0", color = "") +
  theme(legend.position = "bottom", legend.direction = 'horizontal', legend.text = element_text(size = 7))
cowplot::ggsave2(p1, file = "../plots/RNA_detection_freq_exon_rate_bi.pdf", width = 4, height = 4)
