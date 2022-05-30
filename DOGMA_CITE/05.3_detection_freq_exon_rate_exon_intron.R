library(Seurat)
library(BuenColors)
library(data.table)
library(dplyr)

setwd('~/RWorkSpace/CITE-seq/Duerr/DOGMA-seq/DIG_CITE_rerun_1/code/')

load("../data/CITE_data.RData")
cite <- data
rm(data)

counts <- Read10X_h5('../../../Duerr_20210610_CITEseq_1_include_introns/outs/raw_feature_bc_matrix.h5')
rna_counts <- counts$`Gene Expression`

cite_exon_intron <- CreateSeuratObject(
  counts = rna_counts,
  assay = "RNA"
)

cite_exon_intron <- cite_exon_intron[,colnames(cite)]

CITE <- cite@assays$RNA@counts
CITE_exon_intron <- cite_exon_intron@assays$RNA@counts
CITE_exon_intron <- CITE_exon_intron[rownames(CITE),colnames(CITE)]

CITE_exon <- CITE
CITE_intron <- CITE_exon_intron - CITE_exon

save(CITE_exon, CITE_intron, file = '../output/CITE_exon_intron.RData')

exon <- rowSums(CITE_exon)
exon_intron <- rowSums(CITE_exon_intron)
exon[exon > exon_intron] <- exon_intron[exon > exon_intron]
cite_exon_rate <- exon/exon_intron

save(cite_exon_rate, file = '../output/CITE_exon_rate.RData')

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
DIG_intron <- DIG - DIG_exon

save(DIG_exon, DIG_intron, file = '../output/DIG_exon_intron.RData')

load('../output/exon_rate.RData')
load('../output/CITE_exon_rate.RData')
load('../output/DIG_exon_intron.RData')
load('../output/CITE_exon_intron.RData')

DIG <- apply(DIG_exon + DIG_intron, 1, function(x){sum(ifelse(x > 0, 1, 0))/length(x)})
CITE <- apply(CITE_exon + CITE_intron, 1, function(x){sum(ifelse(x > 0, 1, 0))/length(x)})
table(names(DIG) == names(CITE))

df <- as.data.frame(cbind(DIG, CITE))
df$Type <- 'Other Genes'
df$Type[grep('^MT-', rownames(df))] <- 'MT- Genes'
df$Type[grep('^RPL|^RPS', rownames(df))] <- 'RPL/RPS Genes'
df$Type <- factor(df$Type, levels = c('Other Genes', 'RPL/RPS Genes', 'MT- Genes'))
table(rownames(df) == names(exon_rate))
df$DIG_Exon <- exon_rate
df$DIG_Exon_bi <- ifelse(df$DIG_Exon > 0.5, 'Exon-dominated Genes', 'Intron-dominated Genes')
table(rownames(df) == names(cite_exon_rate))
df$CITE_Exon <- cite_exon_rate
df$CITE_Exon_bi <- ifelse(df$CITE_Exon > 0.5, 'Exon-dominated Genes', 'Intron-dominated Genes')

p1 <- ggplot(df, aes(x = CITE, y = DIG, color = Type)) +
  geom_point() +
  pretty_plot(fontsize = 10) + L_border() + 
  scale_color_manual(values = c("darkgrey", "dodgerblue3", "firebrick")) +
  labs(x = "CITE Fraction > 0", y = "DIG Fraction > 0", color = "") +
  theme(legend.position = "bottom", legend.direction = 'horizontal')
cowplot::ggsave2(p1, file = "../plots/RNA_detection_freq_exon_intron.pdf", width = 4, height = 4)

p1 <- ggplot(df, aes(x = CITE, y = DIG, color = DIG_Exon)) +
  geom_point() +
  pretty_plot(fontsize = 10) + L_border() + 
  scale_color_gradient(low = "dodgerblue3", high = "firebrick") +
  labs(x = "CITE Fraction > 0", y = "DIG Fraction > 0", color = "Proportion of Exonic UMIs") +
  theme(legend.position = "bottom", legend.direction = 'horizontal')
cowplot::ggsave2(p1, file = "../plots/RNA_detection_freq_DIG_exon_rate_exon_intron.pdf", width = 4, height = 4)

p1 <- ggplot(df, aes(x = CITE, y = DIG, color = DIG_Exon_bi)) +
  geom_point() +
  pretty_plot(fontsize = 10) + L_border() + 
  scale_color_manual(values = c("dodgerblue3", "firebrick")) +
  labs(x = "CITE Fraction > 0", y = "DIG Fraction > 0", color = "") +
  theme(legend.position = "bottom", legend.direction = 'horizontal', legend.text = element_text(size = 7))
cowplot::ggsave2(p1, file = "../plots/RNA_detection_freq_DIG_exon_rate_bi_exon_intron.pdf", width = 4, height = 4)

p1 <- ggplot(df, aes(x = CITE, y = DIG, color = CITE_Exon)) +
  geom_point() +
  pretty_plot(fontsize = 10) + L_border() + 
  scale_color_gradient(low = "dodgerblue3", high = "firebrick") +
  labs(x = "CITE Fraction > 0", y = "DIG Fraction > 0", color = "Proportion of Exonic UMIs") +
  theme(legend.position = "bottom", legend.direction = 'horizontal')
cowplot::ggsave2(p1, file = "../plots/RNA_detection_freq_CITE_exon_rate_exon_intron.pdf", width = 4, height = 4)

p1 <- ggplot(df, aes(x = CITE, y = DIG, color = CITE_Exon_bi)) +
  geom_point() +
  pretty_plot(fontsize = 10) + L_border() + 
  scale_color_manual(values = c("dodgerblue3", "firebrick")) +
  labs(x = "CITE Fraction > 0", y = "DIG Fraction > 0", color = "") +
  theme(legend.position = "bottom", legend.direction = 'horizontal', legend.text = element_text(size = 7))
cowplot::ggsave2(p1, file = "../plots/RNA_detection_freq_CITE_exon_rate_bi_exon_intron.pdf", width = 4, height = 4)

#################
DIG <- apply(DIG_intron, 1, function(x){sum(ifelse(x > 0, 1, 0))/length(x)})
CITE <- apply(CITE_intron, 1, function(x){sum(ifelse(x > 0, 1, 0))/length(x)})
table(names(DIG) == names(CITE))

df <- as.data.frame(cbind(DIG, CITE))
df$Type <- 'Other Genes'
df$Type[grep('^MT-', rownames(df))] <- 'MT- Genes'
df$Type[grep('^RPL|^RPS', rownames(df))] <- 'RPL/RPS Genes'
df$Type <- factor(df$Type, levels = c('Other Genes', 'RPL/RPS Genes', 'MT- Genes'))
table(rownames(df) == names(exon_rate))
df$DIG_Exon <- exon_rate
df$DIG_Exon_bi <- ifelse(df$DIG_Exon > 0.5, 'Exon-dominated Genes', 'Intron-dominated Genes')
table(rownames(df) == names(cite_exon_rate))
df$CITE_Exon <- cite_exon_rate
df$CITE_Exon_bi <- ifelse(df$CITE_Exon > 0.5, 'Exon-dominated Genes', 'Intron-dominated Genes')

p1 <- ggplot(df, aes(x = CITE, y = DIG, color = Type)) +
  geom_point() +
  pretty_plot(fontsize = 10) + L_border() + 
  scale_color_manual(values = c("darkgrey", "dodgerblue3", "firebrick")) +
  labs(x = "CITE Fraction > 0", y = "DIG Fraction > 0", color = "") +
  theme(legend.position = "bottom", legend.direction = 'horizontal')
cowplot::ggsave2(p1, file = "../plots/RNA_detection_freq_intron_only.pdf", width = 4, height = 4)

p1 <- ggplot(df, aes(x = CITE, y = DIG, color = DIG_Exon)) +
  geom_point() +
  pretty_plot(fontsize = 10) + L_border() + 
  scale_color_gradient(low = "dodgerblue3", high = "firebrick") +
  labs(x = "CITE Fraction > 0", y = "DIG Fraction > 0", color = "Proportion of Exonic UMIs") +
  theme(legend.position = "bottom", legend.direction = 'horizontal')
cowplot::ggsave2(p1, file = "../plots/RNA_detection_freq_DIG_exon_rate_intron_only.pdf", width = 4, height = 4)

p1 <- ggplot(df, aes(x = CITE, y = DIG, color = DIG_Exon_bi)) +
  geom_point() +
  pretty_plot(fontsize = 10) + L_border() + 
  scale_color_manual(values = c("dodgerblue3", "firebrick")) +
  labs(x = "CITE Fraction > 0", y = "DIG Fraction > 0", color = "") +
  theme(legend.position = "bottom", legend.direction = 'horizontal', legend.text = element_text(size = 7))
cowplot::ggsave2(p1, file = "../plots/RNA_detection_freq_DIG_exon_rate_bi_intron_only.pdf", width = 4, height = 4)

p1 <- ggplot(df, aes(x = CITE, y = DIG, color = CITE_Exon)) +
  geom_point() +
  pretty_plot(fontsize = 10) + L_border() + 
  scale_color_gradient(low = "dodgerblue3", high = "firebrick") +
  labs(x = "CITE Fraction > 0", y = "DIG Fraction > 0", color = "Proportion of Exonic UMIs") +
  theme(legend.position = "bottom", legend.direction = 'horizontal')
cowplot::ggsave2(p1, file = "../plots/RNA_detection_freq_CITE_exon_rate_intron_only.pdf", width = 4, height = 4)

p1 <- ggplot(df, aes(x = CITE, y = DIG, color = CITE_Exon_bi)) +
  geom_point() +
  pretty_plot(fontsize = 10) + L_border() + 
  scale_color_manual(values = c("dodgerblue3", "firebrick")) +
  labs(x = "CITE Fraction > 0", y = "DIG Fraction > 0", color = "") +
  theme(legend.position = "bottom", legend.direction = 'horizontal', legend.text = element_text(size = 7))
cowplot::ggsave2(p1, file = "../plots/RNA_detection_freq_CITE_exon_rate_bi_intron_only.pdf", width = 4, height = 4)

table(names(exon_rate) == names(cite_exon_rate))

df <- data.frame(DIG = exon_rate, CITE = cite_exon_rate)
df$Type <- 'Other Genes'
df$Type[grep('^MT-', rownames(df))] <- 'MT- Genes'
df$Type[grep('^RPL|^RPS', rownames(df))] <- 'RPL/RPS Genes'
df$Type <- factor(df$Type, levels = c('Other Genes', 'RPL/RPS Genes', 'MT- Genes'))
p1 <- ggplot(df, aes(x = CITE, y = DIG, color = Type))+
  geom_point()+
  pretty_plot(fontsize = 10) + L_border() + 
  scale_color_manual(values = c("darkgrey", "dodgerblue3", "firebrick")) +
  labs(x = "Proportion of Exonic UMIs (CITE)", y = "Proportion of Exonic UMIs (DIG)", color = "") +
  theme(legend.position = "bottom", legend.direction = 'horizontal')
cowplot::ggsave2(p1, file = "../plots/RNA_detection_freq_exon_intron_cor.pdf", width = 4, height = 4)

