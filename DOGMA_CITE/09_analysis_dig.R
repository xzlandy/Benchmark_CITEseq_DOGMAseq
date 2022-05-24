library(Seurat)
library(Signac)
library(dplyr)
library(data.table)
library(stringr)
library(harmony)

setwd('~/CITE-seq/Duerr/DOGMA-seq/DIG_CITE_rerun_1/code/')

load('../output/dig_harmony_independent.RData')

dig$wnn_clusters <- factor(dig$wnn_clusters, levels = 0:14)
Idents(dig) <- "wnn_clusters"

reference <- readRDS('~/CITE-seq/Seurat/PBMC.RDS')

DefaultAssay(dig) <- 'SCT'

anchors <- FindTransferAnchors(
  reference = reference,
  query = dig,
  normalization.method = "SCT",
  reference.reduction = "spca",
  dims = 1:50
)

dig <- MapQuery(
  anchorset = anchors,
  query = dig,
  reference = reference,
  refdata = list(
    celltype.l1 = "celltype.l1",
    celltype.l2 = "celltype.l2",
    predicted_ADT = "ADT"
  ),
  reference.reduction = "spca", 
  reduction.model = "wnn.umap"
)

p1 <- DimPlot(dig, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 4)
p2 <- DimPlot(dig, reduction = 'wnn.umap', label = TRUE, repel = TRUE, group.by = "predicted.celltype.l1", label.size = 4) 
p3 <- DimPlot(dig, reduction = 'wnn.umap', label = TRUE, repel = TRUE, group.by = "predicted.celltype.l2", label.size = 4) 
pdf('../plots/harmony/dig_predicted.pdf', width = 24, height = 8)
p1 | p2 | p3
dev.off()

p1 <- DimPlot(dig, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 4)
p2 <- DimPlot(dig, reduction = 'wnn.umap', label = TRUE, repel = TRUE, group.by = "predicted.celltype.l1", label.size = 4) 
pdf('../plots/harmony/dig_predicted_2.pdf', width = 16, height = 8)
p1 | p2
dev.off()

library(ggplot2)
data_plot <- table(dig$predicted.celltype.l1, dig$condition)
data_plot <- reshape2::melt(prop.table(data_plot, 2))
colnames(data_plot) <- c('Clusters', 'Conditions', 'Proportions')
data_plot$Clusters <- as.factor(data_plot$Clusters)
data_plot$Proportions <- data_plot$Proportions*100

p1 <- ggplot(data_plot, aes(x = Clusters, y = Proportions, col = Conditions, fill = Conditions))+
  geom_col(position = 'dodge')+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0))
pdf('../plots/harmony/dig_proportion_predicted_slide.pdf', width = 8, height = 4)
p1
dev.off()

data_plot <- table(dig$predicted.celltype.l2, dig$condition)
data_plot <- reshape2::melt(prop.table(data_plot, 2))
colnames(data_plot) <- c('Clusters', 'Conditions', 'Proportions')
data_plot$Clusters <- as.factor(data_plot$Clusters)
data_plot$Proportions <- data_plot$Proportions*100

p1 <- ggplot(data_plot, aes(x = Clusters, y = Proportions, col = Conditions, fill = Conditions))+
  geom_col(position = 'dodge')+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0))
pdf('../plots/harmony/dig_proportion_predicted_l2_slide.pdf', width = 8, height = 4)
p1
dev.off()

write.csv(dig@meta.data, file = '../plots/harmony/dig_predicted.csv')

library(chromVAR)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg38)

# Scan the DNA sequence of each peak for the presence of each motif, and create a Motif object
DefaultAssay(dig) <- "peaks"
pwm_set <- getMatrixSet(x = JASPAR2020, opts = list(species = 9606, all_versions = FALSE))
motif.matrix <- CreateMotifMatrix(features = granges(dig), pwm = pwm_set, genome = 'hg38', use.counts = FALSE)
motif.object <- CreateMotifObject(data = motif.matrix, pwm = pwm_set)
dig <- SetAssayData(dig, assay = 'peaks', slot = 'motifs', new.data = motif.object)

# Note that this step can take 30-60 minutes
dig <- RunChromVAR(
  object = dig,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

p1 <- FeaturePlot(dig, c('sct_IL17A'), reduction = 'wnn.umap', min.cutoff = 'q2', max.cutoff = 'q98') + ggtitle('IL17A (RNA)')
p2 <- FeaturePlot(dig, c('sct_CCR6'), reduction = 'wnn.umap', min.cutoff = 'q2', max.cutoff = 'q98') + ggtitle('CCR6 (RNA)')
p3 <- FeaturePlot(dig, c('sct_RORC'), reduction = 'wnn.umap', min.cutoff = 'q2', max.cutoff = 'q98') + ggtitle('RORC (RNA)')
p4 <- FeaturePlot(dig, c('adt_CD194-A0071'), reduction = 'wnn.umap', cols = c("lightgrey", "darkgreen"), min.cutoff = 'q2', max.cutoff = 'q98') + ggtitle('CCR4 (ADT)')
p5 <- FeaturePlot(dig, c('adt_CD196-A0143'), reduction = 'wnn.umap', cols = c("lightgrey", "darkgreen"), min.cutoff = 'q2', max.cutoff = 'q98') + ggtitle('CCR6 (ADT)')
p6 <- FeaturePlot(dig, c('chromvar_MA1151.1'), reduction = 'wnn.umap', cols = c("lightgrey", "darkred"), min.cutoff = 'q2', max.cutoff = 'q98') + ggtitle('RORC (ATAC)')

pdf('../plots/harmony/dig_Th17.pdf', width = 12, height = 8)
(p1 | p2 | p3)/ (p4 | p5 | p6)
dev.off()

p1 <- FeaturePlot(dig, c('sct_IFNG'), reduction = 'wnn.umap', min.cutoff = 'q2', max.cutoff = 'q98') + ggtitle('IFNG (RNA)')
p2 <- FeaturePlot(dig, c('sct_IL12RB2'), reduction = 'wnn.umap', min.cutoff = 'q2', max.cutoff = 'q98') + ggtitle('IL12RB2 (RNA)')
p3 <- FeaturePlot(dig, c('sct_TBX21'), reduction = 'wnn.umap', min.cutoff = 'q2', max.cutoff = 'q98') + ggtitle('TBX21 (RNA)')
p4 <- FeaturePlot(dig, c('adt_CD183-A0140'), reduction = 'wnn.umap', cols = c("lightgrey", "darkgreen"), min.cutoff = 'q2', max.cutoff = 'q98') + ggtitle('CXCR3 (ADT)')
p5 <- FeaturePlot(dig, c('adt_CD195-A0141'), reduction = 'wnn.umap', cols = c("lightgrey", "darkgreen"), min.cutoff = 'q2', max.cutoff = 'q98') + ggtitle('CCR5 (ADT)')
p6 <- FeaturePlot(dig, c('chromvar_MA0690.1'), reduction = 'wnn.umap', cols = c("lightgrey", "darkred"), min.cutoff = 'q2', max.cutoff = 'q98') + ggtitle('TBX21 (ATAC)')

pdf('../plots/harmony/dig_Th1.pdf', width = 12, height = 8)
(p1 | p2 | p3)/ (p4 | p5 | p6)
dev.off()

save(dig, file = '../output/dig_harmony_independent_predicted_chrom.RData')
