library(Seurat)
library(Signac)
library(dplyr)
library(data.table)
library(stringr)
library(harmony)
library(patchwork)

setwd('~/RWorkSpace/CITE-seq/Duerr/DOGMA-seq/DIG_CITE_rerun_1/code/')

load('../output/cite_harmony_independent.RData')

cite$wnn_clusters <- factor(cite$wnn_clusters, levels = 0:10)
Idents(cite) <- "wnn_clusters"

reference <- readRDS('~/CITE-seq/Seurat/PBMC.RDS')

DefaultAssay(cite) <- 'SCT'

anchors <- FindTransferAnchors(
  reference = reference,
  query = cite,
  normalization.method = "SCT",
  reference.reduction = "spca",
  dims = 1:50
)

cite <- MapQuery(
  anchorset = anchors,
  query = cite,
  reference = reference,
  refdata = list(
    celltype.l1 = "celltype.l1",
    celltype.l2 = "celltype.l2",
    predicted_ADT = "ADT"
  ),
  reference.reduction = "spca", 
  reduction.model = "wnn.umap"
)

p1 <- DimPlot(cite, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 4)
p2 <- DimPlot(cite, reduction = 'wnn.umap', label = TRUE, repel = TRUE, group.by = "predicted.celltype.l1", label.size = 4) 
p3 <- DimPlot(cite, reduction = 'wnn.umap', label = TRUE, repel = TRUE, group.by = "predicted.celltype.l2", label.size = 4) 
pdf('../plots/harmony/cite_predicted.pdf', width = 24, height = 8)
p1 | p2 | p3
dev.off()

p1 <- DimPlot(cite, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 4)
p2 <- DimPlot(cite, reduction = 'wnn.umap', label = TRUE, repel = TRUE, group.by = "predicted.celltype.l1", label.size = 4) 
pdf('../plots/harmony/cite_predicted_2.pdf', width = 16, height = 8)
p1 | p2
dev.off()

library(ggplot2)
data_plot <- table(cite$predicted.celltype.l1, cite$condition)
data_plot <- reshape2::melt(prop.table(data_plot, 2))
colnames(data_plot) <- c('Clusters', 'Conditions', 'Proportions')
data_plot$Clusters <- as.factor(data_plot$Clusters)
data_plot$Proportions <- data_plot$Proportions*100

p1 <- ggplot(data_plot, aes(x = Clusters, y = Proportions, col = Conditions, fill = Conditions))+
  geom_col(position = 'dodge')+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0))
pdf('../plots/harmony/cite_proportion_predicted_slide.pdf', width = 8, height = 4)
p1
dev.off()

data_plot <- table(cite$predicted.celltype.l2, cite$condition)
data_plot <- reshape2::melt(prop.table(data_plot, 2))
colnames(data_plot) <- c('Clusters', 'Conditions', 'Proportions')
data_plot$Clusters <- as.factor(data_plot$Clusters)
data_plot$Proportions <- data_plot$Proportions*100

p1 <- ggplot(data_plot, aes(x = Clusters, y = Proportions, col = Conditions, fill = Conditions))+
  geom_col(position = 'dodge')+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0))
pdf('../plots/harmony/cite_proportion_predicted_l2_slide.pdf', width = 8, height = 4)
p1
dev.off()

write.csv(cite@meta.data, file = '../plots/harmony/cite_predicted.csv')

p1 <- FeaturePlot(cite, c('sct_IL17A'), reduction = 'wnn.umap', min.cutoff = 'q2', max.cutoff = 'q98') + ggtitle('IL17A (RNA)')
p2 <- FeaturePlot(cite, c('sct_CCR6'), reduction = 'wnn.umap', min.cutoff = 'q2', max.cutoff = 'q98') + ggtitle('CCR6 (RNA)')
p3 <- FeaturePlot(cite, c('sct_RORC'), reduction = 'wnn.umap', min.cutoff = 'q2', max.cutoff = 'q98') + ggtitle('RORC (RNA)')
p4 <- FeaturePlot(cite, c('adt_CD194-A0071'), reduction = 'wnn.umap', cols = c("lightgrey", "darkgreen"), min.cutoff = 'q2', max.cutoff = 'q98') + ggtitle('CCR4 (ADT)')
p5 <- FeaturePlot(cite, c('adt_CD196-A0143'), reduction = 'wnn.umap', cols = c("lightgrey", "darkgreen"), min.cutoff = 'q2', max.cutoff = 'q98') + ggtitle('CCR6 (ADT)')

pdf('../plots/harmony/cite_Th17.pdf', width = 12, height = 8)
wrap_plots(p1, p2, p3, p4, p5, ncol = 3)
dev.off()

p1 <- FeaturePlot(cite, c('sct_IFNG'), reduction = 'wnn.umap', min.cutoff = 'q2', max.cutoff = 'q98') + ggtitle('IFNG (RNA)')
p2 <- FeaturePlot(cite, c('sct_IL12RB2'), reduction = 'wnn.umap', min.cutoff = 'q2', max.cutoff = 'q98') + ggtitle('IL12RB2 (RNA)')
p3 <- FeaturePlot(cite, c('sct_TBX21'), reduction = 'wnn.umap', min.cutoff = 'q2', max.cutoff = 'q98') + ggtitle('TBX21 (RNA)')
p4 <- FeaturePlot(cite, c('adt_CD183-A0140'), reduction = 'wnn.umap', cols = c("lightgrey", "darkgreen"), min.cutoff = 'q2', max.cutoff = 'q98') + ggtitle('CXCR3 (ADT)')
p5 <- FeaturePlot(cite, c('adt_CD195-A0141'), reduction = 'wnn.umap', cols = c("lightgrey", "darkgreen"), min.cutoff = 'q2', max.cutoff = 'q98') + ggtitle('CCR5 (ADT)')

pdf('../plots/harmony/cite_Th1.pdf', width = 12, height = 8)
wrap_plots(p1, p2, p3, p4, p5, ncol = 3)
dev.off()

save(cite, file = '../output/cite_harmony_independent_predicted.RData')
