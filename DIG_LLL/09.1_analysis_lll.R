library(Seurat)
library(Signac)
library(dplyr)
library(data.table)
library(stringr)
library(harmony)

setwd('~/RWorkSpace/CITE-seq/Duerr/DOGMA-seq/DIG_LLL/code/')

load('../output/lll_harmony.RData')

DefaultAssay(lll) <- "RNA"
lll <- SCTransform(lll, method = "glmGamPoi")

# ADT analysis
DefaultAssay(lll) <- 'ADT'
# we will use all ADT features for dimensional reduction
# we set a dimensional reduction name to avoid overwriting the 
VariableFeatures(lll) <- rownames(lll[["ADT"]])
lll <- NormalizeData(lll, normalization.method = 'CLR', margin = 2)

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(lll) <- "peaks"
lll <- RunTFIDF(lll)
lll <- FindTopFeatures(lll, min.cutoff = 'q0')
lll <- RunSVD(lll)

lll$wnn_clusters <- factor(lll$wnn_clusters, levels = 0:12)
Idents(lll) <- "wnn_clusters"

load('~/RWorkSpace/CITE-seq/Seurat/reference.RData')

DefaultAssay(lll) <- 'SCT'

anchors <- FindTransferAnchors(
  reference = reference,
  query = lll,
  normalization.method = "SCT",
  reference.reduction = "spca",
  dims = 1:50
)

lll <- MapQuery(
  anchorset = anchors,
  query = lll,
  reference = reference,
  refdata = list(
    celltype.l1 = "celltype.l1",
    celltype.l2 = "celltype.l2",
    predicted_ADT = "ADT"
  ),
  reference.reduction = "spca", 
  reduction.model = "wnn.umap"
)

p1 <- DimPlot(lll, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 4)
p2 <- DimPlot(lll, reduction = 'wnn.umap', label = TRUE, repel = TRUE, group.by = "predicted.celltype.l1", label.size = 4) 
p3 <- DimPlot(lll, reduction = 'wnn.umap', label = TRUE, repel = TRUE, group.by = "predicted.celltype.l2", label.size = 4) 
pdf('../plots/harmony/lll_predicted.pdf', width = 24, height = 8)
p1 | p2 | p3
dev.off()

p1 <- DimPlot(lll, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 4)
p2 <- DimPlot(lll, reduction = 'wnn.umap', label = TRUE, repel = TRUE, group.by = "predicted.celltype.l1", label.size = 4) 
pdf('../plots/harmony/lll_predicted_2.pdf', width = 16, height = 8)
p1 | p2
dev.off()

library(ggplot2)
data_plot <- table(lll$predicted.celltype.l1, lll$condition)
data_plot <- reshape2::melt(prop.table(data_plot, 2))
colnames(data_plot) <- c('Clusters', 'Conditions', 'Proportions')
data_plot$Clusters <- as.factor(data_plot$Clusters)
data_plot$Proportions <- data_plot$Proportions*100

p1 <- ggplot(data_plot, aes(x = Clusters, y = Proportions, col = Conditions, fill = Conditions))+
  geom_col(position = 'dodge')+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0))
pdf('../plots/harmony/lll_proportion_predicted_slide.pdf', width = 8, height = 4)
p1
dev.off()

data_plot <- table(lll$predicted.celltype.l2, lll$condition)
data_plot <- reshape2::melt(prop.table(data_plot, 2))
colnames(data_plot) <- c('Clusters', 'Conditions', 'Proportions')
data_plot$Clusters <- as.factor(data_plot$Clusters)
data_plot$Proportions <- data_plot$Proportions*100

p1 <- ggplot(data_plot, aes(x = Clusters, y = Proportions, col = Conditions, fill = Conditions))+
  geom_col(position = 'dodge')+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0))
pdf('../plots/harmony/lll_proportion_predicted_l2_slide.pdf', width = 8, height = 4)
p1
dev.off()

write.csv(lll@meta.data, file = '../plots/harmony/lll_predicted.csv')

library(chromVAR)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg38)

# Scan the DNA sequence of each peak for the presence of each motif, and create a Motif object
DefaultAssay(lll) <- "peaks"
pwm_set <- getMatrixSet(x = JASPAR2020, opts = list(species = 9606, all_versions = FALSE))
motif.matrix <- CreateMotifMatrix(features = granges(lll), pwm = pwm_set, genome = 'hg38', use.counts = FALSE)
motif.object <- CreateMotifObject(data = motif.matrix, pwm = pwm_set)
lll <- SetAssayData(lll, assay = 'peaks', slot = 'motifs', new.data = motif.object)

# Note that this step can take 30-60 minutes
lll <- RunChromVAR(
  object = lll,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

p1 <- FeaturePlot(lll, c('sct_IL17A'), reduction = 'wnn.umap', min.cutoff = 'q2', max.cutoff = 'q98') + ggtitle('IL17A (RNA)')
p2 <- FeaturePlot(lll, c('sct_CCR6'), reduction = 'wnn.umap', min.cutoff = 'q2', max.cutoff = 'q98') + ggtitle('CCR6 (RNA)')
p3 <- FeaturePlot(lll, c('sct_RORC'), reduction = 'wnn.umap', min.cutoff = 'q2', max.cutoff = 'q98') + ggtitle('RORC (RNA)')
p4 <- FeaturePlot(lll, c('adt_CD194-A0071'), reduction = 'wnn.umap', cols = c("lightgrey", "darkgreen"), min.cutoff = 'q2', max.cutoff = 'q98') + ggtitle('CCR4 (ADT)')
p5 <- FeaturePlot(lll, c('adt_CD196-A0143'), reduction = 'wnn.umap', cols = c("lightgrey", "darkgreen"), min.cutoff = 'q2', max.cutoff = 'q98') + ggtitle('CCR6 (ADT)')
p6 <- FeaturePlot(lll, c('chromvar_MA1151.1'), reduction = 'wnn.umap', cols = c("lightgrey", "darkred"), min.cutoff = 'q2', max.cutoff = 'q98') + ggtitle('RORC (ATAC)')

pdf('../plots/harmony/lll_Th17.pdf', width = 12, height = 8)
(p1 | p2 | p3)/ (p4 | p5 | p6)
dev.off()

p1 <- FeaturePlot(lll, c('sct_IFNG'), reduction = 'wnn.umap', min.cutoff = 'q2', max.cutoff = 'q98') + ggtitle('IFNG (RNA)')
p2 <- FeaturePlot(lll, c('sct_IL12RB2'), reduction = 'wnn.umap', min.cutoff = 'q2', max.cutoff = 'q98') + ggtitle('IL12RB2 (RNA)')
p3 <- FeaturePlot(lll, c('sct_TBX21'), reduction = 'wnn.umap', min.cutoff = 'q2', max.cutoff = 'q98') + ggtitle('TBX21 (RNA)')
p4 <- FeaturePlot(lll, c('adt_CD183-A0140'), reduction = 'wnn.umap', cols = c("lightgrey", "darkgreen"), min.cutoff = 'q2', max.cutoff = 'q98') + ggtitle('CXCR3 (ADT)')
p5 <- FeaturePlot(lll, c('adt_CD195-A0141'), reduction = 'wnn.umap', cols = c("lightgrey", "darkgreen"), min.cutoff = 'q2', max.cutoff = 'q98') + ggtitle('CCR5 (ADT)')
p6 <- FeaturePlot(lll, c('chromvar_MA0690.1'), reduction = 'wnn.umap', cols = c("lightgrey", "darkred"), min.cutoff = 'q2', max.cutoff = 'q98') + ggtitle('TBX21 (ATAC)')

pdf('../plots/harmony/lll_Th1.pdf', width = 12, height = 8)
(p1 | p2 | p3)/ (p4 | p5 | p6)
dev.off()

p1 <- FeaturePlot(lll, c('sct_CCR4'), reduction = 'wnn.umap', min.cutoff = 'q2', max.cutoff = 'q98') + ggtitle('CCR4 (RNA)')
p2 <- FeaturePlot(lll, c('sct_CCR6'), reduction = 'wnn.umap', min.cutoff = 'q2', max.cutoff = 'q98') + ggtitle('CCR6 (RNA)')
p3 <- FeaturePlot(lll, c('sct_RORC'), reduction = 'wnn.umap', min.cutoff = 'q2', max.cutoff = 'q98') + ggtitle('RORC (RNA)')
p4 <- FeaturePlot(lll, c('adt_CD194-A0071'), reduction = 'wnn.umap', cols = c("lightgrey", "darkgreen"), min.cutoff = 'q2', max.cutoff = 'q98') + ggtitle('CCR4 (ADT)')
p5 <- FeaturePlot(lll, c('adt_CD196-A0143'), reduction = 'wnn.umap', cols = c("lightgrey", "darkgreen"), min.cutoff = 'q2', max.cutoff = 'q98') + ggtitle('CCR6 (ADT)')
p6 <- FeaturePlot(lll, c('chromvar_MA1151.1'), reduction = 'wnn.umap', cols = c("lightgrey", "darkred"), min.cutoff = 'q2', max.cutoff = 'q98') + ggtitle('RORC (ATAC)')

pdf('../plots/harmony/lll_Th17_alt.pdf', width = 12, height = 8)
(p1 | p2 | p3)/ (p4 | p5 | p6)
dev.off()

p1 <- FeaturePlot(lll, c('sct_CXCR3'), reduction = 'wnn.umap', min.cutoff = 'q2', max.cutoff = 'q98') + ggtitle('CXCR3 (RNA)')
p2 <- FeaturePlot(lll, c('sct_CCR5'), reduction = 'wnn.umap', min.cutoff = 'q2', max.cutoff = 'q98') + ggtitle('CCR5 (RNA)')
p3 <- FeaturePlot(lll, c('sct_TBX21'), reduction = 'wnn.umap', min.cutoff = 'q2', max.cutoff = 'q98') + ggtitle('TBX21 (RNA)')
p4 <- FeaturePlot(lll, c('adt_CD183-A0140'), reduction = 'wnn.umap', cols = c("lightgrey", "darkgreen"), min.cutoff = 'q2', max.cutoff = 'q98') + ggtitle('CXCR3 (ADT)')
p5 <- FeaturePlot(lll, c('adt_CD195-A0141'), reduction = 'wnn.umap', cols = c("lightgrey", "darkgreen"), min.cutoff = 'q2', max.cutoff = 'q98') + ggtitle('CCR5 (ADT)')
p6 <- FeaturePlot(lll, c('chromvar_MA0690.1'), reduction = 'wnn.umap', cols = c("lightgrey", "darkred"), min.cutoff = 'q2', max.cutoff = 'q98') + ggtitle('TBX21 (ATAC)')

pdf('../plots/harmony/lll_Th1_alt.pdf', width = 12, height = 8)
(p1 | p2 | p3)/ (p4 | p5 | p6)
dev.off()

save(lll, file = '../output/lll_harmony_predicted_chrom.RData')
