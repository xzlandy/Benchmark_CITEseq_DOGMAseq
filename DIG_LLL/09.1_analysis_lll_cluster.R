library(Seurat)
library(Signac)
library(dplyr)
library(data.table)
library(stringr)
library(harmony)
library(patchwork)

setwd('~/RWorkSpace/CITE-seq/Duerr/DOGMA-seq/DIG_LLL/code/')

load('../output/lll_harmony_predicted_chrom.RData')

DefaultAssay(lll) <- "RNA"
lll <- FindNeighbors(lll, dims = 1:50, reduction = 'harmony_RNA')
lll <- FindClusters(lll, resolution = 0.2, graph.name = 'SCT_snn')
lll$rna_clusters <- Idents(lll)

# ADT analysis
DefaultAssay(lll) <- 'ADT'
lll <- FindNeighbors(lll, dims = 1:30, reduction = 'harmony_ADT')
lll <- FindClusters(lll, resolution = 0.2, graph.name = 'ADT_snn')
lll$adt_clusters <- Idents(lll)

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(lll) <- "peaks"
lll <- FindNeighbors(lll, dims = 2:50, reduction = 'harmony_peaks')
lll <- FindClusters(lll, resolution = 0.2, graph.name = 'peaks_snn', algorithm = 3)
lll$atac_clusters <- Idents(lll)

lll <- FindMultiModalNeighbors(object = lll,
                               reduction.list = list("harmony_RNA", "harmony_peaks", "harmony_ADT"),
                               dims.list = list(1:50, 2:50, 1:30))
# lll <- RunUMAP(lll, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_" )
lll <- FindClusters(lll, graph.name = "wsnn", algorithm = 3, resolution = 0.2)

lll$wnn_clusters <- Idents(lll)

Idents(lll) <- 'wnn_clusters'
p1 <- DimPlot(lll, reduction = 'wnn.umap', label = T, 
              repel = TRUE, label.size = 4)
Idents(lll) <- 'rna_clusters'
p2 <- DimPlot(lll, reduction = 'rna.umap', label = T, 
              repel = TRUE, label.size = 4)
Idents(lll) <- 'adt_clusters'
p3 <- DimPlot(lll, reduction = 'adt.umap', label = T, 
              repel = TRUE, label.size = 4)
Idents(lll) <- 'atac_clusters'
p4 <- DimPlot(lll, reduction = 'atac.umap', label = T, 
              repel = TRUE, label.size = 4)
Idents(lll) <- 'predicted.celltype.l1'
lll@active.ident <- factor(lll@active.ident, levels = names(table(lll$predicted.celltype.l1)))
p5 <- DimPlot(lll, reduction = 'wnn.umap', label = T, 
              repel = TRUE, label.size = 4)
p6 <- DimPlot(lll, reduction = 'rna.umap', label = T, 
              repel = TRUE, label.size = 4)
p7 <- DimPlot(lll, reduction = 'adt.umap', label = T, 
              repel = TRUE, label.size = 4)
p8 <- DimPlot(lll, reduction = 'atac.umap', label = T, 
              repel = TRUE, label.size = 4)
pdf('../plots/harmony/lll_all_slide.pdf', width = 32, height = 16)
wrap_plots(p1,p2,p3,p4,p5,p6,p7,p8, ncol = 4)
dev.off()

save(lll, file = '../output/lll_harmony_predicted_chrom_cluster.RData')
