library(Seurat)
library(Signac)
library(dplyr)
library(data.table)
library(stringr)
library(harmony)
library(patchwork)

setwd('~/RWorkSpace/CITE-seq/Duerr/DOGMA-seq/DIG_LLL/code/')

load('../output/dig_harmony_predicted_chrom.RData')

DefaultAssay(dig) <- "RNA"
dig <- FindNeighbors(dig, dims = 1:50, reduction = 'harmony_RNA')
dig <- FindClusters(dig, resolution = 0.2, graph.name = 'SCT_snn')
dig$rna_clusters <- Idents(dig)

# ADT analysis
DefaultAssay(dig) <- 'ADT'
dig <- FindNeighbors(dig, dims = 1:30, reduction = 'harmony_ADT')
dig <- FindClusters(dig, resolution = 0.2, graph.name = 'ADT_snn')
dig$adt_clusters <- Idents(dig)

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(dig) <- "peaks"
dig <- FindNeighbors(dig, dims = 2:50, reduction = 'harmony_peaks')
dig <- FindClusters(dig, resolution = 0.2, graph.name = 'peaks_snn', algorithm = 3)
dig$atac_clusters <- Idents(dig)

dig <- FindMultiModalNeighbors(object = dig,
                               reduction.list = list("harmony_RNA", "harmony_peaks", "harmony_ADT"),
                               dims.list = list(1:50, 2:50, 1:30))
# dig <- RunUMAP(dig, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_" )
dig <- FindClusters(dig, graph.name = "wsnn", algorithm = 3, resolution = 0.2)

dig$wnn_clusters <- Idents(dig)

Idents(dig) <- 'wnn_clusters'
p1 <- DimPlot(dig, reduction = 'wnn.umap', label = T, 
              repel = TRUE, label.size = 4)
Idents(dig) <- 'rna_clusters'
p2 <- DimPlot(dig, reduction = 'rna.umap', label = T, 
              repel = TRUE, label.size = 4)
Idents(dig) <- 'adt_clusters'
p3 <- DimPlot(dig, reduction = 'adt.umap', label = T, 
              repel = TRUE, label.size = 4)
Idents(dig) <- 'atac_clusters'
p4 <- DimPlot(dig, reduction = 'atac.umap', label = T, 
              repel = TRUE, label.size = 4)
Idents(dig) <- 'predicted.celltype.l1'
dig@active.ident <- factor(dig@active.ident, levels = names(table(dig$predicted.celltype.l1)))
p5 <- DimPlot(dig, reduction = 'wnn.umap', label = T, 
              repel = TRUE, label.size = 4)
p6 <- DimPlot(dig, reduction = 'rna.umap', label = T, 
              repel = TRUE, label.size = 4)
p7 <- DimPlot(dig, reduction = 'adt.umap', label = T, 
              repel = TRUE, label.size = 4)
p8 <- DimPlot(dig, reduction = 'atac.umap', label = T, 
              repel = TRUE, label.size = 4)
pdf('../plots/harmony/dig_all_slide.pdf', width = 32, height = 16)
wrap_plots(p1,p2,p3,p4,p5,p6,p7,p8, ncol = 4)
dev.off()

save(dig, file = '../output/dig_harmony_predicted_chrom_cluster.RData')
