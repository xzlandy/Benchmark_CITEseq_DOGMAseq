library(Seurat)
library(Signac)
library(dplyr)
library(data.table)
library(stringr)
library(harmony)
library(celda)
library(SingleCellExperiment)
library(ggplot2)

setwd('~/CITE-seq/Duerr/DOGMA-seq/DIG_CITE_rerun_1/code/')

load('../output/cite_harmony_independent_predicted.RData')

DefaultAssay(cite) <- "RNA"
sce <- as.SingleCellExperiment(cite)

raw <- Read10X_h5('../../../Duerr_20210610_CITEseq_1/outs/raw_feature_bc_matrix.h5')
raw <- CreateSeuratObject(raw$`Gene Expression`)
raw <- as.SingleCellExperiment(raw)
table(rownames(sce) == rownames(raw))
colnames(raw) <- paste0('CITE_', colnames(raw))
table(colnames(sce) %in% colnames(raw))

filtered <- Read10X_h5('../../../Duerr_20210610_CITEseq_1/outs/filtered_feature_bc_matrix.h5')
filtered <- CreateSeuratObject(filtered$`Gene Expression`)
filtered <- as.SingleCellExperiment(filtered)
table(rownames(sce) == rownames(filtered))
colnames(filtered) <- paste0('CITE_', colnames(filtered))
table(colnames(sce) %in% colnames(filtered))

raw_filtered <- raw[,!colnames(raw) %in% colnames(filtered)]
save(raw_filtered, file = '../output/cite_contamination.RData')

sce <- decontX(sce, background = raw_filtered)

umap <- reducedDim(sce, "decontX_UMAP")
plotDimReduceCluster(x = sce$predicted.celltype.l1,
                     dim1 = umap[, 1], dim2 = umap[, 2])
plotDecontXContamination(sce)
table(rownames(colData(sce)) == rownames(cite@meta.data))
cite$contamination <- colData(sce)$decontX_contamination*100

pdf('../plots/harmony/cite_contamination_original.pdf', width = 4, height = 4)
FeaturePlot(cite, 'contamination', reduction = 'wnn.umap') + ggtitle('Contamination Rate (%)') + scale_color_gradient(low = 'lightgrey', high = 'blue', limits = c(0,100), oob = scales::squish)
dev.off()

write.csv(cite@meta.data, file = '../output/cite_harmony_predicted_chrom_contamination_original.csv', quote = F, row.names = T, col.names = T)

load('../output/dig_contamination.RData')
dig_raw_filtered <- raw_filtered
load('../output/cite_contamination.RData')
raw_filtered <- cbind(dig_raw_filtered, raw_filtered)

sce <- decontX(sce, background = raw_filtered)

umap <- reducedDim(sce, "decontX_UMAP")
plotDimReduceCluster(x = sce$predicted.celltype.l1,
                     dim1 = umap[, 1], dim2 = umap[, 2])
plotDecontXContamination(sce)
table(rownames(colData(sce)) == rownames(cite@meta.data))
cite$contamination <- colData(sce)$decontX_contamination*100

pdf('../plots/harmony/cite_contamination.pdf', width = 4, height = 4)
FeaturePlot(cite, 'contamination', reduction = 'wnn.umap') + ggtitle('Contamination Rate (%)') + scale_color_gradient(low = 'lightgrey', high = 'blue', limits = c(0,100), oob = scales::squish)
dev.off()

write.csv(cite@meta.data, file = '../output/cite_harmony_predicted_chrom_contamination.csv', quote = F, row.names = T, col.names = T)

sce <- decontX(sce)

umap <- reducedDim(sce, "decontX_UMAP")
plotDimReduceCluster(x = sce$predicted.celltype.l1,
                     dim1 = umap[, 1], dim2 = umap[, 2])
plotDecontXContamination(sce)
table(rownames(colData(sce)) == rownames(cite@meta.data))
cite$contamination <- colData(sce)$decontX_contamination*100

pdf('../plots/harmony/cite_contamination_null.pdf', width = 4, height = 4)
FeaturePlot(cite, 'contamination', reduction = 'wnn.umap') + ggtitle('Contamination Rate (%)') + scale_color_gradient(low = 'lightgrey', high = 'blue', limits = c(0,100), oob = scales::squish)
dev.off()

write.csv(cite@meta.data, file = '../output/cite_harmony_predicted_chrom_contamination_null.csv', quote = F, row.names = T, col.names = T)
