library(Seurat)
library(Signac)
library(dplyr)
library(data.table)
library(stringr)
library(harmony)
library(celda)

setwd('~/RWorkSpace/CITE-seq/Duerr/DOGMA-seq/DIG_LLL/code/')

load('../output/lll_harmony_predicted_chrom.RData')

DefaultAssay(lll) <- "RNA"
sce <- as.SingleCellExperiment(lll)

raw <- Read10X_h5('../../../Duerr_20210419_DOGMAseq_PFA_LLL/outs/raw_feature_bc_matrix.h5')
raw <- CreateSeuratObject(raw$`Gene Expression`)
raw <- as.SingleCellExperiment(raw)
table(rownames(sce) == rownames(raw))
colnames(raw) <- paste0('LLL_', colnames(raw))
table(colnames(sce) %in% colnames(raw))

filtered <- Read10X_h5('../../../Duerr_20210419_DOGMAseq_PFA_LLL/outs/filtered_feature_bc_matrix.h5')
filtered <- CreateSeuratObject(filtered$`Gene Expression`)
filtered <- as.SingleCellExperiment(filtered)
table(rownames(sce) == rownames(filtered))
colnames(filtered) <- paste0('LLL_', colnames(filtered))
table(colnames(sce) %in% colnames(filtered))

raw_filtered <- raw[,!colnames(raw) %in% colnames(filtered)]
save(raw_filtered, file = '../output/lll_contamination.RData')

sce <- decontX(sce, background = raw_filtered)
decontx <- decontXcounts(sce)
saveRDS(decontx, file = '../output/lll_decontx.RDS')

umap <- reducedDim(sce, "decontX_UMAP")
plotDimReduceCluster(x = sce$predicted.celltype.l1,
                     dim1 = umap[, 1], dim2 = umap[, 2])
plotDecontXContamination(sce)
table(rownames(colData(sce)) == rownames(lll@meta.data))
lll$contamination <- colData(sce)$decontX_contamination*100

pdf('../plots/harmony/lll_contamination_original.pdf', width = 4, height = 4)
FeaturePlot(lll, 'contamination', reduction = 'wnn.umap') + ggtitle('Contamination Rate (%)') + scale_color_gradient(low = 'lightgrey', high = 'blue', limits = c(0,100), oob = scales::squish)
dev.off()

write.csv(lll@meta.data, file = '../output/lll_harmony_predicted_chrom_contamination_original.csv', quote = F, row.names = T, col.names = T)

load('../output/dig_contamination.RData')
dig_raw_filtered <- raw_filtered
load('../output/lll_contamination.RData')
raw_filtered <- cbind(dig_raw_filtered, raw_filtered)

sce <- decontX(sce, background = raw_filtered)

umap <- reducedDim(sce, "decontX_UMAP")
plotDimReduceCluster(x = sce$predicted.celltype.l1,
                     dim1 = umap[, 1], dim2 = umap[, 2])
plotDecontXContamination(sce)
table(rownames(colData(sce)) == rownames(lll@meta.data))
lll$contamination <- colData(sce)$decontX_contamination*100

pdf('../plots/harmony/lll_contamination.pdf', width = 4, height = 4)
FeaturePlot(lll, 'contamination', reduction = 'wnn.umap') + ggtitle('Contamination Rate (%)') + scale_color_gradient(low = 'lightgrey', high = 'blue', limits = c(0,100), oob = scales::squish)
dev.off()

write.csv(lll@meta.data, file = '../output/lll_harmony_predicted_chrom_contamination.csv', quote = F, row.names = T, col.names = T)

sce <- decontX(sce)

umap <- reducedDim(sce, "decontX_UMAP")
plotDimReduceCluster(x = sce$predicted.celltype.l1,
                     dim1 = umap[, 1], dim2 = umap[, 2])
plotDecontXContamination(sce)
table(rownames(colData(sce)) == rownames(lll@meta.data))
lll$contamination <- colData(sce)$decontX_contamination*100

pdf('../plots/harmony/lll_contamination_null.pdf', width = 4, height = 4)
FeaturePlot(lll, 'contamination', reduction = 'wnn.umap') + ggtitle('Contamination Rate (%)') + scale_color_gradient(low = 'lightgrey', high = 'blue', limits = c(0,100), oob = scales::squish)
dev.off()

write.csv(lll@meta.data, file = '../output/lll_harmony_predicted_chrom_contamination_null.csv', quote = F, row.names = T, col.names = T)
