library(Seurat)
library(Signac)
library(dplyr)
library(data.table)
library(stringr)
library(harmony)

setwd('~/RWorkSpace/CITE-seq/Duerr/DOGMA-seq/DIG_LLL/code/')

load('../data/DIG_data.RData')
dig <- data
load('../data/LLL_data.RData')
lll <- data
rm(data)

data <- merge(dig, lll, add.cell.ids = c('DIG', 'LLL'))
rm(dig, lll)
data$tech <- data$orig.ident

DefaultAssay(data) <- "RNA"
data <- SCTransform(data, method = "glmGamPoi") %>% RunPCA(reduction.name = "pca") %>%
  RunHarmony(group.by.vars = c("tech"), reduction = 'pca', assay.use = 'SCT',project.dim = FALSE,  reduction.save = "harmony_RNA")
data <- RunUMAP(data, reduction = "harmony_RNA", dims = 1:50, assay = 'SCT', reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_')

DefaultAssay(data) <- 'ADT'
# we will use all ADT features for dimensional reduction
# we set a dimensional reduction name to avoid overwriting the 
VariableFeatures(data) <- rownames(data[["ADT"]])
data <- NormalizeData(data, normalization.method = 'CLR', margin = 2) %>% 
  ScaleData() %>% RunPCA(reduction.name = 'apca') %>% 
  RunHarmony(group.by.vars = c("tech"), reduction = 'apca', assay.use = 'ADT',project.dim = FALSE,  reduction.save = "harmony_ADT")
data <- RunUMAP(data, reduction = "harmony_ADT", dims = 1:30, assay = 'ADT', reduction.name = 'adt.umap', reduction.key = 'adtUMAP_')

DefaultAssay(data) <- "peaks"
data <- RunTFIDF(data)
data <- FindTopFeatures(data, min.cutoff = 'q0')
data <- RunSVD(data)
data <- RunHarmony(data, group.by.vars = c("tech"), reduction = 'lsi', assay.use = 'peaks',project.dim = FALSE,  reduction.save = "harmony_peaks")
data <- RunUMAP(data, reduction = "harmony_peaks", dims = 2:50, assay = 'peaks', reduction.name = 'atac.umap', reduction.key = 'atacUMAP_')

# Now run multimodal neighbors and embedding
data <- FindMultiModalNeighbors(object = data,
                                reduction.list = list("harmony_RNA", "harmony_peaks", "harmony_ADT"),
                                dims.list = list(1:50, 2:50, 1:30))
data <- RunUMAP(data, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_" )
data <- FindClusters(data, graph.name = "wsnn", algorithm = 3, resolution = 0.2)

data$harmony_wnn_clusters <- Idents(data)

save(data, file = '../output/data_harmony.RData')

pdf('../plots/harmony/UMAP_harmony_weighted_cluster_split.pdf', width = 16, height = 8)
DimPlot(data, reduction = 'wnn.umap', split.by = 'tech', label = T, repel = TRUE, label.size = 4)
dev.off()

pdf('../plots/harmony/UMAP_harmony_RNA_cluster_split.pdf', width = 16, height = 8)
DimPlot(data, reduction = 'rna.umap', split.by = 'tech', label = T, repel = TRUE, label.size = 4)
dev.off()

pdf('../plots/harmony/UMAP_harmony_ADT_cluster_split.pdf', width = 16, height = 8)
DimPlot(data, reduction = 'adt.umap', split.by = 'tech', label = T, repel = TRUE, label.size = 4)
dev.off()

pdf('../plots/harmony/UMAP_harmony_ATAC_cluster_split.pdf', width = 16, height = 8)
DimPlot(data, reduction = 'atac.umap', split.by = 'tech', label = T, repel = TRUE, label.size = 4)
dev.off()

dig <- data[,data$tech == 'DIG']
lll <- data[,data$tech == 'LLL']
rm(data)

save(dig, file = '../output/dig_harmony.RData')
save(lll, file = '../output/lll_harmony.RData')
