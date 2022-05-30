library(Seurat)
library(Signac)
library(dplyr)
library(data.table)
library(stringr)
library(harmony)

setwd('~/RWorkSpace/CITE-seq/Duerr/DOGMA-seq/DIG_CITE_rerun_1/code/')

load('../data/DIG_data.RData')
dig <- data
load('../data/CITE_data.RData')
cite <- data
rm(data)

dig$tech <- 'DIG'
cite$tech <- 'CITE'

dig$sample <- str_split_fixed(dig$condition, '_Act', 2)[,1]
# cite$condition <- cite$orig.ident
cite$sample <- str_split_fixed(cite$condition, '_Act', 2)[,1]

DefaultAssay(dig) <- 'ADT'
DefaultAssay(cite) <- 'ADT'

dig <- FindVariableFeatures(dig, nfeatures = 100)
cite <- FindVariableFeatures(cite, nfeatures = 100)

VariableFeatures_ADT_dig <- VariableFeatures(dig)
VariableFeatures_ADT_cite <- VariableFeatures(cite)

IntegrationFeatures_ADT <- intersect(VariableFeatures_ADT_dig, VariableFeatures_ADT_cite)

# VariableFeatures_ADT_dig <- VariableFeatures(dig)
# tmp <- str_split_fixed(str_split_fixed(VariableFeatures_ADT_dig, '-A0', 2)[,1], '-A1', 2)[,1]
# VariableFeatures_ADT_dig <- data.frame(origin = VariableFeatures_ADT_dig, modified = tmp)
# VariableFeatures_ADT_cite <- VariableFeatures(cite)
# tmp <- str_split_fixed(VariableFeatures_ADT_cite, '\\.1', 2)[,1]
# VariableFeatures_ADT_cite <- data.frame(origin = VariableFeatures_ADT_cite, modified = tmp)
# 
# IntegrationFeatures_ADT <- VariableFeatures_ADT_dig[VariableFeatures_ADT_dig$modified %in% intersect(VariableFeatures_ADT_dig$modified, VariableFeatures_ADT_cite$modified),]$origin
# 
# cite@assays$ADT@counts@Dimnames[[1]] <- rownames(dig)
# rownames(cite@assays$ADT@data) <- rownames(dig)
# rownames(cite@assays$ADT@scale.data) <- rownames(dig)

DefaultAssay(dig) <- 'SCT'
DefaultAssay(cite) <- 'SCT'

VariableFeatures_SCT_dig <- VariableFeatures(dig)
VariableFeatures_SCT_cite <- VariableFeatures(cite)

IntegrationFeatures_SCT <- intersect(VariableFeatures_SCT_dig, VariableFeatures_SCT_cite)

data <- merge(dig, cite, add.cell.ids = c('DIG', 'CITE'))
rm(dig, cite)

DefaultAssay(data) <- "SCT"
VariableFeatures(data) <- IntegrationFeatures_SCT
data <- RunPCA(data, reduction.name = "pca") %>% RunHarmony(group.by.vars = c("tech"), reduction = 'pca', assay.use = 'SCT',project.dim = FALSE,  reduction.save = "harmony_RNA")
data <- RunUMAP(data, reduction = "harmony_RNA", dims = 1:50, assay = 'SCT', reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_')

DefaultAssay(data) <- 'ADT'
# we will use all ADT features for dimensional reduction
# we set a dimensional reduction name to avoid overwriting the 
VariableFeatures(data) <- IntegrationFeatures_ADT
data <- ScaleData(data) %>% RunPCA(reduction.name = 'apca') %>% 
  RunHarmony(group.by.vars = c("tech"), reduction = 'apca', assay.use = 'ADT',project.dim = FALSE,  reduction.save = "harmony_ADT")
data <- RunUMAP(data, reduction = "harmony_ADT", dims = 1:30, assay = 'ADT', reduction.name = 'adt.umap', reduction.key = 'adtUMAP_')

data <- FindMultiModalNeighbors(object = data,
                                reduction.list = list("harmony_RNA", "harmony_ADT"),
                                dims.list = list(1:50, 1:30))
data <- RunUMAP(data, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_" )
data <- FindClusters(data, graph.name = "wsnn", algorithm = 3, resolution = 0.2)

data$harmony_wnn_clusters <- Idents(data)

save(data, file = '../output/data_harmony_independent.RData')

pdf('../plots/harmony/UMAP_harmony_weighted_cluster_split.pdf', width = 16, height = 8)
DimPlot(data, reduction = 'wnn.umap', split.by = 'tech', label = T, repel = TRUE, label.size = 4)
dev.off()

pdf('../plots/harmony/UMAP_harmony_RNA_cluster_split.pdf', width = 16, height = 8)
DimPlot(data, reduction = 'rna.umap', split.by = 'tech', label = T, repel = TRUE, label.size = 4)
dev.off()

pdf('../plots/harmony/UMAP_harmony_ADT_cluster_split.pdf', width = 16, height = 8)
DimPlot(data, reduction = 'adt.umap', split.by = 'tech', label = T, repel = TRUE, label.size = 4)
dev.off()

dig <- data[,data$tech == 'DIG']
cite <- data[,data$tech == 'CITE']
rm(data)

save(dig, file = '../output/dig_harmony_independent.RData')
save(cite, file = '../output/cite_harmony_independent.RData')
