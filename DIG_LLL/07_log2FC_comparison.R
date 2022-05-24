library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggplot2)
library(patchwork)
library(data.table)
library(Matrix)
library(dplyr)
library(future)

plan("multisession")

setwd('~/RWorkSpace/CITE-seq/Duerr/DOGMA-seq/DIG_LLL/code/')

load('../data/DIG_data.RData')

library(openxlsx)
Idents(data) <- 'condition'
rna.markers <- FindAllMarkers(data, assay = 'SCT', min.pct = 0, logfc.threshold = 0)
rna.markers.1 <- FindMarkers(data, assay = 'SCT', ident.1 = 'ACT_IL1B_IL23', 'ACT_IL1B_IL23_PGE2', min.pct = 0, logfc.threshold = 0)
rna.markers.2 <- FindMarkers(data, assay = 'SCT', ident.1 = 'ACT_IL1B_IL23', 'ACT_IL1B_IL23_TGFB', min.pct = 0, logfc.threshold = 0)
rna.markers.3 <- FindMarkers(data, assay = 'SCT', ident.1 = 'ACT_IL1B_IL23', 'ACT_IL1B_IL23_PGE2_TGFB', min.pct = 0, logfc.threshold = 0)
rna.markers.4 <- FindMarkers(data, assay = 'SCT', ident.1 = 'ACT_IL1B_IL23_PGE2', 'ACT_IL1B_IL23_TGFB', min.pct = 0, logfc.threshold = 0)
rna.markers.5 <- FindMarkers(data, assay = 'SCT', ident.1 = 'ACT_IL1B_IL23_PGE2', 'ACT_IL1B_IL23_PGE2_TGFB', min.pct = 0, logfc.threshold = 0)
rna.markers.6 <- FindMarkers(data, assay = 'SCT', ident.1 = 'ACT_IL1B_IL23_TGFB', 'ACT_IL1B_IL23_PGE2_TGFB', min.pct = 0, logfc.threshold = 0)

DIG <- list(rna.markers, rna.markers.1, rna.markers.2, rna.markers.3, rna.markers.4, rna.markers.5, rna.markers.6)
rm(data)

load('../data/LLL_data.RData')

Idents(data) <- 'condition'
rna.markers <- FindAllMarkers(data, assay = 'SCT', min.pct = 0, logfc.threshold = 0)
rna.markers.1 <- FindMarkers(data, assay = 'SCT', ident.1 = 'ACT_IL1B_IL23', 'ACT_IL1B_IL23_PGE2', min.pct = 0, logfc.threshold = 0)
rna.markers.2 <- FindMarkers(data, assay = 'SCT', ident.1 = 'ACT_IL1B_IL23', 'ACT_IL1B_IL23_TGFB', min.pct = 0, logfc.threshold = 0)
rna.markers.3 <- FindMarkers(data, assay = 'SCT', ident.1 = 'ACT_IL1B_IL23', 'ACT_IL1B_IL23_PGE2_TGFB', min.pct = 0, logfc.threshold = 0)
rna.markers.4 <- FindMarkers(data, assay = 'SCT', ident.1 = 'ACT_IL1B_IL23_PGE2', 'ACT_IL1B_IL23_TGFB', min.pct = 0, logfc.threshold = 0)
rna.markers.5 <- FindMarkers(data, assay = 'SCT', ident.1 = 'ACT_IL1B_IL23_PGE2', 'ACT_IL1B_IL23_PGE2_TGFB', min.pct = 0, logfc.threshold = 0)
rna.markers.6 <- FindMarkers(data, assay = 'SCT', ident.1 = 'ACT_IL1B_IL23_TGFB', 'ACT_IL1B_IL23_PGE2_TGFB', min.pct = 0, logfc.threshold = 0)

LLL <- list(rna.markers, rna.markers.1, rna.markers.2, rna.markers.3, rna.markers.4, rna.markers.5, rna.markers.6)
rm(data)

save(DIG, LLL, file = '../output/RNA_markers.RData')

load('../output/RNA_markers.RData')

library(BuenColors)

comparison <- function(index, name){
  tmp1 <- DIG[[index]]
  tmp1$gene <- rownames(tmp1)
  tmp2 <- LLL[[index]]
  tmp2$gene <- rownames(tmp2)
  tmp <- merge(tmp1, tmp2, by = 'gene', suffixes = c('_DIG', '_LLL'))
  
  p1 <- ggplot(tmp, aes(x = avg_log2FC_LLL, y = avg_log2FC_DIG)) +
    geom_point() +
    pretty_plot(fontsize = 10) + L_border() + 
    labs(x = "LLL log2FC", y = "DIG log2FC")
  cowplot::ggsave2(p1, file = paste0("../plots/log2FC_", name, ".pdf"), width = 4, height = 4)
}

comparison(2, 'ACT_IL1B_IL23_vs_ACT_IL1B_IL23_PGE2')
comparison(3, 'ACT_IL1B_IL23_vs_ACT_IL1B_IL23_TGFB')
comparison(4, 'ACT_IL1B_IL23_vs_ACT_IL1B_IL23_PGE2_TGFB')
comparison(5, 'ACT_IL1B_IL23_PGE2_vs_ACT_IL1B_IL23_TGFB')
comparison(6, 'ACT_IL1B_IL23_PGE2_vs_ACT_IL1B_IL23_PGE2_TGFB')
comparison(7, 'ACT_IL1B_IL23_TGFB_vs_ACT_IL1B_IL23_PGE2_TGFB')

comparison_label <- function(index, name){
  tmp1 <- DIG[[index]]
  tmp1$gene <- rownames(tmp1)
  tmp2 <- LLL[[index]]
  tmp2$gene <- rownames(tmp2)
  tmp <- merge(tmp1, tmp2, by = 'gene', suffixes = c('_DIG', '_LLL'))
  
  p1 <- ggplot(tmp, aes(x = avg_log2FC_LLL, y = avg_log2FC_DIG)) +
    geom_point() +
    geom_text_repel(max.overlaps = Inf, aes(label=ifelse(gene %in% c('IL7R', 'PDE3B', 'PLCL1', 'TMSB4X', 'CCL4', 'IFNG', 'LTA', 'GZMB', 'IL2', 'CCL3'),gene,''), color=ifelse(avg_log2FC_DIG > 0, 'positive', 'negative')))+
    pretty_plot(fontsize = 10) + L_border() + 
    labs(x = "LLL log2FC", y = "DIG log2FC") +
    scale_color_manual(values = c("dodgerblue3", "firebrick")) +
    theme(legend.position = "none")
  
  cowplot::ggsave2(p1, file = paste0("../plots/log2FC_", name, "_label.pdf"), width = 4, height = 4)
}

comparison_label(2, 'ACT_IL1B_IL23_vs_ACT_IL1B_IL23_PGE2')
