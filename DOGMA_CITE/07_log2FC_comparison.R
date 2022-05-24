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

setwd('~/CITE-seq/Duerr/DOGMA-seq/DIG_CITE_rerun_1/code/')

load('../data/DIG_data.RData')

library(openxlsx)
Idents(data) <- 'condition'
rna.markers.1 <- FindMarkers(data, assay = 'SCT', ident.1 = 'SB775372_Act_IL1B_IL23', 'SB775372_Act_IL1B_IL23_PGE2', min.pct = 0, logfc.threshold = 0)
rna.markers.2 <- FindMarkers(data, assay = 'SCT', ident.1 = 'SB775372_Act_IL1B_IL23', 'SB775372_Act_IL1B_IL23_TGFB', min.pct = 0, logfc.threshold = 0)
rna.markers.3 <- FindMarkers(data, assay = 'SCT', ident.1 = 'SB775372_Act_IL1B_IL23', 'SB775372_Act_IL1B_IL23_PGE2_TGFB', min.pct = 0, logfc.threshold = 0)
rna.markers.4 <- FindMarkers(data, assay = 'SCT', ident.1 = 'SB775372_Act_IL1B_IL23_PGE2', 'SB775372_Act_IL1B_IL23_TGFB', min.pct = 0, logfc.threshold = 0)
rna.markers.5 <- FindMarkers(data, assay = 'SCT', ident.1 = 'SB775372_Act_IL1B_IL23_PGE2', 'SB775372_Act_IL1B_IL23_PGE2_TGFB', min.pct = 0, logfc.threshold = 0)
rna.markers.6 <- FindMarkers(data, assay = 'SCT', ident.1 = 'SB775372_Act_IL1B_IL23_TGFB', 'SB775372_Act_IL1B_IL23_PGE2_TGFB', min.pct = 0, logfc.threshold = 0)
rna.markers.7 <- FindMarkers(data, assay = 'SCT', ident.1 = 'SB775393_Act_IL1B_IL23', 'SB775393_Act_IL1B_IL23_PGE2', min.pct = 0, logfc.threshold = 0)
rna.markers.8 <- FindMarkers(data, assay = 'SCT', ident.1 = 'SB775393_Act_IL1B_IL23', 'SB775393_Act_IL1B_IL23_TGFB', min.pct = 0, logfc.threshold = 0)
rna.markers.9 <- FindMarkers(data, assay = 'SCT', ident.1 = 'SB775393_Act_IL1B_IL23', 'SB775393_Act_IL1B_IL23_PGE2_TGFB', min.pct = 0, logfc.threshold = 0)
rna.markers.10 <- FindMarkers(data, assay = 'SCT', ident.1 = 'SB775393_Act_IL1B_IL23_PGE2', 'SB775393_Act_IL1B_IL23_TGFB', min.pct = 0, logfc.threshold = 0)
rna.markers.11 <- FindMarkers(data, assay = 'SCT', ident.1 = 'SB775393_Act_IL1B_IL23_PGE2', 'SB775393_Act_IL1B_IL23_PGE2_TGFB', min.pct = 0, logfc.threshold = 0)
rna.markers.12 <- FindMarkers(data, assay = 'SCT', ident.1 = 'SB775393_Act_IL1B_IL23_TGFB', 'SB775393_Act_IL1B_IL23_PGE2_TGFB', min.pct = 0, logfc.threshold = 0)

DIG <- list(rna.markers.1, rna.markers.2, rna.markers.3, rna.markers.4, rna.markers.5, rna.markers.6, rna.markers.7, rna.markers.8, rna.markers.9, rna.markers.10, rna.markers.11, rna.markers.12)
rm(data)

load('../data/CITE_data.RData')

# data$condition <- data$orig.ident
Idents(data) <- 'condition'
rna.markers.1 <- FindMarkers(data, assay = 'SCT', ident.1 = 'SB775372_Act_IL1B_IL23', 'SB775372_Act_IL1B_IL23_PGE2', min.pct = 0, logfc.threshold = 0)
rna.markers.2 <- FindMarkers(data, assay = 'SCT', ident.1 = 'SB775372_Act_IL1B_IL23', 'SB775372_Act_IL1B_IL23_TGFB', min.pct = 0, logfc.threshold = 0)
rna.markers.3 <- FindMarkers(data, assay = 'SCT', ident.1 = 'SB775372_Act_IL1B_IL23', 'SB775372_Act_IL1B_IL23_PGE2_TGFB', min.pct = 0, logfc.threshold = 0)
rna.markers.4 <- FindMarkers(data, assay = 'SCT', ident.1 = 'SB775372_Act_IL1B_IL23_PGE2', 'SB775372_Act_IL1B_IL23_TGFB', min.pct = 0, logfc.threshold = 0)
rna.markers.5 <- FindMarkers(data, assay = 'SCT', ident.1 = 'SB775372_Act_IL1B_IL23_PGE2', 'SB775372_Act_IL1B_IL23_PGE2_TGFB', min.pct = 0, logfc.threshold = 0)
rna.markers.6 <- FindMarkers(data, assay = 'SCT', ident.1 = 'SB775372_Act_IL1B_IL23_TGFB', 'SB775372_Act_IL1B_IL23_PGE2_TGFB', min.pct = 0, logfc.threshold = 0)
rna.markers.7 <- FindMarkers(data, assay = 'SCT', ident.1 = 'SB775393_Act_IL1B_IL23', 'SB775393_Act_IL1B_IL23_PGE2', min.pct = 0, logfc.threshold = 0)
rna.markers.8 <- FindMarkers(data, assay = 'SCT', ident.1 = 'SB775393_Act_IL1B_IL23', 'SB775393_Act_IL1B_IL23_TGFB', min.pct = 0, logfc.threshold = 0)
rna.markers.9 <- FindMarkers(data, assay = 'SCT', ident.1 = 'SB775393_Act_IL1B_IL23', 'SB775393_Act_IL1B_IL23_PGE2_TGFB', min.pct = 0, logfc.threshold = 0)
rna.markers.10 <- FindMarkers(data, assay = 'SCT', ident.1 = 'SB775393_Act_IL1B_IL23_PGE2', 'SB775393_Act_IL1B_IL23_TGFB', min.pct = 0, logfc.threshold = 0)
rna.markers.11 <- FindMarkers(data, assay = 'SCT', ident.1 = 'SB775393_Act_IL1B_IL23_PGE2', 'SB775393_Act_IL1B_IL23_PGE2_TGFB', min.pct = 0, logfc.threshold = 0)
rna.markers.12 <- FindMarkers(data, assay = 'SCT', ident.1 = 'SB775393_Act_IL1B_IL23_TGFB', 'SB775393_Act_IL1B_IL23_PGE2_TGFB', min.pct = 0, logfc.threshold = 0)

CITE <- list(rna.markers.1, rna.markers.2, rna.markers.3, rna.markers.4, rna.markers.5, rna.markers.6, rna.markers.7, rna.markers.8, rna.markers.9, rna.markers.10, rna.markers.11, rna.markers.12)
rm(data)

save(DIG, CITE, file = '../output/RNA_markers.RData')

load('../output/RNA_markers.RData')
load('../output/exon_rate.RData')

library(BuenColors)
library(plotly)
library(ggrepel)

comparison <- function(index, name){
  tmp1 <- DIG[[index]]
  tmp1$gene <- rownames(tmp1)
  tmp2 <- CITE[[index]]
  tmp2$gene <- rownames(tmp2)
  tmp <- merge(tmp1, tmp2, by = 'gene', suffixes = c('_DIG', '_LLL'))
  
  exon_rate_tmp <- exon_rate[tmp$gene]
  tmp$exon <- exon_rate_tmp
  tmp$exon_bi <- ifelse(tmp$exon > 0.5, 'Exon-dominated Genes', 'Intron-dominated Genes')
  
  p1 <- ggplot(tmp, aes(x = avg_log2FC_LLL, y = avg_log2FC_DIG, color = exon)) +
    geom_point() +
    pretty_plot(fontsize = 10) + L_border() + 
    scale_color_gradient(low = "dodgerblue3", high = "firebrick") +
    labs(x = "CITE log2FC", y = "DIG log2FC", color = "Proportion of Exonic UMIs") +
    theme(legend.position = "bottom", legend.direction = 'horizontal')
  cowplot::ggsave2(p1, file = paste0("../plots/log2FC_", name, ".pdf"), width = 4, height = 4)
  
  p1 <- ggplot(tmp, aes(x = avg_log2FC_LLL, y = avg_log2FC_DIG, color = exon_bi)) +
    geom_point() +
    pretty_plot(fontsize = 10) + L_border() + 
    scale_color_manual(values = c("dodgerblue3", "firebrick")) +
    labs(x = "CITE log2FC", y = "DIG log2FC", color = "") +
    theme(legend.position = "bottom", legend.direction = 'horizontal')
  cowplot::ggsave2(p1, file = paste0("../plots/log2FC_", name, "_bi.pdf"), width = 4, height = 4)
}

comparison(1, 'SB775372_Act_IL1B_IL23_vs_Act_IL1B_IL23_PGE2')
comparison(2, 'SB775372_Act_IL1B_IL23_vs_Act_IL1B_IL23_TGFB')
comparison(3, 'SB775372_Act_IL1B_IL23_vs_Act_IL1B_IL23_PGE2_TGFB')
comparison(4, 'SB775372_Act_IL1B_IL23_PGE2_vs_Act_IL1B_IL23_TGFB')
comparison(5, 'SB775372_Act_IL1B_IL23_PGE2_vs_Act_IL1B_IL23_PGE2_TGFB')
comparison(6, 'SB775372_Act_IL1B_IL23_TGFB_vs_Act_IL1B_IL23_PGE2_TGFB')
comparison(7, 'SB775393_Act_IL1B_IL23_vs_Act_IL1B_IL23_PGE2')
comparison(8, 'SB775393_Act_IL1B_IL23_vs_Act_IL1B_IL23_TGFB')
comparison(9, 'SB775393_Act_IL1B_IL23_vs_Act_IL1B_IL23_PGE2_TGFB')
comparison(10, 'SB775393_Act_IL1B_IL23_PGE2_vs_Act_IL1B_IL23_TGFB')
comparison(11, 'SB775393_Act_IL1B_IL23_PGE2_vs_Act_IL1B_IL23_PGE2_TGFB')
comparison(12, 'SB775393_Act_IL1B_IL23_TGFB_vs_Act_IL1B_IL23_PGE2_TGFB')

comparison_label <- function(index, name){
  tmp1 <- DIG[[index]]
  tmp1$gene <- rownames(tmp1)
  tmp2 <- CITE[[index]]
  tmp2$gene <- rownames(tmp2)
  tmp <- merge(tmp1, tmp2, by = 'gene', suffixes = c('_DIG', '_LLL'))
  
  exon_rate_tmp <- exon_rate[tmp$gene]
  tmp$exon <- exon_rate_tmp
  tmp$exon_bi <- ifelse(tmp$exon > 0.5, 'Exon-dominated Genes', 'Intron-dominated Genes')
  
  p1 <- ggplot(tmp, aes(x = avg_log2FC_LLL, y = avg_log2FC_DIG, color = exon_bi)) +
    geom_point() +
    geom_text_repel(max.overlaps = Inf, aes(label=ifelse(gene %in% c('IL7R', 'CXCR4', 'PDE3B', 'PLCL1', 'TMSB4X', 'CCL4', 'IFNG', 'LTA', 'GZMB', 'IL2', 'CCL3'),gene,'')))+
    pretty_plot(fontsize = 10) + L_border() + 
    scale_color_manual(values = c("dodgerblue3", "firebrick")) +
    labs(x = "CITE log2FC", y = "DIG log2FC", color = "") +
    theme(legend.position = "bottom", legend.direction = 'horizontal')
  cowplot::ggsave2(p1, file = paste0("../plots/log2FC_", name, "_bi_label.pdf"), width = 4, height = 4)
}

comparison_label(1, 'SB775372_Act_IL1B_IL23_vs_Act_IL1B_IL23_PGE2')
comparison_label(7, 'SB775393_Act_IL1B_IL23_vs_Act_IL1B_IL23_PGE2')
