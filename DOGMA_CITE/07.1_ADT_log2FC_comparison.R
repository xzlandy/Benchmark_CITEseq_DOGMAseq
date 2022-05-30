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
library(stringr)

plan("multisession")

setwd('~/RWorkSpace/CITE-seq/Duerr/DOGMA-seq/DIG_CITE_rerun_1/code/')

load('../data/DIG_data.RData')

library(openxlsx)
Idents(data) <- 'condition'
rna.markers.1 <- FindMarkers(data, assay = 'ADT', ident.1 = 'SB775372_Act_IL1B_IL23', 'SB775372_Act_IL1B_IL23_PGE2', min.pct = 0, logfc.threshold = 0)
rna.markers.2 <- FindMarkers(data, assay = 'ADT', ident.1 = 'SB775372_Act_IL1B_IL23', 'SB775372_Act_IL1B_IL23_TGFB', min.pct = 0, logfc.threshold = 0)
rna.markers.3 <- FindMarkers(data, assay = 'ADT', ident.1 = 'SB775372_Act_IL1B_IL23', 'SB775372_Act_IL1B_IL23_PGE2_TGFB', min.pct = 0, logfc.threshold = 0)
rna.markers.4 <- FindMarkers(data, assay = 'ADT', ident.1 = 'SB775372_Act_IL1B_IL23_PGE2', 'SB775372_Act_IL1B_IL23_TGFB', min.pct = 0, logfc.threshold = 0)
rna.markers.5 <- FindMarkers(data, assay = 'ADT', ident.1 = 'SB775372_Act_IL1B_IL23_PGE2', 'SB775372_Act_IL1B_IL23_PGE2_TGFB', min.pct = 0, logfc.threshold = 0)
rna.markers.6 <- FindMarkers(data, assay = 'ADT', ident.1 = 'SB775372_Act_IL1B_IL23_TGFB', 'SB775372_Act_IL1B_IL23_PGE2_TGFB', min.pct = 0, logfc.threshold = 0)
rna.markers.7 <- FindMarkers(data, assay = 'ADT', ident.1 = 'SB775393_Act_IL1B_IL23', 'SB775393_Act_IL1B_IL23_PGE2', min.pct = 0, logfc.threshold = 0)
rna.markers.8 <- FindMarkers(data, assay = 'ADT', ident.1 = 'SB775393_Act_IL1B_IL23', 'SB775393_Act_IL1B_IL23_TGFB', min.pct = 0, logfc.threshold = 0)
rna.markers.9 <- FindMarkers(data, assay = 'ADT', ident.1 = 'SB775393_Act_IL1B_IL23', 'SB775393_Act_IL1B_IL23_PGE2_TGFB', min.pct = 0, logfc.threshold = 0)
rna.markers.10 <- FindMarkers(data, assay = 'ADT', ident.1 = 'SB775393_Act_IL1B_IL23_PGE2', 'SB775393_Act_IL1B_IL23_TGFB', min.pct = 0, logfc.threshold = 0)
rna.markers.11 <- FindMarkers(data, assay = 'ADT', ident.1 = 'SB775393_Act_IL1B_IL23_PGE2', 'SB775393_Act_IL1B_IL23_PGE2_TGFB', min.pct = 0, logfc.threshold = 0)
rna.markers.12 <- FindMarkers(data, assay = 'ADT', ident.1 = 'SB775393_Act_IL1B_IL23_TGFB', 'SB775393_Act_IL1B_IL23_PGE2_TGFB', min.pct = 0, logfc.threshold = 0)

DIG <- list(rna.markers.1, rna.markers.2, rna.markers.3, rna.markers.4, rna.markers.5, rna.markers.6, rna.markers.7, rna.markers.8, rna.markers.9, rna.markers.10, rna.markers.11, rna.markers.12)
rm(data)

load('../data/CITE_data.RData')

# data$condition <- data$orig.ident
Idents(data) <- 'condition'
rna.markers.1 <- FindMarkers(data, assay = 'ADT', ident.1 = 'SB775372_Act_IL1B_IL23', 'SB775372_Act_IL1B_IL23_PGE2', min.pct = 0, logfc.threshold = 0)
rna.markers.2 <- FindMarkers(data, assay = 'ADT', ident.1 = 'SB775372_Act_IL1B_IL23', 'SB775372_Act_IL1B_IL23_TGFB', min.pct = 0, logfc.threshold = 0)
rna.markers.3 <- FindMarkers(data, assay = 'ADT', ident.1 = 'SB775372_Act_IL1B_IL23', 'SB775372_Act_IL1B_IL23_PGE2_TGFB', min.pct = 0, logfc.threshold = 0)
rna.markers.4 <- FindMarkers(data, assay = 'ADT', ident.1 = 'SB775372_Act_IL1B_IL23_PGE2', 'SB775372_Act_IL1B_IL23_TGFB', min.pct = 0, logfc.threshold = 0)
rna.markers.5 <- FindMarkers(data, assay = 'ADT', ident.1 = 'SB775372_Act_IL1B_IL23_PGE2', 'SB775372_Act_IL1B_IL23_PGE2_TGFB', min.pct = 0, logfc.threshold = 0)
rna.markers.6 <- FindMarkers(data, assay = 'ADT', ident.1 = 'SB775372_Act_IL1B_IL23_TGFB', 'SB775372_Act_IL1B_IL23_PGE2_TGFB', min.pct = 0, logfc.threshold = 0)
rna.markers.7 <- FindMarkers(data, assay = 'ADT', ident.1 = 'SB775393_Act_IL1B_IL23', 'SB775393_Act_IL1B_IL23_PGE2', min.pct = 0, logfc.threshold = 0)
rna.markers.8 <- FindMarkers(data, assay = 'ADT', ident.1 = 'SB775393_Act_IL1B_IL23', 'SB775393_Act_IL1B_IL23_TGFB', min.pct = 0, logfc.threshold = 0)
rna.markers.9 <- FindMarkers(data, assay = 'ADT', ident.1 = 'SB775393_Act_IL1B_IL23', 'SB775393_Act_IL1B_IL23_PGE2_TGFB', min.pct = 0, logfc.threshold = 0)
rna.markers.10 <- FindMarkers(data, assay = 'ADT', ident.1 = 'SB775393_Act_IL1B_IL23_PGE2', 'SB775393_Act_IL1B_IL23_TGFB', min.pct = 0, logfc.threshold = 0)
rna.markers.11 <- FindMarkers(data, assay = 'ADT', ident.1 = 'SB775393_Act_IL1B_IL23_PGE2', 'SB775393_Act_IL1B_IL23_PGE2_TGFB', min.pct = 0, logfc.threshold = 0)
rna.markers.12 <- FindMarkers(data, assay = 'ADT', ident.1 = 'SB775393_Act_IL1B_IL23_TGFB', 'SB775393_Act_IL1B_IL23_PGE2_TGFB', min.pct = 0, logfc.threshold = 0)

CITE <- list(rna.markers.1, rna.markers.2, rna.markers.3, rna.markers.4, rna.markers.5, rna.markers.6, rna.markers.7, rna.markers.8, rna.markers.9, rna.markers.10, rna.markers.11, rna.markers.12)
rm(data)

save(DIG, CITE, file = '../output/ADT_markers.RData')

load('../output/ADT_markers.RData')

library(BuenColors)
library(plotly)
library(ggrepel)

comparison <- function(index, name){
  tmp1 <- DIG[[index]]
  tmp1$gene <- str_split_fixed(str_split_fixed(rownames(tmp1), '-A0', 2)[,1], '-A1', 2)[,1]
  tmp2 <- CITE[[index]]
  tmp2$gene <- str_split_fixed(rownames(tmp2), '\\.1', 2)[,1]
  tmp <- merge(tmp1, tmp2, by = 'gene', suffixes = c('_DIG', '_LLL'))
  
  p1 <- ggplot(tmp, aes(x = avg_log2FC_LLL, y = avg_log2FC_DIG)) +
    geom_point() +
    pretty_plot(fontsize = 10) + L_border() + 
    labs(x = "CITE log2FC", y = "DIG log2FC")
  cowplot::ggsave2(p1, file = paste0("../plots/ADT_log2FC_", name, ".pdf"), width = 4, height = 4)
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
  tmp1$gene <- str_split_fixed(str_split_fixed(rownames(tmp1), '-A0', 2)[,1], '-A1', 2)[,1]
  tmp2 <- CITE[[index]]
  tmp2$gene <- str_split_fixed(str_split_fixed(rownames(tmp2), '-A0', 2)[,1], '-A1', 2)[,1]
  tmp <- merge(tmp1, tmp2, by = 'gene', suffixes = c('_DIG', '_LLL'))
  
  p1 <- ggplot(tmp, aes(x = avg_log2FC_LLL, y = avg_log2FC_DIG)) +
    geom_point() +
    geom_text_repel(aes(label=ifelse(gene %in% c('CD127', 'CD62L', 'CD63', 'CD25', 'CD71', 'CD69', 'CD54'),gene,''), color=ifelse(avg_log2FC_DIG > 0, 'positive', 'negative')))+
    pretty_plot(fontsize = 10) + L_border() + 
    labs(x = "CITE log2FC", y = "DIG log2FC", color = "") +
    scale_color_manual(values = c("dodgerblue3", "firebrick")) +
    theme(legend.position = "none")
  cowplot::ggsave2(p1, file = paste0("../plots/ADT_log2FC_", name, "_label.pdf"), width = 4, height = 4)
}

comparison_label(1, 'SB775372_Act_IL1B_IL23_vs_Act_IL1B_IL23_PGE2')
comparison_label(7, 'SB775393_Act_IL1B_IL23_vs_Act_IL1B_IL23_PGE2')
