library(Seurat)
library(Signac)
library(dplyr)
library(data.table)
library(stringr)
library(harmony)
library(celda)
library(ROGUE)
library(parallel)
library(cowplot)
library(ggrepel)
library(patchwork)

setwd('~/RWorkSpace/CITE-seq/Duerr/DOGMA-seq/DIG_CITE_rerun_1/code/')

draw_plot <- function(path1, path2){
  dig <- readRDS(path1)
  lll <- readRDS(path2)
  dig$tech <- 'DIG'
  lll$tech <- 'CITE'
  data_plot <- rbind(dig, lll)
  data_plot <- data_plot[!is.na(data_plot$resolution),]
  data_plot$tech <- factor(data_plot$tech, levels = c('CITE', 'DIG'))
  
  p <- ggplot(data_plot, aes(x = resolution, y = ave.rogue, col = tech))+
    geom_point()+
    geom_line()+
    geom_text_repel(max.overlaps = Inf, aes(label = cluster), min.segment.length = 0)+
    scale_color_manual(values = c("dodgerblue3", "firebrick")) +
    theme_cowplot()+
    labs(x = 'Clustering Resolution', y = 'Average ROGUE', col = '')+
    xlim(0,0.3)
  return(p) 
}

rna <- draw_plot('../output/dig_rna_rogue.RDS',
                 '../output/cite_rna_rogue.RDS')+ggtitle('RNA')+theme(legend.position = 'none', plot.title = element_text(hjust = 0.5))
adt <- draw_plot('../output/dig_adt_rogue.RDS',
                 '../output/cite_adt_rogue.RDS')+ggtitle('ADT')+theme(plot.title = element_text(hjust = 0.5))
wnn <- draw_plot('../output/dig_wnn_rogue.RDS',
                 '../output/cite_wnn_rogue.RDS')+ggtitle('3WNN/2WNN')+theme(legend.position = 'none', plot.title = element_text(hjust = 0.5))
pdf('../plots/rogue.pdf', width = 12, height = 4)
wrap_plots(wnn, rna, adt, ncol = 3)
dev.off()

rna <- draw_plot('../output/dig_rna_rogue.RDS',
                 '../output/cite_rna_rogue_wt_intron.RDS')+ggtitle('RNA')+theme(legend.position = 'none', plot.title = element_text(hjust = 0.5))
adt <- draw_plot('../output/dig_adt_rogue.RDS',
                 '../output/cite_adt_rogue_wt_intron.RDS')+ggtitle('ADT')+theme(plot.title = element_text(hjust = 0.5))
wnn <- draw_plot('../output/dig_wnn_rogue.RDS',
                 '../output/cite_wnn_rogue_wt_intron.RDS')+ggtitle('3WNN/2WNN')+theme(legend.position = 'none', plot.title = element_text(hjust = 0.5))
pdf('../plots/rogue_wt_intron.pdf', width = 12, height = 4)
wrap_plots(wnn, rna, adt, ncol = 3)
dev.off()

draw_plot_new <- function(path1, path2){
  dig <- readRDS(path1)
  lll <- readRDS(path2)
  dig$tech <- 'DIG'
  lll$tech <- 'CITE'
  data_plot <- rbind(dig, lll)
  data_plot <- data_plot[!is.na(data_plot$resolution),]
  data_plot <- data_plot[data_plot$resolution <= 0.3,]
  data_plot$tech <- factor(data_plot$tech, levels = c('CITE', 'DIG'))
  
  p <- ggplot(data_plot, aes(x = cluster, y = ave.rogue, col = tech))+
    geom_point()+
    geom_line()+
    # geom_text_repel(max.overlaps = Inf, aes(label = cluster), min.segment.length = 0)+
    scale_color_manual(values = c("dodgerblue3", "firebrick")) +
    theme_cowplot()+
    labs(x = 'Number of Clusters', y = 'Average ROGUE', col = '')+
    scale_x_continuous(breaks=seq(min(data_plot$cluster), max(data_plot$cluster), 1))
  return(p) 
}

rna <- draw_plot_new('../output/dig_rna_rogue.RDS',
                 '../output/cite_rna_rogue.RDS')+ggtitle('RNA')+theme(legend.position = 'none', plot.title = element_text(hjust = 0.5))
adt <- draw_plot_new('../output/dig_adt_rogue.RDS',
                 '../output/cite_adt_rogue.RDS')+ggtitle('ADT')+theme(plot.title = element_text(hjust = 0.5))
wnn <- draw_plot_new('../output/dig_wnn_rogue.RDS',
                 '../output/cite_wnn_rogue.RDS')+ggtitle('3WNN/2WNN')+theme(legend.position = 'none', plot.title = element_text(hjust = 0.5))
pdf('../plots/rogue_cluster.pdf', width = 12, height = 4)
wrap_plots(wnn, rna, adt, ncol = 3)
dev.off()

rna <- draw_plot_new('../output/dig_rna_rogue.RDS',
                 '../output/cite_rna_rogue_wt_intron.RDS')+ggtitle('RNA')+theme(legend.position = 'none', plot.title = element_text(hjust = 0.5))
adt <- draw_plot_new('../output/dig_adt_rogue.RDS',
                 '../output/cite_adt_rogue_wt_intron.RDS')+ggtitle('ADT')+theme(plot.title = element_text(hjust = 0.5))
wnn <- draw_plot_new('../output/dig_wnn_rogue.RDS',
                 '../output/cite_wnn_rogue_wt_intron.RDS')+ggtitle('3WNN/2WNN')+theme(legend.position = 'none', plot.title = element_text(hjust = 0.5))
pdf('../plots/rogue_wt_intron_cluster.pdf', width = 12, height = 4)
wrap_plots(wnn, rna, adt, ncol = 3)
dev.off()
