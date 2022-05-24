library(Seurat)
library(Signac)
library(dplyr)
library(data.table)
library(stringr)
library(harmony)

setwd('~/CITE-seq/Duerr/DOGMA-seq/DIG_CITE_rerun_1/code/')

dig <- read.csv('../plots/harmony/dig_predicted.csv')
cite <- read.csv('../plots/harmony/cite_predicted.csv')

DIG <- table(dig$predicted.celltype.l1)
DIG <- prop.table(DIG)
CITE <- table(cite$predicted.celltype.l1)
CITE <- prop.table(CITE)

data_plot <- as.data.frame(cbind(DIG, CITE))
data_plot$Type <- rownames(data_plot)
data_plot <- reshape2::melt(data_plot)
colnames(data_plot) <- c('Type', 'Tech', 'Proportions')
data_plot$Proportions <- 100 * data_plot$Proportions

library(ggsci)
library(ggpubr)
library(BuenColors)
pdf('../plots/harmony/proportion.pdf', width = 4, height = 4)
ggbarplot(data_plot, x = "Type", y = "Proportions",
          fill = "Tech",
          width = 0.8,
          position = position_dodge(0.8),
          color = "white",
          label = round(data_plot$Proportion, 1),
          lab.vjust = 0.5,
          lab.hjust = 0.5
) + coord_flip() + theme(axis.text.y = element_text(hjust = 0)) +
  scale_fill_manual(values = c("dodgerblue3", "firebrick")) +
  labs(x = "", y = "Proportion (%)", fill = "") 
dev.off()

DIG <- table(dig$predicted.celltype.l1, dig$condition)
DIG <- reshape2::melt(prop.table(DIG, 2))
DIG$Tech <- 'DIG'
CITE <- table(cite$predicted.celltype.l1, cite$condition)
CITE <- reshape2::melt(prop.table(CITE, 2))
CITE$Tech <- 'CITE'

data_plot <- as.data.frame(rbind(DIG, CITE))
colnames(data_plot) <- c('Type', 'Conditions', 'Proportions', 'Tech')
data_plot$Proportions <- 100 * data_plot$Proportions

pdf('../plots/harmony/proportion_split.pdf', width = 8, height = 8)
ggbarplot(data_plot, x = "Type", y = "Proportions",
          fill = "Conditions",
          width = 0.8,
          position = position_dodge(0.8),
          color = "Tech",
          lab.vjust = 0.5,
          lab.hjust = 0.5
) + coord_flip() + theme(axis.text.y = element_text(hjust = 0)) +
  labs(x = "", y = "Proportion (%)", fill = "", color = "") +
  scale_fill_npg() + 
  scale_color_manual(values = c("dodgerblue3", "firebrick"))+
  theme(legend.position = 'right')
dev.off()
