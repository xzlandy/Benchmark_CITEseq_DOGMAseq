library(Seurat)
library(BuenColors)
library(data.table)
library(dplyr)
library(viridis)

setwd('~/RWorkSpace/CITE-seq/Duerr/DOGMA-seq/DIG_LLL/code/')

load("../data/DIG_data.RData")
dig <- data
load("../data/LLL_data.RData")
lll <- data
rm(data)

dig$orig.ident <- 'DIG'
lll$orig.ident <- 'LLL'

p1 <- TSSPlot(dig, group.by = 'orig.ident')
p2 <- TSSPlot(lll, group.by = 'orig.ident')

TSS <- p1$data
TSS <- rbind(TSS, p2$data)

p1 <- ggplot(data = TSS, mapping = aes(x = position, y = norm.value, color = group))+
  geom_line(stat = "identity", size = 0.2)+
  labs(x = "Distance from TSS (bp)", y = "Mean TSS enrichment score", color = "") +
  pretty_plot(fontsize = 7) + L_border() + 
  scale_color_manual(values = c("dodgerblue3", "firebrick")) +
  ggtitle("TSS enrichment") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = 'none')
cowplot::ggsave2(p1, file = "../plots/TSS.pdf", width = 2, height = 2)

data <- merge(dig, lll, add.cell.ids = c('DIG', 'LLL'))
Idents(data) <- 'orig.ident'
DefaultAssay(data) <- 'ADT'

pdf('../plots/ridge.pdf', width = 20, height = 20)
for (i in 1:6){
  print(RidgePlot(data, features = rownames(data)[(i-1)*25 + 1:25], ncol = 5))
}
print(RidgePlot(data, features = rownames(data)[151:163], ncol = 5))
dev.off()

data@assays$ADT@scale.data <- as.matrix(log10(1 + data@assays$ADT@counts))

pdf('../plots/ridge_log10.pdf', width = 20, height = 20)
for (i in 1:6){
  print(RidgePlot(data, features = rownames(data)[(i-1)*25 + 1:25], ncol = 5, slot = 'scale.data'))
}
print(RidgePlot(data, features = rownames(data)[151:163], ncol = 5, slot = 'scale.data'))
dev.off()