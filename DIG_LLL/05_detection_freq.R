library(Seurat)
library(BuenColors)
library(data.table)
library(dplyr)

setwd('~/RWorkSpace/CITE-seq/Duerr/DOGMA-seq/DIG_LLL/code/')

load("../data/DIG_data.RData")
dig <- data
load("../data/LLL_data.RData")
lll <- data
rm(data)

DIG <- apply(dig@assays$RNA@counts, 1, function(x){sum(ifelse(x > 0, 1, 0))/length(x)})
LLL <- apply(lll@assays$RNA@counts, 1, function(x){sum(ifelse(x > 0, 1, 0))/length(x)})
table(names(DIG) == names(LLL))

df <- as.data.frame(cbind(DIG, LLL))
df$Type <- 'Other Genes'
df$Type[grep('^MT-', rownames(df))] <- 'MT- Genes'
df$Type[grep('^RPL|^RPS', rownames(df))] <- 'RPL/RPS Genes'
df$Type <- factor(df$Type, levels = c('Other Genes', 'RPL/RPS Genes', 'MT- Genes'))

p1 <- ggplot(df, aes(x = LLL, y = DIG, color = Type)) +
  geom_point() +
  pretty_plot(fontsize = 10) + L_border() + 
  scale_color_manual(values = c("darkgrey", "dodgerblue3", "firebrick")) +
  labs(x = "LLL Fraction > 0", y = "DIG Fraction > 0", color = "") +
  theme(legend.position = "bottom", legend.direction = 'horizontal')
cowplot::ggsave2(p1, file = "../plots/RNA_detection_freq.pdf", width = 4, height = 4)

DIG <- apply(dig@assays$RNA@counts, 1, function(x){sum(ifelse(x > 1, 1, 0))/length(x)})
LLL <- apply(lll@assays$RNA@counts, 1, function(x){sum(ifelse(x > 1, 1, 0))/length(x)})
table(names(DIG) == names(LLL))

df <- as.data.frame(cbind(DIG, LLL))
df$Type <- 'Other Genes'
df$Type[grep('^MT-', rownames(df))] <- 'MT- Genes'
df$Type[grep('^RPL|^RPS', rownames(df))] <- 'RPL/RPS Genes'
df$Type <- factor(df$Type, levels = c('Other Genes', 'RPL/RPS Genes', 'MT- Genes'))

p1 <- ggplot(df, aes(x = LLL, y = DIG, color = Type)) +
  geom_point() +
  pretty_plot(fontsize = 10) + L_border() + 
  scale_color_manual(values = c("darkgrey", "dodgerblue3", "firebrick")) +
  labs(x = "LLL Fraction > 1", y = "DIG Fraction > 1", color = "") +
  theme(legend.position = "bottom", legend.direction = 'horizontal')
cowplot::ggsave2(p1, file = "../plots/RNA_detection_freq_larger_than_1.pdf", width = 4, height = 4)

DIG <- apply(dig@assays$ADT@counts, 1, function(x){sum(ifelse(x > 0, 1, 0))/length(x)})
LLL <- apply(lll@assays$ADT@counts, 1, function(x){sum(ifelse(x > 0, 1, 0))/length(x)})
table(names(DIG) == names(LLL))

df_adt <- as.data.frame(cbind(DIG, LLL))
df_adt$Type <- 'Other Proteins'
df_adt$Type[grep('Ctrl', rownames(df_adt))] <- 'Ctrl Proteins'
df_adt$Type <- factor(df_adt$Type, levels = c('Other Proteins', 'Ctrl Proteins'))

p1 <- ggplot(df_adt, aes(x = LLL, y = DIG, color = Type)) +
  geom_point() +
  pretty_plot(fontsize = 10) + L_border() + 
  scale_color_manual(values = c("darkgrey", "dodgerblue3")) +
  labs(x = "LLL Fraction > 0", y = "DIG Fraction > 0", color = "") +
  theme(legend.position = "bottom", legend.direction = 'horizontal')
cowplot::ggsave2(p1, file = "../plots/ADT_detection_freq.pdf", width = 4, height = 4)
