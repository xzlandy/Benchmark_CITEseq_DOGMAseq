library(Seurat)
library(BuenColors)
library(data.table)
library(dplyr)
library(stringr)

setwd('~/CITE-seq/Duerr/DOGMA-seq/DIG_CITE_rerun_1/code/')

load("../data/DIG_data.RData")
dig <- data
load("../data/CITE_data.RData")
cite <- data
rm(data)

DIG <- apply(dig@assays$RNA@counts, 1, function(x){sum(ifelse(x > 0, 1, 0))/length(x)})
CITE <- apply(cite@assays$RNA@counts, 1, function(x){sum(ifelse(x > 0, 1, 0))/length(x)})
table(names(DIG) == names(CITE))

df <- as.data.frame(cbind(DIG, CITE))
df$Type <- 'Other Genes'
df$Type[grep('^MT-', rownames(df))] <- 'MT- Genes'
df$Type[grep('^RPL|^RPS', rownames(df))] <- 'RPL/RPS Genes'
df$Type <- factor(df$Type, levels = c('Other Genes', 'RPL/RPS Genes', 'MT- Genes'))

p1 <- ggplot(df, aes(x = CITE, y = DIG, color = Type)) +
  geom_point() +
  pretty_plot(fontsize = 10) + L_border() + 
  scale_color_manual(values = c("darkgrey", "dodgerblue3", "firebrick")) +
  labs(x = "CITE Fraction > 0", y = "DIG Fraction > 0", color = "") +
  theme(legend.position = "bottom", legend.direction = 'horizontal')
cowplot::ggsave2(p1, file = "../plots/RNA_detection_freq.pdf", width = 4, height = 4)

###############
DIG <- apply(dig@assays$RNA@counts, 1, function(x){sum(ifelse(x > 1, 1, 0))/length(x)})
CITE <- apply(cite@assays$RNA@counts, 1, function(x){sum(ifelse(x > 1, 1, 0))/length(x)})
table(names(DIG) == names(CITE))

df <- as.data.frame(cbind(DIG, CITE))
df$Type <- 'Other Genes'
df$Type[grep('^MT-', rownames(df))] <- 'MT- Genes'
df$Type[grep('^RPL|^RPS', rownames(df))] <- 'RPL/RPS Genes'
df$Type <- factor(df$Type, levels = c('Other Genes', 'RPL/RPS Genes', 'MT- Genes'))

p1 <- ggplot(df, aes(x = CITE, y = DIG, color = Type)) +
  geom_point() +
  pretty_plot(fontsize = 10) + L_border() + 
  scale_color_manual(values = c("darkgrey", "dodgerblue3", "firebrick")) +
  labs(x = "CITE Fraction > 1", y = "DIG Fraction > 1", color = "") +
  theme(legend.position = "bottom", legend.direction = 'horizontal')
cowplot::ggsave2(p1, file = "../plots/RNA_detection_freq_larger_than_1.pdf", width = 4, height = 4)

###############
DIG <- apply(dig@assays$ADT@counts, 1, function(x){sum(ifelse(x > 0, 1, 0))/length(x)})
CITE <- apply(cite@assays$ADT@counts, 1, function(x){sum(ifelse(x > 0, 1, 0))/length(x)})
table(names(DIG) == names(CITE))
# names(CITE) <- names(DIG)
# table(names(DIG) == names(CITE))

df_adt <- as.data.frame(cbind(DIG, CITE))
df_adt$Type <- 'Other Proteins'
df_adt$Type[grep('Ctrl', rownames(df_adt))] <- 'Ctrl Proteins'
df_adt$Type <- factor(df_adt$Type, levels = c('Other Proteins', 'Ctrl Proteins'))

p1 <- ggplot(df_adt, aes(x = CITE, y = DIG, color = Type)) +
  geom_point() +
  pretty_plot(fontsize = 10) + L_border() + 
  scale_color_manual(values = c("darkgrey", "dodgerblue3")) +
  labs(x = "CITE Fraction > 0", y = "DIG Fraction > 0", color = "") +
  theme(legend.position = "bottom", legend.direction = 'horizontal')
cowplot::ggsave2(p1, file = "../plots/ADT_detection_freq.pdf", width = 4, height = 4)

###############
DIG <- apply(dig@assays$ADT@counts, 1, function(x){sum(ifelse(x > 1, 1, 0))/length(x)})
CITE <- apply(cite@assays$ADT@counts, 1, function(x){sum(ifelse(x > 1, 1, 0))/length(x)})
table(names(DIG) == names(CITE))
# names(CITE) <- names(DIG)
# table(names(DIG) == names(CITE))

df_adt <- as.data.frame(cbind(DIG, CITE))
df_adt$Type <- 'Other Proteins'
df_adt$Type[grep('Ctrl', rownames(df_adt))] <- 'Ctrl Proteins'
df_adt$Type <- factor(df_adt$Type, levels = c('Other Proteins', 'Ctrl Proteins'))

p1 <- ggplot(df_adt, aes(x = CITE, y = DIG, color = Type)) +
  geom_point() +
  pretty_plot(fontsize = 10) + L_border() + 
  scale_color_manual(values = c("darkgrey", "dodgerblue3")) +
  labs(x = "CITE Fraction > 1", y = "DIG Fraction > 1", color = "") +
  theme(legend.position = "bottom", legend.direction = 'horizontal')
cowplot::ggsave2(p1, file = "../plots/ADT_detection_freq_larger_than_1.pdf", width = 4, height = 4)

###############
dig$sample <- str_split_fixed(dig$condition, '_Act', 2)[,1]
cite$sample <- str_split_fixed(cite$condition, '_Act', 2)[,1]

dig_1 <- dig[,dig$sample == 'SB775372']
dig_2 <- dig[,dig$sample == 'SB775393']
rm(dig)

cite_1 <- cite[,cite$sample == 'SB775372']
cite_2 <- cite[,cite$sample == 'SB775393']
rm(cite)

DIG_1 <- apply(dig_1@assays$RNA@counts, 1, function(x){sum(ifelse(x > 0, 1, 0))/length(x)})
CITE_1 <- apply(cite_1@assays$RNA@counts, 1, function(x){sum(ifelse(x > 0, 1, 0))/length(x)})
table(names(DIG_1) == names(CITE_1))

df <- as.data.frame(cbind(DIG_1, CITE_1))
df$Type <- 'Other Genes'
df$Type[grep('^MT-', rownames(df))] <- 'MT- Genes'
df$Type[grep('^RPL|^RPS', rownames(df))] <- 'RPL/RPS Genes'
df$Type <- factor(df$Type, levels = c('Other Genes', 'RPL/RPS Genes', 'MT- Genes'))

p1 <- ggplot(df, aes(x = CITE_1, y = DIG_1, color = Type)) +
  geom_point() +
  pretty_plot(fontsize = 10) + L_border() + 
  scale_color_manual(values = c("darkgrey", "dodgerblue3", "firebrick")) +
  labs(x = "CITE Fraction > 0", y = "DIG Fraction > 0", color = "") +
  theme(legend.position = "bottom", legend.direction = 'horizontal')
cowplot::ggsave2(p1, file = "../plots/RNA_detection_freq_SB775372.pdf", width = 4, height = 4)

DIG_1 <- apply(dig_1@assays$ADT@counts, 1, function(x){sum(ifelse(x > 0, 1, 0))/length(x)})
CITE_1 <- apply(cite_1@assays$ADT@counts, 1, function(x){sum(ifelse(x > 0, 1, 0))/length(x)})
table(names(DIG_1) == names(CITE_1))
# names(CITE_1) <- names(DIG_1)
# table(names(DIG_1) == names(CITE_1))

df_adt <- as.data.frame(cbind(DIG_1, CITE_1))
df_adt$Type <- 'Other Proteins'
df_adt$Type[grep('Ctrl', rownames(df_adt))] <- 'Ctrl Proteins'
df_adt$Type <- factor(df_adt$Type, levels = c('Other Proteins', 'Ctrl Proteins'))

p1 <- ggplot(df_adt, aes(x = CITE_1, y = DIG_1, color = Type)) +
  geom_point() +
  pretty_plot(fontsize = 10) + L_border() + 
  scale_color_manual(values = c("darkgrey", "dodgerblue3")) +
  labs(x = "CITE Fraction > 0", y = "DIG Fraction > 0", color = "") +
  theme(legend.position = "bottom", legend.direction = 'horizontal')
cowplot::ggsave2(p1, file = "../plots/ADT_detection_freq_SB775372.pdf", width = 4, height = 4)

###################
DIG_2 <- apply(dig_2@assays$RNA@counts, 1, function(x){sum(ifelse(x > 0, 1, 0))/length(x)})
CITE_2 <- apply(cite_2@assays$RNA@counts, 1, function(x){sum(ifelse(x > 0, 1, 0))/length(x)})
table(names(DIG_2) == names(CITE_2))

df <- as.data.frame(cbind(DIG_2, CITE_2))
df$Type <- 'Other Genes'
df$Type[grep('^MT-', rownames(df))] <- 'MT- Genes'
df$Type[grep('^RPL|^RPS', rownames(df))] <- 'RPL/RPS Genes'
df$Type <- factor(df$Type, levels = c('Other Genes', 'RPL/RPS Genes', 'MT- Genes'))

p1 <- ggplot(df, aes(x = CITE_2, y = DIG_2, color = Type)) +
  geom_point() +
  pretty_plot(fontsize = 10) + L_border() + 
  scale_color_manual(values = c("darkgrey", "dodgerblue3", "firebrick")) +
  labs(x = "CITE Fraction > 0", y = "DIG Fraction > 0", color = "") +
  theme(legend.position = "bottom", legend.direction = 'horizontal')
cowplot::ggsave2(p1, file = "../plots/RNA_detection_freq_SB775393.pdf", width = 4, height = 4)

DIG_2 <- apply(dig_2@assays$ADT@counts, 1, function(x){sum(ifelse(x > 0, 1, 0))/length(x)})
CITE_2 <- apply(cite_2@assays$ADT@counts, 1, function(x){sum(ifelse(x > 0, 1, 0))/length(x)})
table(names(DIG_2) == names(CITE_2))
# names(CITE_2) <- names(DIG_2)
# table(names(DIG_2) == names(CITE_2))

df_adt <- as.data.frame(cbind(DIG_2, CITE_2))
df_adt$Type <- 'Other Proteins'
df_adt$Type[grep('Ctrl', rownames(df_adt))] <- 'Ctrl Proteins'
df_adt$Type <- factor(df_adt$Type, levels = c('Other Proteins', 'Ctrl Proteins'))

p1 <- ggplot(df_adt, aes(x = CITE_2, y = DIG_2, color = Type)) +
  geom_point() +
  pretty_plot(fontsize = 10) + L_border() + 
  scale_color_manual(values = c("darkgrey", "dodgerblue3")) +
  labs(x = "CITE Fraction > 0", y = "DIG Fraction > 0", color = "") +
  theme(legend.position = "bottom", legend.direction = 'horizontal')
cowplot::ggsave2(p1, file = "../plots/ADT_detection_freq_SB775393.pdf", width = 4, height = 4)
