library(Seurat)
library(Signac)
library(dplyr)
library(data.table)
library(stringr)
library(harmony)
library(BuenColors)

setwd('~/RWorkSpace/CITE-seq/Duerr/DOGMA-seq/DIG_LLL/code/')

dig <- read.csv('../output/dig_harmony_predicted_chrom_contamination.csv')
lll <- read.csv('../output/lll_harmony_predicted_chrom_contamination.csv')

DIG <- data.frame(contamination = dig$contamination, tech = 'DIG')
LLL <- data.frame(contamination = lll$contamination, tech = 'LLL')

data_plot <- rbind(DIG, LLL)

p <- ggplot(data_plot, aes(x = tech, y = contamination, color = tech)) +
  geom_boxplot(outlier.shape = NA) + 
  pretty_plot(fontsize = 7) + L_border() + 
  scale_color_manual(values = c("dodgerblue3", "firebrick", "darkgrey")) +
  labs(x = "Condition", y = "Contamination Rate (%)", color = "") +
  theme(legend.position = "none")

data_plot %>% group_by(tech) %>%
  summarize(contamination = median(contamination))

cowplot::ggsave2(cowplot::plot_grid(
  p, nrow = 1
), file = "../plots/contamination.pdf", width = 1, height = 1.5)

dig <- read.csv('../output/dig_harmony_predicted_chrom_contamination_original.csv')
lll <- read.csv('../output/lll_harmony_predicted_chrom_contamination_original.csv')

DIG <- data.frame(contamination = dig$contamination, tech = 'DIG')
LLL <- data.frame(contamination = lll$contamination, tech = 'LLL')

data_plot <- rbind(DIG, LLL)

p <- ggplot(data_plot, aes(x = tech, y = contamination, color = tech)) +
  geom_boxplot(outlier.shape = NA) + 
  pretty_plot(fontsize = 7) + L_border() + 
  scale_color_manual(values = c("dodgerblue3", "firebrick", "darkgrey")) +
  labs(x = "Condition", y = "Contamination Rate (%)", color = "") +
  theme(legend.position = "none")

data_plot %>% group_by(tech) %>%
  summarize(contamination = median(contamination))

cowplot::ggsave2(cowplot::plot_grid(
  p, nrow = 1
), file = "../plots/contamination_original.pdf", width = 1, height = 1.5)

dig <- read.csv('../output/dig_harmony_predicted_chrom_contamination_null.csv')
lll <- read.csv('../output/lll_harmony_predicted_chrom_contamination_null.csv')

DIG <- data.frame(contamination = dig$contamination, tech = 'DIG')
LLL <- data.frame(contamination = lll$contamination, tech = 'LLL')

data_plot <- rbind(DIG, LLL)

p <- ggplot(data_plot, aes(x = tech, y = contamination, color = tech)) +
  geom_boxplot(outlier.shape = NA) + 
  pretty_plot(fontsize = 7) + L_border() + 
  scale_color_manual(values = c("dodgerblue3", "firebrick", "darkgrey")) +
  labs(x = "Condition", y = "Contamination Rate (%)", color = "") +
  theme(legend.position = "none")

data_plot %>% group_by(tech) %>%
  summarize(contamination = median(contamination))

cowplot::ggsave2(cowplot::plot_grid(
  p, nrow = 1
), file = "../plots/contamination_null.pdf", width = 1, height = 1.5)

dig <- readRDS('../output/dig_decontx.RDS')
lll <- readRDS('../output/lll_decontx.RDS')

DIG <- apply(dig, 1, function(x){sum(ifelse(x > 0, 1, 0))/length(x)})
LLL <- apply(lll, 1, function(x){sum(ifelse(x > 0, 1, 0))/length(x)})
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
cowplot::ggsave2(p1, file = "../plots/RNA_detection_freq_decontx.pdf", width = 4, height = 4)

DIG <- apply(dig, 1, function(x){sum(ifelse(x > 1, 1, 0))/length(x)})
LLL <- apply(lll, 1, function(x){sum(ifelse(x > 1, 1, 0))/length(x)})
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
cowplot::ggsave2(p1, file = "../plots/RNA_detection_freq_decontx_larger_than_1.pdf", width = 4, height = 4)
