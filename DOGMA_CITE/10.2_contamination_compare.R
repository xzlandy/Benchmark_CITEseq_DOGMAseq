library(Seurat)
library(Signac)
library(dplyr)
library(data.table)
library(stringr)
library(harmony)
library(BuenColors)

setwd('~/CITE-seq/Duerr/DOGMA-seq/DIG_CITE_rerun_1/code/')

dig <- read.csv('../output/dig_harmony_predicted_chrom_contamination.csv')
cite <- read.csv('../output/cite_harmony_predicted_chrom_contamination.csv')

DIG <- data.frame(contamination = dig$contamination, tech = 'DIG')
CITE <- data.frame(contamination = cite$contamination, tech = 'CITE')

data_plot <- rbind(DIG, CITE)

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
cite <- read.csv('../output/cite_harmony_predicted_chrom_contamination_original.csv')

DIG <- data.frame(contamination = dig$contamination, tech = 'DIG')
CITE <- data.frame(contamination = cite$contamination, tech = 'CITE')

data_plot <- rbind(DIG, CITE)

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
cite <- read.csv('../output/cite_harmony_predicted_chrom_contamination_null.csv')

DIG <- data.frame(contamination = dig$contamination, tech = 'DIG')
CITE <- data.frame(contamination = cite$contamination, tech = 'CITE')

data_plot <- rbind(DIG, CITE)

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
