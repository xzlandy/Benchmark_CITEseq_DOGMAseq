library(Seurat)
library(BuenColors)
library(data.table)
library(dplyr)

setwd('~/RWorkSpace/CITE-seq/Duerr/DOGMA-seq/DIG_CITE_rerun_1/code/')

load("../data/DIG_exon_data.RData")
dig <- data
load('../output/exon_rate.RData')

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
table(rownames(df) == names(exon_rate))
df$Exon <- exon_rate
df$Exon_bi <- ifelse(df$Exon > 0.5, 'Exon-dominated Genes', 'Intron-dominated Genes')

p1 <- ggplot(df, aes(x = CITE, y = DIG, color = Type)) +
  geom_point() +
  pretty_plot(fontsize = 10) + L_border() + 
  scale_color_manual(values = c("darkgrey", "dodgerblue3", "firebrick")) +
  labs(x = "CITE Fraction > 0", y = "DIG Fraction > 0", color = "") +
  theme(legend.position = "bottom", legend.direction = 'horizontal')
cowplot::ggsave2(p1, file = "../plots/RNA_detection_freq_exon_only.pdf", width = 4, height = 4)

p1 <- ggplot(df, aes(x = CITE, y = DIG, color = Exon)) +
  geom_point() +
  pretty_plot(fontsize = 10) + L_border() + 
  scale_color_gradient(low = "dodgerblue3", high = "firebrick") +
  labs(x = "CITE Fraction > 0", y = "DIG Fraction > 0", color = "Proportion of Exonic UMIs") +
  theme(legend.position = "bottom", legend.direction = 'horizontal')
cowplot::ggsave2(p1, file = "../plots/RNA_detection_freq_exon_rate_exon_only.pdf", width = 4, height = 4)

p1 <- ggplot(df, aes(x = CITE, y = DIG, color = Exon_bi)) +
  geom_point() +
  pretty_plot(fontsize = 10) + L_border() + 
  scale_color_manual(values = c("dodgerblue3", "firebrick")) +
  labs(x = "CITE Fraction > 0", y = "DIG Fraction > 0", color = "") +
  theme(legend.position = "bottom", legend.direction = 'horizontal', legend.text = element_text(size = 7))
cowplot::ggsave2(p1, file = "../plots/RNA_detection_freq_exon_rate_bi_exon_only.pdf", width = 4, height = 4)
