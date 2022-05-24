library(Seurat)
library(BuenColors)
library(data.table)
library(dplyr)
library(ggsci)

setwd('~/RWorkSpace/CITE-seq/Duerr/DOGMA-seq/DIG_LLL/code/')

load("../data/DIG_data.RData")
dig <- data

counts <- Read10X_h5('../../../Duerr_20210419_DOGMAseq_DIG_gex_exclude_introns/outs/raw_feature_bc_matrix.h5')
rna_counts <- counts$`Gene Expression`
atac_counts <- counts$Peaks
# create a Seurat object containing the RNA data
metadata <- read.csv(
  file = "../../../Duerr_20210419_DOGMAseq_DIG_gex_exclude_introns/outs/per_barcode_metrics.csv",
  header = TRUE,
  row.names = 1
)

dig_exon <- CreateSeuratObject(
  counts = rna_counts,
  assay = "RNA",
  meta.data = metadata
)

dig_exon <- dig_exon[,colnames(dig)]

DIG <- dig@assays$RNA@counts
DIG_exon <- dig_exon@assays$RNA@counts
DIG_exon <- DIG_exon[rownames(DIG),colnames(DIG)]

exon <- rowSums(DIG_exon)
exon_intron <- rowSums(DIG)
exon[exon > exon_intron] <- exon_intron[exon > exon_intron]
DIG_exon_rate <- exon/exon_intron

save(DIG_exon_rate, file = '../output/DIG_exon_rate.RData')

load("../data/LLL_data.RData")
lll <- data

counts <- Read10X_h5('../../../Duerr_20210419_DOGMAseq_PFA_LLL_gex_exclude_introns/outs/raw_feature_bc_matrix.h5')
rna_counts <- counts$`Gene Expression`
atac_counts <- counts$Peaks
# create a Seurat object containing the RNA data
metadata <- read.csv(
  file = "../../../Duerr_20210419_DOGMAseq_PFA_LLL_gex_exclude_introns/outs/per_barcode_metrics.csv",
  header = TRUE,
  row.names = 1
)

lll_exon <- CreateSeuratObject(
  counts = rna_counts,
  assay = "RNA",
  meta.data = metadata
)

lll_exon <- lll_exon[,colnames(lll)]

LLL <- lll@assays$RNA@counts
LLL_exon <- lll_exon@assays$RNA@counts
LLL_exon <- LLL_exon[rownames(LLL),colnames(LLL)]

exon <- rowSums(LLL_exon)
exon_intron <- rowSums(LLL)
exon[exon > exon_intron] <- exon_intron[exon > exon_intron]
LLL_exon_rate <- exon/exon_intron

save(LLL_exon_rate, file = '../output/LLL_exon_rate.RData')

load('../output/DIG_exon_rate.RData')
load('../output/LLL_exon_rate.RData')
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
table(rownames(df) == names(DIG_exon_rate))
table(rownames(df) == names(LLL_exon_rate))
df$DIG_Exon <- DIG_exon_rate
df$LLL_Exon <- LLL_exon_rate
df$DIG_Exon_bi <- ifelse(df$DIG_Exon > 0.5, 'Exon-dominated Genes', 'Intron-dominated Genes')
df$LLL_Exon_bi <- ifelse(df$LLL_Exon > 0.5, 'Exon-dominated Genes', 'Intron-dominated Genes')
df$Exon_bi <- NA
df[df$DIG_Exon_bi %in% 'Exon-dominated Genes' & df$LLL_Exon_bi %in% 'Exon-dominated Genes',]$Exon_bi <- 'Exon-dominated Genes (Both)'
df[df$DIG_Exon_bi %in% 'Exon-dominated Genes' & !df$LLL_Exon_bi %in% 'Exon-dominated Genes',]$Exon_bi <- 'Exon-dominated Genes (DIG)'
df[!df$DIG_Exon_bi %in% 'Exon-dominated Genes' & df$LLL_Exon_bi %in% 'Exon-dominated Genes',]$Exon_bi <- 'Exon-dominated Genes (LLL)'
df[df$DIG_Exon_bi %in% 'Intron-dominated Genes' & df$LLL_Exon_bi %in% 'Intron-dominated Genes',]$Exon_bi <- 'Intron-dominated Genes (Both)'
df$Exon_bi <- factor(df$Exon_bi, levels = c('Exon-dominated Genes (Both)', 'Exon-dominated Genes (DIG)', 'Exon-dominated Genes (LLL)', 'Intron-dominated Genes (Both)'))

dist2d <- function(a,b,c) {
  v1 <- b - c
  v2 <- a - b
  m <- cbind(v1,v2)
  d <- abs(det(m))/sqrt(sum(v1*v1))
  return(d)
} 

dist <- dist2d(c(df$LLL[1], df$DIG[1]), c(0,0), c(1,1))
for(i in 2:nrow(df)){
  dist <- c(dist, dist2d(c(df$LLL[i], df$DIG[i]), c(0,0), c(1,1)))
}
df$dist <- dist
df$diff <- df$LLL - df$DIG
write.xlsx(df, file = '../output/RNA_detection_freq_exon_rate.xlsx', colNames = T, rowNames = T)

p1 <- ggplot(df, aes(x = LLL, y = DIG, color = DIG_Exon_bi)) +
  geom_point() +
  pretty_plot(fontsize = 10) + L_border() + 
  scale_color_manual(values = c("dodgerblue3", "firebrick")) +
  labs(x = "LLL Fraction > 0", y = "DIG Fraction > 0", color = "") +
  theme(legend.position = "bottom", legend.direction = 'horizontal', legend.text = element_text(size = 7))
cowplot::ggsave2(p1, file = "../plots/RNA_detection_freq_DIG_exon_rate.pdf", width = 4, height = 4)

p1 <- ggplot(df, aes(x = LLL, y = DIG, color = LLL_Exon_bi)) +
  geom_point() +
  pretty_plot(fontsize = 10) + L_border() + 
  scale_color_manual(values = c("dodgerblue3", "firebrick")) +
  labs(x = "LLL Fraction > 0", y = "DIG Fraction > 0", color = "") +
  theme(legend.position = "bottom", legend.direction = 'horizontal', legend.text = element_text(size = 7))
cowplot::ggsave2(p1, file = "../plots/RNA_detection_freq_LLL_exon_rate.pdf", width = 4, height = 4)

p1 <- ggplot(df, aes(x = LLL, y = DIG, color = Exon_bi)) +
  geom_point() +
  pretty_plot(fontsize = 10) + L_border() + 
  scale_color_lancet()+
  labs(x = "LLL Fraction > 0", y = "DIG Fraction > 0", color = "") +
  theme(legend.position = "bottom", legend.direction = 'horizontal', legend.text = element_text(size = 7))+
  guides(color=guide_legend(nrow=2,byrow=F))
cowplot::ggsave2(p1, file = "../plots/RNA_detection_freq_exon_rate.pdf", width = 4, height = 4)
