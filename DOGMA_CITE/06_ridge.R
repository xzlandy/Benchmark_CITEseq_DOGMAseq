library(Seurat)
library(BuenColors)
library(data.table)
library(dplyr)
library(viridis)
library(stringr)

setwd('~/RWorkSpace/CITE-seq/Duerr/DOGMA-seq/DIG_CITE_rerun_1/code/')

load("../data/DIG_data.RData")
dig <- data
load("../data/CITE_data.RData")
cite <- data
rm(data)

dig$tech <- 'DIG'
cite$tech <- 'CITE'

dig$sample <- str_split_fixed(dig$condition, '_Act', 2)[,1]
cite$sample <- str_split_fixed(cite$condition, '_Act', 2)[,1]

DefaultAssay(dig) <- 'ADT'
DefaultAssay(cite) <- 'ADT'

table(rownames(cite@assays$ADT@data) == rownames(dig@assays$ADT@data))
# cite@assays$ADT@counts@Dimnames[[1]] <- rownames(dig)
# rownames(cite@assays$ADT@data) <- rownames(dig)
# rownames(cite@assays$ADT@scale.data) <- rownames(dig)

data <- merge(dig, cite, add.cell.ids = c('DIG', 'CITE'))
Idents(data) <- 'tech'
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

data_1 <- data[,data$sample == 'SB775372']
data_2 <- data[,data$sample == 'SB775393']

pdf('../plots/ridge_SB775372.pdf', width = 20, height = 20)
for (i in 1:6){
  print(RidgePlot(data_1, features = rownames(data_1)[(i-1)*25 + 1:25], ncol = 5))
}
print(RidgePlot(data_1, features = rownames(data_1)[151:163], ncol = 5))
dev.off()

pdf('../plots/ridge_SB775393.pdf', width = 20, height = 20)
for (i in 1:6){
  print(RidgePlot(data_2, features = rownames(data_2)[(i-1)*25 + 1:25], ncol = 5))
}
print(RidgePlot(data_2, features = rownames(data_2)[151:163], ncol = 5))
dev.off()

pdf('../plots/ridge_SB775372_log10.pdf', width = 20, height = 20)
for (i in 1:6){
  print(RidgePlot(data_1, features = rownames(data_1)[(i-1)*25 + 1:25], ncol = 5, slot = 'scale.data'))
}
print(RidgePlot(data_1, features = rownames(data_1)[151:163], ncol = 5, slot = 'scale.data'))
dev.off()

pdf('../plots/ridge_SB775393_log10.pdf', width = 20, height = 20)
for (i in 1:6){
  print(RidgePlot(data_2, features = rownames(data_2)[(i-1)*25 + 1:25], ncol = 5, slot = 'scale.data'))
}
print(RidgePlot(data_2, features = rownames(data_2)[151:163], ncol = 5, slot = 'scale.data'))
dev.off()
