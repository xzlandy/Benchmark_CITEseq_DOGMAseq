library(Seurat)
library(Signac)
library(dplyr)
library(data.table)
library(stringr)

setwd('~/CITE-seq/Duerr/DOGMA-seq/DIG_LLL/code/')

load('../data/DIG_data.RData')
dig <- as.character(colnames(data))
load('../data/LLL_data.RData')
lll <- as.character(colnames(data))
rm(data)

save(dig, lll, file = '../data/barcodes.RData')