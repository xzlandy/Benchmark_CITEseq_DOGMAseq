library(Seurat)
library(Signac)
library(dplyr)
library(data.table)
library(stringr)

setwd('~/RWorkSpace/CITE-seq/Duerr/DOGMA-seq/DIG_CITE_rerun_1/code/')

load('../data/DIG_data.RData')
dig <- as.character(colnames(data))
load('../data/CITE_data_wt_intron.RData')
cite <- as.character(colnames(data))
rm(data)

save(dig, cite, file = '../data/barcodes_wt_intron.RData')
