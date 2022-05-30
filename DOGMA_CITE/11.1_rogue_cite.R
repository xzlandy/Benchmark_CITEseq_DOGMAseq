library(Seurat)
library(Signac)
library(dplyr)
library(data.table)
library(stringr)
library(harmony)
library(celda)
library(ROGUE)
library(parallel)
library(cowplot)
library(ggplot2)

setwd('~/RWorkSpace/CITE-seq/Duerr/DOGMA-seq/DIG_CITE_rerun_1/code/')

load('../data/CITE_data.RData')
DefaultAssay(data) <- "RNA"
data <- FindClusters(data, resolution = seq(0.02, 0.4, 0.02), graph.name = 'SCT_snn')

expr <- as.matrix(data@assays$SCT@counts)
ave.rogue.res <- mclapply(1:20, function(i){
  name <- colnames(data@meta.data)[grep('SCT_snn_res.', colnames(data@meta.data))][i]
  rogue.res <- rogue(expr, labels = data@meta.data[,name], samples = data$orig.ident, platform = "UMI", span = 0.6)
  ave <- c(name, mean(as.matrix(rogue.res)[1,]), length(unique(data@meta.data[,name])))
  return(ave)
}, mc.cores = 20)
ave.rogue.res.df <- as.data.frame(matrix(unlist(ave.rogue.res), ncol = 3, byrow = T))
ave.rogue.res.df$resolution <- as.numeric(str_remove(ave.rogue.res.df$V1, 'SCT_snn_res.'))
colnames(ave.rogue.res.df)[2:3] <- c('ave.rogue', 'cluster')
ave.rogue.res.df$ave.rogue <- as.numeric(ave.rogue.res.df$ave.rogue)
ave.rogue.res.df$cluster <- as.numeric(ave.rogue.res.df$cluster)

p1 <- ggplot(ave.rogue.res.df, aes(x = resolution, y = ave.rogue))+
  geom_point()+
  geom_line()+
  theme_cowplot()+
  labs(x = 'Clustering Resolution', y = 'Average ROGUE')

p2 <- ggplot(ave.rogue.res.df, aes(x = cluster, y = ave.rogue))+
  geom_point()+
  geom_line()+
  theme_cowplot()+
  labs(x = 'Number of Clusters', y = 'Average ROGUE')
p1 | p2

saveRDS(ave.rogue.res.df, file = '../output/cite_rna_rogue.RDS')

DefaultAssay(data) <- "ADT"
data <- FindClusters(data, resolution = seq(0.02, 0.4, 0.02), graph.name = 'ADT_snn')

ave.rogue.res <- mclapply(1:20, function(i){
  name <- colnames(data@meta.data)[grep('ADT_snn_res.', colnames(data@meta.data))][i]
  rogue.res <- rogue(expr, labels = data@meta.data[,name], samples = data$orig.ident, platform = "UMI", span = 0.6)
  ave <- c(name, mean(as.matrix(rogue.res)[1,]), length(unique(data@meta.data[,name])))
  return(ave)
}, mc.cores = 20)
ave.rogue.res.df <- as.data.frame(matrix(unlist(ave.rogue.res), ncol = 3, byrow = T))
ave.rogue.res.df$resolution <- as.numeric(str_remove(ave.rogue.res.df$V1, 'ADT_snn_res.'))
colnames(ave.rogue.res.df)[2:3] <- c('ave.rogue', 'cluster')
ave.rogue.res.df$ave.rogue <- as.numeric(ave.rogue.res.df$ave.rogue)
ave.rogue.res.df$cluster <- as.numeric(ave.rogue.res.df$cluster)

p1 <- ggplot(ave.rogue.res.df, aes(x = resolution, y = ave.rogue))+
  geom_point(na.rm = T)+
  geom_line(na.rm = T)+
  theme_cowplot()+
  labs(x = 'Clustering Resolution', y = 'Average ROGUE')

p2 <- ggplot(ave.rogue.res.df, aes(x = cluster, y = ave.rogue))+
  geom_point(na.rm = T)+
  geom_line(na.rm = T)+
  theme_cowplot()+
  labs(x = 'Number of Clusters', y = 'Average ROGUE')
p1 | p2

saveRDS(ave.rogue.res.df, file = '../output/cite_adt_rogue.RDS')

# DefaultAssay(data) <- "peaks"
# data <- FindClusters(data, resolution = seq(0.02, 0.4, 0.02), graph.name = 'peaks_snn', algorithm = 3)
# 
# ave.rogue.res <- mclapply(1:20, function(i){
#   name <- colnames(data@meta.data)[grep('peaks_snn_res.', colnames(data@meta.data))][i]
#   rogue.res <- rogue(expr, labels = data@meta.data[,name], samples = data$orig.ident, platform = "UMI", span = 0.6)
#   ave <- c(name, mean(as.matrix(rogue.res)[1,]), length(unique(data@meta.data[,name])))
#   return(ave)
# }, mc.cores = 20)
# ave.rogue.res.df <- as.data.frame(matrix(unlist(ave.rogue.res), ncol = 3, byrow = T))
# ave.rogue.res.df$resolution <- as.numeric(str_remove(ave.rogue.res.df$V1, 'peaks_snn_res.'))
# colnames(ave.rogue.res.df)[2:3] <- c('ave.rogue', 'cluster')
# ave.rogue.res.df$ave.rogue <- as.numeric(ave.rogue.res.df$ave.rogue)
# ave.rogue.res.df$cluster <- as.numeric(ave.rogue.res.df$cluster)
# 
# p1 <- ggplot(ave.rogue.res.df, aes(x = resolution, y = ave.rogue))+
#   geom_point(na.rm = T)+
#   geom_line(na.rm = T)+
#   theme_cowplot()+
#   labs(x = 'Clustering Resolution', y = 'Average ROGUE')
# 
# p2 <- ggplot(ave.rogue.res.df, aes(x = cluster, y = ave.rogue))+
#   geom_point(na.rm = T)+
#   geom_line(na.rm = T)+
#   theme_cowplot()+
#   labs(x = 'Number of Clusters', y = 'Average ROGUE')
# p1 | p2
# 
# saveRDS(ave.rogue.res.df, file = '../output/cite_atac_rogue.RDS')

data <- FindClusters(data, graph.name = "wsnn", algorithm = 3, resolution = seq(0.02, 0.4, 0.02))

ave.rogue.res <- mclapply(1:20, function(i){
  name <- colnames(data@meta.data)[grep('wsnn_res.', colnames(data@meta.data))][i]
  rogue.res <- rogue(expr, labels = data@meta.data[,name], samples = data$orig.ident, platform = "UMI", span = 0.6)
  ave <- c(name, mean(as.matrix(rogue.res)[1,]), length(unique(data@meta.data[,name])))
  return(ave)
}, mc.cores = 20)
ave.rogue.res.df <- as.data.frame(matrix(unlist(ave.rogue.res), ncol = 3, byrow = T))
ave.rogue.res.df$resolution <- as.numeric(str_remove(ave.rogue.res.df$V1, 'wsnn_res.'))
colnames(ave.rogue.res.df)[2:3] <- c('ave.rogue', 'cluster')
ave.rogue.res.df$ave.rogue <- as.numeric(ave.rogue.res.df$ave.rogue)
ave.rogue.res.df$cluster <- as.numeric(ave.rogue.res.df$cluster)

p1 <- ggplot(ave.rogue.res.df, aes(x = resolution, y = ave.rogue))+
  geom_point(na.rm = T)+
  geom_line(na.rm = T)+
  theme_cowplot()+
  labs(x = 'Clustering Resolution', y = 'Average ROGUE')

p2 <- ggplot(ave.rogue.res.df, aes(x = cluster, y = ave.rogue))+
  geom_point(na.rm = T)+
  geom_line(na.rm = T)+
  theme_cowplot()+
  labs(x = 'Number of Clusters', y = 'Average ROGUE')
p1 | p2

saveRDS(ave.rogue.res.df, file = '../output/cite_wnn_rogue.RDS')
write.csv(data@meta.data, file = '../output/cite_rogue_meta.csv')