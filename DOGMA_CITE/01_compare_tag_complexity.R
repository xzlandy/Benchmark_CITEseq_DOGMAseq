library(Seurat)
library(Signac)
library(dplyr)
library(data.table)
library(stringr)

setwd('~/RWorkSpace/CITE-seq/Duerr/DOGMA-seq/DIG_CITE_rerun_1/code/')

# Assign the bus set to the corresponding 
process_tag_reference <- function(ref_file){
  
  ref_df <- fread(ref_file, header = FALSE)
  ref_df2 <- data.frame(
    seq = ref_df[[1]][c(FALSE,TRUE)],
    name = gsub(">", "", ref_df[[1]][c(TRUE,FALSE)])
  )
  ref_df2$idx <- as.character(1:dim(ref_df2)[1])
  ref_df2$name_short <- str_split_fixed(ref_df2$name, "-", 2)[,1]
  
  # Some of these Abs have an extra -# that is necessary
  doubles <- data.frame(table(ref_df2$name_short)) %>% filter(Freq > 50) %>% pull(Var1) %>% as.character()
  ref_df2$extra_idx <- str_split_fixed(ref_df2$name, "-", 3)[,2]
  ref_df2$use <- ifelse(ref_df2$name_short %in% doubles, paste0(ref_df2$name_short,"_", ref_df2$extra_idx),ref_df2$name_short)
  
  vec <- as.character(ref_df2$use); names(vec) <- as.character(ref_df2$idx)
  vec
}

# Import PBMC references
pbmc_vec <- process_tag_reference("../data/FeaturesMismatch.fa")

source('../../global_functions/estimateLibraryComplexity.R')

# Function to process the per-barcode summary statistics
estimate_tag_rates <- function(bus_file, bcs, vec){
  bcs <- substr(bcs, 1, 16)
  dff <- fread(bus_file, col.names = c("barcode", "umi", "set", "idx")) %>%
    mutate(Ab = vec[as.character(set)], umi_ab = paste0(Ab,"_",umi)) %>%
    filter(barcode %in% bcs) %>% 
    group_by(barcode) %>%
    summarise(count = n_distinct(umi_ab), n_total = n()) %>%
    mutate(n_duplicated = n_total - count) 
  dff$complexity <- sapply(1:(dim(dff)[1]), function(i){
    estimateLibrarySize(dff$n_total[i],dff$count[i])
  })
  dff
}

# Pull HQ barcodes from embedding 
load('../data/barcodes.RData')

# Compute per-barcode statistics
DIG <- estimate_tag_rates("../data/DIG_mo.bustools.tsv", dig, pbmc_vec)
CITE <- estimate_tag_rates("../data/CITE_mo.bustools.tsv", cite, pbmc_vec)

save(DIG, CITE, file = "../output/complexity_matrices.rda")


if(FALSE){

  load("../output/complexity_matrices.rda")

  compute_duplicate_percentage <- function(df){
    df %>% filter(complexity < 100000) %>%
      mutate(percent_duplicates = n_duplicated/n_total*100) %>%
      pull(percent_duplicates) %>% summary
  }
  compute_duplicate_percentage(DIG)
  compute_duplicate_percentage(CITE)

  DIG$assay <- "DIG"
  CITE$assay <- "CITE"

  all_df <- rbind(DIG, CITE)

  library(ggplot2)
  library(BuenColors)
  p1 <- ggplot(rbind(DIG, CITE), aes(x = assay, y = log10(complexity), color = assay)) +
    geom_boxplot(outlier.shape = NA) +
    pretty_plot(fontsize = 7) + L_border() + scale_y_continuous(limits = c(2, 5)) +
    scale_color_manual(values = c("dodgerblue3", "firebrick")) +
    labs(x = "Condition", y = "log10 Tag Complexity", color = "") +
    theme(legend.position = "none")
  cowplot::ggsave2(p1, file = "../plots/tag_complexity_assay_MO.pdf", width = 1, height = 1.5)

  all_df %>% group_by(assay) %>% summarize(median(complexity), count = n())
}

# # A tibble: 2 × 3
# assay `median(complexity)` count
# <chr>                <dbl> <int>
# 1 CITE                  4112 12599
# 2 DIG                   4576 14172

if(FALSE){

  load("../output/complexity_matrices.rda")
  load("../data/DIG_data.RData")
  dig <- data
  load("../data/CITE_data.RData")
  cite <- data
  rm(data)

  compute_duplicate_percentage <- function(df){
    df %>% filter(complexity < 100000) %>%
      mutate(percent_duplicates = n_duplicated/n_total*100) %>%
      pull(percent_duplicates) %>% summary
  }
  compute_duplicate_percentage(DIG)
  compute_duplicate_percentage(CITE)

  DIG$assay <- "DIG"
  CITE$assay <- "CITE"
  DIG$condition <- dig$condition
  CITE$condition <- cite$condition

  all_df <- rbind(DIG, CITE)
  all_df$sample <- str_split_fixed(all_df$condition, '_Act', 2)[,1]

  library(ggplot2)
  library(BuenColors)
  p1 <- ggplot(rbind(DIG, CITE), aes(x = assay, y = log10(complexity), color = condition)) +
    geom_boxplot(outlier.shape = NA) +
    pretty_plot(fontsize = 7) + L_border() + scale_y_continuous(limits = c(2, 5)) +
    labs(x = "Condition", y = "log10 Tag Complexity", color = "")
  cowplot::ggsave2(p1, file = "../plots/tag_complexity_assay_MO_split.pdf", width = 3, height = 1.5)

  p1 <- ggplot(all_df, aes(x = assay, y = log10(complexity), color = sample)) +
    geom_boxplot(outlier.shape = NA) +
    pretty_plot(fontsize = 7) + L_border() + scale_y_continuous(limits = c(2, 5)) +
    scale_color_manual(values = c("dodgerblue3", "firebrick")) +
    labs(x = "Condition", y = "log10 Tag Complexity", color = "") +
    theme(legend.position = "bottom", legend.direction = 'horizontal')
  cowplot::ggsave2(p1, file = "../plots/tag_complexity_assay_MO_split_sample.pdf", width = 2, height = 2)

  all_df %>% group_by(assay, sample) %>% summarize(median(complexity), count = n())
}

# # A tibble: 4 × 4
# # Groups:   assay [2]
# assay sample   `median(complexity)` count
# <chr> <chr>                   <dbl> <int>
# 1 CITE  SB775372                 4245  6902
# 2 CITE  SB775393                 3970  5697
# 3 DIG   SB775372                 4756  7943
# 4 DIG   SB775393                 4381  6229
