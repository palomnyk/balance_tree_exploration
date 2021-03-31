# Author: Aaron Yerke (aaronyerke@gmail.com)
# Script for testing each transformation against sequence depth
# This was helepful: https://github.com/jsilve24/philr/blob/master/vignettes/philr-intro.Rmd#L142

rm(list = ls()) #clear workspace

##-------------------Load Depencencies------------------------------##
if (!requireNamespace("randomForest", quietly = TRUE)){
  install.packages("randomForest")
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer")
}
library(RColorBrewer)
#set color palette
palette( brewer.pal(7,"Accent") )
library(randomForest)
library(compositions)
library(vegan)
library("DESeq2")

##-Establish directory layout---------------------------------------##
home_dir <- file.path('~','git','balance_tree_exploration')
project <- "Vanderbilt"
#home_dir <- file.path('cloud','project')
output_dir <- file.path(home_dir, project, 'output')
setwd(file.path(home_dir))

##-Functions--------------------------------------------------------##
source(file.path(home_dir, "r_libraries", "statistical_functions.R"))
source(file.path(home_dir, "r_libraries", "table_manipulations.R"))

##-Import tables and data preprocessing-----------------------------##
philr_trans_ref <- read.table(file.path(output_dir, "tables", "ref_tree_ps_philr_transform.csv"),
                              sep = ",",
                              header = TRUE)

philr_trans_upgma <- read.table(file.path(output_dir, "tables", "philr_denovo_tree_UPGMA.csv"),
                                sep = ",",
                                header = TRUE)

otu_tab <- read.table(file.path(output_dir, "tables", "genus_otu_table.csv"),
                      sep = ",",
                      header = TRUE)

ln_otu_tab <- read.table(file.path(output_dir, "tables", "lognorm_genus_otu_table.csv"),
                         sep = ",",
                         header = TRUE)

filtered_otu <- philr_tutorial_normalization(otu_tab)

# otu_ratio <- simple_ratio_transform(otu_tab)

asv_table <- read.table(file.path(output_dir, "tables", "ForwardReads_DADA2.txt"),
                        sep = "\t",
                        header = TRUE)
# 
# my_alr <- as.data.frame(compositions::alr(as.matrix(asv_table)))
# 
# my_clr <- as.data.frame(compositions::clr(as.matrix(asv_table)))

total_seqs <- readRDS(file.path(output_dir, "r_objects", "raw_ASV_total_row_seqs.rds"))

metadata <- read.table(file.path(home_dir, project, "patient_metadata.tsv"), 
                       sep="\t", 
                       header=TRUE, 
                       row.names = "Run", 
                       check.names = FALSE,
                       stringsAsFactors=TRUE)
# print(paste("metadata dim:", dim(metadata)))
metadata <- metadata[row.names(metadata) %in% row.names(asv_table), ]
# print(paste("metadata dim:", dim(metadata)))
print(sapply(metadata,class))

# my_datasets <- list(asv_table, ln_otu_tab, my_alr, my_clr, philr_trans_ref, philr_trans_upgma )
# my_ds_names <- c("Pre-filtered ASV", "lognorm OTU", "Pre-filter alr", "Pre-filter clr", "ref tree philR", "UPGMA philr trans")

my_datasets <- list( ln_otu_tab, filtered_otu, philr_trans_ref, philr_trans_upgma )
my_ds_names <- c( "lognorm OTU", "filtered OTU", "ref tree philR", "UPGMA philr trans")

# total_seqs <- rowSums(asv_table)
# saveRDS(total_seqs, file.path(output_dir, "r_objects","raw_ASV_total_row_seqs.rds"))

#Plot shannon diversity against log10(total_seqs)
#Plot other normalization methods: otu log 100 and alr and clr 

pca_total_seqs_plot <- function(df, ts, main_tit, mds){
  ##-Create PCA-------------------------------------------------------##
  my_prcmp <- prcomp(df, 
                     center = TRUE,
                     scale = TRUE)
  ##-Extract PCA matrix and convert to dataframe----------------------##
  myPCA <- data.frame(my_prcmp$x)
  my_spear <- cor.test(log10(ts), myPCA[,mds],
                       method = "kendall")
  plot(log10(ts), myPCA[,2],
       main = paste0( main_tit, " PCA", mds,"\nr_sq: ", my_spear$estimate^2),
       sub = paste0("Kendall Cor: ", my_spear$estimate))
}

pdf(file = file.path(output_dir, "graphics", "explore_seq_depth_artifact_PCA23.pdf"))
for( ds in 1:length(my_datasets)){
  my_table <- as.data.frame(my_datasets[ds])
  # View(my_table)
  pca_total_seqs_plot(df = my_table, 
                      ts = total_seqs, 
                      main_tit = paste(project, my_ds_names[ds]),
                      mds = 2)
  pca_total_seqs_plot(df = my_table, 
                      ts = total_seqs, 
                      main_tit = paste(project, my_ds_names[ds]),
                      mds = 3)
}
dev.off()

# my_veg <- diversity(philr_trans_ref)

t_otu_tab = t(otu_tab)

all(rownames(metadata) == colnames(t_otu_tab))

metadata <- metadata[colnames(t_otu_tab),]

dds <- DESeq2::DESeqDataSetFromMatrix(countData=t_otu_tab,
                                      colData=metadata, 
                                      design= ~ 1)

dds <- estimateSizeFactors(dds)

normlzd_dds <- counts(dds, normalized=T)

normlzd_dds <- t(data.frame(normlzd_dds))

pca_total_seqs_plot(df = t(normlzd_dds), 
                    ts = total_seqs, 
                    main_tit = paste(project, "DESeq2"),
                    mds = 1)

pca_total_seqs_plot(df = t(normlzd_dds), 
                    ts = total_seqs, 
                    main_tit = paste(project, "DESeq2"),
                    mds = 2)

pca_total_seqs_plot(df = t(normlzd_dds), 
                    ts = total_seqs, 
                    main_tit = paste(project, "DESeq2"),
                    mds = 2)





