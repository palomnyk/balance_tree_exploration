# Author: Aaron Yerke (aaronyerke@gmail.com)
# Script for testing each transformation against sequence depth
# This was helepful: https://github.com/jsilve24/philr/blob/master/vignettes/philr-intro.Rmd#L142

rm(list = ls()) #clear workspace

##-functions--------------------------------------------------------##
pca_total_seqs_plot <- function(df, ts, main_tit, mds){
  ##-Create PCA-------------------------------------------------------##
  my_prcmp <- prcomp(df, 
                     center = TRUE)#,
                    # scale = TRUE)
  ##-Extract PCA matrix and convert to dataframe----------------------##
  myPCA <- data.frame(my_prcmp$x)
  my_spear <- cor.test(log10(ts), myPCA[,mds],
                       method = "kendall")
  plot(log10(ts), myPCA[,2],
       main = paste0( main_tit, " PCA", mds,"\nr_sq: ", my_spear$estimate^2),
       sub = paste0("Kendall Cor: ", my_spear$estimate))
  return(my_spear$estimate)
}

##-Load Depencencies------------------------------------------------##
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
asv_table <- read.table(file.path(output_dir, "tables", "ForwardReads_DADA2.txt"),
                        sep = "\t",
                        header = TRUE)

total_seqs <- rowSums(asv_table)

asv_table <- asv_table[total_seqs > 1000,]

philr_trans_ref <- read.table(file.path(output_dir, "tables", "phil_ref_tree_1000_sd_filtered.csv"),
                              sep = ",",
                              header = TRUE)

philr_trans_ref <- philr_trans_ref[row.names(asv_table), ]

philr_trans_upgma <- read.table(file.path(output_dir, "tables", "phil_denovo_tree_UPGMA_1000_sd_filt.csv"),
                                sep = ",",
                                header = TRUE)

philr_trans_upgma <- philr_trans_upgma[row.names(asv_table), ]

# total_seqs <- data.frame(total_seqs, row.names = row.names(asv_table))
# write.table(total_seqs,
#             file = file.path(output_dir, "tables", "total_seq_depth.csv"),
#             sep = ",")

otu_tab <- read.table(file.path(output_dir, "tables", "genus_otu_table.csv"),
                      sep = ",",
                      header = TRUE)

otu_tab <- otu_tab[row.names(asv_table), ]

ln_otu_tab <- lognorm(otu_tab)

# otu_ratio <- simple_ratio_transform(otu_tab)

# my_alr <- as.data.frame(compositions::alr(as.matrix(otu_tab), ivar=ncol(otu_tab)))

my_clr <- as.data.frame(compositions::clr(as.matrix(otu_tab)))

##-DESeq2 transform-------------------------------------------------##
t_otu_tab = t(otu_tab)

metadata <- read.table(file.path(home_dir, project, "patient_metadata.tsv"), 
                       sep="\t", 
                       header=TRUE, 
                       row.names = "Run", 
                       check.names = FALSE,
                       stringsAsFactors=TRUE)

metadata <- metadata[row.names(asv_table), ]

all(rownames(metadata) == colnames(t_otu_tab))

metadata <- metadata[colnames(t_otu_tab),]

dds <- DESeq2::DESeqDataSetFromMatrix(countData=t_otu_tab,
                                      colData=metadata, 
                                      design= ~ 1)

dds <- estimateSizeFactors(dds)

normlzd_dds <- counts(dds, normalized=T)

normlzd_dds <- t(data.frame(normlzd_dds))


# my_datasets <- list(asv_table, ln_otu_tab, my_alr, my_clr, philr_trans_ref, philr_trans_upgma, normlzd_dds )
# my_ds_names <- c("filtered ASV", "lognorm OTU", "filter alr", "filter clr", "ref tree philR", "UPGMA philr trans", "DESeq2")

my_datasets <- list(asv_table, ln_otu_tab, my_clr, philr_trans_ref, philr_trans_upgma, normlzd_dds )
my_ds_names <- c("filtered ASV", "lognorm OTU", "filter clr", "ref tree philR", "UPGMA philr trans", "DESeq2")


# my_datasets <- list( ln_otu_tab, my_alr, my_clr, philr_trans_ref, philr_trans_upgma )
# my_ds_names <- c( "lognorm OTU", "alr", "clr", "ref tree philR", "UPGMA philr trans")


# my_datasets <- list( ln_otu_tab, filtered_otu, philr_trans_ref, philr_trans_upgma )
# my_ds_names <- c( "lognorm OTU", "filtered OTU", "ref tree philR", "UPGMA philr trans")

#Plot shannon diversity against log10(total_seqs)
#Plot other normalization methods: otu log 100 and alr and clr 

pca1_kend <- vector(length = length(my_datasets))
pca2_kend <- vector(length = length(my_datasets))
pca3_kend <- vector(length = length(my_datasets))

pdf(file = file.path(output_dir, "graphics", "explore_seq_depth_artifact_seq_depth_1000_PCA123_c.pdf"))
for( ds in 1:length(my_datasets)){
  my_table <- as.data.frame(my_datasets[ds])
  # View(my_table)
  pca1_kend[ds] <- pca_total_seqs_plot(df = my_table, 
                      ts = total_seqs[total_seqs > 1000], 
                      main_tit = paste(project, my_ds_names[ds]),
                      mds = 1)
  pca2_kend[ds] <- pca_total_seqs_plot(df = my_table, 
                      ts = total_seqs[total_seqs > 1000], 
                      main_tit = paste(project, my_ds_names[ds]),
                      mds = 2)
  pca3_kend[ds] <- pca_total_seqs_plot(df = my_table, 
                      ts = total_seqs[total_seqs > 1000], 
                      main_tit = paste(project, my_ds_names[ds]),
                      mds = 3)
}
dev.off()

pca1_rsq <- pca1_kend^2
pca2_rsq <- pca2_kend^2
pca3_rsq <- pca3_kend^2


my_results <- data.frame(row.names = my_ds_names,
                         pca1_kend, pca1_rsq,
                         pca2_kend, pca2_rsq,
                         pca3_kend, pca3_rsq)

my_kend <- data.frame(row.names = my_ds_names,
                      pca1_kend,
                      pca2_kend,
                      pca3_kend,)

write.table(my_kend, 
            file = file.path(output_dir, "tables", "PCA123_vs_total_seqs_min_1000_kend.csv"),
            sep = ","
            )
