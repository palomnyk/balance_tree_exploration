# Author: Aaron Yerke (aaronyerke@gmail.com)
# Script for testing each transformation against different sequence depth to find the best one
# This was helepful: https://github.com/jsilve24/philr/blob/master/vignettes/philr-intro.Rmd#L142

rm(list = ls()) #clear workspace

##-Load Depencencies------------------------------------------------##
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("vegan", quietly = TRUE)) BiocManager::install("vegan")
library("vegan")
library("ggplot2")

##-Establish directory layout---------------------------------------##
home_dir <- file.path('~','git','balance_tree_exploration')
project <- "Vanderbilt"
#home_dir <- file.path('cloud','project')
output_dir <- file.path(home_dir, project, 'output')
setwd(file.path(home_dir))

##-functions--------------------------------------------------------##
find_step_size <- function(num_steps, max_dep, min_dep) {
  return( (max_dep - min_dep) / num_steps )
}

source(file.path(home_dir, "r_libraries", "statistical_functions.R"))
source(file.path(home_dir, "r_libraries", "table_manipulations.R"))

##-Import tables and data preprocessing-----------------------------##
asv_table <- readRDS(file.path(output_dir, "r_objects", "ForwardReads_DADA2.rds"))

##-Create new dataframe---------------------------------------------##

set.seed(54)

total_seqs <- rowSums(asv_table)
max_depth <- max(total_seqs)
min_depth <- 10#min(total_seqs)
max_row <- asv_table[match(max_depth, total_seqs), ]
num_steps <- 200
step_length <- find_step_size(num_steps, max_depth, min_depth)

new_df <- data.frame(max_row)

num_drop <- floor(ncol(asv_table)/num_steps)

for(s in 1:num_steps){
  new_col <- new_df[,ncol(new_df)]
  non_zeros <- as.integer(which(max_row != 0))
  new_zeros <- sample(non_zeros, num_drop)
  new_col[new_zeros] <- 0
  new_df[,ncol(new_df)+1] <- new_col
}

new_df <- t(new_df)

new_row_sums <- rowSums(new_df)
# new_row_sums <- data.frame(new_row_sums, row.names = row.names(new_df))

min_seq_depths <- c(0, 500, 1000, 5000, 10000, 20000, 40000)
mds_depth <- 5

mds_lev <- vector(mode = "integer", length = length(min_seq_depths) * mds_depth)
seq_depth <- vector(mode = "integer", length = length(min_seq_depths) * mds_depth)
var_exp <- vector(mode = "integer", length(min_seq_depths) * mds_depth)
spear_cor <- vector(mode = "integer", length(min_seq_depths) * mds_depth)
samples_left <- vector(mode = "integer", length(min_seq_depths) * mds_depth)


counter <- 1

for(s in 1:length(min_seq_depths)){
  seq_d <- min_seq_depths[s]#new sequencing depth
  my_table <- new_df[new_row_sums >= seq_d,]#dataset 1
  my_prcmp <- prcomp(my_table, 
                     center = TRUE,
                     rank = mds_depth)
  ##-Extract PCA matrix and convert to dataframe----------------------##
  myPCA <- data.frame(my_prcmp$x)
  my_var_exp <- my_prcmp$sdev^2/sum(my_prcmp$sdev^2)
  for (md in 1:mds_depth){
    mds_lev[counter] <- md
    seq_depth[counter] <- seq_d
    var_exp[counter] <- my_var_exp[md]
    spear_cor[counter] <- cor(new_row_sums[new_row_sums >= seq_d], myPCA[,md], method = "spearman")
    samples_left[counter] <- nrow(my_table)
    counter <- counter + 1
  }
}
result_df <- data.frame(mds_lev, seq_depth, spear_cor, var_exp, samples_left)
print("created resulting DF")

##-Create plots-----------------------------------------------------##
pdf(file = file.path(output_dir, "graphics", "sing_samp_rare_exp.pdf"))
for (i in 1:max(result_df$mds_lev)){
  pca_only <- result_df[mds_lev == i, ]
  g <- ggplot2::ggplot(pca_only, aes(x=seq_depth, y=spear_cor^2)) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::ggtitle(paste0(project, ": PCA",  i, " vs total sequences per sample")) +
    ggplot2::xlab("Min sequence depth per sample") +
    ggplot2::ylab("Spearman R^2") + 
    theme(axis.text.x = element_text(angle = 90)) +
    ggplot2::theme_minimal()
  print(g)
  
  g <- ggplot2::ggplot(pca_only, 
                       aes(x=seq_depth, y=samples_left)) +
    ggplot2::geom_point(aes()) +
    ggplot2::geom_line(aes()) +
    ggplot2::ggtitle(paste0(project, ": PCA",  i, " vs total sequences per sample")) +
    ggplot2::xlab("Min sequence depth per sample") +
    ggplot2::ylab("Samples") + 
    theme(axis.text.x = element_text(angle = 90)) +
    ggplot2::theme_minimal()
  g
  print(g)
}
dev.off()


