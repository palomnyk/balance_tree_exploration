#!/usr/bin/env Rscript
# Author: Aaron Yerke
# Filter out data2 results with less than X seq depth

##-Functions--------------------------------------------------------##
filt_seq_dpth <- function(min_seq_depth, df) {
  df <- df[rowSums(df) > min_seq_depth, ]
}

##-Establish directory layout---------------------------------------##
home_dir <- file.path('~','git','balance_tree_exploration')
project <- "Vanderbilt"
output_dir <- file.path(home_dir, project, 'output')

min_sd <- 1000
asv_table <- asv_table <- read.table(file.path(output_dir, "tables", "ForwardReads_DADA2.txt"),
                                     sep = "\t",
                                     header = TRUE)

filt_1000_seq_dpt_asv <- filt_seq_dpth(min_sd, asv_table)

write.table(filt_1000_seq_dpt_asv, 
            file = file.path(output_dir, "tables", "filt_1000_seq_dpt_asv.tsv"),
            sep = "/t")