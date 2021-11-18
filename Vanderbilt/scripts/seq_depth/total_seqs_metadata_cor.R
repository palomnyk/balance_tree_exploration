# Author: Aaron Yerke (aaronyerke@gmail.com)
# Script for creating metadata/sequencing abundance pvalues.
# This was helepful: https://github.com/jsilve24/philr/blob/master/vignettes/philr-intro.Rmd#L142

rm(list = ls()) #clear workspace

##-functions--------------------------------------------------------##

##-Load Depencencies------------------------------------------------##

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
asv_table <- data.frame(readRDS(file.path(output_dir, "r_objects", "ForwardReads_DADA2.rds")))

metadata <- read.table(file.path(home_dir, project, "patient_metadata.tsv"), 
                       sep="\t", 
                       header=TRUE, 
                       row.names = "Run", 
                       check.names = FALSE,
                       stringsAsFactors=TRUE)
metadata <- metadata[row.names(asv_table),]
metadata <- metadata[vapply(metadata, function(x) length(unique(x)) > 1, logical(1L))]

total_seqs <- rowSums(asv_table)
total_seqs <- data.frame(total_seqs, row.names = row.names(asv_table))

num_cols <- unlist(lapply(metadata, !is.character))# Identify numeric columns
metadata_num <- metadata[ , num_cols] 

##-Metadata comparisions--------------------------------------------##
pvals <- c()
for(meta in 1:ncol(metadata)){
  my_meta <- metadata[,meta]
  my_lm <- lm(unlist(total_seqs) ~ my_meta)
  my_pval <- anova(my_lm)$"Pr(>F)"[1]
  pvals <- c(pvals, my_pval)
}

adj_pvals <- p.adjust(pvals, "BH")

anov_meta <- t(data.frame(pvals, adj_pvals, row.names = names(metadata)))
write.table(formatC(anov_meta, digits = 2, format = "e"), 
            file=file.path(output_dir, "tables", "metaVseqdep_anova.tsv"), 
            sep = "\t")


