#!/usr/bin/env Rscript
# Author: Aaron Yerke (aaronyerke@gmail.com)
# Script for making ratio table of OTUS
# This was helepful: https://github.com/jsilve24/philr/blob/master/vignettes/philr-intro.Rmd#L142

# rm(list = ls()) #clear workspace

if (!requireNamespace("optparse", quietly = TRUE)){
  install.packages("optparse")
}
library("optparse")

option_list <- list(
  optparse::make_option(c("-d", "--homedir"), type="character", 
                        default=file.path('~','git','balance_tree_exploration'), 
                        help="dataset dir path", metavar="character"),
  optparse::make_option(c("-p", "--project"), type="character", default=NULL, 
                        help="project folder", metavar="character"),
  optparse::make_option(c("-t", "--taxonomic_level"), type="numeric", default=NULL,
                        help="taxonimic level for making otu table", metavar="numeric")
); 

opt_parser <- optparse::OptionParser(option_list=option_list);

opt <- optparse::parse_args(opt_parser);

print(opt)

##-load other dependencies------------------------------------------##
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
  BiocManager::install("phyloseq")
}
library("phyloseq")
##-Establish directory layout---------------------------------------##
home_dir <- opt$homedir
project <- opt$project
output_dir <- file.path(home_dir, project, 'output')

##-Import R objects and data preprocessing--------------------------##
con <- gzfile(file.path( output_dir, "r_objects","ref_tree_ps_philr_transform.rds"))
ps <- readRDS(con)
close(con)

##-Create "OTU" table-----------------------------------------------##
no <- data.frame(ps@otu_table)
tax <- data.frame(ps@tax_table)
otu_labels <- vector(length = ncol(no), mode = "character")
#assign otus to the asvs
for (i in 1:ncol(no)){
  asv <- names(no)[i]
  otu_labels[i] <- tolower(as.character(tax[asv, opt$taxonomic_level]))#this line sets the taxonomic level
}
uniq_otus <- unique(otu_labels)
uniq_otus <- uniq_otus[!is.na(uniq_otus)]

# set up main data.frame (read.table allows for prenaming columns)
# otu_tab <- read.table(text = "",
#                           colClasses = rep(x = "numeric", length(uniq_otus)),
#                           col.names = uniq_otus)
# otu_tab <- list()
otu_tab <- data.frame(row.names = row.names(no))

for(otu1 in 1:length(uniq_otus)){
  my_otu_lab <- uniq_otus[otu1]
  my_otus <- otu_labels == my_otu_lab
  my_otus <- replace(my_otus, is.na(my_otus), F)
  #need to add vectors so that they are the same length as nrow(ps.philr)
  my_otus <- as.data.frame(no[,my_otus])
  my_otus <- unlist(rowSums(my_otus, na.rm=T))
  otu_tab[,otu1] <- my_otus
}
colnames(otu_tab) <- uniq_otus

write.table(otu_tab, 
            file = file.path(output_dir, "tables", "otu_table.csv"),
            sep = ",")

saveRDS(otu_tab, file.path(output_dir, "r_objects", "otu_table.rds"))
