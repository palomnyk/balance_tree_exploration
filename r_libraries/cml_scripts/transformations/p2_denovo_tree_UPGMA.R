#!/usr/bin/env Rscript
# Author: Aaron Yerke
# This is a pipeline that was created to build denovo trees.
# Many resources were used to create this pipeline:
#   Processing the sequences through DADA2, making trees with phangorn, and putting them into phyloseq objects:
#     This is by far the biggest source:
#     https://github.com/spholmes/F1000_workflow/blob/master/MicrobiomeWorkflow/MicrobiomeWorkflowII.Rmd
#     or (in article format)
#     https://f1000research.com/articles/5-1492/v2

##-cml argument processing------------------------------------------##
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
  optparse::make_option(c("-m", "--metadata"), type="character", default=NULL,
              help="metadata file path with filename", metavar="character"),
  optparse::make_option(c("-l", "--metadata_delim"), type="character", default=NULL,
              help="metadata file deliminator", metavar="character"),
  optparse::make_option(c("-r", "--metadata_rowname"), type="character", default=NULL,
              help="metadata file row to use for row names", metavar="character")
  ); 

opt_parser <- optparse::OptionParser(option_list=option_list);

opt <- optparse::parse_args(opt_parser);

print(opt)

# ‘ape’, ‘dplyr’, ‘reshape2’, ‘plyr’
# .cran_packages <- c("ggplot2", "gridExtra")
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager", type = "source", repos = "http://archive.linux.duke.edu/cran/")
  if (!requireNamespace("phangorn", quietly = TRUE))
    install.packages("phangorn",type = "source", repos = "http://archive.linux.duke.edu/cran/")
  BiocManager::install("phyloseq")
  BiocManager::install("DECIPHER")
}

library("DECIPHER")
library("phangorn")
library("phyloseq")

##----------------Establish directory layout------------------------##
home_dir <- opt$homedir
project <- opt$project
output_dir = file.path(home_dir, project, 'output')

# setwd(file.path(home_dir))

print("Established directory layout")

##---------------------Import R objects-----------------------------##
con <- gzfile(file.path( output_dir, "r_objects", "ForwardReads_DADA2.rds"))
seqtab = readRDS(con)
close(con)

con <- gzfile(file.path( output_dir, "r_objects", "ForwardReads_DADA2_alignment.rds"))
alignment <- readRDS(con)
close(con)

con <- gzfile(file.path( output_dir, "r_objects", "ForwardReads_DADA2_taxonomy.rds"))
taxTab <- readRDS(con)
close(con)

print("Imported R objects")

##-import tables----------------------------------------------------##
myMeta = read.table(opt$metadata,
                    sep=opt$metadata_delim,
                    header=TRUE,
                    row.names = opt$metadata_rowname,
                    check.names = FALSE,
                    stringsAsFactors=FALSE)

print("Imported tables")

##------------------------Build tree--------------------------------##
phangAlign <- phangorn::phyDat(as(alignment, "matrix"), type="DNA")
dm <- phangorn::dist.ml(phangAlign)#create distance matrix
treeNJ <- phangorn::upgma(dm) #make tree
fit <- phangorn::pml(treeNJ, data=phangAlign)#fit model
# fitGTR <- update(fit, k=4, inv=0.2)#fit model with updated parameters
# fitGTR <- phangorn::optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
#                     rearrangement = "stochastic", control = pml.control(trace = 0))
print("phangorn completed")

ps <- phyloseq::phyloseq(otu_table(seqtab, taxa_are_rows=FALSE),
               sample_data(myMeta),
               tax_table(taxTab),
               phy_tree(treeNJ))
# ps
print("Created ps")

#examine tree
library("ape")

pdf(file.path(output_dir, "graphics", paste0("upgma_denovo","_2", ".pdf")))

plot_tree(ps, "treeonly", nodeplotblank, ladderize="left")

plot_tree(ps, ladderize="left", color="host_phenotype")

dev.off()

saveRDS(ps, file.path(output_dir, "r_objects","denovo_tree_UPGMA_phyloseq_obj.rds"))
