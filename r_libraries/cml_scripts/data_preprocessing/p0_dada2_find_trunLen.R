#!/usr/bin/env Rscript
# This script generates images that can be used to finetune the length of 
# It expects the dir layout to be as though you ran "create new subproject_dirs.R"
# the truncLen parameter of the filterAndTrim DADA2 function.

# rm(list = ls()) #clear workspace

##-cml argument processing------------------------------------------##
if (!requireNamespace("optparse", quietly = TRUE)){
  install.packages("optparse")
}
library("optparse")

option_list <- list(
  make_option(c("-d", "--homedir"), type="character", 
              default=file.path('~','git','balance_tree_exploration'), 
              help="dataset dir path", metavar="character"),
  make_option(c("-p", "--project"), type="character", default=NULL, 
              help="project folder", metavar="character"),
  make_option(c("-f", "--r1_pattern"), type="character", default="_R1.fastq.gz", 
              help="pattern for forward filenames for paired-end reads", metavar="character"),
  make_option(c("-r", "--r2_pattern"), type="character", default="_R2.fastq.gz", 
              help="pattern for reverse filenames for paired-end reads", metavar="character")
); 

opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);

print(opt)

##-Load Depencencies------------------------------------------------##
library(dada2); packageVersion("dada2")

##-Establish directory layout---------------------------------------##
home_dir <- opt$homedir
project <- opt$project
#home_dir <- file.path('cloud','project')
output_dir <- file.path(home_dir, project, 'output')
f_path <- file.path(home_dir, project, "downloaded_seqs") # CHANGE ME to the directory containing the fastq files after unzipping.list.files(f_path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(f_path, pattern=opt$r1_pattern, full.names = TRUE))
fnRs <- sort(list.files(f_path, pattern=opt$r2_pattern, full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

plotQualF <- plotQualityProfile(fnFs[1:2])

plotQualR <- plotQualityProfile(fnRs[1:2])
dir.create(file.path(output_dir, "graphics"), showWarnings = FALSE)
setwd(file.path(output_dir, "graphics"))

png(filename="plotQualF.png")
plot(plotQualF)
dev.off()

png(filename="plotQualR.png")
plot(plotQualR)
dev.off()