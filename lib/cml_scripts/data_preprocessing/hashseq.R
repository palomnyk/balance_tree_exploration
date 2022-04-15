# Author: Aaron Yerke (aaronyerke@gmail.com)
# Script for running hashseq https://github.com/FarnazFouladi/HashSeq

rm(list = ls()) #clear workspace

options(java.parameters="-Xmx1024m")

print("Installing dependencies.")
# install.packages("devtools")
if (!requireNamespace("HashSeq", quietly = TRUE)){
  devtools::install_github("FarnazFouladi/HashSeq")
}
library("HashSeq")
if (!requireNamespace("optparse", quietly = TRUE)){
  install.packages("optparse")
}
library("optparse")

print("Loading cml options")
option_list <- list(
  make_option(c("-i", "--input_dir"), type="character", 
              default=file.path('~','git','balance_tree_exploration'), 
              help="dataset dir path", metavar="home dir"),
  make_option(c("-o", "--output_dir"), type="character", default=NULL, 
              help="output_dir", metavar="project"),
  make_option(c("-t", "--threshhold"), type="character", default=1000,
              help="hashseq threshold")); 

opt_parser <- optparse::OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);

print(opt)

hash_in <- path.expand(opt$input_dir)
file.exists(hash_in)
hash_out <- path.expand(opt$output_dir)
file.exists(hash_out)

print("Running hashseq")

HashSeq::inferTrueSequences(inputDir = hash_in, 
                            outputDir = hash_out,
                            abundanceThreshold = opt$threshold)

print("R Script complete.")

