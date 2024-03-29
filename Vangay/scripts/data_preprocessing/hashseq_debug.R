# Author: Aaron Yerke (aaronyerke@gmail.com)
# Script for running hashseq https://github.com/FarnazFouladi/HashSeq

rm(list = ls()) #clear workspace

options(java.parameters="-Xmx500g")

print("Installing dependencies.")
# install.packages("devtools")
if (!requireNamespace("HashSeq", quietly = TRUE)){
  devtools::install_github("FarnazFouladi/HashSeq")
}

print("Running hashseq")

HashSeq::inferTrueSequences(inputDir = "/users/amyerke/git/balance_tree_exploration/Vangay/downloaded_seqs", 
                            outputDir = "/users/amyerke/git/balance_tree_exploration/Vangay/output/hashseq",
                            abundanceThreshold = 1000)

print("R Script complete.")

