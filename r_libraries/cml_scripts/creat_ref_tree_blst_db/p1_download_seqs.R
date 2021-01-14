# Author: Aaron Yerke
# Setting up silva tree for usage with ps

if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
  BiocManager::install("ape")
}
library('ape')

home_dir <- file.path('~','git','balance_tree_exploration')
output_dir <- file.path(home_dir, "ref_tree_objs")

setwd(file.path(output_dir))

fastaFile <- "treeFasta.fasta"

tree <- read.tree(file.path(output_dir, "silva","viFy10M5J2nvIBpCLM-QMQ_newick.txt"))
for (i in 1:length(tree$tip.label)){
  lab = tree$tip.label[i]
  id = strsplit(lab, "_")[[1]][1]
  myGen = ape:::read.GenBank(id, as.character=T, species.names = FALSE)
  myDna = lapply(myGen, function(x) paste0(x, collapse = ''))
  cat(paste0(">", id, "\n", paste0(myDna), "\n"), file = fastaFile, append=TRUE)
  Sys.sleep(.45)
}

