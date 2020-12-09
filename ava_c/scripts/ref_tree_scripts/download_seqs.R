# Author: Aaron Yerke
# Setting up silva tree for usage with ps

if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
  BiocManager::install("ape")
}
library('ape')

home_dir = file.path('~','git','Western_gut')

# system2(command = "module load blast/2.9.0+")
# 
# newPath = "/apps/pkg/ncbi-blast/2.5.0+/rhel7_u5/gnu/bin"
# 
# myPath = system(command = "echo $PATH",
#                 wait = TRUE,
#                 intern=TRUE)
# 
# system2(command = "export",
#   args= c("PATH=", newPath, ":$PATH"))
# system2(command = "echo, args  $PATH",
#        wait = TRUE)
# 
# system(command = paste("blastn"),
#        wait = TRUE,)

print(getwd())

fastaFile = "treeFasta.fasta"

tree = read_tree(file.path(home_dir,"philr_pipelines",  "taxonomy" , "silva","viFy10M5J2nvIBpCLM-QMQ_newick.txt"))
# plot_tree(tree, nodelabf=nodeplotblank, label.tips="taxa_names", ladderize="left")
# tree_seqs = vector(length = length(tree$tip.label))
for (i in 1:length(tree$tip.label)){
  lab = tree$tip.label[i]
  id = strsplit(lab, "_")[[1]][1]
  myGen = ape:::read.GenBank(id, as.character=T, species.names = FALSE)
  myDna = lapply(myGen, function(x) paste0(x, collapse = ''))
  cat(paste0(">", lab, "\n", paste0(myDna), "\n"), file = fastaFile, append=TRUE)
  Sys.sleep(.4)
}
