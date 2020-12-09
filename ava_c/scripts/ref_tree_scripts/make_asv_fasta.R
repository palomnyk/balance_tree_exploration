#creating fasta for blasting from DADA2 results
home_dir = file.path('~','git','balance_tree_exploration')
#home_dir = file.path('cloud','project')
output_dir = file.path(home_dir, "ava_c", 'output')

con <- gzfile(file.path( output_dir, "r_objects", "ForwardReads_DADA2.rds"))
seqtab = readRDS(con)

setwd(file.path( output_dir, "tree_process_blast"))

for (i in 1:ncol(seqtab)){
  seq = colnames(seqtab)[i]
  cat(paste0(">", i, "\n", seq, "\n"), 
    file = "dada2seqs.fasta", 
    append=TRUE)
}