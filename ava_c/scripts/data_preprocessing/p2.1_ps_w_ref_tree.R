# Author: Aaron Yerke
# Convert output from p1_dada2_rd1.R to silva based taxonomy. Loosely following :
#   http://benjjneb.github.io/dada2/training.html

if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
  BiocManager::install("dada2",type = "source", checkBuilt = TRUE)
  BiocManager::install("DECIPHER")
}

library("DECIPHER")
library("dada2")
library("phyloseq")

home_dir = file.path('~','git','balance_tree_exploration')
#home_dir = file.path('cloud','project')
output_dir = file.path(home_dir, "ava_c", 'output')

setwd(file.path(home_dir, "ava_c"))
myMeta = read.table(file.path("metadata_SRA.tsv"), 
                    sep="\t", 
                    header=TRUE, 
                    row.names = "sample_name", 
                    check.names = FALSE,
                    stringsAsFactors=FALSE)

##---------------------Import R objects-----------------------------##
con <- gzfile(file.path( output_dir, "r_objects", "ForwardReads_DADA2.rds"))
seqtab = readRDS(con)
close(con)

con <- gzfile(file.path( output_dir, "r_objects", "ForwardReads_DADA2_taxonomy.rds"))
taxTab <- readRDS(con)
close(con)

print("Imported R objects")

tree = read_tree(file.path(home_dir, "taxonomy" , "silva","viFy10M5J2nvIBpCLM-QMQ_newick.txt"))

pdf(file = file.path(output_dir, "graphics", "silva_ref_tree.pdf"))
plot_tree(tree, nodelabf=nodeplotblank, label.tips="taxa_names", ladderize="left")
dev.off()

tree_key = read.table(file.path(output_dir, "tree_process_blast", "parsed_output.csv"),
                      sep = ",",
                      row.names = 1,
                      header = T)

match_df <- data.frame(
                      seqs=character(),
                      num_seqs=integer(),
                      stringsAsFactors=FALSE)

unmatched_ids = c()

new_labels = c()

for (i in 1:length(tree$tip.label)){
  old_lab = tree$tip.label[i]
  id = unlist(strsplit(old_lab, "_"))[1]
  if (id %in% tree_key$sseqid){
    
    index = which(tree_key$sseqid == id)[1]
    new_lab = row.names(tree_key)[index]
    tree$tip.label[i] = new_lab
    
    if (id %in% row.names(match_df)){
      print("in if")
      seqs = match_df[id,"seqs"]
      count = match_df[id,"num_seqs"]
      new_seqs = paste(seqs, new_lab)
      match_df[id,"seqs"] = new_seqs
      match_df[id,"num_seqs"] = count + 1
      
    }else{
    
    
    new_row = data.frame(
      seqs=new_lab,
      num_seqs=1,
      stringsAsFactors=FALSE)
    row.names(new_row) = c(id)
    match_df = rbind(match_df, new_row)
    }
    
    

  }else{
    unmatched_ids = c(unmatched_ids, id)
  }
  # index = which(tree_key$sseqid == id)
  # print("not in")
}
print(paste("unmatched ids:", length(unmatched_ids)))
plot_tree(tree)

print(length(tree$tip.label))
  
tree = prune_taxa(row.names(tree_key), tree)

print(length(tree$tip.label))

pdf(file = file.path(output_dir, "graphics", "ava_c_ref_tree.pdf"))
plot_tree(tree, nodelabf=nodeplotblank, label.tips="taxa_names", ladderize="left")
dev.off()

ps <- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE),
               sample_data(myMeta),
               tax_table(taxTab),
               phy_tree(tree)
)

saveRDS(ps, file.path(output_dir, "r_objects","ref_tree_phyloseq_obj.rds"))
