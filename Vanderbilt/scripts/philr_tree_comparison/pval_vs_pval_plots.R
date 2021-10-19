# Author: Aaron Yerke (aaronyerke@gmail.com)
# Script for doing pvalue comparison of philr made from different trees
rm(list = ls()) #clear workspace

##-functions--------------------------------------------------------##

##-Load Depencencies------------------------------------------------##
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("ape", quietly = TRUE)) BiocManager::install("ape")
if (!requireNamespace("philr", quietly = TRUE)) BiocManager::install("philr")
if (!requireNamespace("phyloseq", quietly = TRUE)) BiocManager::install("phyloseq")
if (!requireNamespace("ggplot2", quietly = TRUE)) BiocManager::install("ggplot2")
if (!requireNamespace("geiger", quietly = TRUE)) install.packages("geiger")
library("phyloseq")
library("philr")
library("ggplot2")
library("ape")
library("geiger")


##-Establish directory layout---------------------------------------##
home_dir <- file.path('~','git','balance_tree_exploration')
project <- "Vanderbilt"
#home_dir <- file.path('cloud','project')
output_dir <- file.path(home_dir, project, 'output')
setwd(file.path(home_dir))

##-Functions--------------------------------------------------------##
source(file.path(home_dir, "r_libraries", "statistical_functions.R"))
source(file.path(home_dir, "r_libraries", "table_manipulations.R"))


##-Set up constants-------------------------------------------------##
rf_cols <- 3:7
num_cycles <- 20
if(num_cycles < 3) stop("num_cycles should be 3 or more")

##-Import tables and data preprocessing-----------------------------##
asv_table <- asv_table <- data.frame(readRDS(file.path(output_dir, "r_objects", "ForwardReads_DADA2.rds")))

#clean up otu tables
ref_ps <- readRDS(file.path(output_dir, "r_objects", "ref_tree_phyloseq_obj.rds"))
clean_otu <- data.frame(ref_ps@otu_table@.Data)
clean_otu <- philr_tutorial_normalization(clean_otu)
print(paste("nrow orginal ref:", nrow(ref_ps@otu_table), "nrow clean ref: ", nrow(clean_otu)))

phy_tree(ref_ps) <- ape::makeNodeLabel(phy_tree(ref_ps), method="number", prefix='n')
ref_ps_clean <- phyloseq::phyloseq( otu_table(clean_otu, taxa_are_rows = F), 
                                    phy_tree(ref_ps@phy_tree),
                                    tax_table(ref_ps@tax_table), 
                                    sample_data(ref_ps@sam_data))
denovo_tree_ps <- readRDS(file.path(output_dir, "r_objects", "denovo_tree_UPGMA_phyloseq_obj.rds"))
clean_den_otu <- philr_tutorial_normalization(data.frame(denovo_tree_ps@otu_table@.Data))
print(paste("nrow orginal denovo:", nrow(denovo_tree_ps@otu_table), "nrow clean denovo otu: ", nrow(clean_den_otu)))
cln_denovo_tree_ps <- phyloseq::phyloseq( otu_table(clean_den_otu, taxa_are_rows = F),
                                          phy_tree(ape::makeNodeLabel(phy_tree(denovo_tree_ps@phy_tree))),
                                          tax_table(denovo_tree_ps@tax_table), 
                                          sample_data(denovo_tree_ps@sam_data))
denovo_tree_ps <- transform_sample_counts(denovo_tree_ps, function(x) x + 1 )
# phy_tree(denovo_tree_ps) <- makeNodeLabel(phy_tree(denovo_tree_ps), method="number", prefix='n')
phy_tree(ref_ps_clean) <- makeNodeLabel(phy_tree(ref_ps_clean), method="number", prefix='n')
phy_tree(cln_denovo_tree_ps) <- makeNodeLabel(phy_tree(cln_denovo_tree_ps), method="number", prefix='n')

#make random tree
rand_tree <- rtree(n = length(ref_ps@phy_tree$tip.label), tip.label = ref_ps@phy_tree$tip.label)
#put int in philr
rand_tree_ps <- phyloseq::phyloseq( otu_table(clean_otu, taxa_are_rows = F),
                                    phy_tree(rand_tree),
                                    tax_table(ref_ps@tax_table),
                                    sample_data(ref_ps@sam_data))
phy_tree(cln_denovo_tree_ps) <- makeNodeLabel(phy_tree(cln_denovo_tree_ps), method="number", prefix='n')

##-Calculate pvalues------------------------------------------------##


#Pairwise Distances from a Phylogenetic Tree
test <- ape::dist.nodes(ref_ps_clean@phy_tree)

ape::is.rooted(ref_ps_clean@phy_tree)

data.tree::isr
# https://www.rdocumentation.org/packages/ape/versions/5.5/topics/nodepath

castor::find_root(ref_ps_clean@phy_tree)

# https://rdrr.io/cran/castor/man/get_subtrees_at_nodes.html
ref_philr <- philr::philr(ref_ps_clean@otu_table, ref_ps_clean@phy_tree)

# tips_per_node <- castor::count_tips_per_node(ref_ps_clean@phy_tree)

my_tree <- ref_ps_clean@phy_tree
select.tip.or.node <- function(element, tree) {
  # https://stackoverflow.com/questions/51696837/r-phylo-object-how-to-connect-node-label-and-node-number
  ifelse(element < Ntip(tree)+1, tree$tip.label[element], tree$node.label[element-Ntip(tree)])
}
edge_table <- data.frame(
  "parent" = my_tree$edge[,1],
  "par.name" = sapply(my_tree$edge[,1], select.tip.or.node, tree = my_tree),
  "child" = my_tree$edge[,2],
  "chi.name" = sapply(my_tree$edge[,2], select.tip.or.node, tree = my_tree)
)# 1 through Ntips is tip and Ntips through Ntips + Nnode is node
tree_map <- ape::prop.part(ref_ps_clean@phy_tree)

pvals <- c()

start_edge <- ape::Ntip(my_tree) + 1
end_edge <- ape::Nnode(my_tree) + ape::Ntip(my_tree)

for (node_num in start_edge:end_edge){
  parent_node_indices <- which(edge_table$parent == node_num)
  child_node_edges <- edge_table$child[parent_node_indices]
  tips_1 <- geiger::tips(my_tree, child_node_edges[1])
  tips_2 <- geiger::tips(my_tree, child_node_edges[2])
  # otu_counts <- ref_ps@otu_table[,tips_1]
  my_cor <- cor( ref_ps@otu_table[,tips_1], ref_ps@otu_table[,tips_2])
  my_cor <- cor.test( ref_ps@otu_table[,tips_1], ref_ps@otu_table[,tips_2])
  # my_lm <- lm(ref_ps@otu_table[,tips_1] ~ ref_ps@otu_table[,tips_2])
  # my_lm <- lm(matrix(ref_ps@otu_table[,tips_1]) ~ matrix(ref_ps@otu_table[,tips_2]))
  # pval <- anova(my_lm)$`Pr(>F)`
  # pval <- t.test(ref_ps@otu_table[,tips_1] ~ ref_ps@otu_table[,tips_2])
  pvals <- c(pvals, pval)
}
geiger::tips(ref_ps_clean@phy_tree, 5)

n8 <- "AGAGTTTGATCCTGGCTCAGATTGAACGCTGGCGGCATGCCTTACACATGCAAGTCGAACGGCAGCATGATCTAGCTTGCTAGATTGATGGCGAGTGGCGAACGGGTGAGTAATACATCGGAACGTGCCCTGTAGTGGGGGATAACTAGTCGAAAGATTAGCTAATACCGCATACGACCTGAGGGTGAAAGTGGGGGACCGCAAGGCCTCATGCTATAGGAGCGGCCGATGTCTGATTAG"

my_tree$node.label[151 - Ntip(my_tree)]

77 - 151

for (node_num in 1:ref_ps_clean@phy_tree$Nnode){
  my_tips <- geiger::tips(ref_ps_clean@phy_tree, 8)
  if (my_tips == n8) print(node_num)
}

if (all(startsWith(child_node_names, "n"))){
  print("made it here")
}

comparePhylo(ref_ps_clean@phy_tree, cln_denovo_tree_ps@phy_tree, plot = FALSE, force.rooted = FALSE,
             use.edge.length = FALSE)

comparePhylo(ref_ps_clean@phy_tree, rand_tree_ps@phy_tree, plot = FALSE, force.rooted = FALSE,
             use.edge.length = FALSE)

comparePhylo(cln_denovo_tree_ps@phy_tree, rand_tree_ps@phy_tree, plot = FALSE, force.rooted = FALSE,
             use.edge.length = FALSE)

my_sub_trees <- ape::subtrees(ref_ps_clean@phy_tree)

castor::find_root(my_sub_trees[[1]]) #77 what?

my_sub_trees[[1]]$tip.label

main_root <- castor::find_root(ref_ps_clean@phy_tree)

ape::summary.phylo(my_tree)

ref_ps1 <- readRDS(file.path(output_dir, "r_objects", "ref_tree_phyloseq_obj.rds"))
my_tree1 <- ref_ps1@phy_tree

ape::summary.phylo(my_tree1)

ape::Ntip(my_tree)
ape::Nnode(my_tree)


