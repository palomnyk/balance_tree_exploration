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
library("phyloseq")
library("philr")
library("ggplot2")
library("ape")


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

print("Cleaning Ref tree otu with philr tutorial normalization")
ref_ps <- readRDS(file.path(output_dir, "r_objects", "ref_tree_phyloseq_obj.rds"))

clean_otu <- data.frame(ref_ps@otu_table@.Data)
clean_otu <- philr_tutorial_normalization(clean_otu)
print(paste("nrow orginal ref:", nrow(ref_ps@otu_table), "nrow clean ref: ", nrow(clean_otu)))

phy_tree(ref_ps) <- ape::makeNodeLabel(phy_tree(ref_ps), method="number", prefix='n')
ref_ps_clean <- phyloseq::phyloseq( otu_table(clean_otu, taxa_are_rows = F), 
                                    phy_tree(ref_ps@phy_tree),
                                    tax_table(ref_ps@tax_table), 
                                    sample_data(ref_ps@sam_data))

print("Cleaning UPGMA tree otu with philr tutorial normalization")
denovo_tree_ps <- readRDS(file.path(output_dir, "r_objects", "denovo_tree_UPGMA_phyloseq_obj.rds"))
clean_den_otu <- philr_tutorial_normalization(data.frame(denovo_tree_ps@otu_table@.Data))
print(paste("nrow orginal denovo:", nrow(denovo_tree_ps@otu_table), "nrow clean denovo otu: ", nrow(clean_den_otu)))
cln_denovo_tree_ps <- phyloseq::phyloseq( otu_table(clean_den_otu, taxa_are_rows = F),
                                          phy_tree(ape::makeNodeLabel(phy_tree(denovo_tree_ps@phy_tree))),
                                          tax_table(denovo_tree_ps@tax_table), 
                                          sample_data(denovo_tree_ps@sam_data))
denovo_tree_ps <- transform_sample_counts(denovo_tree_ps, function(x) x + 1 )
phy_tree(ref_ps_clean) <- makeNodeLabel(phy_tree(ref_ps_clean), method="number", prefix='n')
phy_tree(cln_denovo_tree_ps) <- makeNodeLabel(phy_tree(cln_denovo_tree_ps), method="number", prefix='n')

##-Random num seed--------------------------------------------------##
set.seed(36)
print("making random trees")
ref_rand_list <- list()
for (rand in 1:10){
  rand_tree <- ape::rtree(n = length(ref_ps@phy_tree$tip.label), tip.label = ref_ps@phy_tree$tip.label)
  #put int in philr
  rand_tree_ps <- phyloseq::phyloseq( otu_table(clean_otu, taxa_are_rows = F),
                                      phy_tree(rand_tree),
                                      tax_table(ref_ps@tax_table),
                                      sample_data(ref_ps@sam_data))
  phy_tree(rand_tree_ps) <- ape::makeNodeLabel(phy_tree(rand_tree_ps), method="number", prefix='n')
  ref_rand_list[[rand]] <- rand_tree_ps
  names(ref_rand_list)[rand] <- paste0("ref_rand", rand)
}

print("make random trees for cln upgma taxa")
cln_upgma_rand_list <- list()
for (rand in 1:10){
  rand_tree <- ape::rtree(n = length(cln_denovo_tree_ps@phy_tree$tip.label), tip.label = cln_denovo_tree_ps@phy_tree$tip.label)
  #put int in philr
  rand_tree_ps <- phyloseq::phyloseq( otu_table(cln_denovo_tree_ps, taxa_are_rows = F),
                                      phy_tree(rand_tree),
                                      tax_table(cln_denovo_tree_ps@tax_table),
                                      sample_data(cln_denovo_tree_ps@sam_data))
  phy_tree(rand_tree_ps) <- ape::makeNodeLabel(phy_tree(rand_tree_ps), method="number", prefix='n')
  cln_upgma_rand_list[[rand]] <- rand_tree_ps
  names(cln_upgma_rand_list)[rand] <- paste0("cln_upgma_rand", rand)
}

print("make random trees for original upgma taxa")
orig_upgma_rand_list <- list()
for (rand in 1:10){
  rand_tree <- ape::rtree(n = length(denovo_tree_ps@phy_tree$tip.label), tip.label = denovo_tree_ps@phy_tree$tip.label)
  #put int in philr
  rand_tree_ps <- phyloseq::phyloseq( otu_table(denovo_tree_ps, taxa_are_rows = F),
                                      phy_tree(rand_tree),
                                      tax_table(denovo_tree_ps@tax_table),
                                      sample_data(denovo_tree_ps@sam_data))
  phy_tree(rand_tree_ps) <- ape::makeNodeLabel(phy_tree(rand_tree_ps), method="number", prefix='n')
  orig_upgma_rand_list[[rand]] <- rand_tree_ps
  names(orig_upgma_rand_list)[rand] <- paste0("orig_upgma_rand", rand)
}
##-Calculate pvalues------------------------------------------------##
print("making philr tables")

#for making different philr weights
philr_taxa_weights <- c("uniform","gm.counts","anorm","anorm.x.gm.counts","enorm","enorm.x.gm.counts")
philr_ilr_weights <- c("uniform","blw","blw.sqrt","mean.descendants")

# ref_ps_clean_philr <- philr::philr(df = ref_ps_clean@otu_table, 
#                                    tree = ref_ps_clean@phy_tree,
#                                    )

my_philr <- philr::philr(df = ref_ps_clean@otu_table,
                         tree = ref_ps_clean@phy_tree)

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
  child_nodes <- as.character( edge_table$chi.name[parent_node_indices] )
  if (length(child_nodes) == 2 & all(startsWith( child_nodes, "n"))){
    philr_1 <- my_philr[,as.character(child_nodes[1])]
    philr_2 <- my_philr[,as.character(child_nodes[2])]
    pval <- t.test(philr_1, philr_2)$p.value
    pvals <- c(pvals, pval)
  }
}




