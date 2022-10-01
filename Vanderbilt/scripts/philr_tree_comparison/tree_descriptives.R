# Author: Aaron Yerke (aaronyerke@gmail.com)
# Script for making distribution of AUC of ILR balances to compare to philr

rm(list = ls()) #clear workspace

##-functions--------------------------------------------------------##

##-Load Depencencies------------------------------------------------##
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("ape", quietly = TRUE)) BiocManager::install("ape")
if (!requireNamespace("philr", quietly = TRUE)) BiocManager::install("philr")
if (!requireNamespace("randomForest", quietly = TRUE)) BiocManager::install("randomForest")
if (!requireNamespace("ROCR", quietly = TRUE)) BiocManager::install("ROCR")
if (!requireNamespace("ggpubr", quietly = TRUE)) BiocManager::install("ggpubr")
if (!requireNamespace("phyloseq", quietly = TRUE)) BiocManager::install("phyloseq")
if (!requireNamespace("ggplot2", quietly = TRUE)) BiocManager::install("ggplot2")
if (!requireNamespace("rgr", quietly = TRUE)) install.packages("rgr")
library("rgr")
library("phyloseq")
library("ggpubr")
library("ROCR")
library("philr")
library("ggplot2")
library("randomForest")
library("ape")
library("compositions")
print("finished loading libraries")

##-Establish directory layout---------------------------------------##
home_dir <- file.path('~','git','balance_tree_exploration')
project <- "Vanderbilt"
#home_dir <- file.path('cloud','project')
output_dir <- file.path(home_dir, project, 'output')
setwd(file.path(home_dir))

##-Functions--------------------------------------------------------##
source(file.path(home_dir, "lib", "statistical_functions.R"))
source(file.path(home_dir, "lib", "table_manipulations.R"))

##-Set up constants-------------------------------------------------##
rf_cols <- 3:7
num_cycles <- 20
if(num_cycles < 3) stop("num_cycles should be 3 or more")

print("making empty lists for trees and data.frames")
phylo_objects <- list()
phylo_objects_names <- c()
dataframes <- list()
dataframes_names <- c()

##-Import tables and data preprocessing-----------------------------##
asv_table <- data.frame(readRDS(file.path(output_dir, "r_objects", "ForwardReads_DADA2.rds")))

print("Cleaning Ref tree otu with philr tutorial normalization")
ref_ps <- readRDS(file.path(output_dir, "r_objects", "ref_tree_phyloseq_obj.rds"))
phylo_objects <- append(phylo_objects, ref_ps)
phylo_objects_names <- append(phylo_objects_names, "Silva_ref_orig")

clean_otu <- data.frame(ref_ps@otu_table@.Data)
clean_otu <- philr_tutorial_normalization(clean_otu)
print(paste("nrow orginal ref:", nrow(ref_ps@otu_table), "nrow clean ref: ", nrow(clean_otu)))

phy_tree(ref_ps) <- ape::makeNodeLabel(phy_tree(ref_ps), method="number", prefix='n')
ref_ps_clean <- phyloseq::phyloseq( otu_table(clean_otu, taxa_are_rows = F), 
                                    phy_tree(ref_ps@phy_tree),
                                    tax_table(ref_ps@tax_table), 
                                    sample_data(ref_ps@sam_data))
phylo_objects <- append(phylo_objects, ref_ps_clean)
phylo_objects_names <- append(phylo_objects_names, "Silva_ref_cln")

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
phylo_objects <- append(phylo_objects, cln_denovo_tree_ps)
phylo_objects_names <- append(phylo_objects_names, "cln_upgma")
phylo_objects <- append(phylo_objects, denovo_tree_ps)
phylo_objects_names <- append(phylo_objects_names, "orig_upgma")
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
phylo_objects <- append(phylo_objects, ref_rand_list)
phylo_objects <- append(phylo_objects, cln_upgma_rand_list)
phylo_objects <- append(phylo_objects, orig_upgma_rand_list)
phylo_objects_names <- append(phylo_objects_names, 
                              c(names(ref_rand_list), 
                              names(cln_upgma_rand_list), 
                              names(orig_upgma_rand_list)))

print("creating lognorm, ALR and CLR")
if (dir.exists(file.path(output_dir,"r_objects", "lognorm_asv.rds"))) {
  ln_asv_tab <- readRDS(file.path(output_dir,"r_objects", "lognorm_asv.rds"))
}else{
  ln_asv_tab <- lognorm(asv_table)
  saveRDS(ln_asv_tab, file = file.path(output_dir,"r_objects", "lognorm_asv.rds"))
}
my_zeros <- apply(asv_table, 2, function(x) {
  return(sum(x == 0))
})
alr_col <- which(my_zeros == min(my_zeros))[1]
print("creating ALR")
if (file.exists(file.path(output_dir,"r_objects", "alr_asv.rds"))) {
  my_alr <- readRDS(file.path(output_dir,"r_objects", "alr_asv.rds"))
}else{
  my_alr <- as.data.frame(rgr::alr(as.matrix(asv_table + 1), j = as.numeric(alr_col)))
  saveRDS(my_alr, file = file.path(output_dir,"r_objects", "alr_asv.rds"))
}
print("creating CLR")
if (dir.exists(file.path(output_dir,"r_objects", "clr_asv.rds"))) {
  my_clr <- readRDS(file.path(output_dir,"r_objects", "clr_asv.rds"))
}else{
  my_clr <- as.data.frame(rgr::clr(as.matrix(asv_table + 1)))
  saveRDS(my_clr, file = file.path(output_dir,"r_objects", "clr_asv.rds"))
}

print("adding dfs to dataframe list")
dataframes <- append(dataframes, list(asv_table))
dataframes_names <- append(dataframes_names, "raw_data")
dataframes <- append(dataframes, list(ln_asv_tab))
dataframes <- append(dataframes, list(my_alr))
dataframes <- append(dataframes, list(asv_table))
dataframes_names <- append(dataframes_names, "lognorm")
dataframes_names <- append(dataframes_names, "alr")
dataframes_names <- append(dataframes_names, "my_clr")

##-Create tree attribute vectors------------------------------------##
tree_name <- c()
num_nodes <- c()
num_tips <- c()
ave_branch_length <- c()
is_ultrametric <- c()
var_branch_length <- c()
philr_ncol <- c()
# compute.brlen(phy, method = "Grafen", power = 1, ...)
# base.freq()

##-Create non-tree attribute vectors--------------------------------##
table_name <- c()
ncol_tab <- c()
nrow_tab <- c()
has_zeros <- c()
total_cells <- c()

print("updating tree and table descriptives from phylo objects")
for (phy in 1:length(phylo_objects)){
  print(phy)
  #tree
  my_phy <- phylo_objects[[phy]]
  my_tree <- my_phy@phy_tree
  tree_name <- c(tree_name, phylo_objects_names[phy])
  ave_branch_length <- c(ave_branch_length, mean(my_tree$edge.length))
  var_branch_length <- c(var_branch_length, var(my_tree$edge.length))
  num_nodes <- c(num_nodes, my_tree$Nnode)
  num_tips <- c(num_tips, length(my_tree$tip.label))
  is_ultrametric <- c(is_ultrametric, ape::is.ultrametric(my_tree))
  #table
  my_table <- as.data.frame(my_phy@otu_table)
  if ( phy != 4 && phy < 25 ) {
    if (any(my_table == 0)){
      new_table <- my_table + 1
      my_philr <- philr::philr(df = new_table, tree = my_tree)
      philr_ncol <- c(philr_ncol, ncol(my_philr))
    }else{
      my_philr <- philr::philr(df = my_table, tree = my_tree)
      philr_ncol <- c(philr_ncol, ncol(my_philr))
    }
  }else{
    philr_ncol <- c(philr_ncol, NA)
  }
  table_name <- c(table_name, phylo_objects_names[phy])
  ncol_tab <- c(ncol_tab, ncol(my_table))
  nrow_tab <- c(nrow_tab, nrow(my_table))
  has_zeros <- c(has_zeros, any(my_table == 0))
  total_cells <- c(total_cells, ncol(my_table) * nrow(my_table))
}

print("updating table descriptives from datasets")
for (dt in 1:length(dataframes)){
  my_table <- dataframes[[dt]]
  table_name <- c(table_name, dataframes_names[dt])
  ncol_tab <- c(ncol_tab, ncol(my_table))
  nrow_tab <- c(nrow_tab, nrow(my_table))
  has_zeros <- c(has_zeros, any(my_table == 0))
  total_cells <- c(total_cells, ncol(my_table) * nrow(my_table))
}

print("forming dataframes from tree data")
my_trees <- data.frame(tree_name,
                       num_nodes,
                       num_tips,
                       # ave_branch_length,
                       philr_ncol,
                       is_ultrametric,
                       philr_ncol)

my_tables <- data.frame(table_name, ncol_tab, nrow_tab, 
                        has_zeros, total_cells)
write.table(my_trees,
            sep = ",",
            row.names = FALSE,
            file = file.path(output_dir, "tables", paste0(project, "tree_descriptves", ".csv")))
            

write.table(my_tables,
            sep = ",",
            row.names = FALSE,
            file = file.path(output_dir, "tables", paste0(project, "table_descriptves", ".csv")))



