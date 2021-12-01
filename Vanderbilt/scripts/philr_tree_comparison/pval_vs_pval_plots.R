# Author: Aaron Yerke (aaronyerke@gmail.com)
# Script for doing pvalue comparison of philr made from different trees
rm(list = ls()) #clear workspace

##-functions--------------------------------------------------------##
tree_philr_node_pvals <- function(philr_df, phylo_tree) {
  #function for walking down nodes of a tree and comparing their child nodes
  #don't need this function anymore
  #internal functions:
  select.tip.or.node <- function(element, tree) {
    # https://stackoverflow.com/questions/51696837/r-phylo-object-how-to-connect-node-label-and-node-number
    ifelse(element < Ntip(tree)+1, tree$tip.label[element], tree$node.label[element-Ntip(tree)])
  }#end select.tip.or.node
  my_tree <- ref_ps_clean@phy_tree
  select.tip.or.node <- function(element, tree) {
    # https://stackoverflow.com/questions/51696837/r-phylo-object-how-to-connect-node-label-and-node-number
    ifelse(element < Ntip(tree)+1, tree$tip.label[element], tree$node.label[element-Ntip(tree)])
  }
  edge_table <- data.frame(
    "parent" = phylo_tree$edge[,1],
    "par.name" = sapply(phylo_tree$edge[,1], select.tip.or.node, tree = phylo_tree),
    "child" = phylo_tree$edge[,2],
    "chi.name" = sapply(phylo_tree$edge[,2], select.tip.or.node, tree = phylo_tree)
  )# 1 through Ntips is tip and Ntips through Ntips + Nnode is node
  tree_map <- ape::prop.part(ref_ps_clean@phy_tree)
  
  pvals <- c()
  
  start_edge <- ape::Ntip(phylo_tree) + 1
  end_edge <- ape::Nnode(phylo_tree) + ape::Ntip(phylo_tree)
  
  for (node_num in start_edge:end_edge){
    parent_node_indices <- which(edge_table$parent == node_num)
    child_nodes <- as.character( edge_table$chi.name[parent_node_indices] )
    if (length(child_nodes) == 2 & all(startsWith( child_nodes, "n"))){
      philr_1 <- philr_df[,as.character(child_nodes[1])]
      philr_2 <- philr_df[,as.character(child_nodes[2])]
      pval <- t.test(philr_1, philr_2)$p.value
      pvals <- c(pvals, pval)
    }
  }
  return(pvals)
}

log10_pval_pval_plot <- function(pvals_1, label_1, pvals_2, label_2, meta_c, colmn_nam) {
  # function for making pval plots
  max_len <- min(length(pvals_1), length(pvals_2))
  plot_data <- data.frame(pval1 = sort(log10(sample(pvals_1, max_len, replace = F))), 
                          pval2 = sort(log10(sample(pvals_2, max_len, replace = F))))
  my_cor <- cor(plot_data$pval1, plot_data$pval2)
  g <- ggplot2::ggplot(data = plot_data, aes(x = pval1, y = pval2)) +  
    ggplot2::ggtitle(paste("log10(pval) v  log10(pval)", "Rsq:", round(my_cor, 4),"meta:", colmn_nam)) +
    ggplot2::geom_point() +
    ggplot2::geom_abline(intercept = 0, slope = 1) +
    ggplot2::theme_classic() +
    ggplot2::ylab(label_1) +
    ggplot2::xlab(label_2) +
    ggplot2::labs(color = "Part weight")
  return(g)
}

philr_by_meta_pvals <- function(philr_obj, meta_col) {
  library("philr")
  my_philr <- philr::philr(df = philr_obj@otu_table,
                           tree = philr_obj@phy_tree)
  pvals <- apply(my_philr, 2, function(x) {t.test(x ~ meta_col)$p.value})
  return(pvals)
}
##-Load Depencencies------------------------------------------------##
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# if (!requireNamespace("ape", quietly = TRUE)) BiocManager::install("ape")
# if (!requireNamespace("philr", quietly = TRUE)) BiocManager::install("philr")
# if (!requireNamespace("phyloseq", quietly = TRUE)) BiocManager::install("phyloseq")
# if (!requireNamespace("ggplot2", quietly = TRUE)) BiocManager::install("ggplot2")
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
print("Loading raw dada2 count table")
asv_table <- data.frame(readRDS(file.path(output_dir, "r_objects", "ForwardReads_DADA2.rds")))

print("Cleaning Ref tree otu with philr tutorial normalization")
ref_ps <- readRDS(file.path(output_dir, "r_objects", "ref_tree_phyloseq_obj.rds"))

clean_otu <- data.frame(ref_ps@otu_table@.Data)
clean_otu <- philr_tutorial_normalization(clean_otu)
print(paste("nrow orginal ref:", nrow(ref_ps@otu_table), "nrow clean ref: ", nrow(clean_otu)))

print("loading and munging metadata")
metadata <- read.table(file.path(home_dir, project, "patient_metadata.tsv"), 
                       sep="\t", 
                       header=TRUE, 
                       row.names = "Run", 
                       check.names = FALSE,
                       stringsAsFactors=TRUE)
metadata <- metadata[row.names(metadata) %in% row.names(clean_otu), ]
metadata$type <- droplevels(metadata$type)
metadata$type <- factor(metadata$type)

print("'Cleaning' ref_ps")
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
print("making orig ref random trees")
orig_ref_rand_list <- list()
for (rand in 1:10){
  rand_tree <- rtree(n = length(ref_ps@phy_tree$tip.label), tip.label = ref_ps@phy_tree$tip.label)
  #put int in philr
  rand_tree_ps <- phyloseq::phyloseq( otu_table(clean_otu, taxa_are_rows = F),
                                      phy_tree(rand_tree),
                                      tax_table(ref_ps@tax_table),
                                      sample_data(ref_ps@sam_data))
  phy_tree(rand_tree_ps) <- ape::makeNodeLabel(phy_tree(rand_tree_ps), method="number", prefix='n')
  orig_ref_rand_list[[rand]] <- rand_tree_ps
}
print("make random trees for cln upgma taxa")
cln_upgma_rand_list <- list()
for (rand in 1:10){
  rand_tree <- rtree(n = length(cln_denovo_tree_ps@phy_tree$tip.label), tip.label = cln_denovo_tree_ps@phy_tree$tip.label)
  #put int in philr
  rand_tree_ps <- phyloseq::phyloseq( otu_table(cln_denovo_tree_ps, taxa_are_rows = F),
                                      phy_tree(rand_tree),
                                      tax_table(cln_denovo_tree_ps@tax_table),
                                      sample_data(cln_denovo_tree_ps@sam_data))
  phy_tree(rand_tree_ps) <- ape::makeNodeLabel(phy_tree(rand_tree_ps), method="number", prefix='n')
  cln_upgma_rand_list[[rand]] <- rand_tree_ps
}

print("make random trees for clean ref taxa")
cln_ref_rand_list <- list()
for (rand in 1:10){
  rand_tree <- rtree(n = length(ref_ps_clean@phy_tree$tip.label), tip.label = denovo_tree_ps@phy_tree$tip.label)
  #put int in philr
  rand_tree_ps <- phyloseq::phyloseq(otu_table(ref_ps_clean, taxa_are_rows = F),
                                     phy_tree(rand_tree),
                                     tax_table(ref_ps_clean@tax_table),
                                     sample_data(ref_ps_clean@sam_data))
  phy_tree(rand_tree_ps) <- ape::makeNodeLabel(phy_tree(rand_tree_ps), method="number", prefix='n')
  cln_ref_rand_list[[rand]] <- rand_tree_ps
}

# print("make random trees for original upgma taxa")
# orig_upgma_rand_list <- list()
# for (rand in 1:10){
#   rand_tree <- ape::rtree(n = length(denovo_tree_ps@phy_tree$tip.label), tip.label = denovo_tree_ps@phy_tree$tip.label)
#   #put int in philr
#   rand_tree_ps <- phyloseq::phyloseq( otu_table(denovo_tree_ps, taxa_are_rows = F),
#                                       phy_tree(rand_tree),
#                                       tax_table(denovo_tree_ps@tax_table),
#                                       sample_data(denovo_tree_ps@sam_data))
#   phy_tree(rand_tree_ps) <- ape::makeNodeLabel(phy_tree(rand_tree_ps), method="number", prefix='n')
#   orig_upgma_rand_list[[rand]] <- rand_tree_ps
#   names(orig_upgma_rand_list)[rand] <- paste0("orig_upgma_rand", rand)
# }
##-Calculate pvalues------------------------------------------------##
print("making philr tables")

#for making different philr weights
philr_taxa_weights <- c("uniform","gm.counts","anorm","anorm.x.gm.counts","enorm","enorm.x.gm.counts")
philr_ilr_weights <- c("uniform","blw","blw.sqrt","mean.descendants")

metad_num <- 7
metad_col <- metadata[,7]
metad_nam <- colnames(metadata)[metad_num]

# philr_ds <- list(cln_denovo_tree_ps, ref_ps, ref_ps_clean)
# names_philr_ds <- c("cln_upgma", "ref_orig", "ref_clean")


# philr_pvals <- lapply(philr_ds, function(x){return(philr_by_meta_pvals(x,metad_col))})

ref_ps_clean_pvals <- philr_by_meta_pvals(ref_ps_clean, metad_col)

# test_pvals <- apply(ref_ps_clean_philr, 2, function(x) {t.test(x ~ metad_col)$p.value})

cln_denovo_philr <- philr::philr(df = cln_denovo_tree_ps@otu_table,
                                   tree = cln_denovo_tree_ps@phy_tree)
cln_denovo_pvals <- apply(cln_denovo_philr, 2, function(x) {t.test(x ~ metad_col)$p.value})

ref_ps_orig_philr <- philr::philr(df = ref_ps@otu_table + 1,
                            tree = ref_ps@phy_tree)
ref_ps_orig_pvals <- apply(ref_ps_orig_philr, 2, function(x) {t.test(x ~ metad_col)$p.value})

pdf(file = file.path(output_dir, "graphics", "pval_vs_pval_plots.pdf"))
hist(ref_ps_clean_pvals, breaks = 50)
hist(ref_ps_orig_pvals, breaks = 50)
hist(cln_denovo_pvals, breaks = 50)

#verify qqplots by hand (histogram or node by node)
#Are these pvals reasonable?
#check against other metadata columns
#run the other dataset
stats::qqplot(ref_ps_clean_pvals, ref_ps_orig_pvals,
              main = paste("qqplot", "ref_ps_clean_pvals", "ref_ps_orig_philr"))
log10_pval_pval_plot(ref_ps_orig_pvals, "ref_ps_orig_pvals",
                     ref_ps_clean_pvals, "ref_ps_clean_pvals", metad_col, metad_nam)
stats::qqplot(ref_ps_clean_pvals, cln_denovo_pvals,
              main = paste("qqplot", "ref_ps_clean_pvals", "cln_denovo_pvals"))
log10_pval_pval_plot(ref_ps_clean_pvals, "ref_ps_clean_pvals",
                     cln_denovo_pvals, "cln_denovo_pvals", metad_col, metad_nam)
stats::qqplot(ref_ps_orig_pvals, cln_denovo_pvals,
              main = paste("qqplot", "ref_ps_clean_pvals", "cln_denovo_pvals"))
log10_pval_pval_plot(cln_denovo_pvals, "cln_denovo_pvals",
                     ref_ps_clean_pvals, "ref_ps_clean_pvals", metad_col, metad_nam)

print("Generate random tree p values")
ref_cln_rand_pvals <- lapply(cln_ref_rand_list, function(x) {
  philr_by_meta_pvals(x, metad_col)})
cln_upgma_rand_pvals <- lapply(cln_upgma_rand_list, function(x) {
  philr_by_meta_pvals(x, metad_col)})
ref_orig_rand_pvals <- lapply(orig_ref_rand_list, function(x) {
  philr_by_meta_pvals(x, metad_col)})

print("Generate random tree p values")
cln_upgma_vs_ref_cln_rand_cor <- lapply(ref_cln_rand_pvals, function(x){
  len <- min(length(x), length(cln_denovo_pvals))
  cor(sort(sample(cln_denovo_pvals, len)), sort(sample(x, len)))
})

cln_upgma_vs_cln_upgma_rand_cor <- lapply(cln_upgma_rand_pvals, function(x){
  len <- min(length(x), length(cln_denovo_pvals))
  cor(sort(sample(cln_denovo_pvals, len)), sort(sample(x, len)))
})

cln_upgma_vs_orig_ref_rand_cor <- lapply(ref_orig_rand_pvals, function(x){
  len <- min(length(x), length(cln_denovo_pvals))
  cor(sort(sample(cln_denovo_pvals, len)), sort(sample(x, len)))
})

ref_cln_vs_ref_cln_rand_cor <- lapply(ref_cln_rand_pvals, function(x){
  len <- min(length(x), length(ref_ps_clean_pvals))
  cor(sort(sample(ref_ps_clean_pvals, len)), sort(sample(x, len)))
})
#
ref_cln_vs_cln_upgma_rand_cor <- lapply(cln_upgma_rand_pvals, function(x){
  len <- min(length(x), length(ref_ps_clean_pvals))
  cor(sort(sample(ref_ps_clean_pvals, len)), sort(sample(x, len)))
})

ref_cln_vs_ref_orig_rand_cor <- lapply(ref_orig_rand_pvals, function(x){
  len <- min(length(x), length(ref_ps_clean_pvals))
  cor(sort(sample(ref_ps_clean_pvals, len)), sort(sample(x, len)))
})

ref_orig_vs_ref_cln_rand_cor <- lapply(ref_cln_rand_pvals, function(x){
  len <- min(length(x), length(ref_ps_clean_pvals))
  cor(sort(sample(ref_ps_orig_pvals, len)), sort(sample(x, len)))
})

ref_orig_vs_cln_upgma_rand_cor <- lapply(cln_upgma_rand_pvals, function(x){
  len <- min(length(x), length(ref_ps_clean_pvals))
  cor(sort(sample(ref_ps_orig_pvals, len)), sort(sample(x, len)))
})

ref_orig_vs_ref_orig_rand_cor <- lapply(ref_orig_rand_pvals, function(x){
  len <- min(length(x), length(ref_ps_clean_pvals))
  cor(sort(sample(ref_ps_orig_pvals, len)), sort(sample(x, len)))
})

print("Munging output for boxpots")
bp_labels <- c(rep("cln_upgma_vs_ref_cln_rand_cor",length(cln_upgma_vs_ref_cln_rand_cor)),
               rep("cln_upgma_vs_cln_upgma_rand_cor",length(cln_upgma_vs_cln_upgma_rand_cor)),
               rep("cln_upgma_vs_orig_ref_rand_cor",length(cln_upgma_vs_orig_ref_rand_cor)),
               rep("ref_cln_vs_ref_cln_rand_cor",length(ref_cln_vs_ref_cln_rand_cor)),
               rep("ref_cln_vs_cln_upgma_rand_cor",length(ref_cln_vs_cln_upgma_rand_cor)),
               rep("ref_cln_vs_ref_orig_rand_cor",length(ref_cln_vs_ref_orig_rand_cor)),
               rep("ref_orig_vs_ref_cln_rand_cor",length(ref_orig_vs_ref_cln_rand_cor)),
               rep("ref_orig_vs_cln_upgma_rand_cor",length(ref_orig_vs_cln_upgma_rand_cor)),
               rep("ref_orig_vs_ref_orig_rand_cor",length(ref_orig_vs_ref_orig_rand_cor)))

cor_vals <- c(unlist(cln_upgma_vs_ref_cln_rand_cor),
              unlist(cln_upgma_vs_cln_upgma_rand_cor),
              unlist(cln_upgma_vs_orig_ref_rand_cor),
              unlist(ref_cln_vs_ref_cln_rand_cor),
              unlist(ref_cln_vs_cln_upgma_rand_cor),
              unlist(ref_cln_vs_ref_orig_rand_cor),
              unlist(ref_orig_vs_ref_cln_rand_cor),
              unlist(ref_orig_vs_cln_upgma_rand_cor),
              unlist(ref_orig_vs_ref_orig_rand_cor))

plot_data <- data.frame(bp_labels, cor_vals)

g <- ggplot2::ggplot(data = plot_data, aes(y = bp_labels, x = cor_vals)) +  
  ggplot2::geom_boxplot() +
  ggplot2::ggtitle(paste(project, "Random trees vs non-random-trees, metad:", metad_name)) +
  # ggplot2::geom_hline(yintercept = 0) +
  ggplot2::theme(axis.text.x = element_text(angle = 45),
                 axis.text.y = element_text(angle = 45)) +
  ggplot2::theme_classic() +
  # ggplot2::scale_y_discrete(labels = seq(0, 1, by = 0.2)) +
  # ggplot2::ylab("AUC") +
  # ggplot2::xlab("Tree type") +
  ggplot2::labs(color = "Part weight")
print(g)

dev.off()
