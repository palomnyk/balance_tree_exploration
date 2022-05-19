# Author: Aaron Yerke (aaronyerke@gmail.com)
# Script for making tables to run through my 
# Requires UPGMA_Tree, Silva_tree, Iqtree, and metadata

rm(list = ls()) #clear workspace

print("Loading custom functions.")
raw_ps_to_clean_ps <- function(ps) {
  #requires ape, phyloseq and philr_tutorial_normalization 
  clean_otu = data.frame(ps@otu_table@.Data)
  clean_otu = philr_tutorial_normalization(clean_otu)
  ps_clean = phyloseq::phyloseq( otu_table(clean_otu, taxa_are_rows = F), 
                                 phy_tree(ps@phy_tree),
                                 tax_table(ps@tax_table), 
                                 sample_data(ps@sam_data))
  return(ps_clean)
}

##-Load Depencencies------------------------------------------------##
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("ape", quietly = TRUE)) BiocManager::install("ape")
library("ape")
if (!requireNamespace("philr", quietly = TRUE)) BiocManager::install("philr")
library("philr")
if (!requireNamespace("randomForest", quietly = TRUE)) BiocManager::install("randomForest")
library("randomForest")
if (!requireNamespace("ggpubr", quietly = TRUE)) BiocManager::install("ggpubr")
library("ggpubr")
if (!requireNamespace("phyloseq", quietly = TRUE)) BiocManager::install("phyloseq")
library("phyloseq")
if (!requireNamespace("ggplot2", quietly = TRUE)) BiocManager::install("ggplot2")
library("ggplot2")
if (!requireNamespace("rgr", quietly = TRUE)) install.packages("rgr")
library("rgr")
if (!requireNamespace("pROC", quietly = TRUE)) BiocManager::install("pROC")
library("pROC")
if (!requireNamespace("data.table", quietly = TRUE)) BiocManager::install("data.table")
library("data.table")
if (!requireNamespace("optparse", quietly = TRUE)){
  install.packages("optparse")
}
library("optparse")
print("finished loading libraries")

option_list <- list(
  optparse::make_option(c("-d", "--homedir"), type="character", 
                        default=file.path('~','git','balance_tree_exploration'), 
                        help="dataset dir path", metavar="character"),
  optparse::make_option(c("-p", "--project"), type="character", default=NULL, 
                        help="project folder", metavar="character"),
  optparse::make_option(c("-m", "--metadata"), type="character", default=NULL,
                        help="metadata file path with filename", metavar="character"),
  optparse::make_option(c("-l", "--metadata_delim"), type="character", default="\t",
                        help="metadata file deliminator", metavar="character"),
  optparse::make_option(c("-r", "--metadata_rowname"), type="character", default=NULL,
                        help="metadata file row to use for row names", metavar="character"),
  optparse::make_option(c("-n", "--num_cycles"), type="numeric", default=20,
                        help="Number of times to shuffle data and run loop again", 
                        metavar="character")
); 

opt_parser <- optparse::OptionParser(option_list=option_list);

opt <- optparse::parse_args(opt_parser);

print(opt)

print("Establishing directory layout and other constants.")
home_dir <- opt$homedir
project <- opt$project
output_dir <- file.path(home_dir, project, 'output')

##-Functions--------------------------------------------------------##
source(file.path(home_dir, "lib", "statistical_functions.R"))
source(file.path(home_dir, "lib", "table_manipulations.R"))

##-Set up constants-------------------------------------------------##
num_cycles <- opt$num_cycles
if(num_cycles < 3) stop("num_cycles should be 3 or more")
main_output_text <- "auc_rand_v_ref_v_upgma_v_raw_vert_"
main_output_label <- paste0("auc_rand_v_ref_v_upgma_v_raw_vert_", num_cycles)
philr_taxa_weights <- c("uniform","gm.counts","anorm","anorm.x.gm.counts","enorm","enorm.x.gm.counts")
philr_ilr_weights <- c("uniform","blw","blw.sqrt","mean.descendants")
random_seed <- 36

print("Importing and preprocessing tables.")
asv_table <- data.frame(readRDS(file.path(output_dir, "r_objects", "ForwardReads_DADA2.rds")))

# print("Cleaning Ref tree otu with philr tutorial normalization")
# ref_ps <- readRDS(file.path(output_dir, "r_objects", "ref_tree_phyloseq_obj.rds"))
# phy_tree(ref_ps) <- ape::makeNodeLabel(phy_tree(ref_ps), method="number", prefix='n')
# phyloseq::plot_tree(ref_ps, method = "treeonly", nodelabf=nodeplotblank, title = paste0("orig_ref"))

# ref_ps_clean <- raw_ps_to_clean_ps(ref_ps)
# phyloseq::plot_tree(ref_ps_clean, method = "treeonly", nodelabf=nodeplotblank, title = paste0("cln_ref"))
# print("Cleaning UPGMA tree otu with philr tutorial normalization")

# denovo_tree_ps <- readRDS(file.path(output_dir, "r_objects", "denovo_tree_UPGMA_phyloseq_obj.rds"))
# phy_tree(denovo_tree_ps) <- ape::makeNodeLabel(phy_tree(denovo_tree_ps), method="number", prefix='n')
# phyloseq::plot_tree(denovo_tree_ps, method = "treeonly", nodelabf=nodeplotblank, title = paste0("orig_upgma"))

# cln_denovo_tree_ps <- raw_ps_to_clean_ps(denovo_tree_ps)
# denovo_tree_ps <- phyloseq::transform_sample_counts(denovo_tree_ps, function(x) x + 1 )
# phyloseq::plot_tree(cln_denovo_tree_ps, method = "treeonly", nodelabf=nodeplotblank, title = paste0("cln_upgma"))

# iqtree_orig_ps <- readRDS(file.path(output_dir, "r_objects", "denovo_tree_iqtree_phyloseq_obj.rds"))
# phy_tree(iqtree_orig_ps) <- ape::makeNodeLabel(phy_tree(iqtree_orig_ps), method="number", prefix='n')
# phyloseq::plot_tree(iqtree_orig_ps, method = "treeonly", nodelabf=nodeplotblank, title = paste0("orig_iqtree"))

# cln_iqtree_ps <- raw_ps_to_clean_ps(iqtree_orig_ps)
# phyloseq::plot_tree(cln_iqtree_ps, method = "treeonly", nodelabf=nodeplotblank, title = paste0("cln_iqtree"))
# dev.off()

# print("loading and munging metadata")
# metadata <- read.table(opt$metadata,
#                        sep=opt$metadata_delim,
#                        header=TRUE,
#                        row.names = opt$metadata_rowname,
#                        check.names = FALSE,
#                        stringsAsFactors=TRUE)
# print("Attempting to equalize metadata rows to seq data rows")
# # needed_rows <- row.names(data.frame(ref_ps@otu_table@.Data))
# # metadata <- data.frame(metadata[my_rows, ])
# # metadata <- metadata[ order(row.names(metadata)),]#order the rows in alphanumeric order by rowname
# # metadata <- metadata[match(needed_rows, row.names(needed_rows)),]
# rf_cols <- 1:ncol(metadata)#hack so I don't have to fix this in the function

# print("Attempting to read HashSeq count table")
# hashseq <- data.frame(data.table::fread(file = file.path(output_dir,"hashseq", "hashseq.csv"),
#                                         header=TRUE, data.table=FALSE), row.names = 1)
# # hashseq <- hashseq[, colSums(hashseq != 0) > 0.01*nrow(hashseq)]#remove columns that don't have at least 10%
# # print(paste("HashSeq has", ncol(hashseq), "columns after column reduction."))

# # needed_rows <- row.names(data.frame(ref_ps@otu_table@.Data))
# metadata <- data.frame(metadata[row.names(hashseq), ])
# metadata <- metadata[ order(row.names(metadata)),]#order the rows in alphanumeric order by rowname


# if (base::identical(row.names(hashseq), row.names(metadata))){
#   print("Rownames of hashseq and metadata are the same.")
# }else{
#   print("Problem with hashseq rownames - quitting.")
#   quit_due_row_names()
# }

# ##-Random num seed--------------------------------------------------##
# print(paste("Setting random seed to:", random_seed))
# set.seed(random_seed)
# print("making random trees")
# orig_ref_rand_list <- list()
# for (rand in 1:10){
#   rand_tree <- ape::rtree(n = length(ref_ps@phy_tree$tip.label), tip.label = ref_ps@phy_tree$tip.label)
#   #put int in philr
#   rand_tree_ps <- phyloseq::phyloseq( otu_table(ref_ps_clean@otu_table, taxa_are_rows = F),
#                                       phy_tree(rand_tree),
#                                       tax_table(ref_ps@tax_table),
#                                       sample_data(ref_ps@sam_data))
#   phy_tree(rand_tree_ps) <- ape::makeNodeLabel(phy_tree(rand_tree_ps), method="number", prefix='n')
#   phyloseq::plot_tree(rand_tree_ps,  method = "treeonly", nodelabf=nodeplotblank, title = paste0("orig_ref_rand_", rand))
#   orig_ref_rand_list[[rand]] <- rand_tree_ps
# }

# print("make random trees for cln upgma taxa")
# cln_upgma_rand_list <- list()
# for (rand in 1:10){
#   rand_tree <- ape::rtree(n = length(cln_denovo_tree_ps@phy_tree$tip.label), tip.label = cln_denovo_tree_ps@phy_tree$tip.label)
#   #put int in philr
#   rand_tree_ps <- phyloseq::phyloseq( otu_table(cln_denovo_tree_ps@otu_table, taxa_are_rows = F),
#                                       phy_tree(rand_tree),
#                                       tax_table(cln_denovo_tree_ps@tax_table),
#                                       sample_data(cln_denovo_tree_ps@sam_data))
#   phy_tree(rand_tree_ps) <- ape::makeNodeLabel(phy_tree(rand_tree_ps), method="number", prefix='n')
#   phyloseq::plot_tree(rand_tree_ps, title = paste0("cln_upgma_rand_", rand))
#   cln_upgma_rand_list[[rand]] <- rand_tree_ps
# }

# print("make random trees for clean ref taxa")
# cln_ref_rand_list <- list()
# for (rand in 1:10){
#   rand_tree <- ape::rtree(n = length(ref_ps_clean@phy_tree$tip.label), tip.label = ref_ps_clean@phy_tree$tip.label)
#   #put int in philr
#   rand_tree_ps <- phyloseq::phyloseq(otu_table(ref_ps_clean@otu_table, taxa_are_rows = F),
#                                      phy_tree(rand_tree),
#                                      tax_table(ref_ps_clean@tax_table),
#                                      sample_data(ref_ps_clean@sam_data))
#   phy_tree(rand_tree_ps) <- ape::makeNodeLabel(phy_tree(rand_tree_ps), method="number", prefix='n')
#   phyloseq::plot_tree(rand_tree_ps, method = "treeonly", nodelabf=nodeplotblank, title = paste0("cln_ref_rand_", rand))
#   cln_ref_rand_list[[rand]] <- rand_tree_ps
# }
# print("make random trees for iqtree clean")
# iqtree_clean_rand_list <- list()
# for (rand in 1:10){
#   rand_tree <- ape::rtree(n = length(cln_iqtree_ps@phy_tree$tip.label), tip.label = cln_iqtree_ps@phy_tree$tip.label)
#   #put int in philr
#   rand_tree_ps <- phyloseq::phyloseq(otu_table(cln_iqtree_ps@otu_table, taxa_are_rows = F),
#                                      phy_tree(rand_tree),
#                                      tax_table(cln_iqtree_ps@tax_table),
#                                      sample_data(cln_iqtree_ps@sam_data))
#   phy_tree(rand_tree_ps) <- ape::makeNodeLabel(phy_tree(rand_tree_ps), method="number", prefix='n')
#   phyloseq::plot_tree(rand_tree_ps, method = "treeonly", nodelabf=nodeplotblank, title = paste0("iqtree_rand_", rand))
#   iqtree_clean_rand_list[[rand]] <- rand_tree_ps
# }

print("creating lognorm, ALR and CLR")
if (dir.exists(file.path(output_dir,"r_objects", "lognorm_asv.rds"))) {
  ln_asv_tab <- readRDS(file.path(output_dir,"r_objects", "lognorm_asv.rds"))
}else{
  ln_asv_tab <- lognorm(asv_table)
  saveRDS(ln_asv_tab, file = file.path(output_dir,"r_objects", "lognorm_asv.rds"))
  write.csv(ln_asv_tab, file = file.path(output_dir,"tables", "lognorm_dada2.csv"))
}

my_zeros <- apply(asv_table, 2, function(x) {
  return(sum(x == 0))
})
alr_col <- which(my_zeros == min(my_zeros))[1]
# alr_col_num <- grep(alr_col, colnames(asv_table))
print("creating ALR")
if (file.exists(file.path(output_dir,"r_objects", "alr_asv.rds"))) {
  DADA2_alr <- readRDS(file.path(output_dir,"r_objects", "alr_asv.rds"))
}else{
  DADA2_alr <- as.data.frame(rgr::alr(as.matrix(asv_table + 1), j = as.numeric(alr_col)))
  saveRDS(DADA2_alr, file = file.path(output_dir,"r_objects", "alr_asv.rds"))
  write.csv(DADA2_alr, file = file.path(output_dir,"tables", "alr_asv.csv"))
}
print("creating CLR")
if (dir.exists(file.path(output_dir,"r_objects", "clr_asv.rds"))) {
  DADA2_clr <- readRDS(file.path(output_dir,"r_objects", "clr_asv.rds"))
}else{
  DADA2_clr <- as.data.frame(rgr::clr(as.matrix(asv_table + 1)))
  saveRDS(DADA2_clr, file = file.path(output_dir,"r_objects", "clr_asv.rds"))
  write.csv(DADA2_clr, file = file.path(output_dir,"tables", "clr_asv.csv"))
}

hashseq <- data.frame(data.table::fread(file = file.path(output_dir,"hashseq", "hashseq.csv"),
                                        header=TRUE, data.table=FALSE), row.names = 1)

print("creating hashseq lognorm, ALR and CLR")
if (dir.exists(file.path(output_dir,"r_objects", "lognorm_asv.rds"))) {
  ln_hs_tab <- readRDS(file.path(output_dir,"r_objects", "lognorm_HashSeq.rds"))
}else{
	ln_hs_tab <- lognorm(hashseq)
  saveRDS(ln_hs_tab, file = file.path(output_dir,"r_objects", "lognorm_hashseq.rds"))
  write.csv(ln_hs_tab, file = file.path(output_dir,"tables", "lognorm_hashseq.csv"))
}
HashSeq_clr <- as.data.frame(rgr::clr(as.matrix(hashseq + 1)))
saveRDS(HashSeq_clr, file = file.path(output_dir,"r_objects", "clr_hashseq.rds"))
write.csv(HashSeq_clr, file = file.path(output_dir,"r_objects", "clr_hashseq.csv"))
my_zeros <- apply(asv_table, 2, function(x) {
  return(sum(x == 0))
})
alr_col <- which(my_zeros == min(my_zeros))[1]
HashSeq_alr <- as.data.frame(rgr::alr(as.matrix(hashseq + 1), j = as.numeric(alr_col)))
saveRDS(HashSeq_alr, file = file.path(output_dir,"r_objects", "alr_hashseq.rds"))
write.csv(HashSeq_alr, file = file.path(output_dir,"tables", "alr_hashseq.csv"))
