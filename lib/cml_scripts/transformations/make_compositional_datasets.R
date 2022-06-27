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
# hashseq <- data.frame(data.table::fread(file = file.path(output_dir,"hashseq", "hashseq.csv"),
#                                         header="auto", data.table=FALSE), row.names = 1)

print("creating DADA2 lognorm, ALR and CLR")
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
print("creating DADA2 ALR")
if (file.exists(file.path(output_dir,"r_objects", "alr_asv.rds"))) {
  DADA2_alr <- readRDS(file.path(output_dir,"r_objects", "alr_asv.rds"))
}else{
  DADA2_alr <- as.data.frame(rgr::alr(as.matrix(asv_table + 1), j = as.numeric(alr_col)))
  saveRDS(DADA2_alr, file = file.path(output_dir,"r_objects", "alr_asv.rds"))
  write.csv(DADA2_alr, file = file.path(output_dir,"tables", "alr_asv.csv"))
}
print("creating DADA2 CLR")
if (dir.exists(file.path(output_dir,"r_objects", "clr_asv.rds"))) {
  DADA2_clr <- readRDS(file.path(output_dir,"r_objects", "clr_asv.rds"))
}else{
  DADA2_clr <- as.data.frame(rgr::clr(as.matrix(asv_table + 1)))
  saveRDS(DADA2_clr, file = file.path(output_dir,"r_objects", "clr_asv.rds"))
  write.csv(DADA2_clr, file = file.path(output_dir,"tables", "clr_asv.csv"))
}

# print("Creating hashseq lognorm, ALR and CLR.")
# if (dir.exists(file.path(output_dir,"r_objects", "lognorm_hashseq.rds"))) {
#   ln_hs_tab <- readRDS(file.path(output_dir,"r_objects", "lognorm_HashSeq.rds"))
# }else{
# 	ln_hs_tab <- lognorm(hashseq)
#   saveRDS(ln_hs_tab, file = file.path(output_dir,"r_objects", "lognorm_hashseq.rds"))
#   write.csv(ln_hs_tab, file = file.path(output_dir,"tables", "lognorm_hashseq.csv"))
# }
# print("Making HashSeq clr.")
# if (dir.exists(file.path(output_dir,"r_objects", "r_objects", "clr_hashseq.rds"))) {
#   HashSeq_clr <- readRDS(file.path(output_dir,"r_objects", "clr_hashseq.rds"))
# }else{
#   HashSeq_clr <- as.data.frame(rgr::clr(as.matrix(hashseq + 1)))
#   saveRDS(ln_hs_tab, file = file.path(output_dir,"r_objects", "clr_hashseq.rds"))
#   write.csv(ln_hs_tab, file = file.path(output_dir,"tables", "clr_hashseq.csv"))
# }
# print("Making HashSeq alr.")
# my_zeros <- apply(hashseq, 2, function(x) {
#   return(sum(x == 0))
# })
# alr_col <- which(my_zeros == min(my_zeros))[1]
# if (dir.exists(file.path(output_dir,"r_objects", "r_objects", "alr_hashseq.rds"))) {
#   HashSeq_alr <- readRDS(file.path(output_dir,"r_objects", "alr_hashseq.rds"))
# }else{
#   HashSeq_alr <- as.data.frame(rgr::alr(as.matrix(hashseq + 1), j = as.numeric(alr_col)))
#   saveRDS(HashSeq_alr, file = file.path(output_dir,"r_objects", "alr_hashseq.rds"))
#   write.csv(HashSeq_alr, file = file.path(output_dir,"tables", "alr_hashseq.csv"))
# }
print("Making Silva control")
ref_ps <- readRDS(file.path(output_dir, "r_objects", "ref_tree_phyloseq_obj.rds"))
write.csv(ref_ps@otu_table, file.path(output_dir,"tables", "Silva_ref_counts.csv"))

print("Reached end of script.")
