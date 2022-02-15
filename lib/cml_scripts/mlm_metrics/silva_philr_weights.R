# Author: Aaron Yerke (aaronyerke@gmail.com)
# Script for making all weightschemes in philr for Silva Ref Tree

rm(list = ls()) #clear workspace

##-functions--------------------------------------------------------##
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
if (!requireNamespace("phyloseq", quietly = TRUE)) BiocManager::install("phyloseq")
library("phyloseq")
if (!requireNamespace("optparse", quietly = TRUE)){
  install.packages("optparse")
}
library("optparse")

print("Finished loading libraries")

option_list <- list(
  make_option(c("-d", "--homedir"), type="character", 
              default=file.path('~','git','balance_tree_exploration'), 
              help="dataset dir path", metavar="home dir"),
  make_option(c("-p", "--project"), type="character", default=NULL, 
              help="project folder", metavar="project")
); 

opt_parser <- optparse::OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);

print(opt)

print("Establishing directory constants.")
home_dir <- opt$homedir
project <- opt$project
output_dir <- file.path(home_dir, project, 'output')
setwd(file.path(home_dir))

source(file.path(home_dir, "lib", "statistical_functions.R"))
source(file.path(home_dir, "lib", "table_manipulations.R"))

print("Setting up other constants.")
main_output_label <- paste0("philr_weights_shuffled")
philr_taxa_weights <- c("uniform","gm.counts","anorm","anorm.x.gm.counts","enorm","enorm.x.gm.counts")
philr_ilr_weights <- c("uniform","blw","blw.sqrt","mean.descendants")
ref_tree_weights_data_dir_path <- file.path(output_dir, "tables", "silva_philr_weights")
table_name <- "ref_tree_cln"

print("Importing and prepping data.")
ref_ps <- readRDS(file.path(output_dir, "r_objects", "ref_tree_phyloseq_obj.rds"))
phy_tree(ref_ps) <- ape::makeNodeLabel(phy_tree(ref_ps), method="number", prefix='n')
ref_ps_clean <- raw_ps_to_clean_ps(ref_ps)
print("Cleaning UPGMA tree otu with philr tutorial normalization")


if(!file.exists(ref_tree_weights_data_dir_path)){
  print(paste("Creating folder at ", ref_tree_weights_data_dir_path))
  dir.create(ref_tree_weights_data_dir_path)
}

print("Building output files.")
for (ilr_w in 1:length(philr_ilr_weights)){
  iw <- philr_ilr_weights[ilr_w]
  for (tax_w in 1:length(philr_taxa_weights)){
    pw <- philr_taxa_weights[tax_w]
    table_name_full <- paste0(paste(table_name, iw, pw, sep = "_"),".csv")
    my_table <- philr::philr(ref_ps_clean@otu_table, ref_ps_clean@phy_tree,
                             part.weights = philr_taxa_weights[tax_w],
                             ilr.weights = philr_ilr_weights[ilr_w])
    write.csv(my_table, 
              file = file.path(ref_tree_weights_data_dir_path, table_name_full),
              sep = ",", row.names = TRUE)
  }
}

print("script complete")

