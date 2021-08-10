# Author: Aaron Yerke (aaronyerke@gmail.com)
# Script for making distribution of AUC of ILR balances to compare to philr

rm(list = ls()) #clear workspace

##-functions--------------------------------------------------------##

##-Load Depencencies------------------------------------------------##
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("ROCR", quietly = TRUE)) BiocManager::install("ROCR")
library("compositions")
library("philr")
library("ggplot2")
library("randomForest")
library("ROCR")
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

##-Import tables and data preprocessing-----------------------------##
asv_table <- asv_table <- data.frame(readRDS(file.path(output_dir, "r_objects", "ForwardReads_DADA2.rds")))

ref_ps <- readRDS(file.path(output_dir, "r_objects", "ref_tree_phyloseq_obj.rds"))
clean_otu <- data.frame(ref_ps@otu_table@.Data)
clean_otu <- philr_tutorial_normalization(clean_otu)
phy_tree(ref_ps) <- ape::makeNodeLabel(phy_tree(ref_ps), method="number", prefix='n')
ref_ps_clean <- phyloseq::phyloseq( otu_table(clean_otu, taxa_are_rows = F), 
                                    phy_tree(ref_ps@phy_tree),
                                    tax_table(ref_ps@tax_table), 
                                    sample_data(ref_ps@sam_data))

metadata <- read.table(file.path(home_dir, project, "patient_metadata.tsv"), 
                       sep="\t", 
                       header=TRUE, 
                       row.names = "Run", 
                       check.names = FALSE,
                       stringsAsFactors=TRUE)
metadata$type <- droplevels(metadata$type)
metadata$type <- factor(metadata$type)

metadata <- metadata[row.names(metadata) %in% row.names(clean_otu), ]
##-----------------Create training/testing sets---------------------##
set.seed(36)
train_index <- sample(x = nrow(metadata), size = 0.75*nrow(metadata), replace=FALSE)
test_index <- c(1:nrow(metadata))[!(1:nrow(metadata) %in% train_index)]

#to be filled by output from generate random tree data
all_rand_auc <-c()
rand_metadata_col <- c()

# generate random tree data
for (treee in 1:100){
  #make random tree
  rand_tree <- rtree(n = length(ref_ps@phy_tree$tip.label), tip.label = ref_ps@phy_tree$tip.label)
  #make random tree ps object with clean otu data
  rand_ps <- phyloseq::phyloseq( otu_table(clean_otu, taxa_are_rows = F), 
                                 phy_tree(rand_tree),
                                 tax_table(ref_ps@tax_table), 
                                 sample_data(ref_ps@sam_data))
  
  my_table <- philr::philr(rand_ps@otu_table, rand_ps@phy_tree)
  
  my_table_train <- my_table[train_index,]
  my_table_test <- my_table[test_index,]
  
  for(mta in 3:6){
    print(paste("starting metadata col:", mta))
    resp_var_test <- metadata[,mta][test_index]
    resp_var_train <- metadata[,mta][train_index]
    
    #rf requires rownames on resp var
    names(resp_var_test) <- row.names(my_table_test)
    
    print(paste("attempting to run tree", treee, "through meta", mta))
    rf <- randomForest(my_table_train, resp_var_train)
    
    pred <- predict(rf, my_table_test)
    preds <- ROCR::prediction(as.numeric(pred), as.numeric(resp_var_test))
    auc <- ROCR::performance(preds, "auc")@y.values[[1]]
    print(paste("auc: ", auc))
    #update output
    rand_metadata_col <- append(rand_metadata_col, mta)
    all_rand_auc <- append(all_rand_auc, auc)
  }
}

# calculate philr auc 
philr_taxa_weights <- c("uniform","gm.counts","anorm","anorm.x.gm.counts","enorm","enorm.x.gm.counts")
philr_ilr_weights <- c("uniform","blw","blw.sqrt","mean.descendants")

taxa_weight <- c()
ilr_weight <- c()
all_ref_auc <-c()
ref_metadata_col <- c()

# generate random tree data
for (ilr_w in 1:length(philr_ilr_weights)){
  for (tax_w in 1:length(philr_taxa_weights)){
    my_table <- philr::philr(ref_ps_clean@otu_table, ref_ps_clean@phy_tree, 
                             part.weights = philr_taxa_weights[tax_w],
                             ilr.weights = philr_ilr_weights[ilr_w])
    
    my_table_train <- my_table[train_index,]
    my_table_test <- my_table[test_index,]
    
    for(mta in 3:6){
      print(paste("starting metadata col:", mta))
      resp_var_test <- metadata[,mta][test_index]
      resp_var_train <- metadata[,mta][train_index]
      
      #rf requires rownames on resp var
      names(resp_var_test) <- row.names(my_table_test)
      
      rf <- randomForest(my_table_train, resp_var_train)
      
      pred <- predict(rf, my_table_test)
      
      preds <- ROCR::prediction(as.numeric(pred), as.numeric(resp_var_test))
      auc <- ROCR::performance(preds, "auc")@y.values[[1]]
      print(paste("auc: ", auc))
      #update output
      ref_metadata_col <- append(ref_metadata_col, mta)
      all_ref_auc <- append(all_ref_auc, auc)
      taxa_weight <- c(taxa_weight, philr_taxa_weights[tax_w])
      ilr_weight <- c(ilr_weight, philr_ilr_weights[ilr_w])
    }
  }
}


pdf(file = file.path(output_dir, "graphics","auc_rand_vs_ref_philr.pdf"))
hist(all_rand_auc, breaks = 150, xlab = "AUC", 
     main = "Histogram of all Random Tree philR")

hist(all_ref_auc, breaks = 150, xlab = "AUC", 
     main = "Histogram of all Refernce Tree philR")

#https://www.dataanalytics.org.uk/plot-two-overlapping-histograms-on-one-chart-in-r/
ref_hist <- hist(all_ref_auc, breaks = 300, plot = FALSE)

rand_hist <- hist(all_rand_auc, breaks = 300, plot = FALSE)

plot(ref_hist, col = "blue", 
     main = "Random Tree philR (red) vs Reference Tree philR (blue)",
     xlab = "AUC")
plot(rand_hist, col = rgb(255,192,203, max = 255, alpha = 80, names = "lt.pink"), 
     add = TRUE)

dev.off()
