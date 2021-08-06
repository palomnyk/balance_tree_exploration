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
asv_table <- data.frame(readRDS(file.path(output_dir, "r_objects", "ForwardReads_DADA2.rds")))

ref_ps <- readRDS(file.path(output_dir, "r_objects", "ref_tree_phyloseq_obj.rds"))
asv_tax <- data.frame(readRDS(file.path(output_dir, "r_objects", "ForwardReads_DADA2_taxonomy.rds")))

metadata <- read.table(file.path(home_dir, project, "patient_metadata.tsv"), 
                       sep="\t", 
                       header=TRUE, 
                       row.names = "Run", 
                       check.names = FALSE,
                       stringsAsFactors=TRUE)

metadata <- metadata[row.names(metadata) %in% row.names(asv_table), ]

##-----------------Create training/testing sets---------------------##
set.seed(36)
train_index <- sample(x = nrow(metadata), size = 0.75*nrow(metadata), replace=FALSE)
test_index <- c(1:nrow(metadata))[!(1:nrow(metadata) %in% train_index)]

my_ds_ilr_type <- c("None", rep("ilr", 6), rep("philR",4))
my_ds_method <- c("None", 
                  "basic", "balanced","optimal","PBhclust","PBmaxvar","PBangprox", #ilr types
                  "uniform","blw","blw.sqrt","mean.descendants") #philr types

pdf(file = file.path(output_dir, "graphics","auc_ilrVphilr.pdf"))
for(mta in 3:ncol(metadata)){
  print(paste("starting metadata col:", mta))
  resp_var_test <- metadata[,mta][test_index]
  resp_var_train <- metadata[,mta][train_index]
  
  #to be filled by output
  metadata_col <- c()
  auc_val <- c()
  ds_name <- c()
  ds_method <- c()
  
  for( ds in 1:length(my_ds_method)){
    print(paste("starting ds:", ds))
    if (ds == 1) my_table <- asv_table 
    else {
      if (ds <= 7){
        my_table <- compositions::ilr(asv_table, method = my_ds_method[ds])
      }else {
        my_table <- philr::philr(ref_ps@otu_table, ref_ps@phy_tree, 
                                 ilr.weights = my_ds_method[ds])
      }#end else2
    }#end else1
    my_table_train <- my_table[train_index,]
    my_table_test <- my_table[test_index,]
    
    names(resp_var_test) <- row.names(my_table_test)
    
    rf <- randomForest(my_table_train, resp_var_train)
    
    pred <- predict(rf, my_table_test)
    
    preds <- ROCR::prediction(as.numeric(pred), as.numeric(resp_var_test))
    #update output
    metadata_col <- append(metadata_col, mta)
    auc_val <- append(auc_val, ROCR::performance(preds, "auc"))
    ds_name <- append(ds_name, my_ds_ilr_type[ds])
    ds_method <- append(ds_method, my_ds_method[ds])
    
  }#for ds
  print("Writing table")
  write.table(data.frame(metadata_col, auc_val, ds_name, ds_method), 
              row.names=FALSE, sep = ",",
              file = file.path(output_dir, "tables", 
                               paste0(project, "_ilrVphilr_auc_", colnames(metadata)[mta])))
  
  plot(auc_val ~ rep(0, length(auc_val)), 
       xlim = c(0, max(auc_val)),
       axis = FALSE, type = "n",
       xlab = "",
       col = ds_name,
       main = paste(project, "ilr vs philR AUC", colnames(metadata)[mta]),
       ylab = "")
  axis(1, at = auc_val, labels = auc_val)
  legend('topright', 
         legend = ds_method,
         col = ds_name)
}#for mta
dev.off()


# mta <- 3
# 
# train_index = sample(x = nrow(metadata), size = 0.75*nrow(metadata), replace=FALSE)
# test_index = c(1:nrow(metadata))[!(1:nrow(metadata) %in% train_index)]
# 
# resp_var_test <- metadata[,mta][test_index]
# resp_var_train <- metadata[,mta][train_index]
# 
# my_table <- as.data.frame(asv_table)
# my_table_train <- my_table[train_index,]
# my_table_test <- my_table[test_index,]
# 
# names(resp_var_test) <- row.names(my_table_test)
# 
# rf <- randomForest(my_table_train, resp_var_train)
# 
# pred <- predict(rf, my_table_test)
# 
# library("ROCR")
# preds <- ROCR::prediction(as.numeric(pred), as.numeric(resp_var_test), label.ordering = NULL)
# 
# auc <- ROCR::performance(preds, "auc")
# 
# roc <- ROCR::performance(preds, measure = "tpr", "fpr")
# 
# plot(roc)
# 
# auc@y.values[[1]]



#TODO function for creating R & S in random RS
random_RS <- function(choice_vect, #this is the vector of options
                      num_bal, #number of balances (will be ncol of resulting ilr transform)
                      max_comp = length(choice_vect) - 1 ){ #max number of components per balance
  indices_used <- c()
  r_indices <- c()
  s_indices <- c()
  output_name <- c()
  out_put <- c()
  
  #scratch:
  num_bal <- 3
  
  R1 <- c(1,2,3,4)
  S1 <- 5
}

# asv_table1 <- asv_table[1:50, 1:50]
# asv_table2 <- asv_table[5:10, 5:10]
# test1 <- compositions::ilrBase(asv_table1)
# test2 <- compositions::ilrBase(asv_table2)
# 
# test3 <- compositions::ilrBase(asv_table1, method = "balanced")
# test4 <- compositions::ilrBase(asv_table2, method = "balanced")
# 
# test <- compositions::ilr(asv_table1, method = "basic")
# test <- compositions::ilr(asv_table1, method = "balanced")
# test <- compositions::ilr(asv_table1, method = "optimal")
# test <- compositions::ilr(asv_table1, method = "PBhclust")
# test <- compositions::ilr(asv_table1, method = "PBmaxvar")
# test <- compositions::ilr(asv_table1, method = "PBangprox")



