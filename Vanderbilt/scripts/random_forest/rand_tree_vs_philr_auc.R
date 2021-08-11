# Author: Aaron Yerke (aaronyerke@gmail.com)
# Script for making distribution of AUC of ILR balances to compare to philr

rm(list = ls()) #clear workspace

##-functions--------------------------------------------------------##
make_ilr_taxa_auc_df <- function(ps_obj,
                                 metadata = metadata,
                                 metadata_cols,
                                 train_index = train_index,
                                 test_index = test_index,
                                 philr_ilr_weights = philr_ilr_weights,  
                                 philr_taxa_weights = philr_taxa_weights){
  all_auc <- c()
  metadata_col <- c()
  taxa_weight <- c()
  ilr_weight <- c()
  for (ilr_w in 1:length(philr_ilr_weights)){
    for (tax_w in 1:length(philr_taxa_weights)){
      my_table <- philr::philr(ps_obj@otu_table, ps_obj@phy_tree, 
                               part.weights = philr_taxa_weights[tax_w],
                               ilr.weights = philr_ilr_weights[ilr_w])
      
      my_table_train <- my_table[train_index,]
      my_table_test <- my_table[test_index,]
      
      for(mta in metadata_cols){
        print(paste("starting metadata col:", mta))
        resp_var_test <- metadata[,mta][test_index]
        resp_var_train <- metadata[,mta][train_index]
        
        #rf requires rownames on resp var
        names(resp_var_test) <- row.names(my_table_test)
        
        rf <- randomForest::randomForest(my_table_train, resp_var_train)
        
        pred <- predict(rf, my_table_test)
        
        preds <- ROCR::prediction(as.numeric(pred), as.numeric(resp_var_test))
        auc <- ROCR::performance(preds, "auc")@y.values[[1]]
        print(paste("auc: ", auc))
        #update output
        metadata_col <- append(metadata_col, colnames(metadata)[mta])
        all_auc <- append(all_auc, auc)
        taxa_weight <- c(taxa_weight, philr_taxa_weights[tax_w])
        ilr_weight <- c(ilr_weight, philr_ilr_weights[ilr_w])
      }#for mta
    }#taxa
  }#ilr
  return(data.frame(all_auc,
                    metadata_col,
                    taxa_weight,
                    ilr_weight))
}#end function

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

#clean up otu tables
ref_ps <- readRDS(file.path(output_dir, "r_objects", "ref_tree_phyloseq_obj.rds"))
clean_otu <- data.frame(ref_ps@otu_table@.Data)
clean_otu <- philr_tutorial_normalization(clean_otu)
phy_tree(ref_ps) <- ape::makeNodeLabel(phy_tree(ref_ps), method="number", prefix='n')
ref_ps_clean <- phyloseq::phyloseq( otu_table(clean_otu, taxa_are_rows = F), 
                                    phy_tree(ref_ps@phy_tree),
                                    tax_table(ref_ps@tax_table), 
                                    sample_data(ref_ps@sam_data))

denovo_tree_ps <- readRDS(file.path(output_dir, "r_objects", "denovo_tree_UPGMA_phyloseq_obj.rds"))
clean_den_otu <- philr_tutorial_normalization(data.frame(denovo_tree_ps@otu_table@.Data))
denovo_tree_ps <- phyloseq::phyloseq( otu_table(clean_den_otu, taxa_are_rows = F),
                                      phy_tree(ape::makeNodeLabel(phy_tree(denovo_tree_ps@phy_tree))),
                                      tax_table(denovo_tree_ps@tax_table), 
                                      sample_data(denovo_tree_ps@sam_data))

metadata <- read.table(file.path(home_dir, project, "patient_metadata.tsv"), 
                       sep="\t", 
                       header=TRUE, 
                       row.names = "Run", 
                       check.names = FALSE,
                       stringsAsFactors=TRUE)
metadata$type <- droplevels(metadata$type)
metadata$type <- factor(metadata$type)

metadata <- metadata[row.names(metadata) %in% row.names(clean_otu), ]

rf_cols <- 3:6

#for making different philr weights
philr_taxa_weights <- c("uniform","gm.counts","anorm","anorm.x.gm.counts","enorm","enorm.x.gm.counts")
philr_ilr_weights <- c("uniform","blw","blw.sqrt","mean.descendants")

##-----------------Create training/testing sets---------------------##
set.seed(36)
train_index <- sample(x = nrow(metadata), size = 0.75*nrow(metadata), replace=FALSE)
test_index <- c(1:nrow(metadata))[!(1:nrow(metadata) %in% train_index)]

##-Create plot data-------------------------------------------------##
# generate random tree data
rand_all_plot_data <- data.frame(all_auc = c(),
                             metadata_col = c(),
                             taxa_weight = c(),
                             ilr_weight = c())
for (treee in 1:10){
  #make random tree
  rand_tree <- rtree(n = length(ref_ps@phy_tree$tip.label), tip.label = ref_ps@phy_tree$tip.label)
  #put int in philr
  rand_tree_ps <- phyloseq::phyloseq( otu_table(clean_otu, taxa_are_rows = F), 
                                      phy_tree(rand_tree),
                                      tax_table(ref_ps@tax_table), 
                                      sample_data(ref_ps@sam_data))
  rand_plot_data <- make_ilr_taxa_auc_df(ps_obj = rand_tree_ps,
                                        metadata_cols = rf_cols,
                                        metadata = metadata,
                                        train_index = train_index,
                                        test_index = test_index,
                                        philr_ilr_weights = philr_ilr_weights,  
                                        philr_taxa_weights = philr_taxa_weights)
  rand_all_plot_data <- rbind(rand_all_plot_data, rand_plot_data)
}

# calculate ref philr auc 
ref_plot_data <- make_ilr_taxa_auc_df(ps_obj = ref_ps_clean,
                                      metadata_cols = rf_cols,
                                      metadata = metadata,
                                      train_index = train_index,
                                      test_index = test_index,
                                      philr_ilr_weights = philr_ilr_weights,  
                                      philr_taxa_weights = philr_taxa_weights)

# generate denovo tree data
denovo_plot_data <- make_ilr_taxa_auc_df( ps_obj = denovo_tree_ps,
                                          metadata_cols = rf_cols,
                                          metadata = metadata,
                                          train_index = train_index,
                                          test_index = test_index,
                                          philr_ilr_weights = philr_ilr_weights,  
                                          philr_taxa_weights = philr_taxa_weights)

##-Make all the histograms------------------------------------------##
pdf(file = file.path(output_dir, "graphics","auc_rand_vs_ref_philr_all.pdf"))
plot_data <- data.frame(all_auc = c(rand_all_plot_data$all_auc,
                                    ref_plot_data$all_auc,
                                    denovo_plot_data$all_auc),
                        auc_lab = c(rep("Random Tree", nrow(rand_all_plot_data)),
                                    rep("Ref Tree", nrow(ref_plot_data)),
                                    rep("Denovo Tree", nrow(denovo_plot_data))))
data_set <- "All Metadata AUC"
g <- ggplot2::ggplot(plot_data, aes(x=all_auc, fill=auc_lab)) + 
  ggplot2::geom_histogram(bins = 150, colour="black", size = 0.1) +
  ggplot2::labs(fill = "Tree type") +
  ggplot2::ggtitle(paste0(project, ", ", data_set," Rand Forst AUC")) +
  ggplot2::xlab("AUC") +
  ggplot2::ylab("Samples per bin")
g
#looking at individual metadata
for (mta in unique(colnames(metadata)[rf_cols])){
  print(mta)
  plot_data <- data.frame(all_auc = c(rand_all_plot_data$all_auc[rand_all_plot_data$metadata_col == mta],
                                      ref_plot_data$all_auc[ref_plot_data$metadata_col == mta],
                                      denovo_plot_data$all_auc[denovo_plot_data$metadata_col == mta]),
                          auc_lab = c(rep("Random Tree", sum(rand_all_plot_data$metadata_col == mta)),
                                      rep("Ref Tree", sum(ref_plot_data$metadata_col == mta)),
                                      rep("Denovo Tree", sum(denovo_plot_data$metadata_col == mta))))
  data_set <- mta
  g <- ggplot2::ggplot(plot_data, aes(x=all_auc, fill=auc_lab)) + 
    ggplot2::geom_histogram(bins = 150, colour="black", size = 0.1) +
    ggplot2::labs(fill = "Tree type") +
    ggplot2::ggtitle(paste0(project, ", ", data_set, " Rand Forst AUC")) +
    ggplot2::xlab("AUC") +
    ggplot2::ylab("Samples per bin")
  print(g)
}#end for

my_dfs <- list(rand_all_plot_data, ref_plot_data, denovo_plot_data)
my_df_names <-c("Random Tree", "Ref Tree", "Denovo")

for (ds in 1:length(my_dfs)){
  plot_data <- my_dfs[ds][[1]]
  data_set <- my_df_names[ds]
  g <- ggplot2::ggplot(plot_data, aes(x=all_auc, fill=metadata_col)) + 
    ggplot2::geom_histogram(bins = 150, colour="black", size = 0.1) +
    ggplot2::labs(fill = "Metadata") +
    ggplot2::ggtitle(paste0(project, ", ", data_set, " only Rand Forst AUC ")) +
    ggplot2::xlab("AUC") +
    ggplot2::ylab("Samples per bin")
  print(g)
  
  g <- ggplot2::ggplot(plot_data, aes(x=all_auc, 
                                      fill = ilr_weight)) + 
    ggplot2::geom_histogram(bins = 150, colour="black", size = 0.1) +
    ggplot2::ggtitle(paste0(project, ", ", data_set, " only Rand Forst AUC ILR.weight")) +
    ggplot2::xlab("AUC") +
    ggplot2::ylab("Samples per bin")
  g
  
  g <- ggplot2::ggplot(plot_data, aes(x=all_auc, 
                                      fill = taxa_weight)) + 
    ggplot2::geom_histogram(bins = 150, colour="black", size = 0.1) +
    ggplot2::ggtitle(paste0(project, ", ", data_set, " only Rand Forst AUC Taxa.weight")) +
    ggplot2::xlab("AUC") +
    ggplot2::ylab("Samples per bin")
  g
}

dev.off()


