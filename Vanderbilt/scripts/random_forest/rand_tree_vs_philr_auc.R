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
                                 philr_taxa_weights = philr_taxa_weights,
                                 just_otu = FALSE){
  all_auc <- c()
  metadata_col <- c()
  taxa_weight <- c()
  ilr_weight <- c()
  for (ilr_w in 1:length(philr_ilr_weights)){
    for (tax_w in 1:length(philr_taxa_weights)){
      if (just_otu == TRUE){
        my_table <- ps_obj
      }else{
        my_table <- philr::philr(ps_obj@otu_table, ps_obj@phy_tree, 
                                 part.weights = philr_taxa_weights[tax_w],
                                 ilr.weights = philr_ilr_weights[ilr_w])
      }
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
      if (just_otu == TRUE) break
    }#taxa
    if (just_otu == TRUE) break
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


##-Set up constants-------------------------------------------------##
rf_cols <- 3:6
num_cycles <- 100
if(num_cycles < 3) stop("num_cycles should be 3 or more")

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

#for making different philr weights
philr_taxa_weights <- c("uniform","gm.counts","anorm","anorm.x.gm.counts","enorm","enorm.x.gm.counts")
philr_ilr_weights <- c("uniform","blw","blw.sqrt","mean.descendants")

##-Random num seed--------------------------------------------------##
set.seed(36)

##-Create plot data-------------------------------------------------##
all_plot_data <- data.frame(all_auc = c(),
                            metadata_col = c(),
                            taxa_weight = c(),
                            ilr_weight = c(),
                            random_batch = c(),
                            tree_group = c())
for (ran_num in 1:num_cycles){
  ##-Create training/testing sets-------------------------------------##
  train_index <- sample(x = nrow(metadata), size = 0.75*nrow(metadata), replace=FALSE)
  test_index <- c(1:nrow(metadata))[!(1:nrow(metadata) %in% train_index)]
  
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
  rand_plot_data$tree_group <- rep("random", nrow(rand_plot_data))
  all_plot_data <- rbind(all_plot_data, rand_plot_data)
  
  # calculate ref philr auc 
  ref_plot_data <- make_ilr_taxa_auc_df(ps_obj = ref_ps_clean,
                                        metadata_cols = rf_cols,
                                        metadata = metadata,
                                        train_index = train_index,
                                        test_index = test_index,
                                        philr_ilr_weights = philr_ilr_weights,  
                                        philr_taxa_weights = philr_taxa_weights)
  ref_plot_data$tree_group <- rep("Silva_ref", nrow(ref_plot_data))
  all_plot_data <- rbind(all_plot_data, ref_plot_data)
  
  # generate denovo tree data
  denovo_plot_data <- make_ilr_taxa_auc_df( ps_obj = denovo_tree_ps,
                                            metadata_cols = rf_cols,
                                            metadata = metadata,
                                            train_index = train_index,
                                            test_index = test_index,
                                            philr_ilr_weights = philr_ilr_weights,  
                                            philr_taxa_weights = philr_taxa_weights)
  denovo_plot_data$tree_group <- rep("UPGMA", nrow(denovo_plot_data))
  all_plot_data <- rbind(all_plot_data, denovo_plot_data)
  
  # generate "raw data" data
  raw_plot_data <- make_ilr_taxa_auc_df(ps_obj = asv_table,
                                        metadata_cols = rf_cols,
                                        metadata = metadata,
                                        train_index = train_index,
                                        test_index = test_index,
                                        philr_ilr_weights = philr_ilr_weights,  
                                        philr_taxa_weights = philr_taxa_weights,
                                        just_otu = TRUE )
  raw_plot_data$tree_group <- rep("raw_data", nrow(raw_plot_data))
  all_plot_data <- rbind(all_plot_data, raw_plot_data)
}

weight_table <- data.frame(tree_type = c(),
                         metadata = c(),
                         taxa_pval = c(),
                         ilr_pval = c())


##-Make all the histograms------------------------------------------##
pdf(file = file.path(output_dir, "graphics", paste0("auc_rand_v_ref_v_upgma_v_raw_", num_cycles, ".pdf")))
for (mta in 1:length(unique(all_plot_data$metadata_col))){
  my_meta <- as.character(unique(all_plot_data$metadata_col)[mta])
  message(my_meta)
  
  plot_data <- all_plot_data[all_plot_data$metadata_col == my_meta,]
  
  my_shap_pval <- c()
  my_shap_tg <- c()
  taxa_w_pval <- c()
  ilr_w_pval <- c()
  for (tg in unique(plot_data$tree_group)){
    auc <- plot_data$all_auc[plot_data$tree_group ==  tg]
    ilr_w <- plot_data$ilr_weight[plot_data$tree_group ==  tg]
    taxa_w <- plot_data$taxa_weight[plot_data$tree_group ==  tg]
    rand_shap <- shapiro.test(auc)
    my_shap_tg <- c(my_shap_tg, tg)
    my_shap_pval <- c(my_shap_pval, rand_shap$p.value)
    if (tg != "raw_data"){
      t_pval <- anova(lm(auc ~ taxa_w))$"Pr(>F)"[1]
      i_pval <-  anova(lm(auc ~ ilr_w))$"Pr(>F)"[1]
      new_row <- c(tg, my_meta,  t_pval, i_pval)
      weight_table <- rbind(weight_table, new_row)
    }
  }
  # "tree_type" = c(),
  # "metadata" = c(),
  # "taxa_pval" = c(),
  # "ilr_pval" = c())
  
  
  g <- ggplot2::ggplot(plot_data, aes(all_auc, tree_group)) + 
    ggplot2::geom_boxplot(outlier.shape = NA) +
    ggplot2::geom_jitter(aes(color = as.factor(ilr_weight)),width = 0.1, height = 0.1) +
    ggplot2::ggtitle(paste(project, my_meta, "colored by ilr weight")) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::theme(axis.text.x = element_text(angle = 45)) +
    ggplot2::xlab("AUC") +
    ggplot2::ylab("Tree type") +
    ggplot2::labs(color = "ilr weight")  
  print(g)
  
  g <- ggplot2::ggplot(plot_data, aes(all_auc, tree_group)) + 
    ggplot2::geom_boxplot(outlier.shape = NA) +
    ggplot2::geom_jitter(aes(color = as.factor(taxa_weight)),width = 0.1, height = 0.1) +
    ggplot2::ggtitle(paste(project, my_meta, "colored by taxa weight")) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::theme(axis.text.x = element_text(angle = 45)) +
    ggplot2::xlab("AUC") +
    ggplot2::ylab("Tree type") +
    ggplot2::labs(color = "ilr weight")  
  print(g)
}

dev.off()

weight_table$taxa_adj <- p.adjust(weight_table$taxa_pval, method = "BH")
weight_table$ilr_adj <- p.adjust(weight_table$ilr_pval, method = "BH")

write.table(weight_table,
            sep = ",",
            row.names = FALSE,
            file = file.path(output_dir, "tables", paste0("auc_rand_v_ref_v_upgma_v_raw_", num_cycles, ".csv")))

# #normality
# qqnorm(plot_data$all_auc, main = paste0("Random Tree $", my_meta," shap: ", round(rand_shap$p.value, 6)))
# qqline(plot_data$all_auc)
# g <- ggplot2::ggplot(plot_data, aes(x=all_auc)) + 
#   ggplot2::geom_histogram(bins = 150, colour="black", size = 0.1) +
#   # ggplot2::labs(fill = "Metadata") +
#   ggplot2::ggtitle(paste0(project, " only Rand Forst AUC ", my_meta, " shap: ", round(rand_shap$p.value, 4))) +
#   ggplot2::xlab("AUC") +
#   ggplot2::ylab("Samples per bin")
# print(g)












# Notes from meeting with AF on 11 Aug 2021
# Fit random trees to gausian
# Rug
# Finish on imigrant gut


