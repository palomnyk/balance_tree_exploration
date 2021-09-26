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
  #Function for making random forest AUC values
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  if (!requireNamespace("ROCR", quietly = TRUE)) BiocManager::install("ROCR")
  library("ROCR")
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
      my_table_train <- data.frame(my_table[train_index,])
      my_table_test <- data.frame(my_table[test_index,])
      for(mta in metadata_cols){
        tryCatch(
          { 
            print(paste("starting metadata col:", mta, colnames(metadata)[mta]))
            if (any(is.na(my_table))) break
            resp_var_test <- metadata[,mta][test_index]
            resp_var_train <- metadata[,mta][train_index]
            #rf requires rownames on resp var
            names(resp_var_test) <- row.names(my_table_test)
            rf <- randomForest::randomForest(my_table_train, resp_var_train)
            pred <- predict(rf, my_table_test)
            # print(paste("pred:", pred))
            # print(paste("num factors", nlevels(resp_var_test)))
            roc_data <- data.frame(pred = pred, resp_var_test = resp_var_test)
            if (nlevels(resp_var_test) > 2){
              print("multilevels")
            }else{
              preds <- ROCR::prediction(as.numeric(pred), as.numeric(resp_var_test))
              auc <- ROCR::performance(preds, "auc")@y.values[[1]]
            }
            print("ROC made")
            # auc <- pROC::auc(my_roc)
            # print(paste("auc: ")
            #update all output
            metadata_col <- append(metadata_col, colnames(metadata)[mta])
            all_auc <- append(all_auc, auc)
            taxa_weight <- c(taxa_weight, philr_taxa_weights[tax_w])
            ilr_weight <- c(ilr_weight, philr_ilr_weights[ilr_w])
          },
          error=function(cond) {
            print('Opps, an error is thrown')
            message(cond)
            print(message(cond))
          },
          warning=function(cond) {
            print('Opps, a warning is thrown')
            message(cond)
            print(message(cond))
            }
        )
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
if (!requireNamespace("ape", quietly = TRUE)) BiocManager::install("ape")
if (!requireNamespace("philr", quietly = TRUE)) BiocManager::install("philr")
if (!requireNamespace("randomForest", quietly = TRUE)) BiocManager::install("randomForest")
if (!requireNamespace("ROCR", quietly = TRUE)) BiocManager::install("ROCR")
if (!requireNamespace("ggpubr", quietly = TRUE)) BiocManager::install("ggpubr")
if (!requireNamespace("phyloseq", quietly = TRUE)) BiocManager::install("phyloseq")
if (!requireNamespace("ggplot2", quietly = TRUE)) BiocManager::install("ggplot2")
library("phyloseq")
library("ggpubr")
library("ROCR")
library("philr")
library("ggplot2")
library("randomForest")
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
total_seqs <- rowSums(asv_table)
total_seqs <- data.frame("total_seqs"=total_seqs, "duplicate" = total_seqs,
                         row.names = row.names(asv_table))

#clean up otu tables
ref_ps <- readRDS(file.path(output_dir, "r_objects", "ref_tree_phyloseq_obj.rds"))
clean_otu <- data.frame(ref_ps@otu_table@.Data)
clean_otu <- philr_tutorial_normalization(clean_otu)
print(paste("nrow orginal ref:", nrow(ref_ps@otu_table), "nrow clean ref: ", nrow(clean_otu)))

phy_tree(ref_ps) <- ape::makeNodeLabel(phy_tree(ref_ps), method="number", prefix='n')
ref_ps_clean <- phyloseq::phyloseq( otu_table(clean_otu, taxa_are_rows = F), 
                                    phy_tree(ref_ps@phy_tree),
                                    tax_table(ref_ps@tax_table), 
                                    sample_data(ref_ps@sam_data))

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

metadata <- read.table(file.path(home_dir, project, "patient_metadata.tsv"), 
                       sep="\t", 
                       header=TRUE, 
                       row.names = "Run", 
                       check.names = FALSE,
                       stringsAsFactors=TRUE)
metadata <- metadata[row.names(metadata) %in% row.names(clean_otu), ]
metadata$type <- droplevels(metadata$type)
metadata$type <- factor(metadata$type)

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
skips <- 0
counter <- 0
while (counter < num_cycles & skips < 5){
  ##-Create training/testing sets-------------------------------------##
  train_index <- sample(x = nrow(metadata), size = 0.75*nrow(metadata), replace=FALSE)
  test_index <- c(1:nrow(metadata))[!(1:nrow(metadata) %in% train_index)]
  should_break <- FALSE
  for(mta in rf_cols){
    if( length(unique(metadata[,mta][test_index])) != nlevels(metadata[,mta][test_index]) |
        length(unique(metadata[,mta][train_index])) != nlevels(metadata[,mta][train_index])){
      print("levels not equal")
      should_break <- TRUE
      skips = skips + 1
    }
  }
  if (should_break == FALSE){
    #make random tree
    rand_tree <- rtree(n = length(ref_ps@phy_tree$tip.label), tip.label = ref_ps@phy_tree$tip.label)
    #put int in philr
    rand_tree_ps <- phyloseq::phyloseq( otu_table(clean_otu, taxa_are_rows = F),
                                        phy_tree(rand_tree),
                                        tax_table(ref_ps@tax_table),
                                        sample_data(ref_ps@sam_data))
    phy_tree(rand_tree_ps) <- makeNodeLabel(rand_tree_ps, method="number", prefix='n')
    write("making random AUC")
    rand_plot_data <- make_ilr_taxa_auc_df(ps_obj = rand_tree_ps,
                                           metadata_cols = rf_cols,
                                           metadata = metadata,
                                           train_index = train_index,
                                           test_index = test_index,
                                           philr_ilr_weights = philr_ilr_weights,
                                           philr_taxa_weights = philr_taxa_weights)
    rand_plot_data$tree_group <- rep("random", nrow(rand_plot_data))
    all_plot_data <- rbind(all_plot_data, rand_plot_data)

    write("making ref AUC")
    ref_plot_data <- make_ilr_taxa_auc_df(ps_obj = ref_ps_clean,
                                          metadata_cols = rf_cols,
                                          metadata = metadata,
                                          train_index = train_index,
                                          test_index = test_index,
                                          philr_ilr_weights = philr_ilr_weights,
                                          philr_taxa_weights = philr_taxa_weights)
    ref_plot_data$tree_group <- rep("Silva_ref", nrow(ref_plot_data))
    all_plot_data <- rbind(all_plot_data, ref_plot_data)

    # write("making UPGMA AUC")
    # denovo_plot_data <- make_ilr_taxa_auc_df( ps_obj = denovo_tree_ps,
    #                                           metadata_cols = rf_cols,
    #                                           metadata = metadata,
    #                                           train_index = train_index,
    #                                           test_index = test_index,
    #                                           philr_ilr_weights = philr_ilr_weights,
    #                                           philr_taxa_weights = philr_taxa_weights)
    # denovo_plot_data$tree_group <- rep("UPGMA", nrow(denovo_plot_data))
    # all_plot_data <- rbind(all_plot_data, denovo_plot_data)

    write("making cleaned UPGMA AUC")
    denovo_plot_data <- make_ilr_taxa_auc_df( ps_obj = cln_denovo_tree_ps,
                                              metadata_cols = rf_cols,
                                              metadata = metadata,
                                              train_index = train_index,
                                              test_index = test_index,
                                              philr_ilr_weights = philr_ilr_weights,
                                              philr_taxa_weights = philr_taxa_weights)
    denovo_plot_data$tree_group <- rep("clean_UPGMA", nrow(denovo_plot_data))
    all_plot_data <- rbind(all_plot_data, denovo_plot_data)

    write('generate "raw data" data')
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
    
    write('generate "read depth" data')
    raw_plot_data <- make_ilr_taxa_auc_df(ps_obj = data.frame(total_seqs),
                                          metadata_cols = rf_cols,
                                          metadata = metadata,
                                          train_index = train_index,
                                          test_index = test_index,
                                          philr_ilr_weights = philr_ilr_weights,  
                                          philr_taxa_weights = philr_taxa_weights,
                                          just_otu = TRUE )
    raw_plot_data$tree_group <- rep("read_depth", nrow(raw_plot_data))
    all_plot_data <- rbind(all_plot_data, raw_plot_data)
    
    write('generate lognorm data')
    raw_plot_data <- make_ilr_taxa_auc_df(ps_obj = lognorm(asv_table),
                                          metadata_cols = rf_cols,
                                          metadata = metadata,
                                          train_index = train_index,
                                          test_index = test_index,
                                          philr_ilr_weights = philr_ilr_weights,
                                          philr_taxa_weights = philr_taxa_weights,
                                          just_otu = TRUE )
    raw_plot_data$tree_group <- rep("lognorm", nrow(raw_plot_data))
    all_plot_data <- rbind(all_plot_data, raw_plot_data)
    
    counter <- counter + 1
    skips <- 0
  }
}

write.table(all_plot_data,
            file = file.path(output_dir, "tables", paste0("auc_rand_v_ref_v_upgma_v_raw_vert_", num_cycles, ".csv")),
            sep = ",",
            row.names = FALSE)

all_plot_data <- read.table(file = file.path(output_dir, "tables", 
                                             paste0("auc_rand_v_ref_v_upgma_v_raw_vert_", 
                                                    num_cycles, ".csv")),
            sep = ",", header = TRUE)

weight_table <- data.frame(tree_type = c(F),
                           metadata = c(F),
                           taxa_pval = c(F),
                           ilr_pval = c(F))
weight_counter <- 1
##-Make all the boxplots--------------------------------------------##
pdf(file = file.path(output_dir, "graphics", paste0("auc_rand_v_ref_v_upgma_v_raw_vert_", num_cycles, ".pdf")))
g <- ggplot2::ggplot(all_plot_data, aes(y = all_auc, x= tree_group)) + 
  ggplot2::geom_boxplot() +
  ggplot2::geom_jitter(aes(color = as.factor(ilr_weight)),width = 0.2, height = 0.001) +
  ggplot2::ggtitle(paste(project, "all metadata")) +
  # ggplot2::geom_hline(yintercept = 0) +
  ggplot2::theme(axis.text.x = element_text(angle = 45),
                 axis.text = element_text(size = 20)) +
  ggplot2::theme_classic() +
  # ggplot2::scale_y_discrete(labels = seq(0, 1, by = 0.2)) +
  ggplot2::ylab("AUC") +
  ggplot2::xlab("Tree type") +
  ggplot2::labs(color = "ilr weight")
print(g)

g <- ggplot2::ggplot(all_plot_data, aes(y = all_auc, x= tree_group)) + 
  ggplot2::geom_boxplot() +
  ggplot2::geom_jitter(aes(color = as.factor(taxa_weight)),width = 0.2, height = 0.001) +
  ggplot2::ggtitle(paste(project, "all metadata")) +
  # ggplot2::geom_hline(yintercept = 0) +
  ggplot2::theme(axis.text.x = element_text(angle = 45)) +
  ggplot2::theme_classic() +
  # ggplot2::scale_y_discrete(labels = seq(0, 1, by = 0.2)) +
  ggplot2::ylab("AUC") +
  ggplot2::xlab("Tree type") +
  ggplot2::labs(color = "Part weight")
print(g)

for (mta in 1:length(unique(all_plot_data$metadata_col))){
  my_meta <- as.character(unique(all_plot_data$metadata_col)[mta])
  message(my_meta)
  
  plot_data <- all_plot_data[all_plot_data$metadata_col == my_meta,]
  
  my_shap_pval <- c()
  my_shap_tg <- c()
  taxa_w_pval <- c()
  ilr_w_pval <- c()
  # for (tg in unique(plot_data$tree_group)){
  #   auc <- plot_data$all_auc[plot_data$tree_group ==  tg]
  #   ilr_w <- plot_data$ilr_weight[plot_data$tree_group ==  tg]
  #   taxa_w <- plot_data$taxa_weight[plot_data$tree_group ==  tg]
  #   rand_shap <- shapiro.test(auc)
  #   my_shap_tg <- c(my_shap_tg, tg)
  #   my_shap_pval <- c(my_shap_pval, rand_shap$p.value)
  #   if (tg != "raw_data"){
  #     t_pval <- anova(lm(auc ~ taxa_w))$"Pr(>F)"[1]
  #     i_pval <-  anova(lm(auc ~ ilr_w))$"Pr(>F)"[1]
  #     new_row <- c(tg, my_meta,  t_pval, i_pval)
  #     names(new_row) <- c("tree_type", "metadata", "taxa_pval", "ilr_pval")
  #     weight_table[weight_counter,] <- new_row
  #     weight_counter <- weight_counter + 1
  #   }
  # }
  
  betwn_bar_anova <- anova(lm(data = plot_data,all_auc ~ tree_group))
  my_comparisons <- list( c("raw_data", "Silva_ref"), c( "UPGMA", "Silva_ref"), c("Silva_ref", "random"), c("raw_data", "random") )
  
  g <- ggplot2::ggplot(plot_data, aes(y = all_auc, x= tree_group)) + 
    ggplot2::geom_boxplot() +
    ggplot2::geom_jitter(aes(color = as.factor(ilr_weight)),width = 0.2, height = 0.001) +
    ggplot2::ggtitle(paste(project, my_meta, "col by ilr", "anova:", round(betwn_bar_anova$`Pr(>F)`, 5))) +
    # ggplot2::geom_hline(yintercept = 0) +
    ggplot2::theme(axis.text.x = element_text(angle = 45)) +
    ggplot2::theme_classic() +
    # ggplot2::scale_y_discrete(labels = seq(0, 1, by = 0.2)) +
    ggplot2::ylab("AUC") +
    ggplot2::xlab("Tree type") +
    ggplot2::labs(color = "ilr weight") +
    ggpubr::stat_compare_means(comparisons = my_comparisons)
  print(g)
  
  g <- ggplot2::ggplot(plot_data, aes(y = all_auc, x= tree_group)) + 
    ggplot2::geom_boxplot() +
    ggplot2::geom_jitter(aes(color = as.factor(taxa_weight)),width = 0.2, height = 0.001) +
    ggplot2::ggtitle(paste(project, my_meta, "col by part", "anova:", round(betwn_bar_anova$`Pr(>F)`, 5))) +    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::theme(axis.text.x = element_text(angle = 45),
                   axis.text.y = element_text(angle = 45)) +
    ggplot2::theme_classic() +
    # ggplot2::scale_y_discrete(labels = seq(0, 1, by = 0.2)) +
    ggplot2::ylab("AUC") +
    ggplot2::xlab("Tree type") +
    ggplot2::labs(color = "Part weight") +
    ggpubr::stat_compare_means(comparisons = my_comparisons)
  print(g)
}

dev.off()

# weight_table$taxa_adj <- p.adjust(weight_table$taxa_pval, method = "BH")
# weight_table$ilr_adj <- p.adjust(weight_table$ilr_pval, method = "BH")
# 
# write.table(weight_table,
#             sep = ",",
#             row.names = FALSE,
#             file = file.path(output_dir, "tables", paste0("auc_rand_v_ref_v_upgma_v_raw_", num_cycles, ".csv")))

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


