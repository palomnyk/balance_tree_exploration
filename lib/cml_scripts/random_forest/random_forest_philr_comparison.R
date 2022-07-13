# Author: Aaron Yerke (aaronyerke@gmail.com)
# Script for making distribution of AUC of ILR balances to compare to philr
# Requires UPGMA_Tree, Silva_tree, Iqtree, and metadata

rm(list = ls()) #clear workspace

##-functions--------------------------------------------------------##
make_ilr_taxa_auc_df <- function(ps_obj,
																	metadata = metadata,
																	metadata_cols,
																	train_index = train_index,
																	test_index = test_index,
																	philr_ilr_weights = philr_ilr_weights,  
																	philr_taxa_weights = philr_taxa_weights,
																	just_otu = FALSE,
																	output_fpath = main_output_fpath,
																	cycle,
																	transf_label,
																	random_label = "TODO"){
	#Function for making random forest AUC values
	if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
	if (!requireNamespace("pROC", quietly = TRUE)) BiocManager::install("pROC")
	library("pROC")
  if (!requireNamespace("MLmetrics", quietly = TRUE)) BiocManager::install("MLmetrics")
  library("MLmetrics")
   
	for (ilr_w in 1:length(philr_ilr_weights)){
		for (tax_w in 1:length(philr_taxa_weights)){
			tryCatch({
				if (just_otu == TRUE){
					print("In just otu")
					my_table <- ps_obj
				}else{
				  if (any(ps_obj@otu_table == 0)){
				    print("adding pseudocount to before PhILR transform")
				    ps_obj <- transform_sample_counts(ps_obj, function(x) x+1)
				  }
					my_table <- philr::philr(ps_obj@otu_table, ps_obj@phy_tree, 
																		part.weights = philr_taxa_weights[tax_w],
																		ilr.weights = philr_ilr_weights[ilr_w])
				}
				my_table_train <- my_table[row.names(my_table) %in% train_index,]
				my_table_test <- my_table[row.names(my_table) %in% test_index,]
				# print(paste("dim my_table_train:", dim(my_table_test)))
				},
				error=function(cond) {
					print(paste("Opps, an error1 is thrown with", transf_label))
					message(paste(transf_label, cond))
				},
				warning=function(cond) {
					print(paste("Opps, a warning1 is thrown with", transf_label))
				  message(paste(transf_label, cond))
				}
			)
			# print("head(my_table_train)")
			# print(head(my_table_train))
			for(mta in metadata_cols){
				tryCatch(
					{ 
						print(paste("starting metadata col:", mta, colnames(metadata)[mta]))
						if (any(is.na(my_table))) {
							print("There are NA's - breaking loop.")
							break
						}
						resp_var_test <- metadata[row.names(metadata) %in% test_index,mta]
						# print(resp_var_test)
						# print(paste("Length of resp_var_test:", length(resp_var_test)))
						resp_var_train <- metadata[row.names(metadata) %in% train_index,mta]
						# print(resp_var_train)
						print("Unique resp var test/ resp var train")
						# print(paste(unique(resp_var_test)))
						# print(paste(unique(resp_var_train)))
						#rf requires rownames on resp var
						names(resp_var_test) <- row.names(my_table_test)
						rf <- randomForest::randomForest(my_table_train, resp_var_train)
						print("made rf")
						pred <- predict(rf, my_table_test)
						# print(paste("pred:", pred))
						# print(paste("num factors", length(unique(resp_var_test))))
						# print(paste("levels resp_var_train", nlevels(as.factor(resp_var_train))))
						roc_data <- data.frame(pred = pred, resp_var_test = resp_var_test)
						
						score <- MLmetrics::Accuracy(pred, resp_var_test)
						
# 						if (length(unique(unlist(resp_var_test))) > 2){
# 							print("multilevels")
# 							mult_auc <- c()
# 							for (fact in 1:length(unique(resp_var_test))){#only need to test resp_test
# 								try({
# 									my_fact <- as.character(levels(resp_var_test)[fact])
# 									# print(paste("my_fact:", my_fact))
# 									dumb_resp_test <- as.factor(replace(as.character(resp_var_test), as.character(resp_var_test) != my_fact, "dumb_var"))
# 									# print("dumb_resp")
# 									# print(paste(dumb_resp_test))
# 									dumb_pred <- as.factor(replace(as.character(pred), as.character(pred) != my_fact, "dumb_var"))
# 									# print("dumb_pred")
# 									# print(paste(dumb_resp_test))
# 									my_roc <- pROC::roc(as.numeric(dumb_pred), as.numeric(dumb_resp_test))
# 									# print("my_roc")
# 									mult_auc <- c(mult_auc, pROC::auc(my_roc))
# 									print(mult_auc)
# 								})
# 							}
# 							if (length(mult_auc) > 0 ) {
# 								print("in 'if (length(mult_auc) > 0 )' statement")
# 								score <- mean(mult_auc)
# 							}else{
# 								break
# 							}
# 						}else{
#               my_roc <- pROC::roc(as.numeric(pred), as.numeric(resp_var_test))
#               print("ROC made")
#               score <- pROC::auc(my_roc)
# 						}
						# score <- pROC::auc(my_roc)
						# print(paste("score: ")
						my_df <- rf$importance
						maxImp <- max(rf$importance)
						maxRow <- which(rf$importance == maxImp)

						#Check its existence
						if (file.exists(output_fpath)) {
						  print(paste0("Writing output to ", output_fpath, " ."))
						  # main_header <- "all_score,	metadata_col, taxa_weight,	ilr_weight,	rf_imp_se, rf_type, rf_ntree, trans_group, random_batch, cycle"
						  cat(paste(paste0("\n", score), colnames(metadata)[mta], philr_taxa_weights[tax_w],#all_score,	metadata_col, taxa_weight
						            philr_ilr_weights[ilr_w], row.names(my_df)[maxRow], rf$type, #ilr_weight,	rf_imp_se, rf_type,
						            rf$ntree, transf_label, random_label, cycle, #rf_ntree, trans_group, random_batch, cycle
						            sep = ","), 
						      file = output_fpath, 
						      append=TRUE)
						}
					},
					error=function(cond) {
					  print(paste("Opps, an error2 is thrown with", transf_label))
					  message(paste(transf_label, cond))
					},
					warning=function(cond) {
					  print(paste("Opps, a warning2 is thrown with", transf_label))
					  message(paste(transf_label, cond))
					}
				)
			}#for mta
			if (just_otu == TRUE) break
		}#taxa
		if (just_otu == TRUE) break
	}#ilr
	return(TRUE)
}#end function
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
if (!requireNamespace("MLmetrics", quietly = TRUE)) BiocManager::install("MLmetrics")
library("MLmetrics")
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
main_output_text <- "random_forest_score_R_"
main_output_label <- paste0(main_output_text, num_cycles)
main_output_fn <- paste0(main_output_label, ".csv")
main_output_fpath <- file.path(output_dir, "tables", main_output_fn)
# philr_taxa_weights <- c("uniform","gm.counts","anorm","anorm.x.gm.counts","enorm","enorm.x.gm.counts")
# philr_ilr_weights <- c("uniform","blw","blw.sqrt","mean.descendants")
philr_taxa_weights <- c("enorm")
philr_ilr_weights <- c("mean.descendants")
random_seed <- 36
main_header <- "all_score, metadata_col, rf_imp_se, rf_type, rf_ntree, trans_group, cycle"

print(paste("Initilizing", main_output_fpath, "."))
print(paste("File header: ", main_header))
cat(main_header, 
    file = main_output_fpath, 
    append=FALSE)

##-Import tables and data preprocessing-----------------------------##
print("loading and munging metadata")
metadata <- read.table(opt$metadata, 
                       sep=opt$metadata_delim, 
                       header=TRUE, 
                       row.names = opt$metadata_rowname, 
                       check.names = FALSE,
                       stringsAsFactors=TRUE)
print("Removing non-factor columns from metadata to decrease run time.")
metadata <- metadata[,sapply(metadata, is.factor)]
metadata <- metadata[,sapply(metadata, function(x){ length(levels(x)) < 10})]#hack to remove date columns if they get counted as factors
# df[,-which(sapply(df, class) == "factor")]
rf_cols <- 1:ncol(metadata)#hack so I don't have to fix this in the function

# --------------------------------------------------------------------------
print("Defining functions")
# --------------------------------------------------------------------------
add_PhILR_dfs_to_table <- function(lst,
                           root_folder,
                           base_fn,
                           # philr_part_weights = c("anorm","enorm"),
                           # philr_ilr_weights = c("blw.sqrt","mean.descendants"),
                           philr_part_weights = c("anorm"),
                           philr_ilr_weights = c("blw.sqrt"),
                           color = "w"){
  if (!dir.exists(root_folder)){
    print(paste0(root_folder,"does not exist. Use PhILR_random_trees_and_counts_tables.R to create it."))
    break
  }
  for (ilr_w in philr_ilr_weights){
    for (tax_w in philr_taxa_weights){
      my_label <- paste(base_fn, ilr_w, tax_w, sep = "_")
      table_fn <- paste0(my_label, ".csv")
      lst[[length(lst) + 1]] <- c(my_label, file.path(root_folder, table_fn), ',', color)
    }
  }
  return(lst)
}

add_random_tree_PhILRs_to_table  <- function(lst, 
                                    root_folder, 
                                    base_fn, 
                                    # philr_part_weights = c("anorm","enorm"),
                                    # philr_ilr_weights = c("blw.sqrt","mean.descendants"),
                                    philr_part_weights = c("anorm"),
                                    philr_ilr_weights = c("blw.sqrt"),
                                    color = "w", 
                                    num_rand_trees = 10){
  print(paste0("Adding random trees from ", base_fn, "."))
  if (!dir.exists(root_folder)){
    print(paste0(root_folder,"does not exist. Use PhILR_random_trees_and_counts_tables.R to create it."))
    break
  }
  for (rand in 1:num_rand_trees){
    for (ilr_w in philr_ilr_weights){
      for (tax_w in philr_taxa_weights){
        my_label <- paste(paste0("Shuffle", rand, "_PhILR"), base_fn, ilr_w, tax_w, sep = "_")
        table_fn <- paste0(my_label, ".csv")
        lst[[length(lst) + 1]] <- c(my_label, file.path(root_folder, table_fn), ',', color)
      }
    }
  }
  return(lst)
}

df_factory <- function(my_path, my_sep){
  # Use for building df from "tables" list in random forest loop
  df <- read.table(my_path, sep=my_sep, 
                  header=0, row.names=0,
                  check.names = FALSE,
                  stringsAsFactors=TRUE)
  return(df)
}

tables <- list()
tables[[length(tables) + 1]] <- c("DaDa2",file.path(output_dir, "tables", "ForwardReads_DADA2.txt"),"\t","r")
# tables[[length(tables) + 1]] <- c("HashSeq", file.path(output_dir,  "hashseq", "hashseq.csv"),",", "r")
tables[[length(tables) + 1]] <- c("lognorm_DADA2", file.path(output_dir, "tables", "lognorm_dada2.csv"), ",", "y")
# tables[[length(tables) + 1]] <- c("lognorm_HashSeq", file.path(output_dir,"tables", "lognorm_hashseq.csv"), ",", "y")
tables[[length(tables) + 1]] <- c("alr_DADA2", file.path(output_dir, "tables", "alr_asv.csv"), ",", "g")
# tables[[length(tables) + 1]] <- c("alr_HashSeq", file.path(output_dir,"tables", "alr_hashseq.csv"), ",", "g")
tables[[length(tables) + 1]] <- c("clr_DADA2", file.path(output_dir, "tables", "clr_asv.csv"), ",", "g")
# tables[[length(tables) + 1]] <- c("clr_HashSeq", file.path(output_dir,"tables", "clr_hashseq.csv"), ",", "m")
tables[[length(tables) + 1]] <- c("Silva_DADA2", file.path(output_dir,"tables", "Silva_DADA2", "Silva_DADA2.csv"), ",", "#64baeb")
tables <- add_PhILR_dfs_to_table(tables, file.path(output_dir, "tables", "Silva_DADA2"), "Silva_DADA2", color = "#c8e5f5")
tables <- add_random_tree_PhILRs_to_table(tables, file.path(output_dir, "tables", "Silva_DADA2"), "Silva_DADA2", color = "#1474aa", num_rand_trees=5)
tables[[length(tables) + 1]] <- c("Filtered_Silva_DADA2", file.path(output_dir,"tables", "Filtered_Silva_DADA2", "Filtered_Silva_DADA2.csv"), ",", "#f78646")
tables <- add_PhILR_dfs_to_table(tables, file.path(output_dir, "tables", "Filtered_Silva_DADA2"), "Filtered_Silva_DADA2", color = "#f5cbb3")
tables <- add_random_tree_PhILRs_to_table(tables, file.path(output_dir, "tables", "Filtered_Silva_DADA2"), "Filtered_Silva_DADA2", color = "#8c390b", num_rand_trees=5)
tables[[length(tables) + 1]] <- c("Filtered_UPGMA_DADA2", file.path(output_dir,"tables", "Filtered_UPGMA_DADA2", "Filtered_UPGMA_DADA2.csv"), ",", "#0000ff")
tables <- add_PhILR_dfs_to_table(tables, file.path(output_dir, "tables", "Filtered_UPGMA_DADA2"), "Filtered_UPGMA_DADA2", color = "#b7b7f3")
tables <- add_random_tree_PhILRs_to_table(tables, file.path(output_dir, "tables", "Filtered_UPGMA_DADA2"), "Filtered_UPGMA_DADA2", color = "#050598", num_rand_trees=5)
tables[[length(tables) + 1]] <- c("Filtered_IQtree", file.path(output_dir,"tables", "Filtered_IQtree", "Filtered_IQtree.csv"), ",", "orange")
tables <- add_PhILR_dfs_to_table(tables, file.path(output_dir, "tables", "Filtered_IQtree"), "Filtered_IQtree", color = "#f7d8a0")
tables <- add_random_tree_PhILRs_to_table(tables, file.path(output_dir, "tables", "Filtered_IQtree"), "Filtered_IQtree", color = "#e29302", num_rand_trees=5)

print(tables)

skips <- 0
counter <- 0

print(paste("Counter:", counter, " entering main loop"))
for (counter in 1:num_cycles) {
  ##-Create training/testing sets-------------------------------------##
  train_index <- row.names(metadata)[sample(x = nrow(metadata), size = 0.75*nrow(metadata), replace=FALSE)]
  test_index <- row.names(metadata)[c(1:nrow(metadata))[!(1:nrow(metadata) %in% train_index)]]
  for(mta in rf_cols){
    print(paste("starting metadata col:", mta, colnames(metadata)[mta]))
    print(paste("Finding rows where", colnames(metadata)[mta], "has NA or Empty values from training and testing selectors."))
    na_or_empty_index <- which(is.na(metadata[,mta]) | metadata[,mta] == "")
    na_or_empty_rows <- row.names(metadata)[na_or_empty_index]
    print(length(train_index))
    train_index <- setdiff(train_index, na_or_empty_rows)
    print(train_index)
    test_index <- setdiff(test_index, na_or_empty_rows)
    print(length(train_index))
    # tryCatch({
      for (tabl in tables){
        print(tabl)
        transf_label <- tabl[1]
        my_table <- data.frame(data.table::fread(file = tabl[2],
                                                header=TRUE, data.table=FALSE),
                               row.names = 1)
        if(setequal(row.names(metadata), row.names(my_table)) == FALSE){
          print(paste0("These samples are in the metdata, but not ", transf_label,": ", setdiff(row.names(metadata), row.names(my_table))))
          print(paste0("These samples are in the ", transf_label, " but not metadata: ", setdiff(row.names(my_table), row.names(metadata))))
          warning(paste("Metadata dataframe and", transf_label, "dataframe must have the same rows (order is not important)."))   
        }
        my_table_train <- my_table[row.names(my_table) %in% train_index,]
        my_table_test <- my_table[row.names(my_table) %in% test_index,]
        resp_var_train <- metadata[row.names(metadata) %in% train_index,mta]
        resp_var_test <- metadata[row.names(metadata) %in% test_index,mta]
        resp_var_train <- droplevels(resp_var_train, exclude= NA, "")
        resp_var_test <- droplevels(resp_var_test, exclude= NA, "")
        print(levels(resp_var_test))
        print(length(resp_var_train))
        names(resp_var_test) <- row.names(my_table_test)
        rf <- randomForest::randomForest(my_table_train, resp_var_train)
        print("made rf")
        pred <- predict(rf, my_table_test)
        roc_data <- data.frame(pred = pred, resp_var_test = resp_var_test)
        score <- MLmetrics::Accuracy(pred, resp_var_test)
        # 						if (length(unique(unlist(resp_var_test))) > 2){
        # 							print("multilevels")
        # 							mult_auc <- c()
        # 							for (fact in 1:length(unique(resp_var_test))){#only need to test resp_test
        # 								try({
        # 									my_fact <- as.character(levels(resp_var_test)[fact])
        # 									# print(paste("my_fact:", my_fact))
        # 									dumb_resp_test <- as.factor(replace(as.character(resp_var_test), as.character(resp_var_test) != my_fact, "dumb_var"))
        # 									# print("dumb_resp")
        # 									# print(paste(dumb_resp_test))
        # 									dumb_pred <- as.factor(replace(as.character(pred), as.character(pred) != my_fact, "dumb_var"))
        # 									# print("dumb_pred")
        # 									# print(paste(dumb_resp_test))
        # 									my_roc <- pROC::roc(as.numeric(dumb_pred), as.numeric(dumb_resp_test))
        # 									# print("my_roc")
        # 									mult_auc <- c(mult_auc, pROC::auc(my_roc))
        # 									print(mult_auc)
        # 								})
        # 							}
        # 							if (length(mult_auc) > 0 ) {
        # 								print("in 'if (length(mult_auc) > 0 )' statement")
        # 								score <- mean(mult_auc)
        # 							}else{
        # 								break
        # 							}
        # 						}else{
        #               my_roc <- pROC::roc(as.numeric(pred), as.numeric(resp_var_test))
        #               print("ROC made")
        #               score <- pROC::auc(my_roc)
        # 						}
        # score <- pROC::auc(my_roc)
        print(paste("score:", score))
        my_df <- rf$importance
        maxImp <- max(rf$importance)
        maxRow <- which(rf$importance == maxImp)
        
        #Check its existence
        if (file.exists(main_output_fpath)) {
          print(paste0("Writing output to ", main_output_fpath, " ."))
          # main_header <- "all_score,	metadata_col,	rf_imp_se, rf_type, rf_ntree, trans_group, random_batch, cycle"
          cat(paste(paste0("\n", score), colnames(metadata)[mta], row.names(my_df)[maxRow], rf$type, #ilr_weight,	rf_imp_se, rf_type,
                    rf$ntree, transf_label, counter, #rf_ntree, trans_group, cycle
                    sep = ","),
              file = main_output_fpath,
              append=TRUE)
          }#if
        }#for loop
  #     },#try catch
  #     error=function(cond) {
  #       print(paste("Opps, an error2 is thrown with", transf_label))
  #       message(paste(transf_label, cond))
  #     },
  #     warning=function(cond) {
  #       print(paste("Opps, a warning2 is thrown with", transf_label))
  #       message(paste(transf_label, cond))
  #     },
  #     finally=.Call(CfreadCleanup)
  #   )
  }#for mta
  print(paste("completed loop:", counter))
}

print("Reading in data from file.")
all_plot_data <- data.frame(read.table(file = main_output_fpath,
            sep = ",", header = TRUE))

unique(all_plot_data$rf_type)

best_seqs <- table(all_plot_data$rf_imp_seq)

write.table(best_seqs, row.names = FALSE,
            file = file.path(output_dir, "tables", paste0("best_seqs_",main_output_label, ".csv")))

for (tg in all_plot_data$trans_group){
  print(tg)
  my_tg_data <- all_plot_data[all_plot_data$trans_group == tg,]
  print(table(my_tg_data$rf_imp_seq))
}

weight_table <- data.frame(tree_type = c(F),
                           metadata = c(F),
                           taxa_pval = c(F),
                           ilr_pval = c(F))
weight_counter <- 1

print("Make all the boxplots")
##-Make all the boxplots--------------------------------------------##
print(paste("Making empty vectors to fill during plot building"))
pval <- c()
pw_name <- c()
iw_name <- c()
metadata_col <- c()
transformation <- c()
distance_metric <- c()
pg_num <- c()

pdf(file = file.path(output_dir, "graphics", paste0("new_bp_", main_output_label, ".pdf")))
index <- 1
for (mta in 1:length(unique(all_plot_data$metadata_col))){
  my_meta <- as.character(unique(all_plot_data$metadata_col)[mta])
  message(my_meta)
  
  plot_data <- all_plot_data[all_plot_data$metadata_col == my_meta,]
  
  not_uniform_tw <- which(plot_data$taxa_weight != "uniform" )
  not_uniform_iw <- which(plot_data$ilr_weight != "uniform" )
  philr_ds <- unique(plot_data$trans_group[c(not_uniform_tw,not_uniform_iw)] ) #pulls out philr only data
  non_philr_ds <- unique(plot_data$trans_group[ !(plot_data$trans_group %in% philr_ds)])
  non_philr_ds_pd <- subset(plot_data, trans_group %in% non_philr_ds)
  philr_ds_pd <- data.frame(subset(plot_data, trans_group %in% philr_ds))
  rownames(philr_ds_pd) <- seq(length=nrow(philr_ds_pd))
  # my_means <- c()
  # for (tg in 1:length(unique(new_pd$trans_group))){
  #   trans_g <- new_pd$trans_group[tg]
  #   my_vals <- which(new_pd$trans_group == trans_g)
  #   my_means <- c(my_means, mean(new_pd$all_score[my_vals]))
  #   names(my_means)[tg] <- trans_g
  # }
  new_pd <- rbind(non_philr_ds_pd, philr_ds_pd)
  new_pd$trans_group <- factor(new_pd$trans_group, levels = c(non_philr_ds, philr_ds))
  
  for(tw in unique(plot_data$taxa_weight)){
    for(iw in unique(plot_data$ilr_weight)){
      philr_pd_tw_iw <- subset(philr_ds_pd, taxa_weight == tw & ilr_weight == iw)
      jitter_pd <- rbind(non_philr_ds_pd, philr_pd_tw_iw)
      jitter_pd$trans_group <- factor(jitter_pd$trans_group, levels = c(non_philr_ds, philr_ds))
      #need to show means from new_pd, but show jitter of tw and iw
      #or could just show selected points but show overal mean for each tw and iw
      back_ground_points <- subset(philr_ds_pd, taxa_weight != tw & ilr_weight != iw)
      bg_jitter <- rbind(non_philr_ds_pd, back_ground_points)
      bg_jitter$trans_group <- factor(bg_jitter$trans_group, levels = c(non_philr_ds, philr_ds))
      g <- ggplot2::ggplot(new_pd, aes(y = all_score, x = trans_group)) + 
        ggplot2::geom_boxplot(data = jitter_pd, color = "blue", alpha = 0.5) +
        ggplot2::geom_boxplot(data = bg_jitter, color = "red", alpha = 0.5) +
        ggplot2::ggtitle(label = paste("Jones", my_meta, ", part_weight:", tw, ", ilr_weight:", iw)) +
        # ggplot2::ggtitle( label = paste("num_tg:", length(unique(new_pd$trans_group)))) +
        ggplot2::theme_classic() +
        ggplot2::scale_x_discrete(guide = guide_axis(angle = 90)) +
        ggplot2::ylab("AUC") +
        ggplot2::xlab("Tree type")
      print(g)
      #build vectors for table
      for (grp in unique(philr_pd_tw_iw$trans_group)){
        # print(grp)
        my_case <- philr_pd_tw_iw[philr_pd_tw_iw$trans_group == grp, ]$all_score
        my_control <- back_ground_points[back_ground_points$trans_group == grp, ]$all_score
        my_test <- t.test(my_case, my_control)
        my_pval <- my_test$p.value
        pval <- c(pval, my_pval)
        pw_name <- c(pw_name, tw)
        iw_name <- c(iw_name, iw)
        metadata_col <- c(metadata_col, my_meta)
        transformation <- c(transformation, grp)
        pg_num <- c(pg_num, index)
      }
      index <- index + 1
    }#end for iw
  }#end for tw
}

dev.off()

dFrame <- data.frame( pval, pw_name, iw_name, metadata_col, transformation, pg_num)
dFrame$adj_pval <- p.adjust(dFrame$pval, method = "BH" )	
dFrame <- dFrame [order(dFrame$adj_pval),]

write.table(dFrame, file=file.path(output_dir, "tables", paste0("new_bp_", main_output_label, ".csv")), 
            sep=",", 
            row.names=FALSE)

print(paste("completed"))

pdf(file = file.path(output_dir, "graphics", paste0("new_bp_NO_WEIGHT", main_output_label, ".pdf")))
for (mta in 1:length(unique(all_plot_data$metadata_col))){
  my_meta <- as.character(unique(all_plot_data$metadata_col)[mta])
  message(my_meta)
  
  plot_data <- all_plot_data[all_plot_data$metadata_col == my_meta,]
  
  not_uniform_tw <- which(plot_data$taxa_weight != "uniform" )
  not_uniform_iw <- which(plot_data$ilr_weight != "uniform" )
  philr_ds <- unique(plot_data$trans_group[c(not_uniform_tw,not_uniform_iw)] ) #pulls out philr only data
  non_philr_ds <- unique(plot_data$trans_group[ !(plot_data$trans_group %in% philr_ds)])
  non_philr_ds_pd <- subset(plot_data, trans_group %in% non_philr_ds)
  philr_ds_pd <- data.frame(subset(plot_data, trans_group %in% philr_ds))
  rownames(philr_ds_pd) <- seq(length=nrow(philr_ds_pd))

  new_pd <- rbind(non_philr_ds_pd, philr_ds_pd)
  new_pd$trans_group <- factor(new_pd$trans_group, levels = c(non_philr_ds, philr_ds))
  
  for(tw in unique(plot_data$taxa_weight)){
    for(iw in unique(plot_data$ilr_weight)){
      g <- ggplot2::ggplot(new_pd, aes(y = all_score, x = trans_group)) + 
        ggplot2::geom_boxplot( color = "black",) +
        ggplot2::ggtitle(label = paste(project, my_meta)) +
        # ggplot2::ggtitle( label = paste("num_tg:", length(unique(new_pd$trans_group)))) +
        ggplot2::theme_classic() +
        ggplot2::scale_x_discrete(guide = guide_axis(angle = 90)) +
        ggplot2::ylab("AUC") +
        ggplot2::xlab("Tree type")
      print(g)
      #build vectors for table
      for (grp in unique(philr_pd_tw_iw$trans_group)){
        # print(grp)
        my_case <- philr_pd_tw_iw[philr_pd_tw_iw$trans_group == grp, ]$all_score
        my_control <- back_ground_points[back_ground_points$trans_group == grp, ]$all_score
        my_test <- t.test(my_case, my_control)
        my_pval <- my_test$p.value
        pval <- c(pval, my_pval)
        pw_name <- c(pw_name, tw)
        iw_name <- c(iw_name, iw)
        metadata_col <- c(metadata_col, my_meta)
        transformation <- c(transformation, grp)
        pg_num <- c(pg_num, index)
      }
      index <- index + 1
    }#end for iw
  }#end for tw
}

dev.off()

pdf(file = file.path(output_dir, "graphics", paste0("new_bp_all_weights", main_output_label, ".pdf")))
for (mta in 1:length(unique(all_plot_data$metadata_col))){
  my_meta <- as.character(unique(all_plot_data$metadata_col)[mta])
  message(my_meta)
  #select plot data for each metadata cat
  plot_data <- all_plot_data[all_plot_data$metadata_col == my_meta,]
  
  # Organize the boxes so that count tables are first
  my_transforms <- unique(plot_data$trans_group)
  my_raw <- my_transforms[grep("raw", my_transforms, ignore.case = TRUE)]
  my_alr <- my_transforms[grep("alr", my_transforms, ignore.case = TRUE)]
  my_clr <- my_transforms[grep("clr", my_transforms, ignore.case = TRUE)]
  my_lognorm <- my_transforms[grep("lognorm", my_transforms, ignore.case = TRUE)]
  rand_philr <- my_transforms[grep("rand_", my_transforms, ignore.case = TRUE)]
  other_DADA2 <- my_transforms[grep("DADA2_", my_transforms, ignore.case = TRUE)]
  my_factors <- c(my_raw, my_alr, my_clr, my_lognorm, other_DADA2, rand_philr)
  plot_data$trans_group <- factor(plot_data$trans_group, levels = c(my_factors))#put non-philr first
  
  g <- ggplot2::ggplot(plot_data, aes(y = all_score, x = trans_group, group=trans_group)) + 
    ggplot2::geom_boxplot() +
    ggplot2::ggtitle(label = paste(project, my_meta)) +
    ggplot2::geom_hline(yintercept = meta_mean, color="red") +
    # scale_fill_discrete(labels=new_labels) +
    # ggplot2::ggtitle( label = paste("num_tg:", length(unique(new_pd$trans_group)))) +
    ggplot2::theme_classic() +
    ggplot2::scale_x_discrete(guide = guide_axis(angle = 90)) +
    ggplot2::ylab("AUC") +
    ggplot2::xlab("Tree type")
  print(g)
}

dev.off()

print("finished R script")



