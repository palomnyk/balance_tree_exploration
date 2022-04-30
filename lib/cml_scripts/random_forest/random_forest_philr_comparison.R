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
																	just_otu = FALSE){
	#Function for making random forest AUC values
	if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
	if (!requireNamespace("pROC", quietly = TRUE)) BiocManager::install("pROC")
	library("pROC")
	all_auc <- c()
	metadata_col <- c()
	taxa_weight <- c()
	ilr_weight <- c()
	rf_imp_seq = c()
	rf_type = c()
	rf_ntree = c()
	rf_inbag = c()
	for (ilr_w in 1:length(philr_ilr_weights)){
		for (tax_w in 1:length(philr_taxa_weights)){
			tryCatch({
				if (just_otu == TRUE){
					print("In just otu")
					my_table <- ps_obj
				}else{
					my_table <- philr::philr(ps_obj@otu_table, ps_obj@phy_tree, 
																		part.weights = philr_taxa_weights[tax_w],
																		ilr.weights = philr_ilr_weights[ilr_w])
				}
				my_table_train <- data.frame(my_table[train_index,])
				my_table_test <- data.frame(my_table[test_index,])
				},
				error=function(cond) {
					print('Opps, an error1 is thrown')
					message(cond)
					print(message(cond))
				},
				warning=function(cond) {
					print('Opps, a warning1 is thrown')
					message(cond)
					print(message(cond))
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
						resp_var_test <- factor(metadata[,mta][test_index])
						resp_var_train <- factor(metadata[,mta][train_index])
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
						
						if (length(unique(unlist(resp_var_test))) > 2){
							print("multilevels")
							mult_auc <- c()
							for (fact in 1:length(unique(resp_var_test))){#only need to test resp_test
								try({
									my_fact <- as.character(levels(resp_var_test)[fact])
									# print(paste("my_fact:", my_fact))
									dumb_resp_test <- as.factor(replace(as.character(resp_var_test), as.character(resp_var_test) != my_fact, "dumb_var"))
									# print("dumb_resp")
									print(paste(dumb_resp_test))
									dumb_pred <- as.factor(replace(as.character(pred), as.character(pred) != my_fact, "dumb_var"))
									# print("dumb_pred")
									print(paste(dumb_resp_test))
									my_roc <- pROC::roc(as.numeric(dumb_pred), as.numeric(dumb_resp_test))
									# print("my_roc")
									mult_auc <- c(mult_auc, pROC::auc(my_roc))
									print(mult_auc)
								})
							}
							if (length(mult_auc) > 0 ) {
								print("in 'if (length(mult_auc) > 0 )' statement")
								auc <- mean(mult_auc)
							}else{
								break
							}
						}else{
              my_roc <- pROC::roc(as.numeric(pred), as.numeric(resp_var_test))
              print("ROC made")
              auc <- pROC::auc(my_roc)
						}
						# auc <- pROC::auc(my_roc)
						# print(paste("auc: ")
						print("update all output")
						metadata_col <- append(metadata_col, colnames(metadata)[mta])
						all_auc <- append(all_auc, auc)
						taxa_weight <- c(taxa_weight, philr_taxa_weights[tax_w])
						ilr_weight <- c(ilr_weight, philr_ilr_weights[ilr_w])
						##-Section for saving random forest parameters----------------------##
						my_df <- rf$importance
						maxImp <- max(rf$importance)
						maxRow <- which(rf$importance == maxImp)
						rf_imp_seq = c(rf_imp_seq, row.names(my_df)[maxRow])
						rf_type = c(rf_type, rf$type)
						rf_ntree = c(rf_ntree, rf$ntree)
						##------------------------------------------------------------------##
					},
					error=function(cond) {
						print('Opps, an error2 is thrown')
						message(cond)
						print(message(cond))
					},
					warning=function(cond) {
						print('Opps, a warning2 is thrown')
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
										ilr_weight,
										rf_imp_seq,
										rf_type,
										rf_ntree))
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
main_output_label <- paste0("auc_rand_v_ref_v_upgma_v_raw_vert_", num_cycles)
philr_taxa_weights <- c("uniform","gm.counts","anorm","anorm.x.gm.counts","enorm","enorm.x.gm.counts")
philr_ilr_weights <- c("uniform","blw","blw.sqrt","mean.descendants")
random_seed <- 36

##-Import tables and data preprocessing-----------------------------##
asv_table <- data.frame(readRDS(file.path(output_dir, "r_objects", "ForwardReads_DADA2.rds")))
total_seqs <- rowSums(asv_table)
total_seqs <- data.frame("total_seqs"=total_seqs, "duplicate" = total_seqs,
                         row.names = row.names(asv_table))

pdf(file = file.path(output_dir, "graphics", paste0("trees_", main_output_label, ".pdf")))
print("Cleaning Ref tree otu with philr tutorial normalization")
ref_ps <- readRDS(file.path(output_dir, "r_objects", "ref_tree_phyloseq_obj.rds"))
phy_tree(ref_ps) <- ape::makeNodeLabel(phy_tree(ref_ps), method="number", prefix='n')
phyloseq::plot_tree(ref_ps, method = "treeonly", nodelabf=nodeplotblank, title = paste0("orig_ref"))

ref_ps_clean <- raw_ps_to_clean_ps(ref_ps)
phyloseq::plot_tree(ref_ps_clean, method = "treeonly", nodelabf=nodeplotblank, title = paste0("cln_ref"))
print("Cleaning UPGMA tree otu with philr tutorial normalization")

denovo_tree_ps <- readRDS(file.path(output_dir, "r_objects", "denovo_tree_UPGMA_phyloseq_obj.rds"))
phy_tree(denovo_tree_ps) <- ape::makeNodeLabel(phy_tree(denovo_tree_ps), method="number", prefix='n')
phyloseq::plot_tree(denovo_tree_ps, method = "treeonly", nodelabf=nodeplotblank, title = paste0("orig_upgma"))

cln_denovo_tree_ps <- raw_ps_to_clean_ps(denovo_tree_ps)
denovo_tree_ps <- phyloseq::transform_sample_counts(denovo_tree_ps, function(x) x + 1 )
phyloseq::plot_tree(cln_denovo_tree_ps, method = "treeonly", nodelabf=nodeplotblank, title = paste0("cln_upgma"))

iqtree_orig_ps <- readRDS(file.path(output_dir, "r_objects", "denovo_tree_iqtree_phyloseq_obj.rds"))
phy_tree(iqtree_orig_ps) <- ape::makeNodeLabel(phy_tree(iqtree_orig_ps), method="number", prefix='n')
phyloseq::plot_tree(iqtree_orig_ps, method = "treeonly", nodelabf=nodeplotblank, title = paste0("orig_iqtree"))

cln_iqtree_ps <- raw_ps_to_clean_ps(iqtree_orig_ps)
phyloseq::plot_tree(cln_iqtree_ps, method = "treeonly", nodelabf=nodeplotblank, title = paste0("cln_iqtree"))
dev.off()

print("loading and munging metadata")
metadata <- read.table(opt$metadata, 
                       sep=opt$metadata_delim, 
                       header=TRUE, 
                       row.names = opt$metadata_rowname, 
                       check.names = FALSE,
                       stringsAsFactors=TRUE)
print("attempting to equalize metadata rows to seq data rows")
my_rows <- row.names(metadata) %in% row.names(data.frame(ref_ps@otu_table@.Data))
metadata <- data.frame(metadata[my_rows, ])
rf_cols <- 1:ncol(metadata)

print("Attempting to read HashSeq count table")
hashseq <- data.frame(data.table::fread(file = file.path(output_dir,"hashseq", "SvTable.txt"),
                                        header=TRUE, data.table=FALSE), row.names = 1)
hashseq <- hashseq[, colSums(hashseq != 0) > 0.1*nrow(hashseq)]
print("attempting to equalize hashseq rows to seq data rows")
my_names <- sapply(as.character(row.names(hashseq)), function(x) {
  my_val = strsplit(x,"_")[[1]][1]
  return(my_val)
  })
row.names(hashseq) <- my_names
print
hashseq <- hashseq[my_rows,]

if (!base::identical(row.names(hashseq), row.names(metadata))){
  print("Problem with hashseq rownames.")
  quit_due_row_names()
}

##-Random num seed--------------------------------------------------##
print(paste("Setting random seed to:", random_seed))
set.seed(random_seed)
print("making random trees")
orig_ref_rand_list <- list()
for (rand in 1:10){
  rand_tree <- ape::rtree(n = length(ref_ps@phy_tree$tip.label), tip.label = ref_ps@phy_tree$tip.label)
  #put int in philr
  rand_tree_ps <- phyloseq::phyloseq( otu_table(ref_ps_clean@otu_table, taxa_are_rows = F),
                                      phy_tree(rand_tree),
                                      tax_table(ref_ps@tax_table),
                                      sample_data(ref_ps@sam_data))
  phy_tree(rand_tree_ps) <- ape::makeNodeLabel(phy_tree(rand_tree_ps), method="number", prefix='n')
  phyloseq::plot_tree(rand_tree_ps,  method = "treeonly", nodelabf=nodeplotblank, title = paste0("orig_ref_rand_", rand))
  orig_ref_rand_list[[rand]] <- rand_tree_ps
}

print("make random trees for cln upgma taxa")
cln_upgma_rand_list <- list()
for (rand in 1:10){
  rand_tree <- ape::rtree(n = length(cln_denovo_tree_ps@phy_tree$tip.label), tip.label = cln_denovo_tree_ps@phy_tree$tip.label)
  #put int in philr
  rand_tree_ps <- phyloseq::phyloseq( otu_table(cln_denovo_tree_ps@otu_table, taxa_are_rows = F),
                                      phy_tree(rand_tree),
                                      tax_table(cln_denovo_tree_ps@tax_table),
                                      sample_data(cln_denovo_tree_ps@sam_data))
  phy_tree(rand_tree_ps) <- ape::makeNodeLabel(phy_tree(rand_tree_ps), method="number", prefix='n')
  phyloseq::plot_tree(rand_tree_ps, title = paste0("cln_upgma_rand_", rand))
  cln_upgma_rand_list[[rand]] <- rand_tree_ps
}

print("make random trees for clean ref taxa")
cln_ref_rand_list <- list()
for (rand in 1:10){
  rand_tree <- ape::rtree(n = length(ref_ps_clean@phy_tree$tip.label), tip.label = ref_ps_clean@phy_tree$tip.label)
  #put int in philr
  rand_tree_ps <- phyloseq::phyloseq(otu_table(ref_ps_clean@otu_table, taxa_are_rows = F),
                                     phy_tree(rand_tree),
                                     tax_table(ref_ps_clean@tax_table),
                                     sample_data(ref_ps_clean@sam_data))
  phy_tree(rand_tree_ps) <- ape::makeNodeLabel(phy_tree(rand_tree_ps), method="number", prefix='n')
  phyloseq::plot_tree(rand_tree_ps, method = "treeonly", nodelabf=nodeplotblank, title = paste0("cln_ref_rand_", rand))
  cln_ref_rand_list[[rand]] <- rand_tree_ps
}
print("make random trees for iqtree clean")
iqtree_clean_rand_list <- list()
for (rand in 1:10){
  rand_tree <- ape::rtree(n = length(cln_iqtree_ps@phy_tree$tip.label), tip.label = cln_iqtree_ps@phy_tree$tip.label)
  #put int in philr
  rand_tree_ps <- phyloseq::phyloseq(otu_table(cln_iqtree_ps@otu_table, taxa_are_rows = F),
                                     phy_tree(rand_tree),
                                     tax_table(cln_iqtree_ps@tax_table),
                                     sample_data(cln_iqtree_ps@sam_data))
  phy_tree(rand_tree_ps) <- ape::makeNodeLabel(phy_tree(rand_tree_ps), method="number", prefix='n')
  phyloseq::plot_tree(rand_tree_ps, method = "treeonly", nodelabf=nodeplotblank, title = paste0("iqtree_rand_", rand))
  iqtree_clean_rand_list[[rand]] <- rand_tree_ps
}

print("creating lognorm, ALR and CLR")
if (dir.exists(file.path(output_dir,"r_objects", "lognorm_asv.rds"))) {
  ln_asv_tab <- readRDS(file.path(output_dir,"r_objects", "lognorm_asv.rds"))
}else{
  ln_asv_tab <- lognorm(asv_table)
  saveRDS(ln_asv_tab, file = file.path(output_dir,"r_objects", "lognorm_asv.rds"))
}

my_zeros <- apply(asv_table, 2, function(x) {
  return(sum(x == 0))
})
alr_col <- which(my_zeros == min(my_zeros))[1]
# alr_col_num <- grep(alr_col, colnames(asv_table))
print("creating ALR")
if (file.exists(file.path(output_dir,"r_objects", "alr_asv.rds"))) {
  my_alr <- readRDS(file.path(output_dir,"r_objects", "alr_asv.rds"))
}else{
  my_alr <- as.data.frame(rgr::alr(as.matrix(asv_table + 1), j = as.numeric(alr_col)))
  saveRDS(my_alr, file = file.path(output_dir,"r_objects", "alr_asv.rds"))
}
print("creating CLR")
if (dir.exists(file.path(output_dir,"r_objects", "clr_asv.rds"))) {
  my_clr <- readRDS(file.path(output_dir,"r_objects", "clr_asv.rds"))
}else{
  my_clr <- as.data.frame(rgr::clr(as.matrix(asv_table + 1)))
  saveRDS(my_clr, file = file.path(output_dir,"r_objects", "clr_asv.rds"))
}

##-Create plot data-------------------------------------------------##
all_plot_data <- data.frame(all_auc = c(),
                            metadata_col = c(),
                            taxa_weight = c(),
                            ilr_weight = c(),
                            random_batch = c(),
                            trans_group = c(),
                            rf_imp_seq = c(),
                            rf_type = c(),
                            rf_ntree = c())
skips <- 0
counter <- 0

print(paste("counter:", counter, " entering while loop"))
while (counter < num_cycles & skips < 5){
  ##-Create training/testing sets-------------------------------------##
  train_index <- sample(x = nrow(metadata), size = 0.75*nrow(metadata), replace=FALSE)
  test_index <- c(1:nrow(metadata))[!(1:nrow(metadata) %in% train_index)]
  print(paste("counter:", counter, " making ref cln random AUC"))
  for( rand_ps in 1:length(cln_ref_rand_list)){
    rand_tree_ps <- cln_ref_rand_list[[rand_ps]]
    rand_plot_data <- make_ilr_taxa_auc_df(ps_obj = rand_tree_ps,
                                            metadata_cols = rf_cols,
                                            metadata = metadata,
                                            train_index = train_index,
                                            test_index = test_index,
                                            philr_ilr_weights = philr_ilr_weights,
                                            philr_taxa_weights = philr_taxa_weights)
    rand_plot_data$trans_group <- rep(paste0("Culled_filtered_random_Silva_PhILR_", rand_ps), nrow(rand_plot_data))
    rand_plot_data$random_batch <- rep(rand_ps, nrow(rand_plot_data))
    all_plot_data <- rbind(all_plot_data, rand_plot_data)
  }

  print(paste("counter:", counter, " making orig ref random AUC"))
  for( rand_ps in 1:length(orig_ref_rand_list)){
    rand_tree_ps <- orig_ref_rand_list[[rand_ps]]
    rand_plot_data <- make_ilr_taxa_auc_df(ps_obj = rand_tree_ps,
                                            metadata_cols = rf_cols,
                                            metadata = metadata,
                                            train_index = train_index,
                                            test_index = test_index,
                                            philr_ilr_weights = philr_ilr_weights,
                                            philr_taxa_weights = philr_taxa_weights)
    rand_plot_data$trans_group <- rep(paste0("Culled_random_Silva_PhILR_", rand_ps), nrow(rand_plot_data))
    rand_plot_data$random_batch <- rep(rand_ps, nrow(rand_plot_data))
    all_plot_data <- rbind(all_plot_data, rand_plot_data)
  }
  print(paste("counter:", counter, " making random cleaned upgma AUC"))
  for( rand_ps in 1:length(cln_upgma_rand_list)){
    rand_tree_ps <- cln_upgma_rand_list[[rand_ps]]
    rand_plot_data <- make_ilr_taxa_auc_df(ps_obj = rand_tree_ps,
                                            metadata_cols = rf_cols,
                                            metadata = metadata,
                                            train_index = train_index,
                                            test_index = test_index,
                                            philr_ilr_weights = philr_ilr_weights,
                                            philr_taxa_weights = philr_taxa_weights)
    rand_plot_data$trans_group <- rep(paste0("Filtered_random_UPGMA_PhILR_", rand_ps), nrow(rand_plot_data))
    rand_plot_data$random_batch <- rep(rand_ps, nrow(rand_plot_data))
    all_plot_data <- rbind(all_plot_data, rand_plot_data)
  }
  print(paste("counter:", counter, " making random cleaned iqtree AUC"))
  for( rand_ps in 1:length(iqtree_clean_rand_list)){
    rand_tree_ps <- iqtree_clean_rand_list[[rand_ps]]
    rand_plot_data <- make_ilr_taxa_auc_df(ps_obj = rand_tree_ps,
                                            metadata_cols = rf_cols,
                                            metadata = metadata,
                                            train_index = train_index,
                                            test_index = test_index,
                                            philr_ilr_weights = philr_ilr_weights,
                                            philr_taxa_weights = philr_taxa_weights)
    rand_plot_data$trans_group <- rep(paste0("Filtered_random_IQTree_PhILR_", rand_ps), nrow(rand_plot_data))
    rand_plot_data$random_batch <- rep(rand_ps, nrow(rand_plot_data))
    all_plot_data <- rbind(all_plot_data, rand_plot_data)
  }

  print(paste("counter:", counter, " making ref cln tree philr AUC"))
  my_plot_data <- make_ilr_taxa_auc_df(ps_obj = ref_ps_clean,
                                        metadata_cols = rf_cols,
                                        metadata = metadata,
                                        train_index = train_index,
                                        test_index = test_index,
                                        philr_ilr_weights = philr_ilr_weights,
                                        philr_taxa_weights = philr_taxa_weights)
  my_plot_data$trans_group <- rep("Culled_filtered_Silva_PhILR", nrow(my_plot_data))
  my_plot_data$random_batch <- rep("None", nrow(my_plot_data))
  all_plot_data <- rbind(all_plot_data, my_plot_data)

  print(paste("counter:", counter, " making ref orig tree philr AUC"))
  my_plot_data <- make_ilr_taxa_auc_df(ps_obj = ref_ps_clean,
                                        metadata_cols = rf_cols,
                                        metadata = metadata,
                                        train_index = train_index,
                                        test_index = test_index,
                                        philr_ilr_weights = philr_ilr_weights,
                                        philr_taxa_weights = philr_taxa_weights)
  my_plot_data$trans_group <- rep("Culled_Silva_PhILR", nrow(my_plot_data))
  my_plot_data$random_batch <- rep("None", nrow(my_plot_data))
  all_plot_data <- rbind(all_plot_data, my_plot_data)

  print(paste("counter:", counter, " making orig ref no trees AUC"))
  my_plot_data <- make_ilr_taxa_auc_df(ps_obj = as.data.frame(ref_ps@otu_table),
                                        metadata_cols = rf_cols,
                                        metadata = metadata,
                                        train_index = train_index,
                                        test_index = test_index,
                                        philr_ilr_weights = philr_ilr_weights,
                                        philr_taxa_weights = philr_taxa_weights,
                                        just_otu = TRUE)
  my_plot_data$trans_group <- rep("Culled_Silva_counts_table", nrow(my_plot_data))
  my_plot_data$random_batch <- rep("None", nrow(my_plot_data))
  all_plot_data <- rbind(all_plot_data, my_plot_data)

  print(paste("counter:", counter, " making clean ref no trees AUC"))
  my_plot_data <- make_ilr_taxa_auc_df(ps_obj = as.data.frame(ref_ps_clean@otu_table),
                                        metadata_cols = rf_cols,
                                        metadata = metadata,
                                        train_index = train_index,
                                        test_index = test_index,
                                        philr_ilr_weights = philr_ilr_weights,
                                        philr_taxa_weights = philr_taxa_weights,
                                        just_otu = TRUE)
  my_plot_data$trans_group <- rep("Culled_Filtered_Silva_counts_table", nrow(my_plot_data))
  my_plot_data$random_batch <- rep("None", nrow(my_plot_data))
  all_plot_data <- rbind(all_plot_data, my_plot_data)

  print(paste("counter:", counter, " making seq only clean UPGMA (no trees) AUC"))
  my_plot_data <- make_ilr_taxa_auc_df(ps_obj = as.data.frame(cln_denovo_tree_ps@otu_table),
                                        metadata_cols = rf_cols,
                                        metadata = metadata,
                                        train_index = train_index,
                                        test_index = test_index,
                                        philr_ilr_weights = philr_ilr_weights,
                                        philr_taxa_weights = philr_taxa_weights,
                                        just_otu = TRUE)
  my_plot_data$trans_group <- rep("Filtered_UPGMA_counts_table", nrow(my_plot_data))
  my_plot_data$random_batch <- rep("None", nrow(my_plot_data))
  all_plot_data <- rbind(all_plot_data, my_plot_data)

  print(paste("counter:", counter, " making iqtree cln tree philr AUC"))
  my_plot_data <- make_ilr_taxa_auc_df(ps_obj = cln_iqtree_ps,
                                        metadata_cols = rf_cols,
                                        metadata = metadata,
                                        train_index = train_index,
                                        test_index = test_index,
                                        philr_ilr_weights = philr_ilr_weights,
                                        philr_taxa_weights = philr_taxa_weights)
  my_plot_data$trans_group <- rep("Filtered_IQTree_PhILR", nrow(my_plot_data))
  my_plot_data$random_batch <- rep("None", nrow(my_plot_data))
  all_plot_data <- rbind(all_plot_data, my_plot_data)

  print(paste("counter:", counter, " making iqtree orig tree philr AUC"))
  my_plot_data <- make_ilr_taxa_auc_df(ps_obj = iqtree_orig_ps,
                                        metadata_cols = rf_cols,
                                        metadata = metadata,
                                        train_index = train_index,
                                        test_index = test_index,
                                        philr_ilr_weights = philr_ilr_weights,
                                        philr_taxa_weights = philr_taxa_weights)
  my_plot_data$trans_group <- rep("IQTREE_Orig_PhILR", nrow(my_plot_data))
  my_plot_data$random_batch <- rep("None", nrow(my_plot_data))
  all_plot_data <- rbind(all_plot_data, my_plot_data)

  print(paste("counter:", counter, " making orig ref no trees AUC"))
  my_plot_data <- make_ilr_taxa_auc_df(ps_obj = as.data.frame(iqtree_orig_ps@otu_table),
                                        metadata_cols = rf_cols,
                                        metadata = metadata,
                                        train_index = train_index,
                                        test_index = test_index,
                                        philr_ilr_weights = philr_ilr_weights,
                                        philr_taxa_weights = philr_taxa_weights,
                                        just_otu = TRUE)
  my_plot_data$trans_group <- rep("Orig_IQTree_counts_table", nrow(my_plot_data))
  my_plot_data$random_batch <- rep("None", nrow(my_plot_data))
  all_plot_data <- rbind(all_plot_data, my_plot_data)

  print(paste("counter:", counter, " making clean iqtree no trees AUC"))
  my_plot_data <- make_ilr_taxa_auc_df(ps_obj = as.data.frame(cln_iqtree_ps@otu_table),
                                        metadata_cols = rf_cols,
                                        metadata = metadata,
                                        train_index = train_index,
                                        test_index = test_index,
                                        philr_ilr_weights = philr_ilr_weights,
                                        philr_taxa_weights = philr_taxa_weights,
                                        just_otu = TRUE)
  my_plot_data$trans_group <- rep("Filtered_IQTree_counts_table", nrow(my_plot_data))
  my_plot_data$random_batch <- rep("None", nrow(my_plot_data))
  all_plot_data <- rbind(all_plot_data, my_plot_data)

  # print("making UPGMA AUC")
  # my_plot_data <- make_ilr_taxa_auc_df( ps_obj = denovo_tree_ps,
  #                                           metadata_cols = rf_cols,
  #                                           metadata = metadata,
  #                                           train_index = train_index,
  #                                           test_index = test_index,
  #                                           philr_ilr_weights = philr_ilr_weights,
  #                                           philr_taxa_weights = philr_taxa_weights)
  # my_plot_data$trans_group <- rep("UPGMA", nrow(my_plot_data))
  # all_plot_data <- rbind(all_plot_data, my_plot_data)

  print(paste("counter:", counter, "making clean upgma AUC"))
  my_plot_data <- make_ilr_taxa_auc_df( ps_obj = cln_denovo_tree_ps,
                                            metadata_cols = rf_cols,
                                            metadata = metadata,
                                            train_index = train_index,
                                            test_index = test_index,
                                            philr_ilr_weights = philr_ilr_weights,
                                            philr_taxa_weights = philr_taxa_weights)
  my_plot_data$random_batch <- rep("None", nrow(my_plot_data))
  my_plot_data$trans_group <- rep("Filtered_UPGMA_PhILR", nrow(my_plot_data))
  all_plot_data <- rbind(all_plot_data, my_plot_data)

  print(paste("counter:", counter, " making seq only orig ref (no trees) AUC"))
  my_plot_data <- make_ilr_taxa_auc_df( ps_obj = as.data.frame(ref_ps@otu_table),
                                        metadata_cols = rf_cols,
                                        metadata = metadata,
                                        train_index = train_index,
                                        test_index = test_index,
                                        philr_ilr_weights = philr_ilr_weights,
                                        philr_taxa_weights = philr_taxa_weights,
                                        just_otu = TRUE)
  my_plot_data$random_batch <- rep("None", nrow(my_plot_data))
  my_plot_data$trans_group <- rep("Culled_Silva_counts_table", nrow(my_plot_data))
  all_plot_data <- rbind(all_plot_data, my_plot_data)

  print(paste("counter:", counter, " generate ", "'raw data' data"))
  my_plot_data <- make_ilr_taxa_auc_df(ps_obj = asv_table,
                                        metadata_cols = rf_cols,
                                        metadata = metadata,
                                        train_index = train_index,
                                        test_index = test_index,
                                        philr_ilr_weights = philr_ilr_weights,
                                        philr_taxa_weights = philr_taxa_weights,
                                        just_otu = TRUE )
  my_plot_data$random_batch <- rep("None", nrow(my_plot_data))
  my_plot_data$trans_group <- rep("Raw_Dada2_counts_table", nrow(my_plot_data))
  all_plot_data <- rbind(all_plot_data, my_plot_data)

  # print(paste("counter:", counter, " generate ", "read depth data"))
  # my_plot_data <- make_ilr_taxa_auc_df(ps_obj = data.frame(total_seqs),
  #                                       metadata_cols = rf_cols,
  #                                       metadata = metadata,
  #                                       train_index = train_index,
  #                                       test_index = test_index,
  #                                       philr_ilr_weights = philr_ilr_weights,
  #                                       philr_taxa_weights = philr_taxa_weights,
  #                                       just_otu = TRUE )
  # my_plot_data$trans_group <- rep("read_depth", nrow(my_plot_data))
  # all_plot_data <- rbind(all_plot_data, my_plot_data)

  print(paste("counter:", counter, " generate ", "lognorm data"))
  my_plot_data <- make_ilr_taxa_auc_df(ps_obj = ln_asv_tab,
                                        metadata_cols = rf_cols,
                                        metadata = metadata,
                                        train_index = train_index,
                                        test_index = test_index,
                                        philr_ilr_weights = philr_ilr_weights,
                                        philr_taxa_weights = philr_taxa_weights,
                                        just_otu = TRUE )
  my_plot_data$random_batch <- rep("None", nrow(my_plot_data))
  my_plot_data$trans_group <- rep("lognorm", nrow(my_plot_data))
  all_plot_data <- rbind(all_plot_data, my_plot_data)

  print(paste("counter:", counter, " generate alr data"))
  my_plot_data <- make_ilr_taxa_auc_df(ps_obj = my_alr,
                                        metadata_cols = rf_cols,
                                        metadata = metadata,
                                        train_index = train_index,
                                        test_index = test_index,
                                        philr_ilr_weights = philr_ilr_weights,
                                        philr_taxa_weights = philr_taxa_weights,
                                        just_otu = TRUE )
  my_plot_data$random_batch <- rep("None", nrow(my_plot_data))
  my_plot_data$trans_group <- rep("alr", nrow(my_plot_data))
  all_plot_data <- rbind(all_plot_data, my_plot_data)

  print(paste("counter:", counter, " generate clr data"))
  my_plot_data <- make_ilr_taxa_auc_df(ps_obj = my_clr,
                                        metadata_cols = rf_cols,
                                        metadata = metadata,
                                        train_index = train_index,
                                        test_index = test_index,
                                        philr_ilr_weights = philr_ilr_weights,
                                        philr_taxa_weights = philr_taxa_weights,
                                        just_otu = TRUE )
  my_plot_data$random_batch <- rep("None", nrow(my_plot_data))
  my_plot_data$trans_group <- rep("clr", nrow(my_plot_data))
  all_plot_data <- rbind(all_plot_data, my_plot_data)
  
  print(paste("counter:", counter, " generate hashseq data"))
  my_plot_data <- make_ilr_taxa_auc_df(ps_obj = hashseq,
                                       metadata_cols = rf_cols,
                                       metadata = metadata,
                                       train_index = train_index,
                                       test_index = test_index,
                                       philr_ilr_weights = philr_ilr_weights,
                                       philr_taxa_weights = philr_taxa_weights,
                                       just_otu = TRUE )
  my_plot_data$random_batch <- rep("None", nrow(my_plot_data))
  my_plot_data$trans_group <- rep("HashSeq", nrow(my_plot_data))
  all_plot_data <- rbind(all_plot_data, my_plot_data)
  
  print(paste("completed loop:", counter))
  counter <- counter + 1
  skips <- 0
}

print(paste0("saving plot data to hardrive:\n",
             as.character(file.path(output_dir, "tables", paste0(main_output_label, ".csv")))))
write.table(all_plot_data,
            file = file.path(output_dir, "tables", paste0(main_output_label, ".csv")),
            sep = ",",
            row.names = FALSE)

all_plot_data <- data.frame(read.table(file = file.path(output_dir, "tables", 
                                             paste0(main_output_label, ".csv")),
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
pdf(file = file.path(output_dir, "graphics", paste0(main_output_label, ".pdf")))
g <- ggplot2::ggplot(all_plot_data, aes(y = all_auc, x= trans_group)) + 
  ggplot2::geom_boxplot() +
  ggplot2::geom_jitter(aes(color = as.factor(ilr_weight)),width = 0.2, height = 0.001) +
  ggplot2::ggtitle(paste(project, "all metadata")) +
  # ggplot2::geom_hline(yintercept = 0) +
  ggplot2::theme(axis.text.x = element_text(angle = 45),
                 axis.text = element_text(size = 20)) +
  ggplot2::theme_classic() +
  # ggplot2::scale_y_discrete(labels = seq(0, 1, by = 0.2)) +
  ggplot2::ylab("AUC") +
  ggplot2::xlab("Transformation") +
  ggplot2::labs(color = "ilr weight")
print(g)

g <- ggplot2::ggplot(all_plot_data, aes(y = all_auc, x= trans_group)) + 
  ggplot2::geom_boxplot() +
  ggplot2::geom_jitter(aes(color = as.factor(taxa_weight)),width = 0.2, height = 0.001) +
  ggplot2::ggtitle(paste(project, "all metadata")) +
  # ggplot2::geom_hline(yintercept = 0) +
  ggplot2::theme(axis.text.x = element_text(angle = 45)) +
  ggplot2::theme_classic() +
  # ggplot2::scale_y_discrete(labels = seq(0, 1, by = 0.2)) +
  ggplot2::ylab("AUC") +
  ggplot2::xlab("Transformation") +
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
    # for (tg in unique(plot_data$trans_group)){
    #   auc <- plot_data$all_auc[plot_data$trans_group ==  tg]
    #   ilr_w <- plot_data$ilr_weight[plot_data$trans_group ==  tg]
    #   taxa_w <- plot_data$taxa_weight[plot_data$trans_group ==  tg]
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
    
    non_philr_ds <- c("cln_upgma_taxa_only", "Silva_ref_cln_taxa_only",
                      "cln_upgma_taxa_only", "raw_data", "lognorm", "alr", 
                      "clr", "Silva_ref_orig_taxa_only")
    
    not_uniform_tw <- which(plot_data$taxa_weight != "uniform" )
    not_uniform_iw <- which(plot_data$ilr_weight != "uniform" )
    philr_ds <- unique(plot_data$trans_group[c(not_uniform_tw,not_uniform_iw)] ) #pulls out philr only data
    non_philr_ds <- unique(plot_data$trans_group[ !(plot_data$trans_group %in% philr_ds)])
    non_philr_ds_pd <- subset(plot_data, trans_group %in% non_philr_ds)
    philr_ds_pd <- subset(plot_data, trans_group %in% philr_ds)
    
    # new_pd <- rbind(non_philr_ds_pd, philr_ds_pd)
    
    for(tw in unique(plot_data$taxa_weight)){
    my_phil_ds_pd <- philr_ds_pd[philr_ds_pd$taxa_weight == tw,]
    new_pd <- rbind(non_philr_ds_pd, philr_ds_pd)
    g <- ggplot2::ggplot(new_pd, aes(y = all_auc, x = trans_group)) + 
      ggplot2::geom_boxplot() +
      ggplot2::geom_jitter(aes(color = as.factor(ilr_weight)), size = 0.75) +
      ggplot2::ggtitle(paste(project, my_meta, "only taxa weight:", tw, "| col by ilr")) +
      ggplot2::theme_classic() +
      ggplot2::scale_x_discrete(guide = guide_axis(angle = 90)) +
      ggplot2::ylab("AUC") +
      ggplot2::xlab("Tree type") +
      ggplot2::labs(color = "ilr weight") 
      # ggpubr::stat_compare_means(comparisons = my_comparisons)
    print(g)
  }
  for(iw in unique(plot_data$ilr_weight)){
    my_phil_ds_pd <- philr_ds_pd[philr_ds_pd$ilr_weight == iw,]
    new_pd <- rbind(non_philr_ds_pd, my_phil_ds_pd)
    g <- ggplot2::ggplot(new_pd, aes(y = all_auc, x = trans_group)) + 
      ggplot2::geom_boxplot() +
      ggplot2::geom_jitter(aes(color = as.factor(taxa_weight)), size = 0.75) +
      ggplot2::ggtitle(paste(project, my_meta, "only ilr", iw, "| col by taxa weight")) +
      ggplot2::theme_classic() +
      ggplot2::scale_x_discrete(guide = guide_axis(angle = 90)) +
      ggplot2::ylab("AUC") +
      ggplot2::xlab("Tree type") +
      ggplot2::labs(color = "ilr weight") 
      # ggpubr::stat_compare_means(comparisons = my_comparisons)
    print(g)
  }
  
  # betwn_bar_anova <- anova(lm(data = plot_data,all_auc ~ trans_group))
  # my_comparisons <- list( c("raw_data", "Silva_ref"), c( "UPGMA", "Silva_ref"), c("Silva_ref", "random"), c("raw_data", "random") )
}

dev.off()

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
  #   my_means <- c(my_means, mean(new_pd$all_auc[my_vals]))
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
      g <- ggplot2::ggplot(new_pd, aes(y = all_auc, x = trans_group)) + 
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
        my_case <- philr_pd_tw_iw[philr_pd_tw_iw$trans_group == grp, ]$all_auc
        my_control <- back_ground_points[back_ground_points$trans_group == grp, ]$all_auc
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
      g <- ggplot2::ggplot(new_pd, aes(y = all_auc, x = trans_group)) + 
        ggplot2::geom_boxplot( color = "blue",) +
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
        my_case <- philr_pd_tw_iw[philr_pd_tw_iw$trans_group == grp, ]$all_auc
        my_control <- back_ground_points[back_ground_points$trans_group == grp, ]$all_auc
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

print("finished R script")
