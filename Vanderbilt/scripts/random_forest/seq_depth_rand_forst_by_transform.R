# Author: Aaron Yerke (aaronyerke@gmail.com)
# Script for testing each transformation against different sequence depth to find the best one
# This was helepful: https://github.com/jsilve24/philr/blob/master/vignettes/philr-intro.Rmd#L142

rm(list = ls()) #clear workspace

##-functions--------------------------------------------------------##
munge_ref_ps <- function(ps){
  #ps must be a phyloseq object
  #requires phyloseq and ape packages to be loaded in the env
  # ps <- filter_taxa(ps, function(x) sum(x > 3) > (0.2*length(x)), TRUE)
  # ps <- filter_taxa(ps, function(x) sd(x)/mean(x) > 3.0, TRUE)
  ps <- transform_sample_counts(ps, function(x) x+1)
  phy_tree(ps) <- makeNodeLabel(phy_tree(ps), method="number", prefix='n')
  return(ps)
}
approx_auc <- function(true_pos, false_pos) {
  #code straight from: https://stackoverflow.com/questions/4903092/calculate-auc-in-r
  return(mean(sample(true_pos,1000,replace=T) > sample(false_pos,1000,replace=T)))
}

##-Load Depencencies------------------------------------------------##
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("ALDEx2", quietly = TRUE)) BiocManager::install("ALDEx2")
if (!requireNamespace("edgeR", quietly = TRUE)) BiocManager::install("edgeR")
if (!requireNamespace("randomForest", quietly = TRUE)) install.packages("randomForest")
if (!require("RColorBrewer")) install.packages("RColorBrewer")
library("compositions")
library("phyloseq")
library("vegan")
library("DESeq2")
library("philr")
library("ape")
library("ALDEx2")
library("ggplot2")
library("RColorBrewer")
#set color palette
palette( brewer.pal(7,"Accent") )
library("randomForest")

##-Establish directory layout---------------------------------------##
home_dir <- file.path('~','git','balance_tree_exploration')
project <- "Vanderbilt"
#home_dir <- file.path('cloud','project')
output_dir <- file.path(home_dir, project, 'output')
setwd(file.path(home_dir))

##-My functions--------------------------------------------------------##
source(file.path(home_dir, "r_libraries", "statistical_functions.R"))
source(file.path(home_dir, "r_libraries", "table_manipulations.R"))

##-Import tables and data preprocessing-----------------------------##
asv_table <- readRDS(file.path(output_dir, "r_objects", "ForwardReads_DADA2.rds"))

ref_ps <- readRDS(file.path(output_dir, "r_objects", "ref_tree_phyloseq_obj.rds"))
asv_tax <- readRDS(file.path(output_dir, "r_objects", "ForwardReads_DADA2_taxonomy.rds"))

metadata <- read.table(file.path(home_dir, project, "patient_metadata.tsv"), 
                       sep="\t", 
                       header=TRUE, 
                       row.names = "Run", 
                       check.names = FALSE,
                       stringsAsFactors=TRUE)
metadata <- metadata[row.names(asv_table),]
metadata$type <- droplevels(metadata$type)

##-----------------Create training/testing sets---------------------##
set.seed(36)
train_index <- sample(x = nrow(metadata), size = 0.75*nrow(metadata), replace=FALSE)
test_index <- c(1:nrow(metadata))[!(1:nrow(metadata) %in% train_index)]


#Plot shannon diversity against log10(total_seqs)
#Plot other normalization methods: otu log 100 and alr and clr 
# my_ds_names <- c( "raw seqs", "clr(raw seqs)", "lognorm raw seqs", "philr ref", "DESeq2", "ALDEx2.clr")
my_ds_names <- c( "raw seqs", "clr(raw seqs)", "lognorm raw seqs", "DESeq2", "ALDEx2.clr")

#lists to hold ROC data for each ds
raw_ROCs <- list()
clr_ROCs <- list()
logrnorm_ROCs <- list()
philr_ROCs <- list()
DESeq_ROCs <- list()
ald_ROCs <- list()
# all_rocs <- list(raw_ROCs, clr_ROCs, logrnorm_ROCs, philr_ROCs, DESeq_ROCs, ald_ROCs)
all_rocs <- list(raw_ROCs, clr_ROCs, logrnorm_ROCs, DESeq_ROCs, ald_ROCs)

# my_ds_names <- c( "raw seqs", "clr(raw seqs)", "lognorm raw seqs")

min_seq_depths <- c(0, 500, 1000, 5000, 10000, 20000, 40000)
mds_depth <- 5
mta = "Treatment"

total_seqs <- rowSums(asv_table)
total_seqs <- data.frame(total_seqs, row.names = row.names(asv_table))


for(s in 1:length(min_seq_depths)){
  seq_d <- min_seq_depths[s]#new sequencing depth
  sd_filt_asv <- asv_table[total_seqs$total_seqs >= seq_d,]#dataset 1
  
  print(paste("sd_filtered dim:", paste(dim(sd_filt_asv))))
  safe_rns <- intersect(row.names(ref_ps@otu_table), row.names(sd_filt_asv)) #rows for this iterate
  new_metadata <- metadata[safe_rns,]
  train_index <- sample(x = nrow(new_metadata), size = 0.75*nrow(new_metadata), replace=FALSE)
  test_index <- c(1:nrow(new_metadata))[!(1:nrow(new_metadata) %in% train_index)]
  
  
  my_clr <- compositions::clr(sd_filt_asv)#dataset 2
  new_tree <- phyloseq::prune_taxa(colnames(sd_filt_asv), ref_ps@phy_tree)#update tree for new phyloseq obj
  new_ref_ps <- phyloseq::prune_samples(safe_rns, ref_ps) #remove non-safe rows from ps
  new_ref_ps <- munge_ref_ps(new_ref_ps)
  print(paste("new dim ref ps:", dim(data.frame(new_ref_ps@otu_table))))
  #create DESeq2 dtaset from new ref ps
  new_DESeq2 <- phyloseq::phyloseq_to_deseq2(new_ref_ps, design= ~ 1)#dataset 5
  new_DESeq2 <- DESeq2::estimateSizeFactors(new_DESeq2)
  new_DESeq2 <- t(counts(new_DESeq2, normalized=T))

  print(paste("new DSeq:", paste(dim(new_DESeq2))))

  # ref_philr <- philr::philr(new_ref_ps@otu_table, new_ref_ps@phy_tree,
  #                           part.weights='enorm.x.gm.counts',
  #                           ilr.weights='blw.sqrt')#dataset 4
  # print(paste("made new philr", dim(as.data.frame(ref_philr))))

  ln_asv <- lognorm(sd_filt_asv)#dataset 6

  ald <- ALDEx2::aldex.clr(sd_filt_asv, mc.samples=12, denom="all", verbose=F)
  ald <-  data.frame(ald@analysisData)
  print(paste("size of ald:", object.size(ald)))
  print(paste("ald dim:", paste(dim(ald))))
  my_datasets <- list(sd_filt_asv, my_clr,  ln_asv, new_DESeq2, ald)
  
  # my_datasets <- list(sd_filt_asv, my_clr,  ln_asv, ref_philr, new_DESeq2, ald)
  print(paste("finished seq depth filter:", s))
  
  for( ds in 1:length(my_datasets)){
    print(my_ds_names[ds])
    my_table <- as.data.frame(my_datasets[ds])
    
    resp_var_test <- new_metadata[,mta][test_index]
    resp_var_train <- new_metadata[,mta][train_index]
    my_table_train <- my_table[train_index,]
    my_table_test <- my_table[test_index,]
    print(dim(my_table))
    rf <- randomForest::randomForest(
      my_table_train, as.factor(resp_var_train)
    )
    my_roc = roc_axes(test_data = my_table_test,
                      true_resp = resp_var_test,
                      ml_model = rf,
                      error_range = 0)
    all_rocs[[ds]][[s]] = my_roc

  }#for ds
  
}#for sd

pdf(file = file.path(output_dir, "graphics", paste0(mta, "_roc_seq_depth.pdf")))
for( ds in 1:length(my_ds_names)){
  
  par(bg = 'grey96')
  plot(true_pos ~ false_pos,
       data = as.data.frame(all_rocs[[ds]][[1]]),
       type = "l",
       xlab = "False positives",
       ylab = "True positives",
       col = palette()[1],
       main = paste(project, mta,"seq_dept:", my_ds_names[ds]),
       ylim = c(0,1))
  for(i in 2:length(min_seq_depths)){
    lines(true_pos ~ false_pos,
          data = as.data.frame(all_rocs[[ds]][[i]]),
          col = palette()[i],
          type = "l")
  }
  abline(
    a = 0,
    b = 1)
  plot.new()
  legend('bottomright', 
         legend = factor(min_seq_depths), 
         col = palette(), 
         pch = 15,
         cex=0.4)
}

for( ds in 1:length(my_ds_names)){
  auc <- c()
  for(i in 1:length(min_seq_depths)){
    my_table <- as.data.frame(all_rocs[[ds]][[i]])
    auc[i] <- approx_auc(my_table$true_pos, my_table$false_pos)
  }
  my_data <- data.frame(min_seq_depths, auc)
  g <- ggplot2::ggplot(data=my_data, aes(x=factor(min_seq_depths), y=auc)) +
    geom_bar(stat="identity") +
    ggplot2::ggtitle(paste0(project, ", Rnmd Frst AUC,", " Transf: ", my_ds_names[ds])) 
  print(g)
}
dev.off()


saveRDS(all_rocs, file.path(output_dir, "r_objects", paste0(mta, "all_rocs_seq_depth_rand_fores.rds")))

# all_rocs <- readRDS(file.path(output_dir, "r_objects", paste0(mta, "all_rocs_seq_depth_rand_fores.rds")))

