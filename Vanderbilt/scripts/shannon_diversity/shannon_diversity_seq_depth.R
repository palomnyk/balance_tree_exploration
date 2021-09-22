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

##-Load Depencencies------------------------------------------------##
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("ALDEx2", quietly = TRUE)) BiocManager::install("ALDEx2")
if (!requireNamespace("reshape2", quietly = TRUE)) BiocManager::install("reshape2")
library("compositions")
library("phyloseq")
library("vegan")
library("DESeq2")
library("philr")
library("ape")
library("ALDEx2")
library("ggplot2")
library("reshape2")

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
asv_table <- readRDS(file.path(output_dir, "r_objects", "ForwardReads_DADA2.rds"))

ref_ps <- readRDS(file.path(output_dir, "r_objects", "ref_tree_phyloseq_obj.rds"))
asv_tax <- readRDS(file.path(output_dir, "r_objects", "ForwardReads_DADA2_taxonomy.rds"))

metadata <- read.table(file.path(home_dir, project, "patient_metadata.tsv"), 
                       sep="\t", 
                       header=TRUE, 
                       row.names = "Run", 
                       check.names = FALSE,
                       stringsAsFactors=TRUE)
metadata <- metadata[row.names(metadata) %in% row.names(asv_table), ]
metadata$type <- droplevels(metadata$type)
metadata$type <- factor(metadata$type)

metad_cols <- 1:2

#Plot shannon diversity against log10(total_seqs)
#Plot other normalization methods: otu log 100 and alr and clr 
min_seq_depths <- c(0, 500, 1000, 5000, 10000, 20000, 40000)
mds_depth <- 5

total_seqs <- rowSums(asv_table)
total_seqs <- data.frame(total_seqs, row.names = row.names(asv_table))

kend <- vector(mode = "numeric", length = length(min_seq_depths) * mds_depth)
perma_r2 <- vector(mode = "numeric", length = length(min_seq_depths) * mds_depth)
mds_lev <- vector(mode = "integer", length = length(min_seq_depths) * mds_depth)
seq_depth <- vector(mode = "integer", length = length(min_seq_depths) * mds_depth)
var_exp <- vector(mode = "integer", length = length(min_seq_depths) * mds_depth)
PCAi_shannon_div_spear_cor <- vector(mode = "integer", length =length(min_seq_depths) * mds_depth)
samples_left <- vector(mode = "integer", length = length(min_seq_depths) * mds_depth)
taxa_left <- vector(mode = "integer", length = length(min_seq_depths) * mds_depth)
zero_count <- vector(mode = "integer", length = length(min_seq_depths) * mds_depth)

meta_cor <- list()
for (i in metad_cols){
  meta_cor[[i]] <- vector(mode = "integer", length = length(min_seq_depths) * mds_depth)
  names(meta_cor)[i] <- colnames(metadata)[i]
}
meta_cor <- data.frame(meta_cor)

# pdf(file = file.path(output_dir, "graphics", "seq_depth_artifact_PCA12345_scatter.pdf"))
counter <- 1
for(s in 1:length(min_seq_depths)){
  seq_d <- min_seq_depths[s]#new sequencing depth
  sd_filt_asv <- asv_table[total_seqs$total_seqs >= seq_d,]#dataset 1
  print(paste("sd_filtered dim:", paste(dim(sd_filt_asv))))
  safe_rns <- intersect(row.names(ref_ps@otu_table), row.names(sd_filt_asv)) #rows for this iterate
  shan_dirv <- vegan::diversity(sd_filt_asv[safe_rns,])
  
  print(paste("finished seq depth filter:", s))
  
  ##-Create a PCA-----------------------------------------------------##
  my_prcmp <- prcomp(sd_filt_asv, 
                     center = TRUE,
                     rank = mds_depth)#,
  # scale = TRUE)
  ##-Extract PCA matrix and convert to dataframe----------------------##
  myPCA <- data.frame(my_prcmp$x)
  my_var_exp <- my_prcmp$sdev^2/sum(my_prcmp$sdev^2)
  shan_div <- vegan::diversity(sd_filt_asv)
  for (md in 1:mds_depth){
    for (m in 1:length(metad_cols)){
      print(m)
      meta_cor[counter, m] <- cor.test(shan_div, metadata[safe_rns,m],
                                       method = "spearman")$estimate
    }
    mds_lev[counter] <- md
    seq_depth[counter] <- seq_d
    var_exp[counter] <- my_var_exp[md]
    PCAi_shannon_div_spear_cor[counter] <- cor(shan_dirv, myPCA[,md], method = "spearman")
    samples_left[counter] <- nrow(sd_filt_asv)
    taxa_left[counter] <- sum(colSums(sd_filt_asv) > 0)
    zero_count[counter] <- sum(sd_filt_asv == 0)
    counter <- counter + 1
  }
}
# dev.off()
result_df <- data.frame( mds_lev, seq_depth, var_exp, PCAi_shannon_div_spear_cor, samples_left, zero_count, taxa_left)
print("created resulting DF")

write.table(result_df, 
            file = file.path(output_dir, "tables", paste0(project, "_PCA_shandiv_filt_results.csv")),
            sep = ",")

result_df <- read.table(file = file.path(output_dir, "tables", paste0(project, "_PCA_shandiv_filt_results.csv")),
                        sep = ",")

result_df <- cbind(result_df, meta_cor)

melted_results <- reshape2::melt(result_df,
                       variable.name = "Metadata",
                       measure.vars = c(names(meta_cor)[metad_cols], "PCAi_shannon_div_spear_cor"))

pdf(file = file.path(output_dir, "graphics", "shan_div_cor_PCA12345_line.pdf"))
for (i in 1:max(result_df$mds_lev)){
  g <- ggplot2::ggplot(melted_results[melted_results$mds_lev == i, ],
                       aes(x=seq_depth, y=value^2, group = Metadata)) +
    ggplot2::geom_point(aes(color = factor(Metadata))) +
    ggplot2::geom_line(aes(color = factor(Metadata))) +
    ggplot2::ggtitle(paste0("PCA",i," Seq depth vs Rsqd of Shannon Diversity vs metadata")) +
    ggplot2::xlab("Min sequence depth per sample") +
    ggplot2::ylab("Spearman corr") +
    ggplot2::labs(fill = "Metadata") +
    theme(axis.text.x = element_text(angle = 90)) +
    ggplot2::theme_minimal()
  g
  print(g)
}
dev.off()


