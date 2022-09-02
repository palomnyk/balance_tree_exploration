# Author: Aaron Yerke (aaronyerke@gmail.com)
# Script for testing each transformation against different rarefication level
# to find the best one.

rm(list = ls()) #clear workspace

##-functions--------------------------------------------------------##
munge_ref_ps <- function(ps){
  #ps must be a phyloseq object
  #requires phyloseq and ape packages to be loaded in the env
  # ps <- filter_taxa(ps, function(x) sum(x > 3) > (0.2*length(x)), TRUE)
  # ps <- filter_taxa(ps, function(x) sd(x)/mean(x) > 3.0, TRUE)
  ps <- phyloseq::transform_sample_counts(ps, function(x) x+1)
  phy_tree(ps) <- ape::makeNodeLabel(phy_tree(ps), method="number", prefix='n')
  return(ps)
}

print("Loading dependencies")
##-Load Depencencies------------------------------------------------##
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("ALDEx2", quietly = TRUE)) BiocManager::install("ALDEx2")
library("ALDEx2")
if (!requireNamespace("optparse", quietly = TRUE)){ install.packages("optparse") }
library("optparse")
if (!requireNamespace("rgr", quietly = TRUE)) install.packages("rgr")
library("rgr")
if (!requireNamespace("phyloseq", quietly = TRUE)) BiocManager::install("phyloseq")
library("phyloseq")
if (!requireNamespace("vegan", quietly = TRUE)) BiocManager::install("vegan")
library("vegan")
if (!requireNamespace("DESeq2", quietly = TRUE)) BiocManager::install("DESeq2")
library("DESeq2")
if (!requireNamespace("philr", quietly = TRUE)) BiocManager::install("philr")
library("philr")
if (!requireNamespace("ape", quietly = TRUE)) BiocManager::install("ape")
library("ape")
if (!requireNamespace("ggplot2", quietly = TRUE)) BiocManager::install("ggplot2")
library("ggplot2")
if (!requireNamespace("vegan", quietly = TRUE)) BiocManager::install("vegan")
library("vegan")


print("Processing cml options")

option_list <- list(
  optparse::make_option(c("-d", "--homedir"), type="character", 
                        default=file.path('~','git','balance_tree_exploration'), 
                        help="dataset dir path", metavar="character"),
  optparse::make_option(c("-p", "--project"), type="character", default=NULL, 
                        help="project folder", metavar="character")
); 

opt_parser <- optparse::OptionParser(option_list=option_list);

opt <- optparse::parse_args(opt_parser);

print(opt)

print("Establishing directory layout and other constants.")
home_dir <- opt$homedir
project <- opt$project
output_dir <- file.path(home_dir, project, 'output')
setwd(file.path(home_dir))

##-Functions--------------------------------------------------------##
source(file.path(home_dir, "lib", "statistical_functions.R"))
source(file.path(home_dir, "lib", "table_manipulations.R"))

##-Import tables and data preprocessing-----------------------------##
asv_table <- readRDS(file.path(output_dir, "r_objects", "ForwardReads_DADA2.rds"))

ref_ps <- readRDS(file.path(output_dir, "r_objects", "ref_tree_phyloseq_obj.rds"))

my_ds_names <- c( "rarefied counts table", "clr", "alr", "lognorm", "PhILR Silva Tree", "DESeq2", "ALDEx2.clr")

total_seqs <- rowSums(asv_table)#read depth
total_seqs <- data.frame(total_seqs, row.names = row.names(asv_table))

median_depth <- stats::median(total_seqs$total_seqs)

mean_depth <- ceiling(base::mean(total_seqs$total_seqs))

print(paste("mean:", mean_depth, "median:", median_depth))

raref_levels <- c(500, seq(1000,mean_depth,ceiling(max(total_seqs$total_seqs)/7)))
print("min read depths")
print(raref_levels)

mds_depth <- 5

kend <- vector(mode = "numeric", length = length(my_ds_names) * length(raref_levels) * mds_depth)
perma_r2 <- vector(mode = "numeric", length = length(my_ds_names) * length(raref_levels) * mds_depth)
ds_num <- vector(mode = "integer", length = length(my_ds_names) * length(raref_levels) * mds_depth)
ds_nam <- vector(mode = "character", length = length(my_ds_names) * length(raref_levels) * mds_depth)
mds_lev <- vector(mode = "integer", length = length(my_ds_names) * length(raref_levels) * mds_depth)
read_depth <- vector(mode = "integer", length = length(my_ds_names) * length(raref_levels) * mds_depth)
var_exp <- vector(mode = "integer", length = length(my_ds_names) * length(raref_levels) * mds_depth)
spear_cor <- vector(mode = "integer", length = length(my_ds_names) * length(raref_levels) * mds_depth)
samples_left <- vector(mode = "integer", length = length(my_ds_names) * length(raref_levels) * mds_depth)
taxa_left <- vector(mode = "integer", length = length(my_ds_names) * length(raref_levels) * mds_depth)
zero_count <- vector(mode = "integer", length = length(my_ds_names) * length(raref_levels) * mds_depth)

pdf(file = file.path(output_dir, "graphics", "raref_artifact_PCA12345_scatter.pdf"))
counter <- 1
for(s in 1:length(raref_levels)){
  seq_d <- raref_levels[s]#new read depth
  rd_filt_asv <- data.frame(vegan::rrarefy(asv_table, seq_d))#dataset 1 (read depth filtered asv)
  rd_filt_rowSums <- base::rowSums(rd_filt_asv)
  print(paste("rd_filt_asv dim:", paste(dim(rd_filt_asv))))
  if (nrow(rd_filt_asv > 2)){
    safe_rns <- intersect(row.names(ref_ps@otu_table), row.names(rd_filt_asv)) #rows for this iterate
    ts <- rowSums(rd_filt_asv[safe_rns,]) #sample read depths
    my_clr <- as.data.frame(rgr::clr(as.matrix(rd_filt_asv + 1)))#dataset 2
    ##-Find best col for alr denominator--------------------------------##
    my_zeros <- apply(rd_filt_asv, 2, function(x) {
      return(sum(x == 0))
    })
    alr_col <- which(my_zeros == min(my_zeros))[1]
    ##------------------------------------------------------------------##
    my_alr <- as.data.frame(rgr::alr(as.matrix(rd_filt_asv + 1), j = as.numeric(alr_col)))#dataset 3
    ##-Update philr tree------------------------------------------------##
    new_tree <- phyloseq::prune_taxa(colnames(rd_filt_asv), ref_ps@phy_tree)#update tree for new phyloseq obj
    new_ref_ps <- phyloseq::prune_samples(safe_rns, ref_ps) #remove non-safe rows from ps
    new_ref_ps <- munge_ref_ps(new_ref_ps)
    print(paste("new dim ref ps:", dim(data.frame(new_ref_ps@otu_table))))
    ##------------------------------------------------------------------##
    #create DESeq2 dtaset from new ref ps
    new_DESeq2 <- phyloseq::phyloseq_to_deseq2(new_ref_ps, design= ~ 1)#dataset 5
    new_DESeq2 <- DESeq2::estimateSizeFactors(new_DESeq2)
    new_DESeq2 <- t(DESeq2::counts(new_DESeq2, normalized=T))
    
    print(paste("new DSeq:", paste(dim(new_DESeq2))))
    
    ref_philr <- philr::philr(new_ref_ps@otu_table, new_ref_ps@phy_tree,
                              part.weights='enorm.x.gm.counts',
                              ilr.weights='blw.sqrt')#dataset 4
    print(paste("made new philr", dim(as.data.frame(ref_philr))))
    
    ln_asv <- lognorm(rd_filt_asv)#dataset 6
    
    ald <- ALDEx2::aldex.clr(rd_filt_asv, mc.samples=12, denom="all", verbose=F)
    ald <- data.frame(ald@analysisData)#dataset 7
    print(paste("size of ald:", object.size(ald)))
    print(paste("ald dim:", paste(dim(ald))))
    
    my_datasets <- list(rd_filt_asv, my_clr, my_alr, ln_asv, ref_philr, new_DESeq2, ald)
    
    print(paste("finished seq depth filter:", s))
    
    for( ds in 1:length(my_datasets)){
      print(my_ds_names[ds])
      my_table <- as.data.frame(my_datasets[ds])
      zeros <- sum(my_table == 0)
      print(dim(my_table))
      print(paste("Creating PCA for", my_ds_names[ds]))
      my_prcmp <- stats::prcomp(my_table, 
                                center = TRUE,
                                rank = mds_depth)#, (scale causes it to fail)
      # scale = TRUE)
      print(paste("Extracting PCA matrix", my_ds_names[ds]))
      ##-Extract PCA matrix and convert to dataframe----------------------##
      myPCA <- data.frame(my_prcmp$x)
      print(paste("Dim of myPCA:", nrow(myPCA), ncol(myPCA)))
      print("Creating proportion of variance explained")
      my_var_exp <- my_prcmp$sdev^2/sum(my_prcmp$sdev^2)
      for (md in 1:mds_depth){
        print(paste("Updating result holding vectors for MDS", md))
        # kend[counter] <- cor.test(log10(total_seqs[total_seqs > seq_d]), myPCA[,md], method = "kendall")$estimate
        ds_num[counter] <- ds
        ds_nam[counter] <- my_ds_names[ds]
        mds_lev[counter] <- md
        read_depth[counter] <- seq_d
        var_exp[counter] <- my_var_exp[md]
        spear_cor[counter] <- cor(rd_filt_rowSums, myPCA[,md], method = "spearman")
        samples_left[counter] <- nrow(my_table)
        taxa_left[counter] <- ncol(my_table)
        zero_count[counter] <- zeros
        counter <- counter + 1
      }
      print(paste("finished ds:", my_ds_names[ds], "read depth:", seq_d))
    }
  }else{
    print(paste(seq_d, "is too low of a read depth filter, not enough samples left, breaking loop"))
    break
  }
}
dev.off()
result_df <- data.frame(ds_num, ds_nam, perma_r2, mds_lev, read_depth, var_exp, spear_cor, samples_left, zero_count, taxa_left)
print("created resulting DF")

write.table(result_df, 
            file = file.path(output_dir, "tables", paste0(project, "_PCA_raref_results.csv")),
            sep = ",")
result_df <- read.table(file = file.path(output_dir, "tables", paste0(project, "_PCA_raref_results.csv")),
                        sep = ",")

pdf(file = file.path(output_dir, "graphics", "raref_rd_artifact_PCA12345_bar.pdf"))
for (i in 1:max(result_df$mds_lev)){
  pca_only <- result_df[mds_lev == i, ]
  g <- ggplot2::ggplot(pca_only,
                       aes(x=as.factor(read_depth), y=spear_cor^2, fill=ds_nam)) +
    ggplot2::geom_bar(width = 0.8, position=position_dodge(width = 1), stat="identity",) +
    ggplot2::ggtitle(paste0(project, ": PCA",  i, " vs read depth after rarefication")) +
    ggplot2::xlab("Rarefication level") +
    ggplot2::ylab(paste0("Read depth to PCA", i, "Spearman Rsq"))
  print(g)
}
dev.off()
print("made bar charts")

pdf(file = file.path(output_dir, "graphics", "raref_rd_artifact_PCA12345_line.pdf"))
for (i in 1:max(result_df$mds_lev)){
  pca_only <- result_df[result_df$mds_lev == i, ]
  g <- ggplot2::ggplot(pca_only, 
                       aes(x=read_depth, y=spear_cor^2, group = ds_nam)) +
    ggplot2::geom_point(aes(color = factor(ds_nam))) +
    ggplot2::geom_line(aes(color = factor(ds_nam))) +
    ggplot2::annotate("text", x = head(pca_only$read_depth, n = length(my_ds_names)), 
                      y = head(c(pca_only$spear_cor^2), n = length(my_ds_names)), 
                      label = head(pca_only$ds_nam, n = length(my_ds_names)),
                      hjust = -0.1) +
    ggplot2::ggtitle(paste0(project, ": PCA",  i, " vs read depth after rarefication")) +
    ggplot2::xlab("Rarefication level") +
    ggplot2::ylab(paste0("Read depth to PCA", i, "Spearman Rsq")) +
    ggplot2::xlab("Rarefication level") +
    ggplot2::ylab("R Squared") + 
    ggplot2::labs(fill = "Transformations") +
    ggplot2::scale_x_continuous(trans='log10') +
    ggplot2::theme(axis.text.x = element_text(angle = 90)) +
    ggplot2::theme_minimal()
  print(g)
  
  #plot remaining taxa
  g <- ggplot2::ggplot(pca_only, 
                       aes(x=read_depth, y=taxa_left, group = ds_nam)) +
    ggplot2::geom_point(aes(color = factor(ds_nam))) +
    ggplot2::geom_line(aes(color = factor(ds_nam))) +
    ggplot2::ggtitle(paste0(project, ": PCA",  i, " vs read depth after rarefication")) +
    ggplot2::xlab("Rarefication level") +
    ggplot2::ylab("Taxa") + 
    ggplot2::labs(fill = "Transformations") +
    ggplot2::theme(axis.text.x = element_text(angle = 90)) +
    ggplot2::theme_minimal()
  g
  print(g)
  
  #plot remaining samples
  g <- ggplot2::ggplot(pca_only, 
                       aes(x=read_depth, y=samples_left, group = ds_nam)) +
    ggplot2::geom_point(aes(color = factor(ds_nam))) +
    ggplot2::geom_line(aes(color = factor(ds_nam))) +
    ggplot2::ggtitle(paste0(project, ": PCA",  i, " vs read depth after rarefication")) +
    ggplot2::xlab("Rarefication level") +
    ggplot2::ylab("Samples") + 
    ggplot2::labs(fill = "Transformations") +
    ggplot2::theme(axis.text.x = element_text(angle = 90)) +
    ggplot2::theme_minimal()
  g
  print(g)
  
  #plot percentage of zeros
  g <- ggplot2::ggplot(pca_only,
                       aes(x=read_depth, y=zero_count/samples_left*taxa_left, group = ds_nam)) +
    ggplot2::geom_point(aes(color = factor(ds_nam))) +
    ggplot2::geom_line(aes(color = factor(ds_nam))) +
    ggplot2::ggtitle(paste0(project, ": PCA",  i, " vs read depth after rarefication")) +
    ggplot2::xlab("Rarefication level") +
    ggplot2::ylab("Percentage of Zeros") +
    ggplot2::labs(fill = "Transformations") +
    ggplot2::theme(axis.text.x = element_text(angle = 90)) +
    ggplot2::theme_minimal()
  g
  print(g)
}
dev.off()
print("made line chart")

