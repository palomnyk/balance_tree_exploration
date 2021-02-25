# Author: Aaron Yerke (aaronyerke@gmail.com)
# Script for testing each transformation against different sequence depth to find the best one
# This was helepful: https://github.com/jsilve24/philr/blob/master/vignettes/philr-intro.Rmd#L142

rm(list = ls()) #clear workspace

##-functions--------------------------------------------------------##
pca_total_seqs_plot <- function(df, ts, main_tit, mds){
  #Function for updating 
  ##-Create PCA-------------------------------------------------------##
  my_prcmp <- prcomp(df, 
                     center = TRUE)#,
  # scale = TRUE)
  ##-Extract PCA matrix and convert to dataframe----------------------##
  myPCA <- data.frame(my_prcmp$x)
  my_spear <- cor.test(log10(ts), myPCA[,mds],
                       method = "kendall")
  plot(log10(ts), myPCA[,2],
       main = paste0( main_tit, " PCA", mds,"\nr_sq: ", my_spear$estimate^2),
       sub = paste0("Kend: ", my_spear$estimate, " Kend^2: ", ),
       )
  return(my_spear$estimate)
}
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
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ALDEx2")
library("compositions")
library("phyloseq")
library("vegan")
library("DESeq2")
library("philr")
library("ape")
library("ALDEx2")
library("ggplot2")

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
# asv_table <- read.table(file.path(output_dir, "tables", "ForwardReads_DADA2.txt"),
#                         sep = "\t",
#                         header = TRUE)

ref_ps <- readRDS(file.path(output_dir, "r_objects", "ref_tree_phyloseq_obj.rds"))
asv_tax <- readRDS(file.path(output_dir, "r_objects", "ForwardReads_DADA2_taxonomy.rds"))


metadata <- read.table(file.path(home_dir, project, "patient_metadata.tsv"), 
                       sep="\t", 
                       header=TRUE, 
                       row.names = "Run", 
                       check.names = FALSE,
                       stringsAsFactors=TRUE)

#Plot shannon diversity against log10(total_seqs)
#Plot other normalization methods: otu log 100 and alr and clr 
# min_seq_depths <- c(0,100,500,1000,2000,5000,10000,50000)
my_ds_names <- c( "raw seqs", "clr(raw seqs)", "lognorm raw seqs", "philr ref", "DESeq2")
min_seq_depths <- c(0, 500, 1000, 5000, 10000, 50000)
mds_depth <- 5

total_seqs <- rowSums(asv_table)
total_seqs <- data.frame(total_seqs, row.names = row.names(asv_table))

kend <- vector(mode = "numeric", length = length(my_ds_names) * length(min_seq_depths) * mds_depth)
perma_r2 <- vector(mode = "numeric", length = length(my_ds_names) * length(min_seq_depths) * mds_depth)
ds_num <- vector(mode = "integer", length = length(my_ds_names) * length(min_seq_depths) * mds_depth)
ds_nam <- vector(mode = "character", length = length(my_ds_names) * length(min_seq_depths) * mds_depth)
mds_lev <- vector(mode = "integer", length = length(my_ds_names) * length(min_seq_depths) * mds_depth)
seq_depth <- vector(mode = "integer", length = length(my_ds_names) * length(min_seq_depths) * mds_depth)
var_exp <- vector(mode = "integer", length = length(my_ds_names) * length(min_seq_depths) * mds_depth)
spear_cor <- vector(mode = "integer", length = length(my_ds_names) * length(min_seq_depths) * mds_depth)

counter <- 1
pdf(file = file.path(output_dir, "graphics", "explore_seq_depth_artifact_var_seq_depth_PCA12345_4ds.pdf"))
for(s in 1:length(min_seq_depths)){
  seq_d <- min_seq_depths[s]#new sequencing depth
  sd_filt_asv <- asv_table[total_seqs$total_seqs >= seq_d,]#dataset 1
  print(paste("sd_filtered dim:", paste(dim(sd_filt_asv))))
  safe_rns <- intersect(row.names(ref_ps@otu_table), row.names(sd_filt_asv))
  my_clr <- compositions::clr(sd_filt_asv)#dataset 3
  new_tree <- phyloseq::prune_taxa(colnames(sd_filt_asv), ref_ps@phy_tree)#update tree for new phyloseq obj
  new_ref_ps <- prune_samples(safe_rns, ref_ps)
  new_ref_ps <- munge_ref_ps(new_ref_ps)
  new_DESeq2 <- phyloseq::phyloseq_to_deseq2(new_ref_ps, design= ~ 1)#dataset 5
  new_DESeq2 <- estimateSizeFactors(new_DESeq2)
  new_DESeq2 <- t(counts(new_DESeq2, normalized=T))
  print(paste("new DSeq:", paste(dim(new_DESeq2))))
  
  print(paste("new dim ref ps:", dim(data.frame(new_ref_ps@otu_table))))
  print("made new ps")
  ln_ps <- lognorm(new_ref_ps@otu_table)#dataset 7

  ref_philr <- philr::philr(new_ref_ps@otu_table, new_ref_ps@phy_tree,
                            part.weights='enorm.x.gm.counts',
                            ilr.weights='blw.sqrt')#dataset 4
  print(paste("made new philr", dim(as.data.frame(ref_philr))))
  # print(any(is.na(data.frame(ref_philr))))

  ln_asv <- lognorm(sd_filt_asv)#dataset 6
  my_datasets <- list(sd_filt_asv, my_clr,  ln_asv, ref_philr, new_DESeq2)
  
  for( ds in 1:length(my_datasets)){
    print(my_ds_names[ds])
    my_table <- as.data.frame(my_datasets[ds])
    print(dim(my_table))
    my_prcmp <- prcomp(my_table, 
                       center = TRUE)#,
    # scale = TRUE)
    ##-Extract PCA matrix and convert to dataframe----------------------##
    myPCA <- data.frame(my_prcmp$x)
    my_var_exp <- my_prcmp$sdev^2/sum(my_prcmp$sdev^2)
    
    # if ( s == 1){
    #   my_cor <- cor(log10(total_seqs[total_seqs >= seq_d]), myPCA[,1], method = "spearman")
    #   plot(myPCA[,1], log10(total_seqs[total_seqs >= seq_d]),
    #        main = paste(my_ds_names[ds], " vs PCA1"),
    #        sub = paste0("Spearman: ", my_cor$estimate))
    # }
    for (md in 1:mds_depth){
      # kend[counter] <- cor.test(log10(total_seqs[total_seqs > seq_d]), myPCA[,md], method = "kendall")$estimate
      # perma_r2[counter] <- adonis2(log10(total_seqs[total_seqs > seq_d]) ~ myPCA[,md])$R2[1]
      ds_num[counter] <- ds
      ds_nam[counter] <- my_ds_names[ds]
      mds_lev[counter] <- md
      seq_depth[counter] <- seq_d
      var_exp[counter] <- my_var_exp[md]
      spear_cor[counter] <- cor(log10(total_seqs[total_seqs >= seq_d]), myPCA[,md], method = "spearman")
      counter <- counter + 1
    }
  }
}

result_df <- data.frame(ds_num, ds_nam, perma_r2, mds_lev, seq_depth, var_exp, spear_cor)
g <- ggplot(result_df, aes(x = as.factor(seq_depth), y = perma_r2)) +
  geom_text(aes(label=mds_lev, color = factor(ds_nam))) +
  ggtitle("Permanova R^2 filtered seq depths")
  # scale_x_discrete(limits=c("2","1","0.5"))
g
dev.off()

for (i in 1:max(result_df$mds_lev)){
  pca_only <- result_df[mds_lev == i, ]
  g <- ggplot2::ggplot(pca_only, 
                       aes(x=as.factor(seq_depth), y=spear_cor^2, fill=ds_nam)) +
    geom_bar(width = 0.8, position=position_dodge(width = 1), stat="identity",) +
    ggtitle(paste0(project, ": PCA",  i, " vs total sequences per sample")) +
    xlab("Min sequence depth per sample") +
    ylab("Spearman Rsq") +
    labs(fill = "Transformations")
    # geom_text(aes(x=as.factor(seq_depth), y=spear_cor^2, 
    #               label = round(var_exp*100, digits = 0)),
    #           position = position_dodge(width = 1))
  # ylab(y_lab)
  ggsave(
    filename = file.path(output_dir, "graphics", paste0(project, "_PCA",  i, "_5ds_var_exp.pdf")),
    plot = g,
    device = "pdf"
  )
}

pdf(file = file.path(output_dir, "graphics", "seq_depth_artifact_PCA12345_lines.pdf"))
for (i in 1:max(result_df$mds_lev)){
  pca_only <- result_df[mds_lev == i, ]
  g <- ggplot2::ggplot(pca_only, 
                       aes(x=seq_depth, y=spear_cor^2, group = ds_nam)) +
    geom_point(aes(color = factor(ds_nam))) +
    geom_line(aes(color = factor(ds_nam))) +
    ggtitle(paste0(project, ": PCA",  i, " vs total sequences per sample")) +
    xlab("Min sequence depth per sample") +
    ylab("Spearman Rsq") + 
    labs(fill = "Transformations") +
    theme_minimal()
  g
  print(g)
  # ggsave(
  #   filename = file.path(output_dir, "graphics", paste0(project, "_PCA",  i, "_5ds_var_exp_points.pdf")),
  #   plot = g,
  #   device = "pdf"
  # )
}
dev.off()



g <- aldex.clr(asv_table, mc.samples=150, denom="all", verbose=F)

# conda create --name R35 R=3.5 

