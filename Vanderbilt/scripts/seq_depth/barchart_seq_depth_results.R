# Author: Aaron Yerke (aaronyerke@gmail.com)
# Script for testing each transformation against sequence depth
# This was helepful: https://github.com/jsilve24/philr/blob/master/vignettes/philr-intro.Rmd#L142

rm(list = ls()) #clear workspace

##-functions--------------------------------------------------------##
pca_total_seqs_plot <- function(df, ts, main_tit, mds){
  ##-Create PCA-------------------------------------------------------##
  my_prcmp <- prcomp(df, 
                     center = TRUE,
                     scale = TRUE)
  ##-Extract PCA matrix and convert to dataframe----------------------##
  myPCA <- data.frame(my_prcmp$x)
  my_spear <- cor.test(log10(ts), myPCA[,mds],
                       method = "kendall")
  plot(log10(ts), myPCA[,2],
       main = paste0( main_tit, " PCA", mds,"\nr_sq: ", my_spear$estimate^2),
       sub = paste0("Kendall Cor: ", my_spear$estimate))
  return(my_spear$estimate)
}

##-Load Depencencies------------------------------------------------##
if (!requireNamespace("randomForest", quietly = TRUE)){
  install.packages("randomForest")
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer")
}
library(RColorBrewer)
#set color palette
palette( brewer.pal(7,"Accent") )
library(randomForest)
library(compositions)
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

##-Import data------------------------------------------------------##
my_corrs <- read.table(file.path(output_dir, "tables", "PCA123_vs_total_seqs_min_1000.csv"),
                       sep = ",",
                       header = 1,
                       row.names = 1)

my_kend <- read.table(file.path(output_dir, "tables", "PCA123_vs_total_seqs_min_1000_kend.csv"),
                       sep = ",",
                       header = 1,
                       row.names = 1)

my_rsq <- my_kend^2
  
plot_df_boxplot <- function(df, main_lab, y_lab){
  my_df <- df
  my_df$group <- row.names(my_df)
  my_df.m <- reshape2::melt(my_df, id.vars = "group")
  g <- ggplot2::ggplot(my_df.m, aes(group, value)) + geom_boxplot() + 
    ggtitle(main_lab) +
    geom_hline(yintercept = 0) +
    theme(axis.text.x = element_text(angle = 45)) +
    xlab("Transformations") +
    ylab(y_lab)
  return(g)
}

plot_df_boxplot(my_kend, 
                paste(project, "PCA12&3 Kendal Cor with Seq Depth"), 
                "Kendal Corr")

plot_df_boxplot(my_rsq, 
                paste(project, "PCA12&3 RSQ with Seq Depth"), 
                "Kendal Rsq")

plot_df_barplot <- function(df, main_lab, y_lab){
  my_df <- df
  my_df$group <- row.names(my_df)
  my_df.m <- reshape2::melt(my_df, id.vars = "group")
  # View(my_df.m)
  g <- ggplot2::ggplot(my_df.m, aes(fill=variable, x=group, y=value)) + 
    geom_bar(position=position_dodge(width = 0.8), stat="identity",) + 
    facet_grid() +
    geom_point(position=position_dodge(width = 0.8)) +
    ggtitle(main_lab) +
    geom_hline(yintercept = 0) +
    theme(axis.text.x = element_text(angle = 0)) +
    xlab("Transformations") +
    theme_classic() +
    ylab(y_lab)
  return(g)
}

plot_df_barplot(my_kend,
                paste(project, "PCA12&3 Kendal Cor with Seq Depth"),
                "Kendal Corr")

plot_df_barplot(my_rsq, 
                paste(project, "PCA12&3 RSQ with Seq Depth"), 
                "Kendal Rsq")



