# Author: Aaron Yerke (aaronyerke@gmail.com)
# Script for testing each transformation against different read depth to find the best one.
# This scipt is just for making the figures
# requires: 

rm(list = ls()) #clear workspace

##-functions--------------------------------------------------------##

print("Loading dependencies")
##-Load Depencencies------------------------------------------------##
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("optparse", quietly = TRUE)){ install.packages("optparse") }
library("optparse")
if (!requireNamespace("ggplot2", quietly = TRUE)) BiocManager::install("ggplot2")
library("ggplot2")
if (!requireNamespace("pracma", quietly = TRUE)) BiocManager::install("pracma")
library("pracma")

print("Processing cml options")

option_list <- list(
  optparse::make_option(c("-d", "--homedir"), type="character", 
                        default=file.path('~','git','balance_tree_exploration'), 
                        help="dataset dir path", metavar="character"),
  optparse::make_option(c("-p", "--project"), type="character", default="Zeller", 
                        help="project folder", metavar="character")
); 

opt_parser <- optparse::OptionParser(option_list=option_list);

opt <- optparse::parse_args(opt_parser);

print(opt)

print("Establishing directory layout and other constants.")
home_dir <- opt$homedir
project <- opt$project
output_dir <- file.path(home_dir, project, 'output')


result_df <- read.table(file = file.path(output_dir, "tables", paste0(project, "_PCA_seqdep_filt_results.csv")),
                        sep = ",")

pdf(file = file.path(output_dir, "graphics", "read_depth_artifact_PCA12345_line.pdf"))
for (i in 1:max(result_df$mds_lev)){
  pca_only <- result_df[result_df$mds_lev == i, ]
  g <- ggplot2::ggplot(pca_only,
                       aes(x=read_depth, y=spear_cor^2, group = ds_nam)) +
    ggplot2::geom_point(aes(color = factor(ds_nam))) +
    ggplot2::geom_line(aes(color = factor(ds_nam))) +
    ggplot2::annotate("text", x = head(pca_only$read_depth, n = length(pca_only$my_ds_names)),
                      y = head(c(pca_only$spear_cor^2), n = length(pca_only$my_ds_names)),
                      label = head(pca_only$ds_nam, n = length(pca_only$my_ds_names)),
                      hjust = -0.1) +
    ggplot2::ggtitle(paste0(project, ": Correlation PCA",  i, " vs read depth")) +
    ggplot2::xlab("Min read depth per sample") +
    ggplot2::ylab("R Squared") +
    ggplot2::labs(color = "Transformations") +
    ggplot2::scale_x_continuous(trans='log10') +
    ggplot2::theme(axis.text.x = element_text(angle = 90)) +
    ggplot2::theme_minimal()
  print(g)

  #plot remaining taxa
  g <- ggplot2::ggplot(pca_only,
                       aes(x=read_depth, y=taxa_left, group = ds_nam)) +
    ggplot2::geom_point(aes(color = factor(ds_nam))) +
    ggplot2::geom_line(aes(color = factor(ds_nam))) +
    ggplot2::ggtitle(paste0(project, ": taxa left",  i, " vs min reads per sample")) +
    ggplot2::xlab("Min read depth per sample") +
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
    ggplot2::ggtitle(paste0(project, ": samples left vs min reads per sample")) +
    ggplot2::xlab("Min read depth per sample") +
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
    ggplot2::ggtitle(paste0(project, ": zeros vs min reads per sample")) +
    ggplot2::xlab("Min read depth per sample") +
    ggplot2::ylab("Percentage of Zeros") +
    ggplot2::labs(fill = "Transformations") +
    ggplot2::theme(axis.text.x = element_text(angle = 90)) +
    ggplot2::theme_minimal()
  g
  print(g)
}
dev.off()
print("made line charts")

#for summary figure
mds_axis <- c()
Transformation <- c()
AUC <- c()
for (i in 1:max(result_df$mds_lev)){
  pca_only <- result_df[result_df$mds_lev == i, ]
  for (trans in unique(pca_only$ds_nam)) {
    my_r2 <- pca_only[pca_only$ds_nam == trans, "spear_cor"]^2
    my_rd <- pca_only[pca_only$ds_nam == trans, "read_depth"]
    my_auc <- pracma::trapz(my_rd, my_r2)
    
    #filling vectors
    AUC <- c(AUC, my_auc)
    Transformation <- c(Transformation, paste0(trans))
    mds_axis <- c(mds_axis, paste0("PCA",i))
  }
}

df1 <- base::data.frame(mds_axis, Transformation, AUC)

pdf(file = file.path(output_dir, "graphics", "read_depth_artifact_PCA12345_bar.pdf"),
    width = 9, height = 5)
ggplot2::ggplot(df1, aes(x = mds_axis, y=AUC, fill = Transformation)) +
  # geom_bar(position = "dodge", stat = "identity") +
  ggplot2::geom_col(position="dodge", width=0.85, color="black") +
  ggplot2::ggtitle(paste0(project, ": transformation AUC by PCA axis")) +
  ggplot2::xlab("PCA Axes") +
  ggplot2::ylab("Area Under Curve") +
  ggplot2::theme_minimal() +
  ggplot2::theme(axis.line = element_line(color="black"),
                 axis.ticks = element_line(color="black"),
                 panel.border = element_blank()) +
  ggplot2::theme(axis.line = element_line(color="black"),
                 axis.ticks = element_line(color="black"),
                 panel.border = element_blank()) +
  theme(text=element_text(size=15), #change font size of all text
        axis.text=element_text(size=15), #change font size of axis text
        axis.title=element_text(size=17), #change font size of axis titles
        plot.title=element_text(size=17), #change font size of plot title
        legend.text=element_text(size=15), #change font size of legend text
        legend.title=element_text(size=16)) #change font size of legend title 
dev.off()
print("made bar charts")







