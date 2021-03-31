# Author: Aaron Yerke (aaronyerke@gmail.com)
# Script for creating taxo anova pvalues.
# This was helepful: https://github.com/jsilve24/philr/blob/master/vignettes/philr-intro.Rmd#L142

rm(list = ls()) #clear workspace

##-functions--------------------------------------------------------##
makeTaxaTable <- function(no, tax, tax_lev_int){
  #no is otu table (dataframe)
  #tax is taxonomic table (dataframe)
  #tax_lev_int is the columns of the tax table usually 1-6
  
  
  otu_labels = vector(length = ncol(no), mode = "character")
  #assign otus to the asvs
  for (i in 1:ncol(no)){
    asv = names(no)[i]
    otu_labels[i] = tolower(as.character(tax[asv, tax_lev_int]))#this line sets the taxonomic level
  }
  uniq_otus = unique(otu_labels)
  
  otu_tab = data.frame(row.names = row.names(no))
  
  for(otu1 in 1:length(uniq_otus)){
    my_otu_lab = uniq_otus[otu1]
    my_otus = otu_labels == my_otu_lab
    my_otus = replace(my_otus, is.na(my_otus), F)
    #need to add vectors so that they are the same length as nrow(ps.philr)
    my_otus = as.data.frame(no[,my_otus])
    my_otus = unlist(rowSums(my_otus, na.rm=T))
    otu_tab[,otu1] =  my_otus
  }
  colnames(otu_tab) = uniq_otus
  return(otu_tab)
}
make_taxa_pie_chart <- function(tabl, title_info) {
  pie_areas <- data.frame("taxa" = colnames(tabl), 
                          "counts" = colSums(tabl),
                          row.names = 1:ncol(tabl))
  
  g <- ggplot(pie_areas, aes(x="", y=counts, fill=taxa)) +
    geom_bar(stat="identity", width=1) +
    coord_polar("y", start=0) +
    ggplot2::ggtitle(title_info) +
    theme(legend.position = "none") +
    theme_classic()
  return(g)
}


##-Load Depencencies------------------------------------------------##
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
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
asv_table <- data.frame(readRDS(file.path(output_dir, "r_objects", "ForwardReads_DADA2.rds")))

ref_ps <- readRDS(file.path(output_dir, "r_objects", "ref_tree_phyloseq_obj.rds"))
asv_tax <- data.frame(readRDS(file.path(output_dir, "r_objects", "ForwardReads_DADA2_taxonomy.rds")))

total_seqs <- rowSums(asv_table)
total_seqs <- data.frame(total_seqs, row.names = row.names(asv_table))

metadata <- read.table(file.path(home_dir, project, "patient_metadata.tsv"), 
                       sep="\t", 
                       header=TRUE, 
                       row.names = "Run", 
                       check.names = FALSE,
                       stringsAsFactors=TRUE)
metadata <- metadata[row.names(asv_table),]

##-Layout loop constants--------------------------------------------##
min_seq_depths <- c(100, 500, 1000, 5000, 10000, 20000, 40000)
mds_depth <- 1

seq_depth <- c()
taxa_int <- c()
taxa_name <- c()
pval <- c()
taxa_lev <- c()
mds_lev <- c()
var_exp <- c()
spear_cor <- c()
taxa_left <- c()
zero_count <- c()
fi_left_zeros <- c()
fo_left_zeros <- c()

counter <- 1
pdf(file = file.path(output_dir, "graphics", "seq_depth_taxa_pie.pdf"))
##-Loop and build result DF-----------------------------------------##
for(s in 1:length(min_seq_depths)){
  seq_d <- min_seq_depths[s]#new sequencing depth
  sd_filt_asv <- asv_table[total_seqs$total_seqs >= seq_d,]#dataset 1
  filter_out <- asv_table[total_seqs$total_seqs < seq_d,]

  for( taxa_col in 1:ncol(asv_tax)){
    my_table <- makeTaxaTable(sd_filt_asv, asv_tax, taxa_col)
    fo_table <- makeTaxaTable(filter_out, asv_tax, taxa_col)
    #make pie_charts
    print(make_taxa_pie_chart(my_table, paste(project, names(asv_tax)[taxa_col], "filter_in", seq_d)))
    print(make_taxa_pie_chart(fo_table, paste(project, names(asv_tax)[taxa_col], "filter_out", seq_d)))
    
    fi_zeros <- sum(my_table == 0)
    fi_nonzeros <- ncol(my_table)*nrow(my_table) - fi_zeros
    fo_zeros <- sum(fo_table == 0)
    fo_nonzeros <- ncol(fo_table)*nrow(fo_table) - fo_zeros
      
    fisher_data <- matrix(c(fi_zeros,fi_nonzeros,fo_zeros,fo_nonzeros), nrow = 2)
    ft <- fisher.test(fisher_data)
    
    my_zeros <- data.frame("zeros" = c(ncol(fo_table)*nrow(fo_table) - fo_zeros, fo_zeros),
                           "cats" = c("non-zeros", "filter_out_zeros"))
    
    g <- ggplot(my_zeros, aes(x="", y=zeros, fill=cats)) +
      geom_bar(stat="identity", width=1) +
      coord_polar("y", start=0) +
      ggplot2::ggtitle(paste(project, names(asv_tax)[taxa_col], "filter out zeros", seq_d, "Fisher: ", ft$p.value)) +
      theme(legend.position = "none")
    print(g)
    
    my_zeros <- data.frame("zeros" = c(fi_nonzeros , fi_zeros),
                           "cats" = c("non-zeros", "filter_in_zeros"))
    
    g <- ggplot(my_zeros, aes(x="", y=zeros, fill=cats)) +
      geom_bar(stat="identity", width=1) +
      coord_polar("y", start=0) +
      ggplot2::ggtitle(paste(project, names(asv_tax)[taxa_col], "filter in zeros", seq_d, "Fisher: ", ft$p.value)) +
      theme(legend.position = "none")
    print(g)
    
    ##-Create a PCA-----------------------------------------------------##
    my_prcmp <- prcomp(my_table, 
                       center = TRUE,
                       rank = mds_depth)#,
    # scale = TRUE)
    ##-Extract PCA matrix and convert to dataframe----------------------##
    myPCA <- data.frame(my_prcmp$x)
    my_var_exp <- my_prcmp$sdev^2/sum(my_prcmp$sdev^2)
    
    for (md in 1:mds_depth){
      mds_lev[counter] <- md
      seq_depth[counter] <- seq_d
      var_exp[counter] <- my_var_exp[md]
      spear_cor[counter] <- cor(total_seqs[total_seqs >= seq_d], myPCA[,md], method = "spearman")
      taxa_left[counter] <- ncol(my_table)
      taxa_lev[counter] <- taxa_col
      taxa_name[counter] <- names(asv_tax)[taxa_col]
      fi_left_zeros[counter] <- fi_zeros
      fo_left_zeros[counter] <- fo_zeros
      counter <- counter + 1
    }
  }
}  
dev.off()

result_df <- data.frame(taxa_lev, taxa_name, mds_lev, seq_depth, var_exp, spear_cor, fi_left_zeros, fo_left_zeros)

for (i in 1:max(result_df$taxa_lev)){
  pca_only <- result_df[taxa_lev == i, ]
  g <- ggplot2::ggplot(pca_only, 
                       aes(x=seq_depth, y=spear_cor^2)) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::ggtitle(paste0(project, ": Taxa Lev: ",  result_df$taxa_name[i], " PCA1 vs seq depth")) +
    ggplot2::xlab("Min sequence depth per sample") +
    ggplot2::ylab("R ^ 2") + 
    ggplot2::labs(fill = "Transformations") +
    theme(axis.text.x = element_text(angle = 90)) +
    ggplot2::theme_minimal()
  print(g)
  
  # g <- ggplot2::ggplot(pca_only, 
  #                      aes(x=seq_depth, y=fi_left_zeros, fill = "blue")) +
  #   ggplot2::geom_point(aes(x=seq_depth, y=fo_left_zeros)) +
  #   ggplot2::geom_line(aes(x=seq_depth, y=fo_left_zeros)) +
  #   ggplot2::ggtitle(paste0(project, ": Taxa Level:",  i, "Zeros vs seq depth")) +
  #   ggplot2::xlab("Min sequence depth per sample") +
  #   ggplot2::ylab("italic(R) ^ 2") + 
  #   ggplot2::labs(fill = "Transformations") +
  #   theme(axis.text.x = element_text(angle = 90)) +
  #   ggplot2::theme_minimal()
  # print(g)
  
}
