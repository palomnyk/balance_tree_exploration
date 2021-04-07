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
  pie_areas <- data.frame("label" = colnames(tabl), 
                          "counts" = colSums(tabl),
                          row.names = 1:ncol(tabl))
  g <- ggplot(pie_areas, aes(x="", y=counts, fill=factor(label))) +
    geom_bar(stat="identity", width=1) +
    coord_polar("y", start=0) +
    ggplot2::ggtitle(title_info) +
    theme(legend.position = "none",
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),) +
    guides(fill=FALSE) +
  return(g)
}
total_cells <- function(df){
  return(dim(df)[1]*dim(df)[2])
}


##-Load Depencencies------------------------------------------------##
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("gridExtra", quietly = TRUE)) install.packages("gridExtra")
library("ggplot2")
library("gridExtra")

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

#columns for first result df
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
taxa_pval <- c()
zero_pval <- c()

counter <- 1
##-Loop and build result DF-----------------------------------------##
for( taxa_col in 1:ncol(asv_tax)){
  tname <- names(asv_tax)[taxa_col]
  pdf(file = file.path(output_dir, "graphics", paste0(tname,"_seq_depth_taxa_pie.pdf")))
  taxa_plots_fi <- list()
  taxa_plots_fo <- list()
  zero_plots_fi <- list()
  zero_plots_fo <- list()
  for(s in 1:length(min_seq_depths)){
    seq_d <- min_seq_depths[s]#new sequencing depth
    sd_filt_asv <- asv_table[total_seqs$total_seqs >= seq_d,]#dataset 1
    filter_out <- asv_table[total_seqs$total_seqs < seq_d,]

    my_table <- makeTaxaTable(sd_filt_asv, asv_tax, taxa_col)
    fo_table <- makeTaxaTable(filter_out, asv_tax, taxa_col)
    #make taxa pie_charts
    taxa_char <- substr(names(asv_tax)[taxa_col],1,1)
    fi_plot <- make_taxa_pie_chart(my_table, paste( taxa_char, "fi t", seq_d))
    fo_plot <- make_taxa_pie_chart(fo_table, paste(taxa_char, "fo t", seq_d))
    taxa_plots_fi[[s]] <- fi_plot
    taxa_plots_fo[[s]] <- fo_plot
    
    #make taxa pie_charts
    fi_zeros <- sum(my_table == 0)/total_cells(my_table)
    fi_nonzeros <- (ncol(my_table)*nrow(my_table) - fi_zeros)/total_cells(my_table)
    fo_zeros <- sum(fo_table == 0)/total_cells(fo_table)
    fo_nonzeros <- (ncol(fo_table)*nrow(fo_table) - fo_zeros)/total_cells(fo_table)

    my_zeros <- data.frame("non-zeros" = c(fo_nonzeros),
                           "zeros" = c(fo_zeros))
    zero_plots_fo[[s]] <- make_taxa_pie_chart(my_zeros,
                                              paste(taxa_char, "fo z", seq_d))
    my_zeros <- data.frame("non-zeros" = c(fi_nonzeros),
                           "zeros" = c(fi_zeros))
    zero_plots_fi[[s]] <- make_taxa_pie_chart(my_zeros,
                                              paste(taxa_char, "fi z", seq_d))

    #Do stats
    fisher_data <- matrix(c(fi_zeros,fi_nonzeros,fo_zeros,fo_nonzeros), nrow = 2)
    zero_ft <- chisq.test(fisher_data, simulate.p.value = T)
    
    fi_taxa <- data.frame("fi" = colSums(my_table),
                         row.names = replace( colnames(my_table), is.na(colnames(my_table)), "UNKNOWN"))
    fo_taxa <- data.frame("fo" = colSums(fo_table),
                          row.names = replace( colnames(fo_table), is.na(colnames(fo_table)), "UNKNOWN"))
    
    comb_tax <- merge(fi_taxa, by = "row.names", fo_taxa, all = T)
    comb_tax <- comb_tax[-c(1)]

    comb_tax[is.na(comb_tax)] <- 0
    print(apply(comb_tax,2,max))
    print(comb_tax[1,])
    col_max <- c(apply(comb_tax,2,max))
    comb_tax <- sweep(comb_tax, 2, apply(comb_tax,2,max), FUN = "/")
    print(comb_tax[1,])
    taxa_chi <- chisq.test(comb_tax$fo, comb_tax$fi, simulate.p.value = T )
    
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
      taxa_name[counter] <- tname
      fi_left_zeros[counter] <- fi_zeros
      fo_left_zeros[counter] <- fo_zeros
      taxa_pval[counter] <- taxa_chi$p.value
      zero_pval[counter] <- zero_ft$p.value
      counter <- counter + 1
    }
  }
  n <- length(taxa_plots_fi)
  # nCol <- floor(sqrt(n))
  do.call("grid.arrange", c(append(taxa_plots_fi, taxa_plots_fo), 
                            ncol=n))
  do.call("grid.arrange", c(append(zero_plots_fi, zero_plots_fo), 
                            ncol=n))
  dev.off()
}  

result_df <- data.frame(taxa_lev, taxa_name, mds_lev, seq_depth, 
                        var_exp, spear_cor, fi_left_zeros, fo_left_zeros,
                        taxa_pval, zero_pval)

# for (i in 2:max(result_df$taxa_lev)){
#   pca_only <- result_df[taxa_lev == i, ]
#   g <- ggplot2::ggplot(pca_only,
#                        aes(x=seq_depth)) +
#     ggplot2::geom_point(aes(y=fo_left_zeros)) +
#     ggplot2::geom_line(aes(y=fo_left_zeros)) +
#     ggplot2::geom_point(aes(y=fi_left_zeros)) +
#     ggplot2::geom_line(aes(y=fi_left_zeros)) +
#     ggplot2::ggtitle(paste0(project, ": Taxa Level:",  pca_only$taxa_name[1], " Zeros vs seq depth")) +
#     ggplot2::xlab("Min sequence depth per sample") +
#     ggplot2::ylab("Zeros") +
#     theme(axis.text.x = element_text(angle = 90)) +
#     ggplot2::theme_minimal()
#   print(g)
# }

g <- ggplot2::ggplot(result_df,
                     aes(x=seq_depth, y=taxa_pval, group = factor(taxa_name))) +
  ggplot2::geom_point(aes(color = factor(taxa_name))) +
  ggplot2::geom_line(aes(color = factor(taxa_name))) +
  ggplot2::ggtitle(paste0(project, ': Normalized filter-in vs filter-out taxa chi sq')) +
  ggplot2::xlab("Min sequence depth per sample") +
  ggplot2::ylab("Taxa vs taxa chi sqr pvalue") 
print(g)

g <- ggplot2::ggplot(result_df,
                     aes(x=seq_depth, y=zero_pval, group = factor(taxa_name))) +
  ggplot2::geom_point(aes(color = factor(taxa_name))) +
  ggplot2::geom_line(aes(color = factor(taxa_name))) +
  ggplot2::ggtitle(paste0(project, ': Normalized filter-in vs filter-out zeros Fisher test')) +
  ggplot2::xlab("Min sequence depth per sample") +
  ggplot2::ylab("Pvalue") 
print(g)

g <- ggplot2::ggplot(result_df,
                     aes(x=seq_depth, y=spear_cor^2, group = factor(taxa_name))) +
  ggplot2::geom_point(aes(color = factor(taxa_name))) +
  ggplot2::geom_line(aes(color = factor(taxa_name))) +
  ggplot2::ggtitle(paste0(project, ": Taxa Lev: ","PCA1 vs seq depth")) +
  ggplot2::xlab("Min sequence depth per sample") +
  ggplot2::ylab("R ^ 2") +
  ggplot2::labs(fill = "Transformations") +
  theme(axis.text.x = element_text(angle = 90)) +
  ggplot2::theme_minimal()
print(g)
