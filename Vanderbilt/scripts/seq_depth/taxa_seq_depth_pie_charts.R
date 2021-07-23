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
    guides(fill=FALSE)
  return(g)
}
total_cells <- function(df){
  return(dim(df)[1]*dim(df)[2])
}
add_row_to_master <- function(single_row_df, master_df, row_name){
  #adds for merging single row of single taxa counts df to master df 
  #for scatterpie
  #both dfs should have rows as count totals and columns as taxa
  # new_names <- setdiff(names(master_df), names(single_row_df))
  new_names <- names(master_df)[! colnames(master_df) %in% colnames(single_row_df)]
  print(paste("new_names length:", length(new_names)))
  # print(paste(colnames(single_row_df)))
  temp_df <- data.frame(matrix(ncol = length(new_names), nrow = 1))
  names(temp_df) <- new_names
  temp_df <- cbind(single_row_df, temp_df)
  temp_df[is.na(temp_df)] <- 0
  row.names(temp_df) <- c(row_name)
  return(rbind(master_df, temp_df))
}

##-Load Depencencies------------------------------------------------##
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("scatterpie", quietly = TRUE)) BiocManager::install("scatterpie")
library("ggplot2")
library("scatterpie")

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

asv_tax_names <- as.character(tolower(unique(unlist(asv_tax))))

asv_na <- which(is.na(asv_tax_names))

asv_tax_names[asv_na] <- "UNKNOWN"

any(is.na(asv_tax_names))

#create data frame with 0 rows and 5 columns
df <- data.frame(matrix(ncol = 5, nrow = 0))

#provide column names
fi_taxa_master <- data.frame(matrix(ncol = length(asv_tax_names), nrow = 0))
colnames(fi_taxa_master) <- asv_tax_names
fo_taxa_master <- data.frame(matrix(ncol = length(asv_tax_names), nrow = 0))
colnames(fo_taxa_master) <- asv_tax_names

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
spear_cor <- c() #cor(total_seqs[total_seqs >= seq_d], myPCA[,md], method = "spearman")
taxa_left <- c() #number of taxa in the kept data
taxa_out <- c() #number of taxa in the data that is filtered away
num_otu_in <- c() #number of otus in the kept data
num_otu_out <- c() #number of otus in the throwaway data
zero_count <- c()
fi_left_zeros <- c()
fo_left_zeros <- c()
taxa_pval <- c()
zero_pval <- c()
fi_zero_percent <- c()
fo_zero_percent <- c()
fi_nonzero_percent <- c()
fo_nonzero_percent <- c()

counter <- 1
##-Loop and build result DF-----------------------------------------##
for( taxa_col in 1:ncol(asv_tax)){
  
  tname <- names(asv_tax)[taxa_col] #taxa name
  taxa_plots_fi <- list()
  taxa_plots_fo <- list()
  zero_plots_fi <- list()
  zero_plots_fo <- list()
  for(s in 1:length(min_seq_depths)){
    seq_d <- min_seq_depths[s]#new sequencing depth
    sd_filt_asv <- asv_table[total_seqs$total_seqs >= seq_d,]#dataset 1
    filter_out <- asv_table[total_seqs$total_seqs < seq_d,]
    my_table <- makeTaxaTable(sd_filt_asv, asv_tax, taxa_col)#filter in
    fo_table <- makeTaxaTable(filter_out, asv_tax, taxa_col)#filter out
    #make taxa pie_charts
    taxa_char <- substr(names(asv_tax)[taxa_col],1,1)
    # fi_plot <- make_taxa_pie_chart(my_table, paste( taxa_char, "fi t", seq_d))
    # fo_plot <- make_taxa_pie_chart(fo_table, paste(taxa_char, "fo t", seq_d))
    # taxa_plots_fi[[s]] <- fi_plot
    # taxa_plots_fo[[s]] <- fo_plot
    
    #make taxa pie_charts
    fi_zeros <- sum(my_table == 0)#/total_cells(my_table)
    fi_nonzeros <- (total_cells(my_table) - fi_zeros)#/total_cells(my_table)
    fo_zeros <- sum(fo_table == 0)#/total_cells(fo_table)
    fo_nonzeros <- (total_cells(fo_table) - fo_zeros)#/total_cells(fo_table)
    
    fo_all_zeros <- data.frame("non-zeros" = c(fo_nonzeros),
                           "zeros" = c(fo_zeros))
    zero_plots_fo[[s]] <- make_taxa_pie_chart(fo_all_zeros,
                                              paste(taxa_char, "fo z", seq_d))
    fi_all_zeros <- data.frame("non-zeros" = c(fi_nonzeros),
                           "zeros" = c(fi_zeros))
    zero_plots_fi[[s]] <- make_taxa_pie_chart(fi_all_zeros,
                                              paste(taxa_char, "fi z", seq_d))

    #Do stats - zeros first
    fisher_data <- rbind( fi_all_zeros, fo_all_zeros, deparse.level = 1 )
    print(paste0("Fisher zero data taxa_lev: ", tname, "seq depth: ", seq_d))
    print(fisher_data)
    zero_ft <- fisher.test(fisher_data, simulate.p.value = T)
    zero_chi <- chisq.test(c(fi_zeros,fi_nonzeros), c(fo_zeros,fo_nonzeros), simulate.p.value = T, rescale.p = T)
    
    #Taxas
    fi_taxa <- data.frame("fi" = colSums(my_table),
                         row.names = replace( colnames(my_table), is.na(colnames(my_table)), "UNKNOWN"))
    fo_taxa <- data.frame("fo" = colSums(fo_table),
                          row.names = replace( colnames(fo_table), is.na(colnames(fo_table)), "UNKNOWN"))
    comb_tax <- merge(fi_taxa, by = "row.names", fo_taxa, all = T)
    comb_tax <- comb_tax[-c(1)]

    comb_tax[is.na(comb_tax)] <- 0
    taxa_chi <- chisq.test(comb_tax$fo, comb_tax$fi, simulate.p.value = T)
    
    condensed_fo_table <- fo_table[, colSums(fo_table != 0) > 0]
    condensed_fi_table <- my_table[, colSums(my_table != 0) > 0]
    
    # fo_taxa_master <- merge(fo_taxa_master, fo_table, all.y = T)
    # fi_taxa_master <- merge(fi_taxa_master, fo_table, all.y = T)
    
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
      taxa_left[counter] <- sum(colSums(my_table != 0) > 0)
      taxa_out[counter] <- sum(colSums(fo_table != 0) > 0)
      taxa_left[counter] <- ncol(my_table)
      taxa_lev[counter] <- taxa_col
      taxa_name[counter] <- tname
      fi_left_zeros[counter] <- fi_zeros
      fo_left_zeros[counter] <- fo_zeros
      taxa_pval[counter] <- taxa_chi$p.value
      zero_pval[counter] <- zero_ft$p.value
      fi_zero_percent[counter] <- sum(my_table == 0)/total_cells(my_table)
      fo_zero_percent[counter] <- sum(fo_table == 0)/total_cells(fo_table)
      fi_nonzero_percent[counter] <- 1 - fi_zero_percent[counter]
      fo_nonzero_percent[counter] <- 1 - fi_zero_percent[counter]
      fo_taxa_master <- add_row_to_master(t(as.data.frame(fo_taxa)), fo_taxa_master, counter)
      fi_taxa_master <- add_row_to_master(t(as.data.frame(fi_taxa)), fi_taxa_master, counter)
      num_otu_in <- sum(fi_taxa) #number of otus in the kept data
      num_otu_out <- sum(fo_taxa) #number of otus in the throwaway data

      
      counter <- counter + 1
    }
    saved_df_name <- paste0("fi_asv_raw_", tname, "_seq_dep_", seq_d)
    write.table(my_table, file=file.path(output_dir, "tables", "seq_dept_taxa_tables", saved_df_name))
    saved_df_name <- paste0("fo_asv_raw_", tname, "_seq_dep_", seq_d)
    write.table(fo_table, file=file.path(output_dir, "tables", "seq_dept_taxa_tables", saved_df_name))
    
  }
}

result_df <- data.frame(taxa_lev, taxa_name, mds_lev, seq_depth, 
                        var_exp, spear_cor, fi_left_zeros, fo_left_zeros,
                        taxa_pval, zero_pval, fi_zero_percent, fo_zero_percent,
                        fi_nonzero_percent, fo_nonzero_percent, taxa_out, taxa_left,
                        num_otu_in, num_otu_out)

write.table(result_df, 
            file = file.path(output_dir, "tables", paste0(project, "taxa_seq_depth_pie_charts_results.csv")),
            sep = ",")

result_df <- read.table(file = file.path(output_dir, "tables", paste0(project, "taxa_seq_depth_pie_charts_results.csv")),
                        sep = ",")


pdf(file = file.path(output_dir, "graphics", paste0(tname,"_seq_depth_taxa_pie.pdf")))
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
                     aes(x=seq_depth, y=log(taxa_pval), group = factor(taxa_name))) +
  ggplot2::geom_point(aes(color = factor(taxa_name))) +
  ggplot2::geom_line(aes(color = factor(taxa_name))) +
  ggplot2::ggtitle(paste0(project, ': Non-normalized filter-in vs filter-out taxa chi sq')) +
  ggplot2::xlab("Min sequence depth per sample") +
  ggplot2::ylab("log(pvalue)") 
print(g)

g <- ggplot2::ggplot(result_df,
                     aes(x=seq_depth, y=log(zero_pval), group = factor(taxa_name),
                         label=round(fo_zero_percent, 2)*100)) +
  ggplot2::geom_point(aes(color = factor(taxa_name))) +
  ggplot2::geom_line(aes(color = factor(taxa_name))) +
  ggplot2::geom_text() +
  ggplot2::ggtitle(paste0(project, ': Non-normalized filter-in vs filter-out zeros Fisher test')) +
  ggplot2::xlab("Min sequence depth per sample") +
  ggplot2::ylab("Log(Pvalue)") 
print(g)

# useful for scatter pie charts
# https://cran.r-project.org/web/packages/scatterpie/vignettes/scatterpie.html
ggplot2::ggplot() +
  geom_scatterpie(aes(x=log10(seq_depth)*100, y=log10(zero_pval), group = factor(taxa_name), r=4),
                  data=result_df, 
                  cols=c("fo_nonzero_percent", "fo_zero_percent")) + 
  ggplot2::geom_line(aes(x=log10(seq_depth)*100, y=log10(zero_pval), color = factor(taxa_name))) +
  # ggplot2::theme(legend.position = "none") +
  ggplot2::coord_equal() +
  ggplot2::ggtitle(paste0(project, ": Filter-out zeros")) +
  ggplot2::ylab("log(Fisher test(filter-in, filter-out zeros))") 
# ggsave(g, filename = file.path(output_dir, "graphics", paste0("test","_seq_depth_taxa_pie.pdf")), height = 5, width = 40)

                           # cols=c("fi_nonzero_percent", "fi_zero_percent")) + coord_equal()

#Zero scatterpie
ggplot2::ggplot() +
  geom_scatterpie(aes(x=log10(seq_depth)*100, y=log10(zero_pval), group = factor(taxa_name), r=4),
                  data=result_df, 
                  cols=c("fi_nonzero_percent", "fi_zero_percent")) + 
  ggplot2::geom_line(aes(x=log10(seq_depth)*100, y=log10(zero_pval), color = factor(taxa_name))) +
  # ggplot2::theme(legend.position = "none") +
  ggplot2::ylab("log(Fisher test(filter-in, filter-out zeros))") +
  ggplot2::coord_equal() +
  ggplot2::ggtitle(paste0(project, ": Filter-in zeros"))

ggplot2::ggplot() +
  geom_scatterpie(aes(x=log10(seq_depth)*100, y=log10(zero_pval), group = factor(taxa_name), r=4),
                  data=result_df, 
                  cols=c("fi_nonzero_percent", "fi_zero_percent")) + 
  ggplot2::geom_line(aes(x=log10(seq_depth)*100, y=log10(zero_pval), color = factor(taxa_name))) +
  # ggplot2::theme(legend.position = "none") +
  ggplot2::ylab("log(Fisher test(filter-in, filter-out zeros))") +
  ggplot2::coord_equal() +
  ggplot2::ggtitle(paste0(project, ": Filter-in zeros"))

#taxa scatterpie
fo_result <- cbind(result_df, fo_taxa_master)
g <- ggplot2::ggplot() +
  geom_scatterpie(aes(x=log(seq_depth)*10, y=taxa_pval*100, group = factor(taxa_name), r=4),
                  data=fo_result,
                  cols=asv_tax_names) + 
  ggplot2::geom_line(aes(x=log(seq_depth)*10, y=taxa_pval*100, color = factor(taxa_name))) +
  ggplot2::theme(legend.position = "none") +
  ggplot2::ylab("chi sqr pval * 100") +
  ggplot2::xlab("log10(seq_depth)*10") +
  ggplot2::coord_equal() +
  ggplot2::ggtitle(paste0(project, ": Filter-out taxa"))
ggsave(g, filename = file.path(output_dir, "graphics", paste0("fo_taxa","_seq_depth_taxa_pie.pdf")), height = 8, width = 12)

fi_result <- cbind(result_df, fi_taxa_master)
g <- ggplot2::ggplot() +
  geom_scatterpie(aes(x=log10(seq_depth)*25, y=log10(taxa_pval)*25, group = factor(taxa_name), r=4),
                  data=fi_result,
                  cols=asv_tax_names) + 
  ggplot2::geom_line(aes(x=log10(seq_depth)*25, y=log10(taxa_pval)*25, color = factor(taxa_name))) +
  ggplot2::theme(legend.position = "none") +
  ggplot2::ylab("log10(chi sqr pval)*25") +
  ggplot2::xlab("log10(chi sqr pval)*25") +
  ggplot2::coord_equal() +
  ggplot2::ggtitle(paste0(project, ": Filter-in taxa"))
ggsave(g, filename = file.path(output_dir, "graphics", paste0("fi_taxa","_seq_depth_taxa_pie.pdf")), height = 8, width = 12)

genus_all <- tolower(unique(unlist(asv_tax$Genus)))
genus_all <- genus_all[!is.na(genus_all)]
fi_master_genus <- fi_taxa_master[, c(genus_all)]
fo_master_genus <- fo_taxa_master[, c(genus_all)]
fo_result_gen <- cbind(result_df, fo_master_genus)[taxa_lev == 6,]
fi_result_gen <- cbind(result_df, fi_master_genus)[taxa_lev == 6,]

fi_result_gen <- subset(fi_result, taxa_lev == 6, 
                        select = c(colnames(result_df), genus_all))

g <- ggplot2::ggplot() +
  geom_scatterpie(aes(x=log10(fi_result_gen$seq_depth)*25, 
                      y=log10(fi_result_gen$taxa_pval)*25),
                  data=fi_result_gen,
                  cols=c(17:ncol(fi_result_gen))) + 
  ggplot2::geom_line(aes(x=log10(seq_depth)*25, y=log10(taxa_pval)*25)) +
  ggplot2::theme(legend.position = "none") +
  ggplot2::ylab("log10(chi sqr pval)*25") +
  ggplot2::xlab("log10(chi sqr pval)*25") +
  ggplot2::coord_equal() +
  ggplot2::ggtitle(paste0(project, ": Filter-in GENUS ONLY"))
# ggsave(g, filename = file.path(output_dir, "graphics", paste0("fi_genus_only","_seq_depth_taxa_pie.pdf")), height = 8, width = 12)
print(g)


setdiff(unlist(colnames(fo_master_genus)),  unlist(colnames(fo_result_gen)) )

# fo_result_gen <- data.frame(fo_result[taxa_lev == 6,])
p <- ggplot2::ggplot() +
  geom_scatterpie(aes(x=log(seq_depth, base = 5), y=log10(taxa_pval)*25, 
                      color = factor(taxa_name)),
                  data=subset(fo_result_gen, taxa_lev == 6), 
                  cols=c("fi_nonzero_percent", "fi_zero_percent"))+
  ggplot2::geom_line(aes(color = factor(taxa_name))) +
  ggplot2::theme(legend.position = "none") +
  ggplot2::ylab("log10(chi sqr pval)*25") +
  ggplot2::xlab("log10(chi sqr pval)*25") +
  ggplot2::coord_equal() +
  ggplot2::ggtitle(paste0(project, ": Filter-out GENUS ONLY"))
# ggsave(g, filename = file.path(output_dir, "graphics", paste0("fo_genus_only","_seq_depth_taxa_pie.pdf")), height = 8, width = 12)
print(p)


g <- ggplot2::ggplot(result_df,
                     aes(x=log(seq_depth), y=fo_zero_percent, 
                         group = factor(taxa_name),
                         )) +
  ggplot2::geom_point(aes(color = factor(taxa_name))) +
  ggplot2::geom_line(aes(color = factor(taxa_name))) +
  # ggplot2::geom_text() +
  ggplot2::ggtitle(paste0(project, ": Discarded data: zero percent")) +
  ggplot2::xlab("Min sequence depth per sample") +
  ggplot2::ylab("Percent") + 
  ggplot2::theme(text = element_text(size = 20))
print(g)



g <- ggplot2::ggplot(result_df,
                     aes(x=log(seq_depth), y=taxa_left, 
                         group = factor(taxa_name),
                     )) +
  ggplot2::geom_point(aes(color = factor(taxa_name))) +
  ggplot2::geom_line(aes(color = factor(taxa_name))) +
  # ggplot2::geom_text() +
  ggplot2::ggtitle(paste0(project, ": Number of taxa in data kept vs log of seq depth")) +
  ggplot2::xlab("Log of min sequence depth per sample") +
  ggplot2::ylab("Number of taxa") 
print(g)

g <- ggplot2::ggplot(result_df,
                     aes(x=log(seq_depth), y=taxa_out, 
                         group = factor(taxa_name),
                     )) +
  ggplot2::geom_point(aes(color = factor(taxa_name))) +
  ggplot2::geom_line(aes(color = factor(taxa_name))) +
  # ggplot2::geom_text() +
  ggplot2::ggtitle(paste0(project, ": Number of taxa in removed data vs log of seq depth")) +
  ggplot2::xlab("Log of min sequence depth per sample") +
  ggplot2::ylab("Number of taxa") 
print(g)

g <- ggplot2::ggplot(result_df,
                     aes(x=log(seq_depth), y=num_otu_in, 
                         group = factor(taxa_name),
                     )) +
  ggplot2::geom_point(aes(color = factor(taxa_name))) +
  ggplot2::geom_line(aes(color = factor(taxa_name))) +
  # ggplot2::geom_text() +
  ggplot2::ggtitle(paste0(project, ": Number of OTUs in data kept vs log of seq depth")) +
  ggplot2::xlab("Log of min sequence depth per sample") +
  ggplot2::ylab("Number of taxa") 
print(g)

g <- ggplot2::ggplot(result_df,
                     aes(x=log(seq_depth), y=num_otu_out, 
                         group = factor(taxa_name),
                     )) +
  ggplot2::geom_point(aes(color = factor(taxa_name))) +
  ggplot2::geom_line(aes(color = factor(taxa_name))) +
  # ggplot2::geom_text() +
  ggplot2::ggtitle(paste0(project, ": Number of OTUs in data tossed vs log of seq depth")) +
  ggplot2::xlab("Log of min sequence depth per sample") +
  ggplot2::ylab("Number of taxa") 
print(g)



g <- ggplot2::ggplot(result_df,
                     aes(x=log(seq_depth), y=fi_zero_percent, 
                         group = factor(taxa_name),
                     )) +
  ggplot2::geom_point(aes(color = factor(taxa_name))) +
  ggplot2::geom_line(aes(color = factor(taxa_name))) +
  # ggplot2::geom_text() +
  ggplot2::ggtitle(paste0(project, ": Kept data: zero percent")) +
  ggplot2::xlab("Min sequence depth per sample") +
  ggplot2::ylab("Percent") 
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

dev.off()

fo_no_zeros <- fo_taxa_master [, colSums(fo_taxa_master != 0) > 0]

t_pval <- vector(mode = "numeric", length = ncol(fo_no_zeros))
tax_name <- vector("character", length = ncol(fo_no_zeros))
taxa_level <- vector("character", length = ncol(fo_no_zeros))

for(co in 1:ncol(fo_no_zeros)){
  taxa <- colnames(fo_no_zeros)[co]
  for(tax_co in 1:ncol(asv_tax)){
    if(taxa %in% tolower(unique(as.vector(asv_tax[,tax_co])))){
      my_tl <- as.character(colnames(asv_tax)[tax_co])
      break
    }#end if(taxa...
  }#end for(tax_co in 1:ncol(asv_tax))...
  fo_col <- subset( fo_result, taxa_name==my_tl)[,taxa]
  fi_col <- subset( fi_result, taxa_name==my_tl)[,taxa]
  # print(taxa)
  # print(fo_col)
  # print(fi_col)
  tax_name[co] <- taxa
  t_pval[co] <- t.test(fo_col, fi_col)$p.value
  taxa_level[co] <- my_tl
}

dFrame <- data.frame(taxa_level, tax_name, t_pval)
dFrame <- dFrame [order(dFrame$t_pval),]
dFrame$adj_pval <- p.adjust( dFrame$t_pval, method = "BH" )	
# write.table(dFrame, file=file.path(output_dir, "tables", "Vanderbilt_pValuesUnivariate_taxaVmetadata.tsv"), 
#             sep="\t", 
#             row.names=FALSE)

pdf(file = file.path(output_dir, "graphics", paste0(project,"_all_seq_depths_filters_scatter.pdf")))
for (r in row.names(dFrame)){
  taxa <- as.character(dFrame[r, "tax_name"])
  # pval <- 
  for(tax_co in 1:ncol(asv_tax)){
    if(taxa %in% tolower(unique(as.vector(asv_tax[,tax_co])))){
      my_tl <- as.character(colnames(asv_tax)[tax_co])
      break
    }#end if(taxa...
  }#end for(tax_co in 1:ncol(asv_tax))...
  fo_col <- subset( fo_result, taxa_name==my_tl)[,c(taxa, "seq_depth")]
  fi_col <- subset( fi_result, taxa_name==my_tl)[,taxa]
  # dat <- data.frame("Counts" = c(fo_col, fi_col),
  #                   "Filter_groups" = c(rep("Filter out", length(fo_col)), 
  #                   rep("Filter in", length(fo_col))))
  # g <- ggplot2::ggplot(aes(x=subset(fo_result, taxa_name==my_tl)[,"seq_depth"],
  #                          y=subset(fo_result, taxa_name==my_tl)[,taxa])) +
    # ggplot2::geom_boxplot(data=dat, 
    #                       aes(x=Filter_groups, y=Counts, 
    #                           fill=Filter_groups),
    #   outlier.colour="black", outlier.shape=16,
    #                     outlier.size=2, notch=FALSE) +
    # ggplot2::geom_point() +
    # ggplot2::geom_line() +
    # ggplot2::ggtitle(paste0(project, 
    #                         " : Lev: ", my_tl," tax: ", taxa, " adj p: ", dFrame[r, "adj_pval"]))
  # print(g)
  plot(x=subset(fo_result, taxa_name==my_tl)[,"seq_depth"],
       y=subset(fo_result, taxa_name==my_tl)[,taxa],
       type = "b",
       main = paste0(project, " : Lev: ", my_tl," tax: ", taxa, " adj p: ", dFrame[r, "adj_pval"]),
       ylab = taxa,
       xlab = "seq depth threshold")
  
  
}
dev.off()


zeros_col_sum <- rowSums(asv_table == 0)
cor(zeros_col_sum, total_seqs, method = "kendall")
