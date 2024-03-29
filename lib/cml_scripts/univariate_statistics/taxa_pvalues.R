# Author: Aaron Yerke (aaronyerke@gmail.com)
# Script for creating taxa anova pvalues.
# This was helepful: https://github.com/jsilve24/philr/blob/master/vignettes/philr-intro.Rmd#L142

rm(list = ls()) #clear workspace

# --------------------------------------------------------------------------
print("Defining functions")
# --------------------------------------------------------------------------
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

# --------------------------------------------------------------------------
print("Loading dependencies")
# --------------------------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("optparse", quietly = TRUE)) install.packages("optparse")
library("optparse")

# --------------------------------------------------------------------------
print("Reading cml arguments")
# --------------------------------------------------------------------------
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
                        help="metadata file row to use for row names", metavar="character")
); 

opt_parser <- optparse::OptionParser(option_list=option_list);

opt <- optparse::parse_args(opt_parser);

print(opt)

# --------------------------------------------------------------------------
print("Establishing directory layout and other constants.")
# --------------------------------------------------------------------------
home_dir <- opt$homedir
project <- opt$project
output_dir <- file.path(home_dir, project, 'output')
setwd(file.path(home_dir))

##-Functions--------------------------------------------------------##
source(file.path(home_dir, "lib", "statistical_functions.R"))
source(file.path(home_dir, "lib", "table_manipulations.R"))

# --------------------------------------------------------------------------
print("Importing and preprocessing data.")
# --------------------------------------------------------------------------
asv_table <- data.frame(readRDS(file.path(output_dir, "r_objects", "ForwardReads_DADA2.rds")))

ref_ps <- readRDS(file.path(output_dir, "r_objects", "ref_tree_phyloseq_obj.rds"))
asv_tax <- data.frame(readRDS(file.path(output_dir, "r_objects", "ForwardReads_DADA2_taxonomy.rds")))

metadata <- read.table(opt$metadata, 
                       sep=opt$metadata_delim, 
                       header=TRUE, 
                       row.names = opt$metadata_rowname, 
                       check.names = FALSE,
                       stringsAsFactors=TRUE)
metadata <- metadata[row.names(asv_table),]
# print("Removing non-factor columns from metadata to decrease run time.")
# metadata <- metadata[,sapply(metadata, is.factor)]
# metadata <- Filter(function(x) length(unique(na.omit(x))) < 6, metadata)#hack to remove date columns if they get counted as factors
# print("Metadata columns:")
# print(colnames(metadata))

# --------------------------------------------------------------------------
print("Main loop.")
# --------------------------------------------------------------------------
pdf(file = file.path(output_dir, "graphics", paste0("univariate_pval_tax_", project,".pdf")))
taxa_int <- c()
taxa_name <- c()
meta_name <- c()
meta_col <- c()
pval <- c()
taxa_lev <- c()
for( taxa_col in 1:ncol(asv_tax)){
  my_table <- makeTaxaTable(asv_table, asv_tax, taxa_col)
  for( tx in 1:ncol(my_table)){
    my_taxa <- names(my_table)[tx]
    for(meta in 1:ncol(metadata)){
      my_meta <- metadata[,meta]
      my_lm <- lm(asv_table[,tx] ~ my_meta)
      my_pval <- anova(my_lm)$"Pr(>F)"[1]
      pval <- c(pval, my_pval)
      taxa_int <- c(taxa_int, taxa_col)
      taxa_lev <- c(taxa_lev, names(asv_tax)[taxa_col])
      taxa_name <- c(taxa_name, my_taxa)
      meta_name <- c(meta_name, names(metadata)[meta])
      meta_col <- c(meta_col, meta)
      graphics::boxplot(asv_table[,tx] ~ my_meta,
                        main=paste(project, my_meta, base::formatC(my_pval,format="e", digits=6)),
                        sub=my_taxa, las=2, xlab = my_taxa, ylab = "sequence:")
    }
  }
}
dev.off()

adj_pval <- p.adjust(pval, method = "BH")

dFrame <- data.frame(taxa_lev, taxa_name, taxa_int, meta_name, meta_col, pval, adj_pval)
dFrame <- dFrame [order(dFrame$pval),]
dFrame$adj_pval <- p.adjust( dFrame$pval, method = "BH" )	
write.table(dFrame, file=file.path(output_dir, "tables", paste0(project, "_pValuesUnivariate_taxaVmetadata.csv")), 
            sep=",", row.names = F)

# --------------------------------------------------------------------------
print("Reached end of R script.")
# --------------------------------------------------------------------------






