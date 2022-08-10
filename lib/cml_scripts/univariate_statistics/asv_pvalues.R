# Author: Aaron Yerke (aaronyerke@gmail.com)
# Script for creating sequence anova pvalues.

rm(list = ls()) #clear workspace

# --------------------------------------------------------------------------
print("Defining functions")
# --------------------------------------------------------------------------

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

# --------------------------------------------------------------------------
print("Main loop.")
# --------------------------------------------------------------------------
asv_name <- c()
meta_name <- c()
meta_col <- c()
pval <- c()
my_table <- asv_table
pdf(file = file.path(output_dir, "graphics", "silva_ref_tree.pdf"))
for( tx in 1:ncol(my_table)){
  my_taxa <- names(my_table)[tx]
  for(meta in 1:ncol(metadata)){
    my_meta <- metadata[,meta]
    my_lm <- lm(asv_table[,tx] ~ my_meta)
    my_pval <- anova(my_lm)$"Pr(>F)"[1]
    pval <- c(pval, my_pval)
    asv_name <- c(asv_name, my_taxa)
    meta_name <- c(meta_name, names(metadata)[meta])
    meta_col <- c(meta_col, meta)
  }
}

adj_pval <- p.adjust(pval, method = "BH")

dFrame <- data.frame(asv_name, meta_name, meta_col, pval, adj_pval)
dFrame <- dFrame [order(dFrame$pval),]
dFrame$adj_pval <- p.adjust( dFrame$pval, method = "BH" )	
write.table(dFrame, file=file.path(output_dir, "tables", paste0(project, "_pValuesUnivariate_sequenceVmetadata.csv")), 
            sep=",", row.names = F)

# --------------------------------------------------------------------------
print("Making boxplots ordered by pval.")
# --------------------------------------------------------------------------
dFrame <- read.csv(file=file.path(output_dir, "tables", paste0(project, "_pValuesUnivariate_sequenceVmetadata.csv")), 
                   sep=",")

pdf(file = file.path(output_dir, "graphics", paste0("univariate_pval_seq_", project,".pdf")))
for( rw in 1:nrow(dFrame)){
  taxon <- dFrame$asv_name[rw]
  my_meta <- metadata[,dFrame$meta_name[rw]]
  boxplot(asv_table[,taxon] ~ my_meta,
  main=paste(project, dFrame$meta_name[rw], base::formatC(dFrame$pval[rw],format="e", digits=6)),
  sub=taxon, las=2, xlab = taxon, ylab = "sequence:")
}
dev.off()

pdf(file = file.path(output_dir, "graphics", paste0("top_100_univariate_pval_seq_", project,".pdf")))
for( rw in 1:100){
  taxon <- dFrame$asv_name[rw]
  my_meta <- metadata[,dFrame$meta_name[rw]]
  boxplot(asv_table[,taxon] ~ my_meta,
          main=paste(project, dFrame$meta_name[rw], base::formatC(dFrame$pval[rw],format="e", digits=6)),
          sub=taxon, las=2, xlab = taxon, ylab = "sequence:")
}
dev.off()

# --------------------------------------------------------------------------
print("Reached end of R script.")
# --------------------------------------------------------------------------


