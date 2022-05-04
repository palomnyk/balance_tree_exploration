# Author: Aaron Yerke (aaronyerke@gmail.com)
# Script for running mundging hashseq output https://github.com/FarnazFouladi/HashSeq
print("Clearing workspace.")
rm(list = ls()) #clear workspace

print("Install dependencies.")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
library("optparse")

print("Loading cml options")
option_list <- list(
  make_option(c("-d", "--homedir"), type="character", 
              default=file.path('~','git','balance_tree_exploration'), 
              help="dataset dir path", metavar="home dir"),
  make_option(c("-p", "--project"), type="character", default=NULL, 
              help="project folder", metavar="project"),
); 

opt_parser <- optparse::OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);

print(opt)

print("Establish directory layout.")
home_dir <- opt$homedir
project <- opt$project
output_dir <- file.path(home_dir, project, 'output')

print("Read data.")
hashseq <- data.frame(data.table::fread(file = file.path(output_dir,"hashseq", "SvTable.txt"),
                                        header=TRUE, data.table=FALSE), row.names = 1)
# hashseq <- hashseq[, colSums(hashseq != 0) > 0.01*nrow(hashseq)]#remove columns that don't have at least 10%
# print(paste("HashSeq has", ncol(hashseq), "columns after column reduction."))
print("Removing '_' and anything after it.")
row.names(hashseq) <- sapply(as.character(row.names(hashseq)), function(x) {
  strsplit(x,"_")[[1]][1]
  })
print("Fixed formating of row names.")
print("First 10 row names.")
print(paste(row.names(hashseq)[1:10], collapse = " "))
print("First 10 row names.")
print(paste(tail(row.names(hashseq), 10), collapse = " "))
print(paste("number of elements:", nrow(hashseq), collapse = " "))

print("Attempting to reorder rownames based on ASV table from dada2\n")
hashseq <- hashseq[ order(row.names(hashseq)),]
print("First 10 row names.")
print(paste(row.names(hashseq)[1:10], collapse = " "))
print("First 10 row names.")
print(paste(tail(row.names(hashseq), 10), collapse = " "))
print(paste("Number of elements:", nrow(hashseq), collapse = " "))
print("Writing munged hashseq table.\n")

write.table(hashseq, file = file.path(output_dir,"hashseq", "hashseq.csv"),
            sep = ",")

print("R script completed.")



