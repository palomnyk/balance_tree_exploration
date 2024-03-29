#!/usr/bin/env Rscript
# script for parsing the blast results to build a tree

# Columns follow this pattern:
# qseqid sseqid pident length evalue bitscore score ppos
# 
# 1.	qseqid	 query (e.g., unknown gene) sequence id
# 2.	sseqid	 subject (e.g., reference genome) sequence id
# 3.	pident	 percentage of identical matches
# 4.	length	 alignment length (sequence overlap)
# 5.	evalue	 expect value
# 6.	bitscore	 bit score
# 7.  score     Raw score
# 8.  ppos      Percentage of positive-scoring matches
# Comparison stategy: compare ppos (if tie, go with biggest alignment length)

print("Begining cml argument processing.")
if (!requireNamespace("optparse", quietly = TRUE)){
  install.packages("optparse")
}
library("optparse")

option_list <- list(
  make_option(c("-d", "--homedir"), type="character", 
              default=file.path('~','git','balance_tree_exploration'), 
              help="dataset dir path", metavar="character"),
  make_option(c("-p", "--project"), type="character", default=NULL, 
              help="project folder name in homedir", metavar="character")
); 

opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);

print(opt)

print("Establishing directory layout.")
home_dir <- opt$homedir
project <- opt$project
output_dir <- file.path(home_dir, project, 'output')

setwd(file.path(output_dir, "tree_process_blast"))

input_file <- "output.txt"

df <- data.frame(#qseqid=character(),
                 sseqid=character(),
                 bitscore=integer(),
                 count_matched=integer(),
                 stringsAsFactors=FALSE)

print("Parsing blast output.")
con <- file(input_file, "r")
while ( TRUE ) {
  line = readLines(con, n = 1)
  if ( length(line) == 0 ) {
    break
  }
  result = unlist(strsplit(line, '\t'))
  qseq = result[1]
  sseq = unlist(strsplit(result[2],"\\|"))[2] #dbj|AB064923|
  evalue = as.numeric(result[5])
  bitsc = as.integer(result[6])
  # print(paste(sseq, bitsc))
  if (evalue < 10^-10){
    if ( qseq %in% row.names(df)){
      
      if (bitsc > df[qseq, "bitscore"]){
        count = df[qseq,count_matched]
        df[qseq,] = list(sseq, bitsc, count + 1)
      }
      
    }else{
      newRow <- data.frame(sseqid=sseq,
                           bitscore=bitsc,
                           count_matched=1) 
      row.names(newRow) = qseq
      df = rbind(df, newRow)
    }

  }

  # break
}
close(con)

print(paste("original nrow:", nrow(df)))

print(paste("number of unique rows:", length(unique(rownames(df)))))

myT <- table(df[,"sseqid"])

print(paste("ave seq/node:", mean(myT), "\nmax seq/node:", max(myT)))

pdf(file = file.path(output_dir, "graphics", "p3_parse_blast.pdf"))
hist(myT, breaks = 150, xlab = "Sequences per node tip", main = "Histogram of seqs per node tip")
barplot(myT, las = 2, xlab = "Sequences per node tip", main = "Histogram of seqs per node tip")
dev.off()

# df = df[!duplicated(df),]
# 
# print(paste("deduplicated nrow:", nrow(df)))
# 
# myT = table(df[,"sseqid"])
# 
# print(paste("ave seq/node:", mean(myT), "\nmax seq/node:", max(myT)))
# 
# hist(myT, breaks = 150, xlab = "Sequences per node tip", main = "Histogram of seqs per node tip")
# barplot(myT, las = 2, xlab = "Sequences per node tip", main = "Histogram of seqs per node tip")

write.csv(df, file = "parsed_output.csv")

con <- gzfile(file.path( output_dir, "r_objects", "ForwardReads_DADA2_taxonomy.rds"))
taxTab <- readRDS(con)
close(con)

print(nrow(taxTab))