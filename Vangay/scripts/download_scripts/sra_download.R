# Author: Aaron Yerke
# Script to download from the SRA
# This was helepful: https://github.com/jsilve24/philr/blob/master/vignettes/philr-intro.Rmd#L142

rm(list = ls()) #clear workspace

##-Load Depencencies------------------------------------------------##
install.packages("SRAdb")
remotes::install_version("SDMTools", "1.1-221")
library(SRAdb)
srafile = getSRAdbFile()
con = dbConnect('SQLite',srafile)
##-Establish directory layout---------------------------------------##
home_dir <- file.path('~','git','balance_tree_exploration')
project <- "Vangay"
download_dir <- file.path(home_dir, project, "downloaded_seqs")
##-Download SRA files ----------------------------------------------##
sra_run_table <- read.table(file.path(home_dir, project, "SraRunTable.txt"),
                            sep = ",",
                            header = TRUE)

my_rows <- sra_run_table$Assay.Type == "AMPLICON"
my_accessions <- sra_run_table$Run[my_rows]

setwd(download_dir)
for (run in my_accessions) {
  getSRAfile(run,con,fileType='sra')
}  
  

