# Author: Aaron Yerke
# This is a pipeline that was created to explore the philr transform.
# Many resources were used to create this pipeline:
#   Processing the sequences through DADA2, making trees with phangorn, and putting them into phyloseq objects:
#     This is by far the biggest source:
#     https://github.com/spholmes/F1000_workflow/blob/master/MicrobiomeWorkflow/MicrobiomeWorkflowII.Rmd
#     or (in article format)
#     https://f1000research.com/articles/5-1492/v2
#     
#   DADA2:
#     https://benjjneb.github.io/dada2/tutorial.html
#   
#   Using philr:
#     https://bioconductor.org/packages/release/bioc/html/philr.html
#     Code example:
#       https://bioconductor.org/packages/release/bioc/vignettes/philr/inst/doc/philr-intro.R

rm(list = ls()) #clear workspace

# ‘ape’, ‘dplyr’, ‘reshape2’, ‘plyr’
# .cran_packages <- c("ggplot2", "gridExtra")
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager", type = "source", repos = "http://archive.linux.duke.edu/cran/")
  if (!requireNamespace("phangorn", quietly = TRUE))
    install.packages("phangorn",type = "source", repos = "http://archive.linux.duke.edu/cran/")
  BiocManager::install("phyloseq")
  BiocManager::install("DECIPHER")
}

library("DECIPHER")
library("phangorn")
library("phyloseq")

##----------------Establish directory layout------------------------##
home_dir = file.path('~','git','balance_tree_exploration')
#home_dir = file.path('cloud','project')
output_dir = file.path(home_dir, "ava_c", 'output')
project = "philr_western_gut"
f_path <- file.path(home_dir, "sequences") # CHANGE ME to the directory containing the fastq files after unzipping.
# list.files(f_path)

setwd(file.path(home_dir))

print("Established directory layout")

##---------------------Import R objects-----------------------------##
con <- gzfile(file.path( output_dir, "r_objects", "ForwardReads_DADA2.rds"))
seqtab = readRDS(con)
close(con)

con <- gzfile(file.path( output_dir, "r_objects", "ForwardReads_DADA2_alignment.rds"))
alignment <- readRDS(con)
close(con)

con <- gzfile(file.path( output_dir, "r_objects", "ForwardReads_DADA2_taxonomy.rds"))
taxTab <- readRDS(con)
close(con)

print("Imported R objects")

##------------------------Build tree--------------------------------##
phangAlign <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phangAlign)#create distance matrix
treeNJ <- upgma(dm) #make tree
fit = pml(treeNJ, data=phangAlign)#fit model
fitGTR <- update(fit, k=4, inv=0.2)#fit model with updated parameters
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
print("phangorn completed")

myMeta = read.table(file.path(home_dir, "ava_c", "fullMetadata.tsv"), 
                    sep="\t", 
                    header=TRUE, 
                    check.names = FALSE,
                    stringsAsFactors=FALSE)

ps <- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE),
               sample_data(myMeta),
               tax_table(taxTab),
               phy_tree(fitGTR$tree))
# ps
print("Created ps")

#examine tree
library("ape")

pdf(paste0(project,"_2", ".pdf"))

plot_tree(ps, "treeonly", nodeplotblank, ladderize="left")

plot_tree(ps, ladderize="left", color="Religion")

dev.off()

# setwd(file.path(home_dir, "philr_pipelines"))
saveRDS(ps, file.path(output_dir, "r_objects","denovo_tree_UPGMA_phyloseq_obj.rds"))
