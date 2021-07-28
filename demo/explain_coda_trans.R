# Author: Aaron Yerke (aaronyerke@gmail.com)
# Script for making a script to demonstrate compsitional data transformations
# Fasta file was run through https://www.ebi.ac.uk/Tools/msa/clustalo/ using the phylip option
# This was helepful: https://github.com/jsilve24/philr/blob/master/vignettes/philr-intro.Rmd#L142

rm(list = ls()) #clear workspace

##-functions--------------------------------------------------------##


##-Load Depencencies------------------------------------------------##
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("ALDEx2", quietly = TRUE)) BiocManager::install("ALDEx2")
if (!require("RColorBrewer")) install.packages("RColorBrewer")
library("RColorBrewer")
library("compositions")
library("phyloseq")
library("vegan")
library("DESeq2")
library("philr")
library("ape")
library("ALDEx2")
library("ggplot2")
library("ape")

##-Establish directory layout---------------------------------------##
home_dir <- file.path('~','git','balance_tree_exploration')
project <- "demo"
#home_dir <- file.path('cloud','project')
output_dir <- file.path(home_dir, project, 'output')
setwd(file.path(home_dir))

##-Functions--------------------------------------------------------##
source(file.path(home_dir, "r_libraries", "statistical_functions.R"))
source(file.path(home_dir, "r_libraries", "table_manipulations.R"))

##-Import tables and data preprocessing-----------------------------##

# The palette with black:
cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#0000000")


##-Create demo data-------------------------------------------------##
set.seed(77)
samp_size <- 8
raw_reads <- base::sample(0:100, size=samp_size, replace = T)
names(raw_reads) <- 1:samp_size
alr_reads <- compositions::alr(raw_reads)
names(alr_reads) <- names(raw_reads)[1:samp_size-1]
print(log((raw_reads[1]/raw_reads[8])))

clr_reads <- compositions::clr(raw_reads)
names(clr_reads) <- names(raw_reads)

##-Test the testing data----------------------------------------------##

for(i in 1:length(raw_reads)){
  print(paste( log((raw_reads[i]/compositions::geometricmean(raw_reads))), 
               raw_reads[i] ))
}

ilr_reads <- compositions::ilr(raw_reads)
test <- compositions::ilr(1:3)
test <- compositions::alr(raw_reads[1:2])


##-Make some graphs-------------------------------------------------##
palette( RColorBrewer::brewer.pal(8,"Dark2") )


dat <- data.frame("data" = as.numeric(raw_reads), 
                  "colrs" = names(raw_reads),
                  "x_ax" = rep(1,8))
g <- ggplot2::ggplot(data = dat,
                     aes(x=x_ax, y=data, col = as.factor(colrs),
                         label=as.factor(colrs)))  +
  ggplot2::geom_point(show.legend = FALSE, size = 5) +
  ggplot2::scale_color_brewer(palette="Dark2") +
  ggplot2::ggtitle(paste0("raw", "_Explain CODA")) +
  ggplot2::scale_x_continuous(breaks = NULL) +
  ggplot2::theme(axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank()) +
  ggplot2::theme_bw()
ggplot2::ggsave(filename = file.path(output_dir, "graphics", 
                                     paste0(project, "_raw_data_demo.pdf")),
                height = 4, width = 1, plot = g)
print(g)

dat <- data.frame("data" = as.numeric(alr_reads), 
                  "colrs" = names(alr_reads),
                  "x_ax" = rep(1,8))
g <- ggplot2::ggplot(data = dat,
                     aes(x="x_ax", y=data, col = as.factor(colrs),
                         label=as.factor(colrs))) +
  ggplot2::geom_point(show.legend = FALSE, size = 5) +
  ggplot2::scale_color_brewer(palette="Dark2") +
  ggplot2::ggtitle(paste0("ALR", "_Explain CODA" )) +
  ggplot2::geom_text(aes(color = "black")) +
  ggplot2::theme(axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank()) +
  ggplot2::theme_bw()
ggplot2::ggsave(filename = file.path(output_dir, "graphics", 
                                     paste0(project, "_ralr_demo.pdf")),
                height = 4, width = 1, plot = g)
print(g)

dat <- data.frame("data" = as.numeric(clr_reads), 
                  "colrs" = names(clr_reads),
                  "x_ax" = rep(1,8))
g <- ggplot2::ggplot(data = dat,
                     aes(x=x_ax, y=data, col = as.factor(colrs),
                         label=as.factor(colrs))) +
  ggplot2::geom_point(show.legend = FALSE, size = 5) +
  ggplot2::scale_color_brewer(palette="Dark2") +
  ggplot2::ggtitle(paste0("CLR", "_Explain CODA" )) +
  # ggplot2::geom_text(aes(color = "black")) +
  theme(axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank()) +
  theme_bw()
ggplot2::ggsave(filename = file.path(output_dir, "graphics", 
                                     paste0(project, "_clr_demo.pdf")),
                height = 4, width = 1, plot = g)
print(g)

##-Import tables and data preprocessing-----------------------------##

demo_tree_fasta <- "
>1
TGCTAAGGCCGATGGCGACCGGCGCA
>2
CGCTAAGTTTGATGGCGACCGGCGCA
>3
TGCTAAGTTTGATGGCGACCGGCGCA
>4
TGCTAAGGCTGATGGCGACCGGCGCA
>5
TGCTTTCTCTTGCTGGCGACCGGCGCA
>6
TTGCAAACCAAAGCTGGCGACCAGCG
>7
TCCTTCGGGACTGATTATTTTGTGAC
>8
TGCTAAGGCTGATGGCGACCGGCGCT
"

demo_tree <- ape::read.tree(file.path(home_dir, "demo", "clustalo-p2m-phylip.dnd"))
demo_taxa <- read.table(file.path(home_dir, "demo", "demo_tax.csv"), sep=",")


ape::plot.phylo(demo_tree)


##-Test for reqs----------------------------------------------------##
print("rooted tree?")
is.rooted(demo_tree) # Is the tree Rooted?
print('All multichotomies resolved?')
is.binary.tree(demo_tree) # All multichotomies resolved?

## ---- message=FALSE, warning=FALSE-----------------------------------------
demo_tree <- ape::makeNodeLabel(demo_tree, method="number", prefix='n')

philr::name.balance(phy_tree(ps), tax_table(ps), 'n1')

##-philr transform--------------------------------------------------##
ps.philr <- philr(ps@otu_table, ps@phy_tree,
                  part.weights='enorm.x.gm.counts',
                  ilr.weights='blw.sqrt')




