# Author: Aaron Yerke (aaronyerke@gmail.com)
# Script for making distribution of AUC of ILR balances to compare to philr

rm(list = ls()) #clear workspace

##-functions--------------------------------------------------------##


##-Load Depencencies------------------------------------------------##
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("ape", quietly = TRUE)) BiocManager::install("ape")
if (!requireNamespace("philr", quietly = TRUE)) BiocManager::install("philr")
if (!requireNamespace("phyloseq", quietly = TRUE)) BiocManager::install("phyloseq")
if (!requireNamespace("ggplot2", quietly = TRUE)) BiocManager::install("ggplot2")
if (!requireNamespace("rgr", quietly = TRUE)) install.packages("rgr")
if (!requireNamespace("vegan", quietly = TRUE)){ install.packages("vegan")}
library("vegan")
library("rgr")
library("phyloseq")
library("philr")
library("ggplot2")
library("ape")
print("finished loading libraries")

##-Establish directory layout---------------------------------------##
home_dir <- file.path('~','git','balance_tree_exploration')
project <- "Vanderbilt"
output_dir <- file.path(home_dir, project, 'output')
setwd(file.path(home_dir))

##-Functions--------------------------------------------------------##
source(file.path(home_dir, "r_libraries", "statistical_functions.R"))
source(file.path(home_dir, "r_libraries", "table_manipulations.R"))

##-Set up constants-------------------------------------------------##
rf_cols <- 3:7
main_output_label <- paste0(project, "_", "PCOA")
dist_metric <- "euclidean"

##-Import tables and data preprocessing-----------------------------##
asv_table <- data.frame(readRDS(file.path(output_dir, "r_objects", "ForwardReads_DADA2.rds")))
total_seqs <- rowSums(asv_table)
total_seqs <- data.frame("total_seqs"=total_seqs, "duplicate" = total_seqs,
                         row.names = row.names(asv_table))

print("Cleaning Ref tree otu with philr tutorial normalization")
ref_ps <- readRDS(file.path(output_dir, "r_objects", "ref_tree_phyloseq_obj.rds"))
clean_otu <- data.frame(ref_ps@otu_table@.Data)
clean_otu <- philr_tutorial_normalization(clean_otu)
print(paste("nrow orginal ref:", nrow(ref_ps@otu_table), "nrow clean ref: ", nrow(clean_otu)))

phy_tree(ref_ps) <- ape::makeNodeLabel(phy_tree(ref_ps), method="number", prefix='n')
ref_ps_clean <- phyloseq::phyloseq( otu_table(clean_otu, taxa_are_rows = F), 
                                    phy_tree(ref_ps@phy_tree),
                                    tax_table(ref_ps@tax_table), 
                                    sample_data(ref_ps@sam_data))
print("Cleaning UPGMA tree otu with philr tutorial normalization")
denovo_tree_ps <- readRDS(file.path(output_dir, "r_objects", "denovo_tree_UPGMA_phyloseq_obj.rds"))
clean_den_otu <- philr_tutorial_normalization(data.frame(denovo_tree_ps@otu_table@.Data))
print(paste("nrow orginal denovo:", nrow(denovo_tree_ps@otu_table), "nrow clean denovo otu: ", nrow(clean_den_otu)))
cln_denovo_tree_ps <- phyloseq::phyloseq( otu_table(clean_den_otu, taxa_are_rows = F),
                                          phy_tree(ape::makeNodeLabel(phy_tree(denovo_tree_ps@phy_tree))),
                                          tax_table(denovo_tree_ps@tax_table), 
                                          sample_data(denovo_tree_ps@sam_data))
denovo_tree_ps <- phyloseq::transform_sample_counts(denovo_tree_ps, function(x) x + 1 )

phy_tree(ref_ps_clean) <- makeNodeLabel(phy_tree(ref_ps_clean), method="number", prefix='n')
phy_tree(cln_denovo_tree_ps) <- makeNodeLabel(phy_tree(cln_denovo_tree_ps), method="number", prefix='n')
ref_ps <- phyloseq::transform_sample_counts(ref_ps, function(x) x + 1 )
phy_tree(ref_ps) <- philr::makeNodeLabel(phy_tree(ref_ps), method="number", prefix='n')
##-Random num seed--------------------------------------------------##
set.seed(36)
print("making random trees")
orig_ref_rand_list <- list()
for (rand in 1:10){
  rand_tree <- rtree(n = length(ref_ps@phy_tree$tip.label), tip.label = ref_ps@phy_tree$tip.label)
  #put int in philr
  rand_tree_ps <- phyloseq::phyloseq( otu_table(clean_otu, taxa_are_rows = F),
                                      phy_tree(rand_tree),
                                      tax_table(ref_ps@tax_table),
                                      sample_data(ref_ps@sam_data))
  phy_tree(rand_tree_ps) <- ape::makeNodeLabel(phy_tree(rand_tree_ps), method="number", prefix='n')
  orig_ref_rand_list[[rand]] <- rand_tree_ps
}

print("make random trees for cln upgma taxa")
cln_upgma_rand_list <- list()
for (rand in 1:10){
  rand_tree <- rtree(n = length(cln_denovo_tree_ps@phy_tree$tip.label), tip.label = cln_denovo_tree_ps@phy_tree$tip.label)
  #put int in philr
  rand_tree_ps <- phyloseq::phyloseq( otu_table(cln_denovo_tree_ps, taxa_are_rows = F),
                                      phy_tree(rand_tree),
                                      tax_table(cln_denovo_tree_ps@tax_table),
                                      sample_data(cln_denovo_tree_ps@sam_data))
  phy_tree(rand_tree_ps) <- ape::makeNodeLabel(phy_tree(rand_tree_ps), method="number", prefix='n')
  cln_upgma_rand_list[[rand]] <- rand_tree_ps
}

print("make random trees for clean ref taxa")
cln_ref_rand_list <- list()
for (rand in 1:10){
  rand_tree <- rtree(n = length(ref_ps_clean@phy_tree$tip.label), tip.label = denovo_tree_ps@phy_tree$tip.label)
  #put int in philr
  rand_tree_ps <- phyloseq::phyloseq(otu_table(ref_ps_clean, taxa_are_rows = F),
                                     phy_tree(rand_tree),
                                     tax_table(ref_ps_clean@tax_table),
                                     sample_data(ref_ps_clean@sam_data))
  phy_tree(rand_tree_ps) <- ape::makeNodeLabel(phy_tree(rand_tree_ps), method="number", prefix='n')
  cln_ref_rand_list[[rand]] <- rand_tree_ps
}

print("creating lognorm, ALR and CLR")
if (dir.exists(file.path(output_dir,"r_objects", "lognorm_asv.rds"))) {
  ln_asv_tab <- readRDS(file.path(output_dir,"r_objects", "lognorm_asv.rds"))
}else{
  ln_asv_tab <- lognorm(asv_table)
  saveRDS(ln_asv_tab, file = file.path(output_dir,"r_objects", "lognorm_asv.rds"))
}

my_zeros <- apply(asv_table, 2, function(x) {
  return(sum(x == 0))
})
alr_col <- which(my_zeros == min(my_zeros))[1]
# alr_col_num <- grep(alr_col, colnames(asv_table))
print("creating ALR")
if (file.exists(file.path(output_dir,"r_objects", "alr_asv.rds"))) {
  my_alr <- readRDS(file.path(output_dir,"r_objects", "alr_asv.rds"))
}else{
  my_alr <- as.data.frame(rgr::alr(as.matrix(asv_table + 1), j = as.numeric(alr_col)))
  saveRDS(my_alr, file = file.path(output_dir,"r_objects", "alr_asv.rds"))
}
print("creating CLR")
if (dir.exists(file.path(output_dir,"r_objects", "clr_asv.rds"))) {
  my_clr <- readRDS(file.path(output_dir,"r_objects", "clr_asv.rds"))
}else{
  my_clr <- as.data.frame(rgr::clr(as.matrix(asv_table + 1)))
  saveRDS(my_clr, file = file.path(output_dir,"r_objects", "clr_asv.rds"))
}
print("loading and munging metadata")
metadata <- read.table(file.path(home_dir, project, "patient_metadata.tsv"), 
                       sep="\t", 
                       header=TRUE, 
                       row.names = "Run", 
                       check.names = FALSE,
                       stringsAsFactors=TRUE)
metadata <- metadata[row.names(metadata) %in% row.names(clean_otu), ]
metadata$type <- droplevels(metadata$type)
metadata$type <- factor(metadata$type)

#for making different philr weights
philr_taxa_weights <- c("uniform","gm.counts","anorm","anorm.x.gm.counts","enorm","enorm.x.gm.counts")
philr_ilr_weights <- c("uniform","blw","blw.sqrt","mean.descendants")

##-Create plot data-------------------------------------------------##
print("Create datasets list")

nonphilr_ds <- list(asv_table, my_alr, my_clr, ln_asv_tab, 
                        as.data.frame(cln_denovo_tree_ps@otu_table), 
                        as.data.frame(ref_ps@otu_table), 
                        as.data.frame(ref_ps_clean@otu_table))
names_nonphilr_ds <- c("asv_table", "alr", "clr", "lognorm", "cln_upgma_count_table",
                           "ref_orig_count_table", "ref_clean_count_table")

# print("removing zeros from re_ps")
# clean_otu <- data.frame(ref_ps@otu_table@.Data)
# clean_otu <- clean_otu + 1
# ref_ps <- phyloseq::phyloseq( otu_table(clean_otu, taxa_are_rows = F), 
#                                     phy_tree(ref_ps@phy_tree),
#                                     tax_table(ref_ps@tax_table), 
#                                     sample_data(ref_ps@sam_data))

philr_ds <- list(cln_denovo_tree_ps, ref_ps, ref_ps_clean)
names_philr_ds <- c("cln_upgma", "ref_orig", "ref_clean")

print("Initializing empty vectors for permanova table")
perma_pval <- c()
perma_r2 <- c()
part_weight <- c()
ilr_weight <- c()
metadata_col <- c()
transformation <- c()
distance_metric <- c()
nonphilr_pg_num <- c()
philr_pg_num <- c()

index <- 1
pdf(file = file.path(output_dir, "graphics", paste0( main_output_label,"_", dist_metric, "_nonphilr.pdf")))
for (ds in 1:length(nonphilr_ds)) {
  my_table <- nonphilr_ds[[ds]]
  my_ds_name <- names_nonphilr_ds[ds]
  ##---------------------PCOA transform and plots---------------------##
  # my_dist <- stats::dist(my_table, method = dist_metric)
  my_pcoa <- phyloseq::ordinate(otu_table(my_table, taxa_are_rows = F), 
                                method = 'PCoA', distance = dist_metric)
  my_MDS <- data.frame(my_pcoa$vectors)[,1:2]

  # my_pcoa1 <- vegan::capscale(my_table~1,distance = dist_metric)
  # my_MDS1 <- data.frame(my_pcoa1$pCCA)
  
  # my_MDS1 <- vegan::eigenvals(my_pcoa1)
  
  # plot(my_MDS$MDS1, my_MDS$MDS2)
  
  for (mta in rf_cols){
    my_meta <- names(metadata)[mta]
    print(my_meta)
    
    perm <- vegan::adonis2(formula = my_MDS ~ metadata[,mta], method = dist_metric)
    
    perma_pval <- c(perma_pval, perm$`Pr(>F)`[1])
    perma_r2 <- c(perma_r2, perm$R2[1])
    nonphilr_pg_num <- c(nonphilr_pg_num, index)
    philr_pg_num <- c(philr_pg_num, 0)
    part_weight <- c(part_weight, 0)
    ilr_weight <- c(ilr_weight, 0)
    metadata_col <- c(metadata_col, my_meta)
    transformation <- c(transformation, my_ds_name)
    distance_metric <- c(distance_metric, dist_metric)
    
    index <- index + 1
    
    g <- ggplot(my_MDS, aes(x = Axis.1, y = Axis.2)) +
      geom_point(aes(col = metadata[,my_meta])) +
      ggplot2::ggtitle(label = paste(project, my_ds_name, my_meta, perm$`Pr(>F)`)) +
      ggplot2::theme_classic() +
      ggplot2::scale_x_discrete(guide = guide_axis(angle = 90))
    print(g)
  }
}  
dev.off()

pdf(file = file.path(output_dir, "graphics", paste0( main_output_label,"_", dist_metric, "_philr.pdf")))
index <- 1
for(tw in philr_taxa_weights){
  for(iw in philr_ilr_weights){
    for (ds in 1:length(philr_ds)) {
      ps_obj <- philr_ds[[ds]]
      my_table <- philr::philr(df = ps_obj@otu_table, tree = ps_obj@phy_tree, 
                               part.weights = tw,
                               ilr.weights = iw)
      my_ds_name <- names_philr_ds[ds]
      ##---------------------PCOA transform and plots---------------------##
      # my_dist <- stats::dist(my_table, method=dist_metric)
      my_pcoa <- phyloseq::ordinate(otu_table(my_table, taxa_are_rows = F), 
                                    method = 'PCoA', distance = dist_metric)
      my_MDS <- data.frame(my_pcoa$vectors)[,1:2]
      
      for (mta in rf_cols){
        my_meta <- names(metadata)[mta]
        print(my_meta)
        
        perm <- vegan::adonis2(formula = my_MDS ~ metadata[,mta], method = dist_metric)
        
        perma_pval <- c(perma_pval, perm$`Pr(>F)`[1])
        perma_r2 <- c(perma_r2, perm$R2[1])
        part_weight <- c(part_weight, tw)
        ilr_weight <- c(ilr_weight, iw)
        metadata_col <- c(metadata_col, my_meta)
        transformation <- c(transformation, my_ds_name)
        distance_metric <- c(distance_metric, dist_metric)
        philr_pg_num <- c(philr_pg_num, index)
        nonphilr_pg_num <- c(nonphilr_pg_num, 0)
        
        index <- index + 1
        
        g <- ggplot(my_MDS, aes(x = Axis.1, y = Axis.2)) +
          geom_point(aes(col = metadata[,my_meta])) +
          ggplot2::ggtitle(label = paste(project, my_ds_name, my_meta, "pw:", tw, "iw:", iw, perm$`Pr(>F)`)) +
          ggplot2::theme_classic() +
          ggplot2::scale_x_discrete(guide = guide_axis(angle = 90))
        print(g)
      }
    }  
  }
}
dev.off()

dFrame <- data.frame(perma_pval, perma_r2, part_weight, ilr_weight, 
                     metadata_col, transformation, distance_metric,
                     philr_pg_num, nonphilr_pg_num)
dFrame$adj_pval <- p.adjust(dFrame$perma_pval, method = "BH" )	
dFrame <- dFrame [order(dFrame$adj_pval),]

write.table(dFrame, file=file.path(output_dir, "tables", paste0(main_output_label, "_", dist_metric, ".csv")), 
            sep=",", 
            row.names=FALSE)

