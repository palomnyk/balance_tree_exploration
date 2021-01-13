# Author: Aaron Yerke
# Library of table manipulating functions

##-Load Dependencies------------------------------------------------##

##_Establish constants----------------------------------------------##
# home_dir = file.path('~','git','Western_gut')
# #home_dir = file.path('cloud','project')
# output_dir = file.path(home_dir, 'output')
# # setwd(file.path(home_dir))

##_Functions--------------------------------------------------------##
equal_num_columns_top_abund <- function(dfs, percent_abund = 0.70) {
  # equalizes all dfs to least num col taking the columns with highest abundance
  # dfs: list of dfs
  # percent_abund: minimun percent of non-zero cells in column
  # INTERNAL FUNCTIONS:
  percent_abund_filt <- function(abund_percent, df) {
    # finds columns that are above X% populated by nonzeros
    # returns a bool vector same length as # cols
    or_abund = apply(df, 2, function(c){
      sum(c!=0) >= abund_percent*nrow(df)})
    return(df[,or_abund])
  }
  
  top_x_abund_filt <- function(top_x, df) {
    #returns a df with highest abundance top_x number of columns
    df_sums = apply(df, 2, sum)
    df_rank = rank(df_sums, ties.method = "random")
    return(df[,df_rank <= top_x])
  }
  # MAIN METHOD OF FUNCTION
  dfs_abund = lapply(dfs, percent_abund_filt, abund_percent = percent_abund)
  lowest_ncol = min(unlist(lapply(dfs, ncol)))
  dfs_top_x_abund = lapply(dfs_abund, top_x_abund_filt, top_x = lowest_ncol)

  return(dfs_top_x_abund)
}

philr_tutorial_normalization <- function(df) {
  # The philr tutorial https://bioconductor.org/packages/release/bioc/vignettes/philr/inst/doc/philr-intro.R uses this block for normalizations: 
  # ps <-  filter_taxa(ps, function(x) sum(x > 3) > (0.2*length(x)), TRUE)
  # ps <-  filter_taxa(ps, function(x) sd(x)/mean(x) > 3.0, TRUE)
  # ps <- transform_sample_counts(ps, function(x) x+1)
  df <- df[apply(df, 2, function(x) sum(x > 3) > (0.2*length(x)))]
  df <-  df[apply(df, 2, function(x) sd(x)/mean(x) > 3.0)]
  df <- df + 1
  return(df)
}
simplifiy_meta_western_gut <- function(meta_df) {
  drop = c("SampleID","Sample.Date", "Subject.ID","Old.Participant.ID","sample_accession",
           "secondary_sample_accession","experiment_accession",             
           "tax_id","scientific_name","instrument_model","library_layout",
           "experiment_title","sample_title", "study_title", "run_alias",
           "fastq_ftp" ,"fastq_galaxy","submitted_ftp"                    
           ,"submitted_galaxy","sra_ftp", "sra_galaxy", "study_accession",
           "Notes.Samples", "Sub.Study", "Notes.Participants", "Birth.Year")
  meta_df[meta_df==""] <- NA
  for (c in 1:ncol(meta_df)){
    if (any(is.na(meta_df[,c]))){
      drop = c(drop, colnames(meta_df)[c])
    }
  }
  drop = unique(drop)
  meta_df = meta_df[ , !(names(meta_df) %in% drop)]
  return(meta_df)
}

philr_node_annot <- function(philr_df, phylo_obj, only_one = TRUE) {
  # gives philr nodes phylogenetic taxonomy
  # philr_df: a philr transformed df
  # phylo_obj: a phyloseq object from which we will get the taxonmy and the phylogenetic tree for name.balance
  # FUNCTION DEPENDENCIES:
  if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
    BiocManager::install("philr")
  }
  library(philr); packageVersion("philr")
  # INTERNAL FUNCTIONS:
  onlyOne = function(namedVector){
    # internal function for finding annotations with only one option (meaning perfect seperation)
    if(length(namedVector) == 1){
      return(names(namedVector)[1])
    }else {
      return("---")
    }
  }
  # MAIN METHOD OF FUNCTION
  le = ncol(philr_df)
  philrAnnot = vector(length = le)
  phylum = vector(length = le)
  class = vector(length = le)
  order = vector(length = le)
  family = vector(length = le)
  genus = vector(length = le)
  
  for( m in 1:ncol(philr_df)){
    votes <- name.balance(phylo_obj@phy_tree, phylo_obj@tax_table, names(philr_df)[m], return.votes = c('up', 'down'))
    philrAnnot[m] = name.balance(ps@phy_tree, ps@tax_table, names(philr_df)[m])
    if (only_one == TRUE){
      genus[m] = paste0(onlyOne(votes$up.votes$Genus), "/", onlyOne(votes$down.votes$Genus))
      family[m] = paste0(onlyOne(votes$up.votes$Family), "/", onlyOne(votes$down.votes$Family))
      order[m] = paste0(onlyOne(votes$up.votes$Order), "/", onlyOne(votes$down.votes$Order))
      class[m] = paste0(onlyOne(votes$up.votes$Class), "/", onlyOne(votes$down.votes$Class))
      phylum[m] = paste0(onlyOne(votes$up.votes$Phylum), "/", onlyOne(votes$down.votes$Phylum))
    }else{
      genus[m] = paste0(names(which.max(votes$up.votes$Genus)), "/", names(which.max(votes$down.votes$Genus)))
      family[m] = paste0(names(which.max(votes$up.votes$Family)), "/", names(which.max(votes$down.votes$Family)))
      order[m] = paste0(names(which.max(votes$up.votes$Order)), "/", names(which.max(votes$down.votes$Order)))
      class[m] = paste0(names(which.max(votes$up.votes$Class)), "/", names(which.max(votes$down.votes$Class)))
      phylum[m] = paste0(names(which.max(votes$up.votes$Phylum)), "/", names(which.max(votes$down.votes$Phylum)))
    }
  }#end for( m in 1:ncol(philr_df)){
  dFrame <- data.frame(philr_node = 1:le, node_names = names(philr_df), philrAnnot, genus, family, order, class, phylum)
  detach(package:philr)
  return(dFrame)
}#end philr_node_annot

simple_ratio_transform <- function(otu_df) {
  # Function for making ratio table of OTUS for ml purposes
  # The resulting table will look like: otu1/otu2, otu2/otu1
  otuRatios = data.frame(row.names = row.names(otu_df))
  for ( rn in 1:nrow(otu_df)){
    ratios = vector(length = ncol(otu_df)^2 - ncol(otu_df))
    index = 1
    for (cn1 in 1:ncol(otu_df)){
      for (cn2 in 1:ncol(otu_df)){
        if (cn2 != cn1){
          ratios[index] = otu_df[rn, cn1] / otu_df[rn, cn2]
          names(ratios)[index] = paste0(colnames(otu_df)[cn1], "/", colnames(otu_df)[cn2])
          index = index + 1
        }# end if cn2 != cn1
      }#end cn2
    }#end cn1
    ratioDf = t(data.frame(ratios))
    colnames(ratioDf) = names(ratios)
    row.names(ratioDf) = row.names(otu_df)[rn]
    
    otuRatios = rbind(otuRatios, ratioDf)
    # otuRatios = rbind( otuRatios, ratios)
  }
  return(otuRatios)
}