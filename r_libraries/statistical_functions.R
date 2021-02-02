# Author: Aaron Yerke
# Library of statistical functions

##-Load Dependencies------------------------------------------------##

##_Establish constants----------------------------------------------##
# home_dir = file.path('~','git','Western_gut')
# #home_dir = file.path('cloud','project')
# output_dir = file.path(home_dir, 'output')
# # setwd(file.path(home_dir))

##_Functions--------------------------------------------------------##
philr_node_annot_anova_pval <- function(philr_df, phylo_obj, meta_df, only_one = TRUE) {
  # gives philr nodes phylogenetic taxonomy and compares the nodes to metadata
  # for ANOVA pval
  # philr_df: a philr transformed df
  # phylo_obj: a phyloseq object from which we will get the taxonmy and the phylogenetic tree for name.balance
  # meta_df: metadata df that you want to iterate over for pvalues
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
  le = ncol(philr_df) * ncol(meta_df)#number of rows in df
  pValues <-vector(length = le)
  philrColumns <- vector(length = le)
  metaNames <- vector(length = le)
  metaIndex <- vector(length = le)
  philrIndex <- vector(length = le)
  
  #taxonomy
  philrAnnot = vector(length = le)
  phylum = vector(length = le)
  class = vector(length = le)
  order = vector(length = le)
  family = vector(length = le)
  genus = vector(length = le)
  
  index <- 1
  
  for( m in 1:ncol(philr_df))
  {
    votes <- name.balance(phylo_obj@phy_tree, phylo_obj@tax_table, names(philr_df)[m], return.votes = c('up', 'down'))
    
    for( p in 1:ncol(meta_df))
    {
      aLm <- lm( philr_df[,m] ~  meta_df[,p])
      pValues[index]  <- 1
      
      try( pValues[index] <- anova(aLm)$"Pr(>F)"[1])
      metaNames[index] <- names(meta_df)[p]
      philrColumns[index] <- names(philr_df)[m]
      metaIndex[index] <- p
      philrIndex[index] <- m
      philrAnnot[index] = name.balance(phylo_obj@phy_tree, phylo_obj@tax_table, names(philr_df)[m])
      
      if (only_one == TRUE){
        genus[index] = paste0(onlyOne(votes$up.votes$Genus), "/", onlyOne(votes$down.votes$Genus))
        family[index] = paste0(onlyOne(votes$up.votes$Family), "/", onlyOne(votes$down.votes$Family))
        order[index] = paste0(onlyOne(votes$up.votes$Order), "/", onlyOne(votes$down.votes$Order))
        class[index] = paste0(onlyOne(votes$up.votes$Class), "/", onlyOne(votes$down.votes$Class))
        phylum[index] = paste0(onlyOne(votes$up.votes$Phylum), "/", onlyOne(votes$down.votes$Phylum))
      }else{
        genus[index] = paste0(names(which.max(votes$up.votes$Genus)), "/", names(which.max(votes$down.votes$Genus)))
        family[index] = paste0(names(which.max(votes$up.votes$Family)), "/", names(which.max(votes$down.votes$Family)))
        order[index] = paste0(names(which.max(votes$up.votes$Order)), "/", names(which.max(votes$down.votes$Order)))
        class[index] = paste0(names(which.max(votes$up.votes$Class)), "/", names(which.max(votes$down.votes$Class)))
        phylum[index] = paste0(names(which.max(votes$up.votes$Phylum)), "/", names(which.max(votes$down.votes$Phylum)))
      }
      index <- index + 1
    }
  }
  
  pValAdj <- p.adjust( pValues, method = "BH" )
  dFrame <- data.frame(pValues,pValAdj,philrColumns,philrIndex, philrAnnot,metaNames,metaIndex, genus, family, order, class, phylum)
  dFrame <- dFrame [order(dFrame$pValAdj),]
  
  detach(package:philr)
  
  return(dFrame)
}#end philr_node_annot_anova_pval

otu_anova_pval <- function(otu_df, meta_df){
  #gives pvalues and adjusted pvalues of otus and metadata
  le = ncol(otu_df) * ncol(meta_df)#number of rows in output df
  pValues <-vector(length = le)
  otuColumns <- vector(length = le)
  metaNames <- vector(length = le)
  metaIndex <- vector(length = le)
  otuIndex <- vector(length = le)
  
  index <- 1
  
  for( m in 1:ncol(otu_df))
  {
    for( p in 1:ncol(meta_df))
    {
      aLm <- lm( otu_df[,m] ~  meta_df[,p])
      
      try( pValues[index] <- anova(aLm)$"Pr(>F)"[1])
      metaNames[index] <- names(meta_df)[p]
      otuColumns[index] <- names(otu_df)[m]
      metaIndex[index] <- p
      otuIndex[index] <- m
      index <- index + 1
    }
  }
  pValAdj <- p.adjust( pValues, method = "BH" )
  dFrame <- data.frame(pValAdj,otuColumns,otuIndex,metaNames,metaIndex)
  dFrame <- dFrame [order(dFrame$pValAdj),]
  return(dFrame)
}

kendall_corr <- function(df1, df2){
  #gives pvalues and adjusted pvalues of otus and metadata
  le = ncol(df1) * ncol(df2)#number of rows in output df
  pValues <- vector(length = le)
  corr_stat <- vector(length = le)
  est <- vector(length = le)
  df1_names <- vector(length = le)
  df1_index <- vector(length = le)
  df2_names <- vector(length = le)
  df2_index <- vector(length = le)
  
  index <- 1
  
  for( m in 1:ncol(df1))
  {
    for( p in 1:ncol(df2))
    {
      a_corr <- cor.test( df1[,m], df2[,p],
                          alternative = "two.sided",
                          method = "kendall")
      
      pValues[index] <- a_corr$p.value
      corr_stat[index] <- a_corr$statistic
      est[index] <- a_corr$estimate
      df2_names[index] <- names(df2)[p]
      df1_names[index] <- names(df1)[m]
      df2_index[index] <- p
      df1_index[index] <- m
      index <- index + 1
    }
  }
  pValAdj <- p.adjust( pValues, method = "BH" )
  dFrame <- data.frame(pValAdj,pValues,corr_stat, est, df1_names,df1_index,df2_names,df2_index)
  dFrame <- dFrame [order(dFrame$pValAdj),]
  return(dFrame)
}

kendall_corr_plots <- function(df1, #i.e. philr_transform
                         df2, #i.e otu ratio
                         df1_annot, #philr_taxonomic annotation vector
                         annot_col, #which level
                         query_stat, #corr stat of interest (pval, estimate, ect)
                         stat_threshold, #what threshold
                         FUN = function(x,y){x > y},#for if statement
                         pdf_path #path for plot output
                         ){
  #gives related stats of kendall correlation
  #annotation table
  le = ncol(df1) * ncol(df2)#number of rows in output df
  pValues <- vector(length = le)
  corr_stat <- vector(length = le)
  est <- vector(length = le)
  df1_names <- vector(length = le)
  df1_index <- vector(length = le)
  df2_names <- vector(length = le)
  df2_index <- vector(length = le)
  df1_tax <- vector(length = le)
  taxon_match <- vector(length = le)
  index <- 1
  
  pdf(file = pdf_path)
  for( m in 1:ncol(df1))
  {
    my_node_taxon <- as.character(df1_annot[m])
    my_node_split <- strsplit(tolower(my_node_taxon), "/")[[1]]
    for( p in 1:ncol(df2))
    {
      a_corr <- cor.test( df1[,m], df2[,p],
                          alternative = "two.sided",
                          method = "kendall")
      my_otu <- names(df2)[p]
      
      #remove the _1 or _2 from some of the OTUs
      my_otu_gsub <- gsub("\\d", "", my_otu)
      my_otu_gsub <- gsub('_','',my_otu_gsub)
      
      pValues[index] <- a_corr$p.value
      corr_stat[index] <- a_corr$statistic
      est[index] <- a_corr$estimate
      df2_names[index] <- my_otu
      df1_names[index] <- names(df1)[m]
      df2_index[index] <- p
      df1_index[index] <- m
      df1_tax[index] <- my_node_taxon
      
      my_intersect <- unlist(intersect(my_node_split, strsplit(tolower(my_otu_gsub), "/")[[1]]))
      # print(my_intersect)
      if (length(my_intersect) > 0 ){
        taxon_match[index] <- TRUE
        a_title <- paste0("node: ", my_node_taxon, "/n", "OTUratio: ", names(df2[p]), "\n", query_stat, ": ", unlist(a_corr["estimate"]))
        plot(x = df1[,m],
             y = df2[,p],
             main = a_title,
             xlab = my_node_taxon,
             ylab = names(df2[p]),
        )
      }else {# end if..
        taxon_match[index] <- FALSE
      }
      index <- index + 1
    }#end for p
  }
  dev.off()
  pValAdj <- p.adjust( pValues, method = "BH" )
  # print(paste(length(taxon_match), length(df1_tax), length(pValAdj),  length(df2_index)))
  dFrame <- data.frame(pValAdj,pValues,corr_stat, est, df1_names,df1_index,df2_names,df2_index,df1_tax, taxon_match)
  dFrame <- dFrame [order(dFrame$pValAdj),]
  print("Completed Kendall Correlation with Plots")
  return(dFrame)
}#end kendall Corr Plot

roc_axes <- function(test_data, 
                    true_resp = resp_var_test,
                    ml_model = rf,
                    error_range = 0.10){
  # This function takes predictor variables (test_data) and tests it 
  # against the correct data (true_resp) and returns a dataframe of the 
  # true postives and true positives as ratios between 0 and 1
  # test_data is the testing variable
  # the true_resp is the correct variable
  # the error range provides an acceptable amount of error for when 
  # the numbers are close, but not exactly the same.
  groups = unique(true_resp)
  true_pos = c(0)
  false_pos = c(0)
  true_neg = c(0)
  # for categorical data
  for (grp in groups){
    if (class(true_resp) %in% c("double", "integer")){
      upper_lim_grp = grp + grp * error_range
      lower_lim_grp = grp - grp * error_range
    }
    for (rw in 1:nrow(test_data)){
      #for when data are not categoric
      if (class(true_resp) %in% c("double", "integer")){
        pred = predict(ml_model, newdata=test_data[rw,])
        upper_lim_tr = true_resp[rw] + true_resp[rw] * error_range
        lower_lim_tr = true_resp[rw] - true_resp[rw] * error_range
        decision = upper_lim_grp > pred & pred > lower_lim_grp
        truth = upper_lim_tr > pred & pred > lower_lim_tr
      }else{
        #for categoric data
        pred = predict(ml_model, newdata=test_data[rw,])
        decision = pred == grp
        truth = true_resp[rw] == grp
      }
      if (decision == truth & truth == T){
        true_pos = c(true_pos, tail(true_pos)+1)
        true_neg = c(true_neg, tail(true_neg))
        false_pos = c(false_pos, tail(false_pos))
      }
      else if (decision == truth & truth == F){
        true_neg = c(true_neg, tail(true_neg)+1)
        true_pos = c(true_pos, tail(true_pos))
        false_pos = c(false_pos, tail(false_pos))
      }
      else{
        true_neg = c(true_neg, tail(true_neg))
        true_pos = c(true_pos, tail(true_pos))
        false_pos = c(false_pos, tail(false_pos)+1)
      }
    }#for (rw in 1:nrow(test_data)){
  }#for (grp in groups){
  print(paste("max(false_pos):", max(false_pos), "max(true_pos):", max(true_pos), "any na:", any(is.na(true_pos))))
  
  if( max(true_neg) != 0){
    true_neg = true_neg/max(true_neg)
  }
  if( max(true_pos) != 0){
    true_pos = true_pos/max(true_pos)
  }
  if( max(false_pos) != 0){
    false_pos = false_pos/max(false_pos)
  }
  return(data.frame(true_pos, false_pos))
}#end roc_data

filt_seq_dpth <- function(min_seq_depth, df) {
  df <- df[rowSums(df) > min_seq_depth, ]
}




