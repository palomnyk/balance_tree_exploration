# Author: Aaron Yerke
# Script for ROC of a all the spreadsheets
# This was helepful: https://github.com/jsilve24/philr/blob/master/vignettes/philr-intro.Rmd#L142

rm(list = ls()) #clear workspace

##-------------------Load Depencencies------------------------------##
if (!requireNamespace("randomForest", quietly = TRUE)){
  install.packages("randomForest")
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer")
}
library(RColorBrewer)
#set color palette
palette( brewer.pal(7,"Accent") )
library(randomForest)
library(compositions)
##----------------Establish directory layout------------------------##
home_dir = file.path('~','git','Western_gut')
#home_dir = file.path('cloud','project')
output_dir = file.path(home_dir, 'output')
setwd(file.path(home_dir))

##---------------------------Functions------------------------------##
# This function takes predictor variables (test_data) and tests it 
# against the correct data (true_resp) and returns a dataframe of the 
# true positives and true positives as ratios between 0 and 1
roc_axes = function(groups = unique(resp_var_test),
                    test_data, 
                    true_resp = resp_var_test,
                    ml_model = rf,
                    error_range = 0.10){
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

##------------Import tables    and data preprocessing---------------##
otu_ratio = read.table(file.path(home_dir, "philr_pipelines", "tables", "otu_ratio_table.csv"),
                        sep = ",",
                        header = TRUE)

ref_tree_philr_trans = read.table(file.path(home_dir, "philr_pipelines", "tables", "ref_tree_ps_philr_transform.csv"),
                      sep = ",",
                      header = TRUE)

upgma_philr_trans <- read.table(file.path(home_dir, "philr_pipelines", "tables", "philr_denovo_tree_UPGMA.csv"),
             sep = ",",
             header = TRUE)

otu_tab = read.table(file.path(home_dir, "philr_pipelines", "tables", "otu_table.csv"),
                     sep = ",",
                     header = TRUE)

asv_table = read.table(file.path(home_dir, "philr_pipelines", "tables", "asv_table.csv"),
                       sep = ",",
                       header = TRUE)

my_alr <- as.data.frame(compositions::alr(as.matrix(asv_table)))

my_clr <- as.data.frame(compositions::clr(as.matrix(asv_table)))

metadata = read.table(file.path("philr_pipelines", "tables", "ps_sample_data.csv"), 
                    sep=",", 
                    header=TRUE)
drop = c("SampleID","Sample.Date", "Subject.ID","Old.Participant.ID","sample_accession",
         "secondary_sample_accession","experiment_accession",             
         "tax_id","scientific_name","instrument_model","library_layout",
         "experiment_title","sample_title", "study_title", "run_alias",
         "fastq_ftp" ,"fastq_galaxy","submitted_ftp"                    
         ,"submitted_galaxy","sra_ftp", "sra_galaxy", "study_accession",
         "Notes.Samples", "Sub.Study", "Notes.Participants", "Birth.Year", "Religion")
metadata[metadata==""] <- NA
for (c in 1:ncol(metadata)){
  if (any(is.na(metadata[,c]))){
    drop = c(drop, colnames(metadata)[c])
  }
}
drop = unique(drop)
metadata = metadata[ , !(names(metadata) %in% drop)]
# lapply(metadata, class)

my_datasets = list(asv_table, otu_tab, otu_ratio, my_alr, my_clr, ref_tree_philr_trans, upgma_philr_trans, cbind(otu_ratio, ref_tree_philr_trans), cbind(otu_tab, otu_ratio))
my_ds_names = c("ASV", "OTU", "OTU ratio", "alr", "clr", "ref tree philR", "UPGMA philr trans", "OTU + ref tree philR", "OTU + OTU ratio")

##-----------------Create training/testing sets---------------------##
set.seed(36)
train_index = sample(x = nrow(metadata), size = 0.75*nrow(metadata), replace=FALSE)
test_index = c(1:nrow(metadata))[!(1:nrow(metadata) %in% train_index)]

pdf(file = file.path(output_dir, "roc_all_metadata1.pdf"))
my_rocs = list()
for(mta in 1:ncol(metadata)){
  resp_var_test = metadata[,mta][test_index]
  resp_var_train = metadata[,mta][train_index]
  for( ds in 1:length(my_datasets)){
    my_table = as.data.frame(my_datasets[ds])
    my_table_train = my_table[train_index,]
    my_table_test = my_table[test_index,]
    
    rf <- randomForest(
      my_table_train, resp_var_train
    )
    my_roc = roc_axes(test_data = my_table_test,
                           true_resp = resp_var_test
    )
    print(paste(ds, mta))
    my_rocs[[ds]] = my_roc
  }#for ds
  par(bg = 'grey96')
  plot(true_pos ~ false_pos,
       data = as.data.frame(my_rocs[1]),
       type = "l",
       xlab = "False positives",
       ylab = "True positives",
       col = palette()[1],
       main = paste("Western Gut", colnames(metadata)[mta]),
       ylim = c(0,1))
  for(i in 2:length(my_ds_names)){
    lines(true_pos ~ false_pos,
          data = as.data.frame(my_rocs[i]),
          col = palette()[i],
          type = "l")
  }
  abline(
    a = 0,
    b = 1)
  legend('bottomright', 
         legend = my_ds_names, 
         col = palette(), 
         pch = 15,
         cex=0.5)
}#for mta
dev.off()
