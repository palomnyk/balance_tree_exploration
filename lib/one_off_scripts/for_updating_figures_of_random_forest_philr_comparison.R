# Author: Aaron Yerke (aaronyerke@gmail.com)
# Script for testing each transformation against different sequence depth to find the best one
# This was helepful: https://github.com/jsilve24/philr/blob/master/vignettes/philr-intro.Rmd#L142

rm(list = ls()) #clear workspace

##-functions--------------------------------------------------------##
munge_ref_ps <- function(ps){
  #ps must be a phyloseq object
  #requires phyloseq and ape packages to be loaded in the env
  # ps <- filter_taxa(ps, function(x) sum(x > 3) > (0.2*length(x)), TRUE)
  # ps <- filter_taxa(ps, function(x) sd(x)/mean(x) > 3.0, TRUE)
  ps <- transform_sample_counts(ps, function(x) x+1)
  phy_tree(ps) <- makeNodeLabel(phy_tree(ps), method="number", prefix='n')
  return(ps)
}

##-Load Depencencies------------------------------------------------##
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# if (!requireNamespace("ALDEx2", quietly = TRUE)) BiocManager::install("ALDEx2")
# library("ALDEx2")
# if (!requireNamespace("optparse", quietly = TRUE)){ install.packages("optparse") }
# library("optparse")
# if (!requireNamespace("rgr", quietly = TRUE)) install.packages("rgr")
# library("rgr")
# if (!requireNamespace("phyloseq", quietly = TRUE)) BiocManager::install("phyloseq")
# library("phyloseq")
# if (!requireNamespace("vegan", quietly = TRUE)) BiocManager::install("vegan")
# library("vegan")
# if (!requireNamespace("DESeq2", quietly = TRUE)) BiocManager::install("DESeq2")
# library("DESeq2")
# if (!requireNamespace("philr", quietly = TRUE)) BiocManager::install("philr")
# library("philr")
# if (!requireNamespace("ape", quietly = TRUE)) BiocManager::install("ape")
# library("ape")
if (!requireNamespace("ggplot2", quietly = TRUE)) BiocManager::install("ggplot2")
library("ggplot2")

##-Establish directory layout---------------------------------------##
home_dir <- file.path('~','git','balance_tree_exploration')
project <- "Vanderbilt"
output_dir <- file.path(home_dir, project, 'output')

##-Set up constants-------------------------------------------------##
num_cycles <- 20
if(num_cycles < 3) stop("num_cycles should be 3 or more")
main_output_text <- "random_forest_auc_R"
main_output_label <- "random_forest_auc_R_20.csv"
# philr_taxa_weights <- c("uniform","gm.counts","anorm","anorm.x.gm.counts","enorm","enorm.x.gm.counts")
# philr_ilr_weights <- c("uniform","blw","blw.sqrt","mean.descendants")
philr_taxa_weights <- c("enorm")
philr_ilr_weights <- c("blw")
random_seed <- 36

all_plot_data <- read.table(file = file.path(output_dir, "tables", main_output_label),
            sep = ",", header = TRUE)

best_seqs <- table(all_plot_data$rf_imp_seq)

write.table(best_seqs, row.names = FALSE,
            file = file.path(output_dir, "tables", paste0("best_seqs_",main_output_label, ".csv")))

# for (tg in all_plot_data$trans_group){
#   print(tg)
#   my_tg_data <- all_plot_data[all_plot_data$trans_group == tg,]
#   print(table(my_tg_data$rf_imp_seq))
# }

weight_table <- data.frame(tree_type = c(F),
                           metadata = c(F),
                           taxa_pval = c(F),
                           ilr_pval = c(F))
weight_counter <- 1

print("Make all the boxplots")
##-Make all the boxplots--------------------------------------------##
# print(paste("Making empty vectors to fill during plot building"))
# pval <- c()
# pw_name <- c()
# iw_name <- c()
# metadata_col <- c()
# transformation <- c()
# distance_metric <- c()
# pg_num <- c()
# 
# pdf(file = file.path(output_dir, "graphics", paste0("new_bp_", main_output_label, ".pdf")))
# index <- 1
# for (mta in 1:length(unique(all_plot_data$metadata_col))){
#   my_meta <- as.character(unique(all_plot_data$metadata_col)[mta])
#   message(my_meta)
#   
#   plot_data <- all_plot_data[all_plot_data$metadata_col == my_meta,]
#   
#   not_uniform_tw <- which(plot_data$taxa_weight != "uniform" )
#   not_uniform_iw <- which(plot_data$ilr_weight != "uniform" )
#   philr_ds <- unique(plot_data$trans_group[c(not_uniform_tw,not_uniform_iw)] ) #pulls out philr only data
#   non_philr_ds <- unique(plot_data$trans_group[ !(plot_data$trans_group %in% philr_ds)])
#   non_philr_ds_pd <- subset(plot_data, trans_group %in% non_philr_ds)
#   philr_ds_pd <- data.frame(subset(plot_data, trans_group %in% philr_ds))
#   rownames(philr_ds_pd) <- seq(length=nrow(philr_ds_pd))
#   # my_means <- c()
#   # for (tg in 1:length(unique(new_pd$trans_group))){
#   #   trans_g <- new_pd$trans_group[tg]
#   #   my_vals <- which(new_pd$trans_group == trans_g)
#   #   my_means <- c(my_means, mean(new_pd$all_auc[my_vals]))
#   #   names(my_means)[tg] <- trans_g
#   # }
#   new_pd <- rbind(non_philr_ds_pd, philr_ds_pd)
#   new_pd$trans_group <- factor(new_pd$trans_group, levels = c(non_philr_ds, philr_ds))
#   
#   
#   for(tw in unique(plot_data$taxa_weight)){
#     for(iw in unique(plot_data$ilr_weight)){
#       philr_pd_tw_iw <- subset(philr_ds_pd, taxa_weight == tw & ilr_weight == iw)
#       jitter_pd <- rbind(non_philr_ds_pd, philr_pd_tw_iw)
#       jitter_pd$trans_group <- factor(jitter_pd$trans_group, levels = c(non_philr_ds, philr_ds))
#       #need to show means from new_pd, but show jitter of tw and iw
#       #or could just show selected points but show overal mean for each tw and iw
#       back_ground_points <- subset(philr_ds_pd, taxa_weight != tw & ilr_weight != iw)
#       bg_jitter <- rbind(non_philr_ds_pd, back_ground_points)
#       bg_jitter$trans_group <- factor(bg_jitter$trans_group, levels = c(non_philr_ds, philr_ds))
#       g <- ggplot2::ggplot(new_pd, aes(y = all_auc, x = trans_group)) + 
#         ggplot2::geom_boxplot(data = jitter_pd, color = "blue", alpha = 0.5) +
#         ggplot2::geom_boxplot(data = bg_jitter, color = "red", alpha = 0.5) +
#         ggplot2::ggtitle(label = paste("Jones", my_meta, ", part_weight:", tw, ", ilr_weight:", iw)) +
#         # ggplot2::ggtitle( label = paste("num_tg:", length(unique(new_pd$trans_group)))) +
#         ggplot2::theme_classic() +
#         ggplot2::scale_x_discrete(guide = guide_axis(angle = 90)) +
#         ggplot2::ylab("AUC") +
#         ggplot2::xlab("Tree type")
#       print(g)
#       #build vectors for table
#       for (grp in unique(philr_pd_tw_iw$trans_group)){
#         # print(grp)
#         my_case <- philr_pd_tw_iw[philr_pd_tw_iw$trans_group == grp, ]$all_auc
#         my_control <- back_ground_points[back_ground_points$trans_group == grp, ]$all_auc
#         my_test <- t.test(my_case, my_control)
#         my_pval <- my_test$p.value
#         pval <- c(pval, my_pval)
#         pw_name <- c(pw_name, tw)
#         iw_name <- c(iw_name, iw)
#         metadata_col <- c(metadata_col, my_meta)
#         transformation <- c(transformation, grp)
#         pg_num <- c(pg_num, index)
#       }
#       index <- index + 1
#     }#end for iw
#   }#end for tw
# }
# 
# dev.off()
# 
# dFrame <- data.frame( pval, pw_name, iw_name, metadata_col, transformation, pg_num)
# dFrame$adj_pval <- p.adjust(dFrame$pval, method = "BH" )	
# dFrame <- dFrame [order(dFrame$adj_pval),]
# 
# write.table(dFrame, file=file.path(output_dir, "tables", paste0("new_bp_", main_output_label, ".csv")), 
#             sep=",", 
#             row.names=FALSE)
# 
# print(paste("completed"))

pdf(file = file.path(output_dir, "graphics", paste0("new_bp_NO_WEIGHT_whole_meta", main_output_label, ".pdf")))
for (mta in 1:length(unique(all_plot_data$metadata_col))){
  my_meta <- as.character(unique(all_plot_data$metadata_col)[mta])
  message(my_meta)
  #select plot data for each metadata cat
  plot_data <- all_plot_data[all_plot_data$metadata_col == my_meta,]

  # Organize the boxes so that count tables are first - these will not 
  # have weighting schemes, therefor select those which do not have 
  # "uniform"
  not_uniform_tw <- which(plot_data$taxa_weight != philr_taxa_weights[1] )
  not_uniform_iw <- which(plot_data$ilr_weight != philr_ilr_weights[1] )
  philr_ds <- unique(plot_data$trans_group[c(not_uniform_tw,not_uniform_iw)] ) #pulls out philr only data
  non_philr_ds <- unique(plot_data$trans_group[ !(plot_data$trans_group %in% philr_ds)])
  non_philr_ds_pd <- subset(plot_data, trans_group %in% non_philr_ds)
  philr_ds_pd <- data.frame(subset(plot_data, trans_group %in% philr_ds))#philr plot data
  rownames(philr_ds_pd) <- seq(length=nrow(philr_ds_pd))

  new_pd <- rbind(non_philr_ds_pd, philr_ds_pd)#reorganize
  new_pd$trans_group <- factor(new_pd$trans_group, levels = c(non_philr_ds, philr_ds))#put non-philr first
  meta_mean <- mean(new_pd$all_auc)
  #add number to the end of each boxplot label
  my_labels <- as.character(levels(new_pd$trans_group))
  # my_new_labels <- lapply()
  # print(my_labels)
  counter <- 0
  random_tree <- FALSE
  new_labels <- sapply(my_labels, function(x) {
    counter <<- counter + 1
    return(paste0(x, "_", counter))
  })
  
  # print(length(new_labels))
  
  print(meta_mean)
	g <- ggplot2::ggplot(new_pd, aes(y = all_auc, x = trans_group, group=trans_group)) + 
		ggplot2::geom_boxplot() +
		ggplot2::ggtitle(label = paste(project, my_meta)) +
		ggplot2::geom_hline(yintercept = meta_mean, color="red") +
	  # scale_fill_discrete(labels=new_labels) +
		# ggplot2::ggtitle( label = paste("num_tg:", length(unique(new_pd$trans_group)))) +
		ggplot2::theme_classic() +
		ggplot2::scale_x_discrete(guide = guide_axis(angle = 90), labels=new_labels) +
		ggplot2::ylab("AUC") +
		ggplot2::xlab("Tree type")
	print(g)
}

dev.off()

print("finished R script")
