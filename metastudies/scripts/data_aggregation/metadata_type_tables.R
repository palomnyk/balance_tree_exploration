# Author: Aaron Yerke (aaronyerke@gmail.com)
# Script for making table of metadata characterics

rm(list = ls()) #clear workspace

##-Load Depencencies------------------------------------------------##
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
print("finished loading libraries")

##-Establish directory layout---------------------------------------##
home_dir <- file.path('~','git','balance_tree_exploration')
projects <- c("Jones", "Zeller","Vangay", "Noguera-Julian")
proj_metaname <- c("patient_metadata.tsv", "patient_metadata.csv", "patient_metadata.tsv", "patient_metadata.tsv")
proj_meta_delim <- c("\t",",","\t","\t")
proj_meta_colname <- c("Run", "Run", "run_accession", "Run")
final_output_dir <- file.path(home_dir, "metastudies", "output")
##-Functions--------------------------------------------------------##

##-Create tree attribute vectors------------------------------------##
binary_meta <- c()
multi_meta <- c()
continous_meta <- c()

for (proj in 1:length(projects)) {
	project <- projects[proj]
	print(project)
	metadata <- data.frame(data.table::fread(file = file.path(home_dir, project, proj_metaname[proj]),
																					header=TRUE, data.table=FALSE), row.names = proj_meta_colname[proj])
	# meta_str <- str(metadata)
	# print(meta_str)[1]
	meta_class <- sapply(metadata,class)
	# print(meta_class[1])
	for (colmn in 1:ncol(metadata)){
		if (meta_class[colmn] == "integer" || meta_class[colmn] == "numeric"){
			continous_meta <- c(continous_meta, colnames(metadata)[colmn])
		}else{
		  # print(unique(metadata)[,colmn])
			if (length(unique(metadata[,colmn], incomparables = c("",NA, NaN))) > 2){
				multi_meta <- c(multi_meta, colnames(metadata)[colmn])
			}else{
				binary_meta <- c(binary_meta, colnames(metadata)[colmn])
			}
		}
	}
}
print(paste("Bi:",length(binary_meta)))
print(paste("Multi:",length(multi_meta)))
print(paste("Cont:",length(continous_meta)))