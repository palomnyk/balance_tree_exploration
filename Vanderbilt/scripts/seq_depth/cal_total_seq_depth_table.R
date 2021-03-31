# Author: Aaron Yerke (aaronyerke@gmail.com)
# Script for making table with total sequencing depth

asv_table <- read.table(file.path(output_dir, "tables", "ForwardReads_DADA2.txt"),
                        sep = "\t",
                        header = TRUE)

total_seqs <- rowSums(asv_table)

total_seqs <- data.frame(total_seqs, row.names = row.names(asv_table))
write.table(total_seqs,
            file = file.path(output_dir, "tables", "total_seq_depth.csv"),
            sep = ",")