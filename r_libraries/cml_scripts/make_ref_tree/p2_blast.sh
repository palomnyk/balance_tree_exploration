module load blast/2.9.0+
#Author: Aaron Yerke (aaronyerke@gmail.com)
home_dir=$1 #first comandline argument is the project name
project=$2 

db_path=${home_dir}/ref_tree_objs/db
cd ~/git/balance_tree_exploration/${project}/output/tree_process_blast

blastn -query dada2seqs.fasta \
  -db ${db_path} \
  -out output.txt \
  -outfmt "6 qseqid sseqid pident length evalue bitscore score ppos"

 # 1.	 qseqid	 query (e.g., unknown gene) sequence id
 # 2.	 sseqid	 subject (e.g., reference genome) sequence id
 # 3.	 pident	 percentage of identical matches
 # 4.	 length	 alignment length (sequence overlap)
 # 5.	 evalue	 expect value
 # 6.	 bitscore	 bit score
 # 7.  score     Raw score
 # 8.  ppos      Percentage of positive-scoring matches