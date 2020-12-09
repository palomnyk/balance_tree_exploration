module load blast/2.9.0+

cd ~/git/balance_tree_exploration/ava_c/output/tree_process_blast

blastn -query dada2seqs.fasta \
  -db ~/git/balance_tree_exploration/taxonomy/silva/db/tree \
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