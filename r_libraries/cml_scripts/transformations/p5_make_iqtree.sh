#Script from creating IQTree on cluster - requires CLUSTALW and IQTREE modules

#Getting fancy with color in print statements
RED='\033[0;31m' #for print red
NC='\033[0m' # No Color

home_dir=$1 #first comandline argument is the project name
project=$2

echo "Found arguments ${home_dir} and ${project}."

#my_fasta is made from 
my_fasta=$home_dir/$project/output/tree_process_blast/dada2seqs.fasta
clust_out=$home_dir/$project/output/tree_process_blast/aligned_dada2seqs.fasta

#check for fastq file that is our starting point
if [ -f "$my_fasta" ]; then
    echo "$my_fasta exists."
else 
    printf "${RED}$my_fasta does not exist,\nrun your sequences through Dada2 and create the fastq with p1 and p2${NC}"
fi

#align the sequences with clustal omega
module load clustal-omega
clustalo --in=$my_fasta --out=$clust_out



