#Script from creating IQTree on cluster - requires CLUSTALW and IQTREE modules

#Getting fancy with color in print statements
RED='\033[0;31m' #for print red
NC='\033[0m' # No Color

home_dir=$1 #first comandline argument is the project name
project=$2

echo "Found arguments ${home_dir} and ${project}."

#my_fasta is made from 
aligned=$home_dir/$project/output/trees/ForwardReads_DADA2_taxonomy.aln

#check for fastq file that is our starting point
if [ -f "$aligned" ]; then
    echo "$aligned exists."
else 
    printf "${RED}$aligned does not exist,\nrun your sequences through Dada2 and create the fastq with p1 and p2${NC}"
fi

cd $home_dir/$project/output/trees/

module load iqtree

iqtree2 -s $aligned -T AUTO

printf "${RED}Reached end of script.${NC}"

