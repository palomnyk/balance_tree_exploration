#!/bin/bash

#SBATCH --partition=Orion                                    		# Partition/queue requested on server
#SBATCH --job-name=run_p2_denovo_tree                                   		# Job name
#SBATCH --time=24:00:00                                      		# Time limit (hrs:min:sec)
#SBATCH --nodes=1                                         			# Number of nodes requested
#SBATCH --ntasks-per-node=1                          			# Number of CPUs (processor cores/tasks)
#SBATCH --mem=50gb                                          		# Memory limit
#SBATCH --mail-type=BEGIN,END,FAIL                              	# Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=amyerke@uncc.edu                            	# Specified email address
#SBATCH --output=/users/amyerke/slurmLogs/%x.%j.out     # Set directory for standard output
#SBATCH --error=/users/amyerke/slurmLogs/%x.%j.err      # Set directory for error log
##SBATCH --uid=amyerke
#SBATCH --get-user-env

### Display the job context
echo Job: $SLURM_JOB_NAME with ID $SLURM_JOB_ID
echo Running on host: `hostname`
echo Using $SLURM_NTASKS processors across $SLURM_NNODES nodes

module load R

srun Rscript ~/git/balance_tree_exploration/r_libraries/cml_scripts/transformations/p2_denovo_tree_UPGMA.R \
  -d ~/git/balance_tree_exploration \
  -p Fodor \
  -m ~/git/balance_tree_exploration/Fodor/SraRunTable.txt \
  -l , \
  -r Run \
  -s _1000_sd_filtered

# Options:
#   -d CHARACTER, --homedir=CHARACTER
# dataset dir path
# 
# -p CHARACTER, --project=CHARACTER
# project folder
# 
# -m CHARACTER, --metadata=CHARACTER
# metadata file path with filename
# 
# -l CHARACTER, --metadata_delim=CHARACTER
# metadata file deliminator
# 
# -r CHARACTER, --metadata_rowname=CHARACTER
# metadata file row to use for row names
# 
# -s CHARACTER, --outputfilesuffix=CHARACTER
# output_file_suffix
# 
# -f NUMERIC, --filter_level=NUMERIC
# taxonimic level for making otu table 1-6
# 
# -h, --help
# Show this help message and exit

echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "======================================================"
