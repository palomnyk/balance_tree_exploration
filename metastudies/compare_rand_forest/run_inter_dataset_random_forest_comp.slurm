#!/bin/bash

#SBATCH --partition=Pis                                    		# Partition/queue requested on server
#SBATCH --job-name=anova_taxa_pvals_all                                   		# Job name
#SBATCH --time=20:00:00                                      		# Time limit (hrs:min:sec)
#SBATCH --nodes=4                                         			# Number of nodes requested
#SBATCH --ntasks-per-node=1                          			# Number of CPUs (processor cores/tasks)
#SBATCH --mem=200gb                                          		# Memory limit
#SBATCH --mail-type=BEGIN,END,FAIL                              	# Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=amyerke@uncc.edu                            	# Specified email address
#SBATCH --output=/users/amyerke/slurmLogs/%x.%j.log     # Set directory for standard output
#SBATCH --error=/users/amyerke/slurmLogs/%x.%j.log      # Set directory for error log
#SBATCH --get-user-env

### Display the job context
echo Job: $SLURM_JOB_NAME with ID $SLURM_JOB_ID
echo "Submit Dir : $SLURM_SUBMIT_DIR"
echo Running on host: `hostname`
echo Using $SLURM_NTASKS processors across $SLURM_NNODES nodes

module load R

echo launching Vanderbilt

python ~/git/balance_tree_exploration/metastudies/compare_rand_forest/intra_dataset_random_forest_comp.py\
  -d ~/git/balance_tree_exploration

echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "======================================================"

