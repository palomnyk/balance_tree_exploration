#!/bin/bash

#SBATCH --partition=Orion                                    		# Partition/queue requested on server
#SBATCH --job-name=rand_tree_vs_philr_auc                                   		# Job name
#SBATCH --time=50:00:00                                      		# Time limit (hrs:min:sec)
#SBATCH --nodes=1                                         			# Number of nodes requested
#SBATCH --ntasks-per-node=1                          			# Number of CPUs (processor cores/tasks)
#SBATCH --mem=50gb                                          		# Memory limit
#SBATCH --mail-type=BEGIN,END,FAIL                              	# Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=amyerke@uncc.edu                            	# Specified email address
#SBATCH --output=/users/amyerke/slurmLogs/vand_%x.%j.out     # Set directory for standard output
#SBATCH --error=/users/amyerke/slurmLogs/vand_%x.%j.err      # Set directory for error log
##SBATCH --uid=amyerke
#SBATCH --get-user-env

### Display the job context
echo Job: $SLURM_JOB_NAME with ID $SLURM_JOB_ID
echo Running on host: `hostname`
echo Using $SLURM_NTASKS processors across $SLURM_NNODES nodes

module load R

srun Rscript ~/git/balance_tree_exploration/module load R

srun Rscript ~/git/balance_tree_exploration/lib/cml_scripts/random_forest/random_forest_philr_comparison.R \
  -d ~/git/balance_tree_exploration \
  -p Vanderbilt \
  -m ~/git/balance_tree_exploration/Zeller/patient_metadata.tsv \
  -r Run \
  -n 20

echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "======================================================"

