#!/bin/bash

#SBATCH --partition=Pisces                                    		# Partition/queue requested on server
#SBATCH --job-name=rf_all                                   		# Job name
#SBATCH --time=30:00:00                                      		# Time limit (hrs:min:sec)
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

echo launching Jones

srun --nodes=1Rscript ~/git/balance_tree_exploration/lib/cml_scripts/random_forest/random_forest_philr_comparison.R \
  -d ~/git/balance_tree_exploration \
  -p Jones \
  -m ~/git/balance_tree_exploration/Jones/patient_metadata.tsv \
  -r Run \
  -n 20

echo launching Vangay

srun --nodes=1 Rscript ~/git/balance_tree_exploration/lib/cml_scripts/random_forest/random_forest_philr_comparison.R \
  -d ~/git/balance_tree_exploration \
  -p Vangay \
  -m ~/git/balance_tree_exploration/Vangay/patient_metadata.tsv \
  -r run_accession \
  -n 20

# echo launching Zeller

# srun --nodes=1 Rscript ~/git/balance_tree_exploration/lib/cml_scripts/random_forest/random_forest_philr_comparison.R \
#   -d ~/git/balance_tree_exploration \
#   -p Zeller \
#   -m ~/git/balance_tree_exploration/Zeller/patient_metadata.csv \
#   -l , \
#   -r Run \
#   -n 20

echo launching Noguera-Julian

srun --nodes=1 Rscript ~/git/balance_tree_exploration/lib/cml_scripts/random_forest/random_forest_philr_comparison.R \
  -d ~/git/balance_tree_exploration \
  -p Noguera-Julian \
  -m ~/git/balance_tree_exploration/Noguera-Julian/patient_metadata.tsv \
  -r Run \
  -n 20

echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "======================================================"

