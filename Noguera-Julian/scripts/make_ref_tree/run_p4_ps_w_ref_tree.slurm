#!/bin/bash

#SBATCH --partition=Orion                                    		# Partition/queue requested on server
#SBATCH --job-name=run_p4_ps_w_ref_tree.R                                   		# Job name
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

srun Rscript ~/git/balance_tree_exploration/r_libraries/cml_scripts/make_ref_tree/p4_ps_w_ref_tree.R -d ~/git/balance_tree_exploration \
  -p Noguera-Julian \
  -m ~/git/balance_tree_exploration/Noguera-Julian/SraRunTable.txt \
  -l , \
  -r Run

echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "======================================================"
