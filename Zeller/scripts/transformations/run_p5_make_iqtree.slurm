#!/bin/bash

#SBATCH --partition=Orion                                    		# Partition/queue requested on server
#SBATCH --job-name=Zeller_clustalo_iqtree                                   		# Job name
#SBATCH --time=48:00:00                                      		# Time limit (hrs:min:sec)
#SBATCH --nodes=6                                         			# Number of nodes requested
#SBATCH --ntasks-per-node=1                          			# Number of CPUs (processor cores/tasks)
#SBATCH --mem=50gb                                          		# Memory limit
#SBATCH --mail-type=BEGIN,END,FAIL                              	# Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=amyerke@uncc.edu                            	# Specified email address
#SBATCH --output=/users/amyerke/slurmLogs/%x.%j.out     # Set directory for standard output
#SBATCH --error=/users/amyerke/slurmLogs/%x.%j.err      # Set directory for error log
#SBATCH --get-user-env

### Display the job context
echo Job: $SLURM_JOB_NAME with ID $SLURM_JOB_ID
echo Running on host: `hostname`
echo Using $SLURM_NTASKS processors across $SLURM_NNODES nodes

srun bash ~/git/balance_tree_exploration/lib/cml_scripts/transformations/p5_make_iqtree.sh ~/git/balance_tree_exploration Zeller

echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "======================================================"

