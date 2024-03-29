#!/bin/bash

#SBATCH --partition=Orion                                    # Partition/queue requested on server
#SBATCH --job-name=Zeller_dada2                                   # Job name
#SBATCH --time=100:00:00                                      # Time limit (hrs:min:sec)
#SBATCH --nodes=1                                         # Number of nodes requested
#SBATCH --ntasks-per-node=1                          # Number of CPUs (processor cores/tasks)
#SBATCH --mem=30gb                                          # Memory limit
#SBATCH --mail-type=BEGIN,END,FAIL                        # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=amyerke@uncc.edu                     # Specified email address
#SBATCH --output=/users/amyerke/slurmLogs/%x.%j.out     # Set directory for standard output
#SBATCH --error=/users/amyerke/slurmLogs/%x.%j.err      # Set directory for error log
##SBATCH --uid=amyerke
#SBATCH --get-user-env

### Display the job context
echo Job: $SLURM_JOB_NAME with ID $SLURM_JOB_ID
echo Running on host: `hostname`
echo Using $SLURM_NTASKS processors across $SLURM_NNODES nodes

module load R

srun Rscript ~/git/balance_tree_exploration/lib/cml_scripts/data_preprocessing/p1_dada2_rd1.R \
-d ~/git/balance_tree_exploration \
-p Zeller \
-f _1.fastq \
-r _2.fastq \
-t 190

echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "======================================================"
