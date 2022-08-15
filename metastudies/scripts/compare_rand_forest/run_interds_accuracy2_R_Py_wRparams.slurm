#!/bin/bash

#SBATCH --partition=Pisces                                    		# Partition/queue requested on server
#SBATCH --job-name=inter_dataset_rand                                   		# Job name
#SBATCH --time=2:00:00                                      		# Time limit (hrs:min:sec)
#SBATCH --nodes=1                                         			# Number of nodes requested
#SBATCH --ntasks-per-node=1                          			# Number of CPUs (processor cores/tasks)
#SBATCH --mem=50gb                                          		# Memory limit
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

module load anaconda3/2020.11
conda activate testing_pyCaret
python ~/git/balance_tree_exploration/lib/cml_scripts/metastudy_comparisons/inter_ds_accuracy_v_accuracy.py\
  -d ~/git/balance_tree_exploration \
	-x wide_random_forest_score_R_20.csv \
	-y sklrn_rand_forst_manual_0.75train_ntree500_m_try0.333.csv \
	-v R_randomForest \
	-w Sklearn_ntree500_mtry0.333 \
	-a trans_group \
	-b dataset \
	-p Vanderbilt Zeller

conda deactivate
module unload anaconda3/2020.11

echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "======================================================"

