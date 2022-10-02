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
  --homedir ~/git/balance_tree_exploration \
	--y_axis_fname wide_random_forest_score_R_20.csv \
	--x_axis_fname sklrn_randmforst_manual_0.75train_ntree500_m_try0.333.csv \
	--y_axis_label R_randomForest \
	--x_axis_label Sklearn_ntree500_mtry0.333 \
	--y_axis_col_name trans_group \
	--x_axis_col_name dataset \
	--y_axis_metadata_col_name metadata_col \
	--x_axis_metadata_col_name metadata \
	--project_list Jones Zeller Vangay Noguera-Julian

conda deactivate
module unload anaconda3/2020.11

echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "======================================================"

