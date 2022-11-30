#!/bin/bash

#SBATCH --partition=Pisces                                    		# Partition/queue requested on server
#SBATCH --job-name=summary_accuracy_average                    # Job name
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
python ~/git/balance_tree_exploration/lib/cml_scripts/metastudy_comparisons/summary_boxplot_accuracy_average.py \
  --homedir ~/git/balance_tree_exploration \
  --input_file wide_random_forest_score_R_20.csv \
  --output_tag R_20reps \
  --transfomration_column_name trans_group \
  --comparison_list True

conda deactivate
module unload anaconda3/2020.11

echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "======================================================"

parser.add_argument("-n", "--input_file",
                  default="sklearn_random_forest_manual_75train.csv",
                  help="Should be one file with this name in the output/tables/ dir in each project", 
                  dest="input_file", metavar="input_file")
parser.add_argument("-t", "--output_tag",
                  default="",
                  help="tag to help keep track of output", 
                  dest="output_tag", metavar="output_tag")
parser.add_argument("-c", "--transfomration_column_name",
                  default="dataset",
                  help="column that has the transformation", 
                  dest="transfomration_column_name", metavar="tcn")