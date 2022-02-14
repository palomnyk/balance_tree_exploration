#!/bin/bash

#SBATCH --partition=Orion                                    		# Partition/queue requested on server
#SBATCH --job-name=model_metrics                                   		# Job name
#SBATCH --time=0:10:00                                      		# Time limit (hrs:min:sec)
#SBATCH --nodes=1                                         			# Number of nodes requested
#SBATCH --ntasks-per-node=1                          			# Number of CPUs (processor cores/tasks)
#SBATCH --mem=10gb                                          		# Memory limit
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

module load anaconda3/2020.11
conda activate testing_pyCaret
srun python ~/git/balance_tree_exploration/lib/cml_scripts/mlm_metrics/sklearn_model_metrics.py\
  -m 6,7\
  -p Noguera-Julian\
  -a True

conda deactivate
module unload anaconda3/2020.11

# parser.add_option("-m", "--metadata_cols",
#                   action="store_false", dest="meta_col", default="int",
#                   help="don't print status messages to stdout")
# parser.add_option("-d", "--homedir",
#                   default=os.path.expanduser(os.path.join("~", "git", "balance_tree_exploration")), 
#                   help="dataset dir path", dest="homedir", metavar="homedir")
# parser.add_option("-p", "--project", default="string",
#                   help="project folder", metavar="project")
# parser.add_option("-a", "--use_all_meta", default=False,
#                   help="use all metadata", metavar="use_all_meta")

echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "======================================================"

