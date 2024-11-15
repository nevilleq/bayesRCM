#!/bin/bash -l
#SBATCH -A zhangl4
#SBATCH -p agsmall
#SBATCH --time=0:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=50
#SBATCH --mem-per-cpu=1g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nevil066@umn.edu
#SBATCH -e run_bayes_unif2_error
#SBATCH -o run_bayes_unif2_out
#SBATCH --array=1-10


#Get to the right directory
#cd /panfs/jay/groups/1/zhangl4/nevil066/bayes/bayesRCM/sim/
cd ~/dissertation/bayesRCM/

#Load R
module load R/4.3.0-openblas-rocky8

#Do the dang thing
SETTING=2
# Print the variables to verify they are set
echo "SLURM_ARRAY_TASK_ID: $SLURM_ARRAY_TASK_ID"
echo "SETTING: $SETTING"

#FSL/Manual ROI
Rscript ./sim/run_sim_bayes_unif.R $SLURM_ARRAY_TASK_ID $SETTING
