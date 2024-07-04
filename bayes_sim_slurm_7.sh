#!/bin/bash -l
#SBATCH -A zhangl4
#SBATCH -p agsmall
#SBATCH --time=7:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=50
#SBATCH --mem-per-cpu=2500mb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nevil066@umn.edu
#SBATCH -e run_bayes_error_7
#SBATCH -o run_bayes_out_7


#Get to the right directory
cd /panfs/jay/groups/1/zhangl4/nevil066/bayes/bayesRCM/sim/

#Load R
module load R/4.3.0-openblas-rocky8

#Do the dang thing

#FSL/Manual ROI
Rscript run_sim_bayes_7.R
