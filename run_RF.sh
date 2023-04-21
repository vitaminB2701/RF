#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --job-name=run_RF
#SBATCH --output=run_RF.out
#SBATCH --mem-per-cpu=2000

module load MATLAB
matlab -nodisplay < run_RF.m