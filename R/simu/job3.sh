#!/bin/bash
#SBATCH --account=co_biostat
#SBATCH --job-name=tevim_test
#SBATCH --partition=savio2
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --time=36:00:00

#SBATCH --mail-user=haodong_li@berkeley.edu
#SBATCH --mail-type=ALL
module load r
Rscript simu_vte.R
