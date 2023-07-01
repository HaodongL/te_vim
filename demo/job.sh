#!/bin/bash
#SBATCH --account=co_biostat
#SBATCH --job-name=tevim_demo
#SBATCH --partition=savio2
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --time=72:00:00

#SBATCH --mail-user=haodong_li@berkeley.edu
#SBATCH --mail-type=ALL
module load r
Rscript demo_run.R
