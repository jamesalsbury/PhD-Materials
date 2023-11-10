#!/bin/bash
#SBATCH --mem=32G              # Request 4GB of memory
#SBATCH --mail-user=smp21js@sheffield.ac.uk
#SBATCH --ntasks-per-node=64
#SBATCH --time=15:00:00

module load rjags/4-13-foss-2022a-R-4.2.1

R CMD BATCH simReproduce.R simReproduce.R.o$JOB_ID