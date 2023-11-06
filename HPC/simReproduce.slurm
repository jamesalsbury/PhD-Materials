#!/bin/bash
#SBATCH --mem=8G              # Request 4GB of memory
#SBATCH --mail-user=jsalsbury1@sheffield.ac.uk  # Replace with your email address
#SBATCH --ntasks-per-node=30

module load rjags/4-13-foss-2022a-R-4.2.1

R CMD BATCH simReproduce.R simReproduce.R.o$JOB_ID