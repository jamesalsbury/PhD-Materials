#!/bin/bash
#SBATCH --mem=16G              # Request 4GB of memory
#SBATCH --mail-user=jsalsbury1@sheffield.ac.uk  # Replace with your email address
#SBATCH --ntasks-per-node=4

module load R/3.6.2-foss-2019b  # Recommended to load a specific version of R

R CMD BATCH simReproduce.R simReproduce.R.o$JOB_ID

