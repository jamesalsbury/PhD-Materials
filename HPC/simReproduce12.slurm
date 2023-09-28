#!/bin/bash
#SBATCH --mem=8G              # Request 4GB of memory
#SBATCH --mail-user=jsalsbury1@sheffield.ac.uk  # Replace with your email address
#SBATCH --ntasks-per-node=32

module load R/4.0.5-foss-2020b  # Recommended to load a specific version of R

R CMD BATCH simReproduce12.R simReproduce12.R.o$JOB_ID