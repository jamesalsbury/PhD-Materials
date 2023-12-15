#!/bin/bash
#SBATCH --mem=8G              # Request 4GB of memory
#SBATCH --mail-user=jsalsbury1@sheffield.ac.uk  # Replace with your email address
#SBATCH --ntasks-per-node=8

module load R/4.0.0-foss-2020a  # Recommended to load a specific version of R
module load JAGS/4.3.0-foss-2020a

R CMD BATCH simReproduceLambdacPiecewise.R simReproduceLambdacPiecewise.R.o$JOB_ID