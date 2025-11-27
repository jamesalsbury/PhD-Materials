#!/bin/bash
#SBATCH --job-name=dte_array
#SBATCH --mem=4G
#SBATCH --array=1-1
#SBATCH --output=array_%A_%a.out
#SBATCH --error=array_%A_%a.err

# Load modules (R first, then rjags)
module load rjags/4-13-foss-2022a-R-4.2.1

# Number of simulations per batch
NSIMS=100000

# Array index becomes the batch seed
SEED=$SLURM_ARRAY_TASK_ID

echo "Running batch $SEED with $NSIMS sims"

Rscript Scenario3_Design1.R $NSIMS $SEED