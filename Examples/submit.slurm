#!/bin/bash
#SBATCH --job-name=dte_array
#SBATCH --time=02:00:00
#SBATCH --mem=4G
#SBATCH --array=1-5
#SBATCH --output=array_%A_%a.out
#SBATCH --error=array_%A_%a.err

module load rjags/4-13-foss-2022a-R-4.2.1

# Number of simulations per batch
NSIMS=500

# Array index becomes the batch seed
SEED=$SLURM_ARRAY_TASK_ID

echo "Running batch $SEED with $NSIMS sims"

Rscript run_assurance_batch.R $NSIMS $SEED
