#!/bin/bash

#SBATCH --partition=broadwl
#SBATCH --account=pi-mstephens
#SBATCH --cpus-per-task=10
#SBATCH --mem=12G
#SBATCH --time=18:00:00

# This script allocates computing resources (CPUs, memory), loads R,
# and runs prefit_poisson_nmf_purified_pbmc.R. See fit_all_models.sh for
# examples illustrating how to use this script.

# Get the command-line arguments.
K=${1}
OUTFILE=${2}

# For reproducibility, I add version numbers to the "module load" calls.
module load R/3.5.1

# Run the R script.
echo "Running prefit_poisson_nmf_purified_pbmc.R"
export OPENBLAS_NUM_THREADS=1
./prefit_poisson_nmf_purified_pbmc.R -k ${K} --nc 10 \
  -o ${OUTFILE}.rds > ${OUTFILE}.Rout
