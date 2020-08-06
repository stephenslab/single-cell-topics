#!/bin/bash

#SBATCH --partition=mstephens
#SBATCH --account=pi-mstephens
#SBATCH --cpus-per-task=10
#SBATCH --mem=12G
#SBATCH --time=22:00:00

# This script allocates computing resources (CPUs, memory), loads R,
# and runs fit_poisson_nmf_purified_pbmc.R. See fit_all_models.sh for
# examples illustrating how to use this script.

# Get the command-line arguments.
K=${1}
METHOD=${2}
if [[ ${3} == 'yes' ]]
then
  EXTRAPOLATE="--extrapolate"
else
  EXTRAPOLATE=""
fi
OUTFILE=${4}
if [[ $# -lt 5 ]]
then
  NUMITER=""
else
  NUMITER="--numiter=${5}"
fi
if [[ $# -lt 6 ]]
then
  PREFIT=""
else
  PREFIT="--prefit=${6}"
fi

# For reproducibility, I add version numbers to the "module load" calls.
module load R/3.5.1

# Run the R script.
echo "Running fit_poisson_nmf_purified_pbmc.R"
export OPENBLAS_NUM_THREADS=1
./fit_poisson_nmf_purified_pbmc.R -k ${K} --method ${METHOD} --nc 10 \
  ${EXTRAPOLATE} ${NUMITER} ${PREFIT} -o ${OUTFILE}.rds > ${OUTFILE}.Rout
