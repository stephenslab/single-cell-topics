#!/bin/bash

#SBATCH --partition=broadwl
#SBATCH --account=pi-mstephens
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=32:00:00

# This script allocates computing resources (CPUs, memory), loads R,
# and runs fit_poisson_nmf.R. See fit_all_models.sh for examples
# illustrating how to use this script.

# Get the command-line arguments.
DAT=${1}
K=${2}
METHOD=${3}
NUMITER=${4}
if [[ ${5} == "yes" ]]
then
  EXTRAPOLATE="--extrapolate"
else
  EXTRAPOLATE=""
fi
OUTFILE=${6}
echo dat=${DAT}
echo k=${K}
echo method=${METHOD}
echo numiter=${NUMITER}
echo extrapolate=${5}
echo outfile=${OUTFILE}

DATFILE=../data/${DAT}.RData
PREFIT=${DAT//"_"/"-"}
PREFIT=../output/${PREFIT}/rds/prefit-${PREFIT}-k=${K}.rds

# For reproducibility, I add version numbers to the "module load" calls.
module load R/3.5.1

# Run the R script.
echo "Running fit_poisson_nmf.R"
export OPENBLAS_NUM_THREADS=1
./fit_poisson_nmf.R --counts ${DATFILE} --prefit ${PREFIT} \
  --method ${METHOD} -n ${NUMITER} ${EXTRAPOLATE} --nc 8 \
  -o ${OUTFILE}.rds > ${OUTFILE}.Rout
