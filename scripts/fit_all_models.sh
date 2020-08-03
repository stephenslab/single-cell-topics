#!/bin/bash

# The shell commands below will submit Slurm jobs to perform the
# Poisson NMF model fitting for all single-cell RNA-seq data sets, and
# for different choices of the model parameters and optimization
# settings.
SCRIPT_PREFIT=prefit_poisson_nmf.sbatch
SCRIPT_FIT=fit_poisson_nmf_purified_pbmc.sbatch

# "Pre-fit" the models.
#                       data                 k    n outfile
sbatch ${SCRIPT_PREFIT} pbmc_purified.RData 13   40 prefit-pbmc-purified-k=13
sbatch ${SCRIPT_PREFIT} pbmc_68k.RData      13   40 prefit-pbmc-68k-k=13
sbatch ${SCRIPT_PREFIT} droplet.RData       13   40 prefit-droplet-k=13
sbatch ${SCRIPT_PREFIT} pulseseq.RData      13   40 prefit-pulseseq-k=13

# Fit rank-2 factorizations, with and without extrapolation.
#                    k method  ex outfile
sbatch ${SCRIPT_FIT} 2 em      no fit-em-k=2
sbatch ${SCRIPT_FIT} 2 ccd     no fit-ccd-k=2
sbatch ${SCRIPT_FIT} 2 scd     no fit-scd-k=2
sbatch ${SCRIPT_FIT} 2 em     yes fit-em-ex-k=2
sbatch ${SCRIPT_FIT} 2 ccd    yes fit-ccd-ex-k=2
sbatch ${SCRIPT_FIT} 2 scd    yes fit-scd-ex-k=2

# Fit rank-3 factorizations, with and without extrapolation.
#                    k method  ex outfile
sbatch ${SCRIPT_FIT} 3 em      no fit-em-k=3
sbatch ${SCRIPT_FIT} 3 ccd     no fit-ccd-k=3
sbatch ${SCRIPT_FIT} 3 scd     no fit-scd-k=3
sbatch ${SCRIPT_FIT} 3 em     yes fit-em-ex-k=3
sbatch ${SCRIPT_FIT} 3 ccd    yes fit-ccd-ex-k=3
sbatch ${SCRIPT_FIT} 3 scd    yes fit-scd-ex-k=3

# Fit rank-4 factorizations, with and without extrapolation.
#                    k method  ex outfile
sbatch ${SCRIPT_FIT} 4 em      no fit-em-k=4
sbatch ${SCRIPT_FIT} 4 ccd     no fit-ccd-k=4
sbatch ${SCRIPT_FIT} 4 scd     no fit-scd-k=4
sbatch ${SCRIPT_FIT} 4 em     yes fit-em-ex-k=4
sbatch ${SCRIPT_FIT} 4 ccd    yes fit-ccd-ex-k=4
sbatch ${SCRIPT_FIT} 4 scd    yes fit-scd-ex-k=4

# Fit rank-5 factorizations, with and without extrapolation.
#                    k method  ex outfile
sbatch ${SCRIPT_FIT} 5 em      no fit-em-k=5
sbatch ${SCRIPT_FIT} 5 ccd     no fit-ccd-k=5
sbatch ${SCRIPT_FIT} 5 scd     no fit-scd-k=5
sbatch ${SCRIPT_FIT} 5 em     yes fit-em-ex-k=5
sbatch ${SCRIPT_FIT} 5 ccd    yes fit-ccd-ex-k=5
sbatch ${SCRIPT_FIT} 5 scd    yes fit-scd-ex-k=5

# Fit rank-6 factorizations, with and without extrapolation.
#                    k method  ex outfile
sbatch ${SCRIPT_FIT} 6 em      no fit-em-k=6
sbatch ${SCRIPT_FIT} 6 ccd     no fit-ccd-k=6
sbatch ${SCRIPT_FIT} 6 scd     no fit-scd-k=6
sbatch ${SCRIPT_FIT} 6 em     yes fit-em-ex-k=6
sbatch ${SCRIPT_FIT} 6 ccd    yes fit-ccd-ex-k=6
sbatch ${SCRIPT_FIT} 6 scd    yes fit-scd-ex-k=6

# Fit rank-7 factorizations, with and without extrapolation.
#                    k method  ex outfile
sbatch ${SCRIPT_FIT} 7 em      no fit-em-k=7
sbatch ${SCRIPT_FIT} 7 ccd     no fit-ccd-k=7
sbatch ${SCRIPT_FIT} 7 scd     no fit-scd-k=7
sbatch ${SCRIPT_FIT} 7 em     yes fit-em-ex-k=7
sbatch ${SCRIPT_FIT} 7 ccd    yes fit-ccd-ex-k=7
sbatch ${SCRIPT_FIT} 7 scd    yes fit-scd-ex-k=7

# Fit rank-8 factorizations, with and without extrapolation.
#                    k method  ex outfile
sbatch ${SCRIPT_FIT} 8 em      no fit-em-k=8
sbatch ${SCRIPT_FIT} 8 ccd     no fit-ccd-k=8
sbatch ${SCRIPT_FIT} 8 scd     no fit-scd-k=8
sbatch ${SCRIPT_FIT} 8 em     yes fit-em-ex-k=8
sbatch ${SCRIPT_FIT} 8 ccd    yes fit-ccd-ex-k=8
sbatch ${SCRIPT_FIT} 8 scd    yes fit-scd-ex-k=8

# Fit rank-9 factorizations, with and without extrapolation.
#                    k method  ex outfile
sbatch ${SCRIPT_FIT} 9 em      no fit-em-k=9
sbatch ${SCRIPT_FIT} 9 ccd     no fit-ccd-k=9
sbatch ${SCRIPT_FIT} 9 scd     no fit-scd-k=9
sbatch ${SCRIPT_FIT} 9 em     yes fit-em-ex-k=9
sbatch ${SCRIPT_FIT} 9 ccd    yes fit-ccd-ex-k=9
sbatch ${SCRIPT_FIT} 9 scd    yes fit-scd-ex-k=9

# Fit rank-10 factorizations, with and without extrapolation.
#                     k method  ex outfile
sbatch ${SCRIPT_FIT} 10 em      no fit-em-k=10
sbatch ${SCRIPT_FIT} 10 ccd     no fit-ccd-k=10
sbatch ${SCRIPT_FIT} 10 scd     no fit-scd-k=10
sbatch ${SCRIPT_FIT} 10 em     yes fit-em-ex-k=10
sbatch ${SCRIPT_FIT} 10 ccd    yes fit-ccd-ex-k=10
sbatch ${SCRIPT_FIT} 10 scd    yes fit-scd-ex-k=10

# Fit rank-11 factorizations, with and without extrapolation.
#                     k method  ex outfile
sbatch ${SCRIPT_FIT} 11 em      no fit-em-k=11
sbatch ${SCRIPT_FIT} 11 ccd     no fit-ccd-k=11
sbatch ${SCRIPT_FIT} 11 scd     no fit-scd-k=11
sbatch ${SCRIPT_FIT} 11 em     yes fit-em-ex-k=11
sbatch ${SCRIPT_FIT} 11 ccd    yes fit-ccd-ex-k=11
sbatch ${SCRIPT_FIT} 11 scd    yes fit-scd-ex-k=11

# Fit rank-12 factorizations, with and without extrapolation.
#                     k method  ex outfile
sbatch ${SCRIPT_FIT} 12 em      no fit-em-k=12
sbatch ${SCRIPT_FIT} 12 ccd     no fit-ccd-k=12
sbatch ${SCRIPT_FIT} 12 scd     no fit-scd-k=12
sbatch ${SCRIPT_FIT} 12 em     yes fit-em-ex-k=12
sbatch ${SCRIPT_FIT} 12 ccd    yes fit-ccd-ex-k=12
sbatch ${SCRIPT_FIT} 12 scd    yes fit-scd-ex-k=12
