#!/bin/bash

# The shell commands below will submit Slurm jobs to perform the
# Poisson NMF model fitting for all single-cell RNA-seq data sets, and
# for different choices of the model parameters and optimization
# settings.
SCRIPT_PREFIT=prefit_poisson_nmf.sbatch
SCRIPT_FIT=fit_poisson_nmf.sbatch

# "Pre-fit" factorizations to the droplet data.
#                       data           k    n outfile
sbatch ${SCRIPT_PREFIT} droplet.RData  2 1000 prefit-droplet-k=2
sbatch ${SCRIPT_PREFIT} droplet.RData  3 1000 prefit-droplet-k=3
sbatch ${SCRIPT_PREFIT} droplet.RData  4 1000 prefit-droplet-k=4
sbatch ${SCRIPT_PREFIT} droplet.RData  5 1000 prefit-droplet-k=5
sbatch ${SCRIPT_PREFIT} droplet.RData  6 1000 prefit-droplet-k=6
sbatch ${SCRIPT_PREFIT} droplet.RData  7 1000 prefit-droplet-k=7
sbatch ${SCRIPT_PREFIT} droplet.RData  8 1000 prefit-droplet-k=8
sbatch ${SCRIPT_PREFIT} droplet.RData  9 1000 prefit-droplet-k=9
sbatch ${SCRIPT_PREFIT} droplet.RData 10 1000 prefit-droplet-k=10
sbatch ${SCRIPT_PREFIT} droplet.RData 11 1000 prefit-droplet-k=11
sbatch ${SCRIPT_PREFIT} droplet.RData 12 1000 prefit-droplet-k=12
sbatch ${SCRIPT_PREFIT} droplet.RData 13 1000 prefit-droplet-k=13

# "Pre-fit" factorizations to the 68k PBMC data.
#                       data            k    n outfile
sbatch ${SCRIPT_PREFIT} pbmc_68k.RData  2 1000 prefit-pbmc-68k-k=2
sbatch ${SCRIPT_PREFIT} pbmc_68k.RData  3 1000 prefit-pbmc-68k-k=3
sbatch ${SCRIPT_PREFIT} pbmc_68k.RData  4 1000 prefit-pbmc-68k-k=4
sbatch ${SCRIPT_PREFIT} pbmc_68k.RData  5 1000 prefit-pbmc-68k-k=5
sbatch ${SCRIPT_PREFIT} pbmc_68k.RData  6 1000 prefit-pbmc-68k-k=6
sbatch ${SCRIPT_PREFIT} pbmc_68k.RData  7 1000 prefit-pbmc-68k-k=7
sbatch ${SCRIPT_PREFIT} pbmc_68k.RData  8 1000 prefit-pbmc-68k-k=8
sbatch ${SCRIPT_PREFIT} pbmc_68k.RData  9 1000 prefit-pbmc-68k-k=9
sbatch ${SCRIPT_PREFIT} pbmc_68k.RData 10 1000 prefit-pbmc-68k-k=10
sbatch ${SCRIPT_PREFIT} pbmc_68k.RData 11 1000 prefit-pbmc-68k-k=11
sbatch ${SCRIPT_PREFIT} pbmc_68k.RData 12 1000 prefit-pbmc-68k-k=12
sbatch ${SCRIPT_PREFIT} pbmc_68k.RData 13 1000 prefit-pbmc-68k-k=13

# "Pre-fit" factorizations to the purified PBMC data.
#                       data                 k    n outfile
sbatch ${SCRIPT_PREFIT} pbmc_purified.RData  2 1000 prefit-pbmc-purified-k=2
sbatch ${SCRIPT_PREFIT} pbmc_purified.RData  3 1000 prefit-pbmc-purified-k=3
sbatch ${SCRIPT_PREFIT} pbmc_purified.RData  4 1000 prefit-pbmc-purified-k=4
sbatch ${SCRIPT_PREFIT} pbmc_purified.RData  5 1000 prefit-pbmc-purified-k=5
sbatch ${SCRIPT_PREFIT} pbmc_purified.RData  6 1000 prefit-pbmc-purified-k=6
sbatch ${SCRIPT_PREFIT} pbmc_purified.RData  7 1000 prefit-pbmc-purified-k=7
sbatch ${SCRIPT_PREFIT} pbmc_purified.RData  8 1000 prefit-pbmc-purified-k=8
sbatch ${SCRIPT_PREFIT} pbmc_purified.RData  9 1000 prefit-pbmc-purified-k=9
sbatch ${SCRIPT_PREFIT} pbmc_purified.RData 10 1000 prefit-pbmc-purified-k=10
sbatch ${SCRIPT_PREFIT} pbmc_purified.RData 11 1000 prefit-pbmc-purified-k=11
sbatch ${SCRIPT_PREFIT} pbmc_purified.RData 12 1000 prefit-pbmc-purified-k=12
sbatch ${SCRIPT_PREFIT} pbmc_purified.RData 13 1000 prefit-pbmc-purified-k=13

# "Pre-fit" factorizations to the pulse-seq data.
#                       data            k    n outfile
sbatch ${SCRIPT_PREFIT} pulseseq.RData  2 1000 prefit-pulseseq-k=2
sbatch ${SCRIPT_PREFIT} pulseseq.RData  3 1000 prefit-pulseseq-k=3
sbatch ${SCRIPT_PREFIT} pulseseq.RData  4 1000 prefit-pulseseq-k=4
sbatch ${SCRIPT_PREFIT} pulseseq.RData  5 1000 prefit-pulseseq-k=5
sbatch ${SCRIPT_PREFIT} pulseseq.RData  6 1000 prefit-pulseseq-k=6
sbatch ${SCRIPT_PREFIT} pulseseq.RData  7 1000 prefit-pulseseq-k=7
sbatch ${SCRIPT_PREFIT} pulseseq.RData  8 1000 prefit-pulseseq-k=8
sbatch ${SCRIPT_PREFIT} pulseseq.RData  9 1000 prefit-pulseseq-k=9
sbatch ${SCRIPT_PREFIT} pulseseq.RData 10 1000 prefit-pulseseq-k=10
sbatch ${SCRIPT_PREFIT} pulseseq.RData 11 1000 prefit-pulseseq-k=11
sbatch ${SCRIPT_PREFIT} pulseseq.RData 12 1000 prefit-pulseseq-k=12
sbatch ${SCRIPT_PREFIT} pulseseq.RData 13 1000 prefit-pulseseq-k=13

# Fit factorizations to droplet data, with and without extrapolation.
#                    data     k method    n  ex outfile
sbatch ${SCRIPT_FIT} droplet  2 em     1000  no fit-droplet-em-k=2
sbatch ${SCRIPT_FIT} droplet  2 ccd    1000  no fit-droplet-ccd-k=2
sbatch ${SCRIPT_FIT} droplet  2 scd    1000  no fit-droplet-scd-k=2
sbatch ${SCRIPT_FIT} droplet  2 em     1000 yes fit-droplet-em-ex-k=2
sbatch ${SCRIPT_FIT} droplet  2 ccd    1000 yes fit-droplet-ccd-ex-k=2
sbatch ${SCRIPT_FIT} droplet  2 scd    1000 yes fit-droplet-scd-ex-k=2

#                    data           k method    n  ex outfile
sbatch ${SCRIPT_FIT} pbmc_purified 13 scd      40 yes fit-pbmc-purified-scd-ex-k-13
sbatch ${SCRIPT_FIT} pulseseq      13 scd      40 yes fit-pulseseq-scd-ex-k-13
