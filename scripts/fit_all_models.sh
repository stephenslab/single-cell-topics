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

sbatch ${SCRIPT_FIT} droplet  3 em     1000  no fit-droplet-em-k=3
sbatch ${SCRIPT_FIT} droplet  3 ccd    1000  no fit-droplet-ccd-k=3
sbatch ${SCRIPT_FIT} droplet  3 scd    1000  no fit-droplet-scd-k=3
sbatch ${SCRIPT_FIT} droplet  3 em     1000 yes fit-droplet-em-ex-k=3
sbatch ${SCRIPT_FIT} droplet  3 ccd    1000 yes fit-droplet-ccd-ex-k=3
sbatch ${SCRIPT_FIT} droplet  3 scd    1000 yes fit-droplet-scd-ex-k=3

sbatch ${SCRIPT_FIT} droplet  4 em     1000  no fit-droplet-em-k=4
sbatch ${SCRIPT_FIT} droplet  4 ccd    1000  no fit-droplet-ccd-k=4
sbatch ${SCRIPT_FIT} droplet  4 scd    1000  no fit-droplet-scd-k=4
sbatch ${SCRIPT_FIT} droplet  4 em     1000 yes fit-droplet-em-ex-k=4
sbatch ${SCRIPT_FIT} droplet  4 ccd    1000 yes fit-droplet-ccd-ex-k=4
sbatch ${SCRIPT_FIT} droplet  4 scd    1000 yes fit-droplet-scd-ex-k=4

sbatch ${SCRIPT_FIT} droplet  5 em     1000  no fit-droplet-em-k=5
sbatch ${SCRIPT_FIT} droplet  5 ccd    1000  no fit-droplet-ccd-k=5
sbatch ${SCRIPT_FIT} droplet  5 scd    1000  no fit-droplet-scd-k=5
sbatch ${SCRIPT_FIT} droplet  5 em     1000 yes fit-droplet-em-ex-k=5
sbatch ${SCRIPT_FIT} droplet  5 ccd    1000 yes fit-droplet-ccd-ex-k=5
sbatch ${SCRIPT_FIT} droplet  5 scd    1000 yes fit-droplet-scd-ex-k=5

sbatch ${SCRIPT_FIT} droplet  6 em     1000  no fit-droplet-em-k=6
sbatch ${SCRIPT_FIT} droplet  6 ccd    1000  no fit-droplet-ccd-k=6
sbatch ${SCRIPT_FIT} droplet  6 scd    1000  no fit-droplet-scd-k=6
sbatch ${SCRIPT_FIT} droplet  6 em     1000 yes fit-droplet-em-ex-k=6
sbatch ${SCRIPT_FIT} droplet  6 ccd    1000 yes fit-droplet-ccd-ex-k=6
sbatch ${SCRIPT_FIT} droplet  6 scd    1000 yes fit-droplet-scd-ex-k=6

sbatch ${SCRIPT_FIT} droplet  7 em     1000  no fit-droplet-em-k=7
sbatch ${SCRIPT_FIT} droplet  7 ccd    1000  no fit-droplet-ccd-k=7
sbatch ${SCRIPT_FIT} droplet  7 scd    1000  no fit-droplet-scd-k=7
sbatch ${SCRIPT_FIT} droplet  7 em     1000 yes fit-droplet-em-ex-k=7
sbatch ${SCRIPT_FIT} droplet  7 ccd    1000 yes fit-droplet-ccd-ex-k=7
sbatch ${SCRIPT_FIT} droplet  7 scd    1000 yes fit-droplet-scd-ex-k=7

sbatch ${SCRIPT_FIT} droplet  8 em     1000  no fit-droplet-em-k=8
sbatch ${SCRIPT_FIT} droplet  8 ccd    1000  no fit-droplet-ccd-k=8
sbatch ${SCRIPT_FIT} droplet  8 scd    1000  no fit-droplet-scd-k=8
sbatch ${SCRIPT_FIT} droplet  8 em     1000 yes fit-droplet-em-ex-k=8
sbatch ${SCRIPT_FIT} droplet  8 ccd    1000 yes fit-droplet-ccd-ex-k=8
sbatch ${SCRIPT_FIT} droplet  8 scd    1000 yes fit-droplet-scd-ex-k=8

sbatch ${SCRIPT_FIT} droplet  9 em     1000  no fit-droplet-em-k=9
sbatch ${SCRIPT_FIT} droplet  9 ccd    1000  no fit-droplet-ccd-k=9
sbatch ${SCRIPT_FIT} droplet  9 scd    1000  no fit-droplet-scd-k=9
sbatch ${SCRIPT_FIT} droplet  9 em     1000 yes fit-droplet-em-ex-k=9
sbatch ${SCRIPT_FIT} droplet  9 ccd    1000 yes fit-droplet-ccd-ex-k=9
sbatch ${SCRIPT_FIT} droplet  9 scd    1000 yes fit-droplet-scd-ex-k=9

sbatch ${SCRIPT_FIT} droplet  10 em     1000  no fit-droplet-em-k=10
sbatch ${SCRIPT_FIT} droplet  10 ccd    1000  no fit-droplet-ccd-k=10
sbatch ${SCRIPT_FIT} droplet  10 scd    1000  no fit-droplet-scd-k=10
sbatch ${SCRIPT_FIT} droplet  10 em     1000 yes fit-droplet-em-ex-k=10
sbatch ${SCRIPT_FIT} droplet  10 ccd    1000 yes fit-droplet-ccd-ex-k=10
sbatch ${SCRIPT_FIT} droplet  10 scd    1000 yes fit-droplet-scd-ex-k=10

sbatch ${SCRIPT_FIT} droplet  11 em     1000  no fit-droplet-em-k=11
sbatch ${SCRIPT_FIT} droplet  11 ccd    1000  no fit-droplet-ccd-k=11
sbatch ${SCRIPT_FIT} droplet  11 scd    1000  no fit-droplet-scd-k=11
sbatch ${SCRIPT_FIT} droplet  11 em     1000 yes fit-droplet-em-ex-k=11
sbatch ${SCRIPT_FIT} droplet  11 ccd    1000 yes fit-droplet-ccd-ex-k=11
sbatch ${SCRIPT_FIT} droplet  11 scd    1000 yes fit-droplet-scd-ex-k=11

sbatch ${SCRIPT_FIT} droplet  12 em     1000  no fit-droplet-em-k=12
sbatch ${SCRIPT_FIT} droplet  12 ccd    1000  no fit-droplet-ccd-k=12
sbatch ${SCRIPT_FIT} droplet  12 scd    1000  no fit-droplet-scd-k=12
sbatch ${SCRIPT_FIT} droplet  12 em     1000 yes fit-droplet-em-ex-k=12
sbatch ${SCRIPT_FIT} droplet  12 ccd    1000 yes fit-droplet-ccd-ex-k=12
sbatch ${SCRIPT_FIT} droplet  12 scd    1000 yes fit-droplet-scd-ex-k=12

sbatch ${SCRIPT_FIT} droplet  13 em     1000  no fit-droplet-em-k=13
sbatch ${SCRIPT_FIT} droplet  13 ccd    1000  no fit-droplet-ccd-k=13
sbatch ${SCRIPT_FIT} droplet  13 scd    1000  no fit-droplet-scd-k=13
sbatch ${SCRIPT_FIT} droplet  13 em     1000 yes fit-droplet-em-ex-k=13
sbatch ${SCRIPT_FIT} droplet  13 ccd    1000 yes fit-droplet-ccd-ex-k=13
sbatch ${SCRIPT_FIT} droplet  13 scd    1000 yes fit-droplet-scd-ex-k=13

# Fit factorizations to 68k PBMC data, with and without extrapolation.
#                    data      k method    n  ex outfile
sbatch ${SCRIPT_FIT} pbmc_68k  2 em     1000  no fit-pbmc-68k-em-k=2
sbatch ${SCRIPT_FIT} pbmc_68k  2 ccd    1000  no fit-pbmc-68k-ccd-k=2
sbatch ${SCRIPT_FIT} pbmc_68k  2 scd    1000  no fit-pbmc-68k-scd-k=2
sbatch ${SCRIPT_FIT} pbmc_68k  2 em     1000 yes fit-pbmc-68k-em-ex-k=2
sbatch ${SCRIPT_FIT} pbmc_68k  2 ccd    1000 yes fit-pbmc-68k-ccd-ex-k=2
sbatch ${SCRIPT_FIT} pbmc_68k  2 scd    1000 yes fit-pbmc-68k-scd-ex-k=2

sbatch ${SCRIPT_FIT} pbmc_68k  3 em     1000  no fit-pbmc-68k-em-k=3
sbatch ${SCRIPT_FIT} pbmc_68k  3 ccd    1000  no fit-pbmc-68k-ccd-k=3
sbatch ${SCRIPT_FIT} pbmc_68k  3 scd    1000  no fit-pbmc-68k-scd-k=3
sbatch ${SCRIPT_FIT} pbmc_68k  3 em     1000 yes fit-pbmc-68k-em-ex-k=3
sbatch ${SCRIPT_FIT} pbmc_68k  3 ccd    1000 yes fit-pbmc-68k-ccd-ex-k=3
sbatch ${SCRIPT_FIT} pbmc_68k  3 scd    1000 yes fit-pbmc-68k-scd-ex-k=3

sbatch ${SCRIPT_FIT} pbmc_68k  4 em     1000  no fit-pbmc-68k-em-k=4
sbatch ${SCRIPT_FIT} pbmc_68k  4 ccd    1000  no fit-pbmc-68k-ccd-k=4
sbatch ${SCRIPT_FIT} pbmc_68k  4 scd    1000  no fit-pbmc-68k-scd-k=4
sbatch ${SCRIPT_FIT} pbmc_68k  4 em     1000 yes fit-pbmc-68k-em-ex-k=4
sbatch ${SCRIPT_FIT} pbmc_68k  4 ccd    1000 yes fit-pbmc-68k-ccd-ex-k=4
sbatch ${SCRIPT_FIT} pbmc_68k  4 scd    1000 yes fit-pbmc-68k-scd-ex-k=4

sbatch ${SCRIPT_FIT} pbmc_68k  5 em     1000  no fit-pbmc-68k-em-k=5
sbatch ${SCRIPT_FIT} pbmc_68k  5 ccd    1000  no fit-pbmc-68k-ccd-k=5
sbatch ${SCRIPT_FIT} pbmc_68k  5 scd    1000  no fit-pbmc-68k-scd-k=5
sbatch ${SCRIPT_FIT} pbmc_68k  5 em     1000 yes fit-pbmc-68k-em-ex-k=5
sbatch ${SCRIPT_FIT} pbmc_68k  5 ccd    1000 yes fit-pbmc-68k-ccd-ex-k=5
sbatch ${SCRIPT_FIT} pbmc_68k  5 scd    1000 yes fit-pbmc-68k-scd-ex-k=5

sbatch ${SCRIPT_FIT} pbmc_68k  6 em     1000  no fit-pbmc-68k-em-k=6
sbatch ${SCRIPT_FIT} pbmc_68k  6 ccd    1000  no fit-pbmc-68k-ccd-k=6
sbatch ${SCRIPT_FIT} pbmc_68k  6 scd    1000  no fit-pbmc-68k-scd-k=6
sbatch ${SCRIPT_FIT} pbmc_68k  6 em     1000 yes fit-pbmc-68k-em-ex-k=6
sbatch ${SCRIPT_FIT} pbmc_68k  6 ccd    1000 yes fit-pbmc-68k-ccd-ex-k=6
sbatch ${SCRIPT_FIT} pbmc_68k  6 scd    1000 yes fit-pbmc-68k-scd-ex-k=6

sbatch ${SCRIPT_FIT} pbmc_68k  7 em     1000  no fit-pbmc-68k-em-k=7
sbatch ${SCRIPT_FIT} pbmc_68k  7 ccd    1000  no fit-pbmc-68k-ccd-k=7
sbatch ${SCRIPT_FIT} pbmc_68k  7 scd    1000  no fit-pbmc-68k-scd-k=7
sbatch ${SCRIPT_FIT} pbmc_68k  7 em     1000 yes fit-pbmc-68k-em-ex-k=7
sbatch ${SCRIPT_FIT} pbmc_68k  7 ccd    1000 yes fit-pbmc-68k-ccd-ex-k=7
sbatch ${SCRIPT_FIT} pbmc_68k  7 scd    1000 yes fit-pbmc-68k-scd-ex-k=7

sbatch ${SCRIPT_FIT} pbmc_68k  8 em     1000  no fit-pbmc-68k-em-k=8
sbatch ${SCRIPT_FIT} pbmc_68k  8 ccd    1000  no fit-pbmc-68k-ccd-k=8
sbatch ${SCRIPT_FIT} pbmc_68k  8 scd    1000  no fit-pbmc-68k-scd-k=8
sbatch ${SCRIPT_FIT} pbmc_68k  8 em     1000 yes fit-pbmc-68k-em-ex-k=8
sbatch ${SCRIPT_FIT} pbmc_68k  8 ccd    1000 yes fit-pbmc-68k-ccd-ex-k=8
sbatch ${SCRIPT_FIT} pbmc_68k  8 scd    1000 yes fit-pbmc-68k-scd-ex-k=8

sbatch ${SCRIPT_FIT} pbmc_68k  9 em     1000  no fit-pbmc-68k-em-k=9
sbatch ${SCRIPT_FIT} pbmc_68k  9 ccd    1000  no fit-pbmc-68k-ccd-k=9
sbatch ${SCRIPT_FIT} pbmc_68k  9 scd    1000  no fit-pbmc-68k-scd-k=9
sbatch ${SCRIPT_FIT} pbmc_68k  9 em     1000 yes fit-pbmc-68k-em-ex-k=9
sbatch ${SCRIPT_FIT} pbmc_68k  9 ccd    1000 yes fit-pbmc-68k-ccd-ex-k=9
sbatch ${SCRIPT_FIT} pbmc_68k  9 scd    1000 yes fit-pbmc-68k-scd-ex-k=9

sbatch ${SCRIPT_FIT} pbmc_68k  10 em     1000  no fit-pbmc-68k-em-k=10
sbatch ${SCRIPT_FIT} pbmc_68k  10 ccd    1000  no fit-pbmc-68k-ccd-k=10
sbatch ${SCRIPT_FIT} pbmc_68k  10 scd    1000  no fit-pbmc-68k-scd-k=10
sbatch ${SCRIPT_FIT} pbmc_68k  10 em     1000 yes fit-pbmc-68k-em-ex-k=10
sbatch ${SCRIPT_FIT} pbmc_68k  10 ccd    1000 yes fit-pbmc-68k-ccd-ex-k=10
sbatch ${SCRIPT_FIT} pbmc_68k  10 scd    1000 yes fit-pbmc-68k-scd-ex-k=10

sbatch ${SCRIPT_FIT} pbmc_68k  11 em     1000  no fit-pbmc-68k-em-k=11
sbatch ${SCRIPT_FIT} pbmc_68k  11 ccd    1000  no fit-pbmc-68k-ccd-k=11
sbatch ${SCRIPT_FIT} pbmc_68k  11 scd    1000  no fit-pbmc-68k-scd-k=11
sbatch ${SCRIPT_FIT} pbmc_68k  11 em     1000 yes fit-pbmc-68k-em-ex-k=11
sbatch ${SCRIPT_FIT} pbmc_68k  11 ccd    1000 yes fit-pbmc-68k-ccd-ex-k=11
sbatch ${SCRIPT_FIT} pbmc_68k  11 scd    1000 yes fit-pbmc-68k-scd-ex-k=11

sbatch ${SCRIPT_FIT} pbmc_68k  12 em     1000  no fit-pbmc-68k-em-k=12
sbatch ${SCRIPT_FIT} pbmc_68k  12 ccd    1000  no fit-pbmc-68k-ccd-k=12
sbatch ${SCRIPT_FIT} pbmc_68k  12 scd    1000  no fit-pbmc-68k-scd-k=12
sbatch ${SCRIPT_FIT} pbmc_68k  12 em     1000 yes fit-pbmc-68k-em-ex-k=12
sbatch ${SCRIPT_FIT} pbmc_68k  12 ccd    1000 yes fit-pbmc-68k-ccd-ex-k=12
sbatch ${SCRIPT_FIT} pbmc_68k  12 scd    1000 yes fit-pbmc-68k-scd-ex-k=12

sbatch ${SCRIPT_FIT} pbmc_68k  13 em     1000  no fit-pbmc-68k-em-k=13
sbatch ${SCRIPT_FIT} pbmc_68k  13 ccd    1000  no fit-pbmc-68k-ccd-k=13
sbatch ${SCRIPT_FIT} pbmc_68k  13 scd    1000  no fit-pbmc-68k-scd-k=13
sbatch ${SCRIPT_FIT} pbmc_68k  13 em     1000 yes fit-pbmc-68k-em-ex-k=13
sbatch ${SCRIPT_FIT} pbmc_68k  13 ccd    1000 yes fit-pbmc-68k-ccd-ex-k=13
sbatch ${SCRIPT_FIT} pbmc_68k  13 scd    1000 yes fit-pbmc-68k-scd-ex-k=13
