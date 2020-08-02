#! /usr/bin/env Rscript
#
# Pre-fit a Poisson non-negative factorization to the purified PBMC
# single-cell RNA-seq data of Zheng et al (2017).
#
# This script is intended to be run from the command-line shell, with
# options that are processed with the optparse package. For example,
# to fit a rank-4 Poisson non-negative matrix factorization by running
# SCD updatees with extrapolation, 2 threads, and with results saved
# to test.rds, run this command:
#
#   ./prefit_poisson_nmf_purified_pbmc.R -k 4 --nc 2 -o test.rds
#
# Running the script without specifying any options will pre-fit a
# rank-3 Poisson non-negative matrix factorization by running EM
# updates (without extrapolation), without multithreading, and with
# results saved to out.rds.
#

# Load a few packages.
library(Matrix)
library(optparse)
library(fastTopics)

# Process the command-line arguments.
parser <- OptionParser()
parser <- add_option(parser,c("--out","-o"),type="character",default="out.rds")
parser <- add_option(parser,c("--k","-k"),type = "integer",default = 3)
parser <- add_option(parser,"--nc",type = "integer",default = 1)
out    <- parse_args(parser)
outfile     <- out$out
k           <- out$k
nc          <- out$nc
rm(parser,out)

# Initialize the sequence of pseudorandom numbers.
set.seed(1)

# LOAD DATA
# ---------
# Load the previously prepared purified PBMC data.
cat("Loading PBMC data.\n")
load("../data/pbmc_purified.RData")
cat(sprintf("Loaded %d x %d counts matrix.\n",nrow(counts),ncol(counts)))

# PRE-FIT POISSON NON-NEGATIVE MATRIX FACTORIZATION
# ---------------------------------------------
# The aim here is to run enough EM updates so that it is difficult for
# the other algorithms to "escape" this local maximum of the
# likelihood surface.
cat("Running 1,000 EM updates to identify a good initialization.\n")
timing <- system.time(
  fit <- fit_poisson_nmf(counts,k = k,numiter = 1000,method = "em",
                         control = list(numiter = 4,nc = nc)))
cat(sprintf("Computation took %0.2f seconds.\n",timing["elapsed"]))

# SAVE RESULTS
# ------------
cat("Saving results.\n")
saveRDS(list(k = k,fit = fit),file = outfile)

# SESSION INFO
# ------------
print(sessionInfo())
