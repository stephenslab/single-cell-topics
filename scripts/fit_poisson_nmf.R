#! /usr/bin/env Rscript
#
# TO DO: Update this description.
#
# Fit a Poisson non-negative factorization to one of the single-cell
# RNA-seq data sets.
#
# This script is intended to be run from the command-line shell, with
# options that are processed with the optparse package. For example,
# to fit a rank-4 Poisson non-negative matrix factorization by running
# 500 SCD updates with extrapolation, 2 threads, with estimates
# initialized from fit.rds, and with results saved to test.rds, run
# this command:
#
#   ./fit_poisson_nmf_purified_pbmc.R -k 4 --method scd --nc 2 \
#     --extrapolate --numiter 500 --prefit fit.rds -o test.rds
#
# Poisson non-negative matrix factorization by running 1,000 EM
# updates without extrapolation, without multithreading, and with
# results saved to out.rds. The estimates of the factors and loadings
# are initialized from ../output/prefit-k=3.rds by default when k = 3.
#

# Load a few packages.
library(Matrix)
library(optparse)
library(fastTopics)

# Process the command-line arguments.
parser <- OptionParser()
parser <- add_option(parser,"--counts",type="character",default="counts.RData")
parser <- add_option(parser,"--prefit",type = "character",default="prefix.rds")
parser <- add_option(parser,c("--out","-o"),type="character",default="out.rds")
parser <- add_option(parser,"--method",type = "character",default = "em")
parser <- add_option(parser,c("--numiter","-n"),type="integer",default=1000)
parser <- add_option(parser,"--extrapolate",action = "store_true")
parser <- add_option(parser,"--nc",type = "integer",default = 1)
out    <- parse_args(parser)
countsfile  <- out$counts
prefitfile  <- out$prefit
outfile     <- out$out
method      <- out$method
numiter     <- out$numiter
extrapolate <- !is.null(out$extrapolate)
nc          <- out$nc
rm(parser,out)

print(countsfile)
print(prefitfile)
print(outfile)
print(method)
print(numiter)
print(extrapolate)
print(nc)

# Initialize the sequence of pseudorandom numbers.
set.seed(1)

# LOAD DATA
# ---------
# Load the previously prepared count data.
cat(sprintf("Loading data from %s.\n",countsfile))
load(countsfile)
cat(sprintf("Loaded %d x %d counts matrix.\n",nrow(counts),ncol(counts)))

# LOAD PRE-FITTED MODEL
# ---------------------
cat(sprintf("Loading pre-fitted model from %s\n",prefitfile))
fit0 <- readRDS(prefitfile)$fit

# FIT POISSON NON-NEGATIVE MATRIX FACTORIZATION
# ---------------------------------------------
# Now we are ready to perform the main model-fitting step.
cat("Fitting Poisson NMF to count data.\n")
control <- list(extrapolate = extrapolate,nc = nc,
                numiter = ifelse(method == "ccd",1,4))
timing <- system.time({
  fit <- fit_poisson_nmf(counts,fit0 = fit0,numiter = numiter,
                         method = method,control = control)
})
cat(sprintf("Computation took %0.2f seconds.\n",timing["elapsed"]))

# SAVE RESULTS
# ------------
cat("Saving results.\n")
saveRDS(list(method = method,control = control,fit = fit),
        file = outfile)

# SESSION INFO
# ------------
print(sessionInfo())
