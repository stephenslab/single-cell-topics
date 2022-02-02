# A short script used to perform the differential expression (DE)
# analysis using the multinomial topic model fitted to the pulseseq
# data, with k = 11 topics. These were the steps taken to load R and
# allocate computing resources for this analysis:
#
#   sinteractive -p broadwl -c 20 --mem=16G --time=60:00:00
#   module load R/3.5.1
#

# Load a few packages.
library(tools)
library(Matrix)
library(fastTopics)

# Initialize the sequence of pseudorandom numbers.
set.seed(1)

# Load the count data.
load("../data/pulseseq.RData")

# Load the K = 11 multinomial topic model fit.
fit <- readRDS("../output/pulseseq/rds/fit-pulseseq-scd-ex-k=11.rds")$fit
fit <- poisson2multinom(fit)

# Perform the DE analysis.
t0 <- proc.time()
de <- de_analysis(fit,counts,pseudocount = 0.1,
                  control = list(ns = 1000,nc = 20,nsplit = 1000))
t1 <- proc.time()
timing <- t1 - t0
cat(sprintf("Computation took %0.2f seconds.\n",timing["elapsed"]))

# Save the results.
save(list = c("seed","genes","de"),
     file = "de-pulseseq.RData")
resaveRdaFiles("de-pulseseq.RData")
