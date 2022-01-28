# TO DO: Explain here what this script does, and how to use it.

# Load a few packages.
library(tools)
library(Matrix)
library(fastTopics)

# Initialize the sequence of pseudorandom numbers.
set.seed(1)

# Load the count data.
load("../data/droplet.RData")

# Load the K = 7 Poisson NMF model fit.
fit <- readRDS("../output/droplet/rds/fit-droplet-scd-ex-k=7.rds")$fit
fit <- poisson2multinom(fit)

# Perform the DE analysis.
t0 <- proc.time()
de1 <- de_analysis(fit,counts,pseudocount = 0.1,
                   control = list(ns = 1000,nc = 4,nsplit = 100))
t1 <- proc.time()
timing <- t1 - t0
cat(sprintf("Computation took %0.2f seconds.\n",timing["elapsed"]))

# Perform a second DE analysis after merging topics X and Y.
t0 <- proc.time()
fit_merged <- merge_topics(fit,c("k5","k7"))
de_merged <- de_analysis(fit_merged,counts,pseudocount = 0.1,
                         control = list(ns = 1000,nc = 4,nsplit = 100))
t1 <- proc.time()
timing <- t1 - t0
cat(sprintf("Computation took %0.2f seconds.\n",timing["elapsed"]))

# Save the results.
save(list = c("genes","de","de_merged"),
     file = "de-droplet.RData")
resaveRdaFiles("de-droplet.RData")
