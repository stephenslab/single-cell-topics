# A short script used to perform the differential expression (DE)
# analysis using the multinomial topic model fitted to the droplet
# data, with k = 7 topics. These were the steps taken to load R and
# allocate computing resources for this analysis:
#
#   sinteractive -p broadwl -c 20 --mem=16G --time=10:00:00
#   module load R/3.5.1
#

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
de <- de_analysis(fit,counts,pseudocount = 0.1,
                  control = list(ns = 1e5,nc = 20,nsplit = 1000))
t1 <- proc.time()
timing <- t1 - t0
cat(sprintf("Computation took %0.2f seconds.\n",timing["elapsed"]))

# Perform a second DE analysis after merging topics X and Y.
t0 <- proc.time()
fit_merged <- merge_topics(fit,c("k5","k7"))
de_merged <- de_analysis(fit_merged,counts,pseudocount = 0.1,
                         control = list(ns = 1e5,nc = 20,nsplit = 1000))
t1 <- proc.time()
timing <- t1 - t0
cat(sprintf("Computation took %0.2f seconds.\n",timing["elapsed"]))

# Save the results.
save(list = c("de","de_merged"),
     file = "de-droplet.RData")
resaveRdaFiles("de-droplet.RData")
