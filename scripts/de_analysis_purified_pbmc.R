# A short script used to perform the differential expression (DE)
# analysis using the multinomial topic model fitted to the mixture of
# purified PBMC data, with k = 6 topics. Computing the log-fold change
# posterior estimates by simulating a Markov chain with 100,000 states
# takes about 56. These were the steps taken to load R and allocate
# computing resources for this analysis:
#
#   sinteractive -p mstephens --account=pi-mstephens -c 20 \
#     --mem=16G --time=60:00:00
#   module load R/3.5.1
#

# Load a few packages.
library(tools)
library(Matrix)
library(fastTopics)

# Initialize the sequence of pseudorandom numbers.
seed <- 1
outfile <- "de-pbmc-purified-seed=1.RData"
print(seed)
print(outfile)
set.seed(seed)

# Load the count data.
load("../data/pbmc_purified.RData")

# Load the K = 6 multinomial topic model fit.
fit <- readRDS(file.path("../output/pbmc-purified/rds",
                         "fit-pbmc-purified-scd-ex-k=6.rds"))$fit
fit <- poisson2multinom(fit)

# Perform the DE analysis.
t0 <- proc.time()
de <- de_analysis(fit,counts,pseudocount = 0.1,
                  control = list(ns = 1e5,nc = 20,nsplit = 1000))
t1 <- proc.time()
timing <- t1 - t0
cat(sprintf("Computation took %0.2f seconds.\n",timing["elapsed"]))

# Save the results.
save(list = c("seed","genes","de"),
     file = outfile)
resaveRdaFiles(outfile)
