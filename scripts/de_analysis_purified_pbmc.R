# TO DO: Explain here what this script is for, and how to use it.
# sinteractive -p mstephens --account=pi-mstephens -c 20 --mem=16G \
#   --time=60:00:00
# Computation takes about 56 h.

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
i   <- match(rownames(counts),rownames(fit$L))
fit <- select_loadings(fit,i)

# Perform the DE analysis.
set.seed(1)
t0 <- proc.time()
de <- de_analysis(fit,counts,pseudocount = 0.1,
                  control = list(ns = 1000,nc = 20,nsplit = 1000))
t1 <- proc.time()
timing <- t1 - t0
cat(sprintf("Computation took %0.2f seconds.\n",timing["elapsed"]))

# Save the results.
save(list = c("seed","genes","de"),
     file = outfile)
