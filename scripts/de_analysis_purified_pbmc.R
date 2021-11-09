# TO DO: Explain here what this script is for, and how to use it.
# sinteractive -p mstephens --account=pi-mstephens -c 20 --mem=16G \
#   --time=24:00:00

# Load a few packages.
library(tools)
library(Matrix)
library(fastTopics)

# Initialize the sequence of pseudorandom numbers.
set.seed(1)

# Load the count data.
load("../data/pbmc_purified.RData")

# Load the K = 6 multinomial topic model fit.
fit <- readRDS(file.path("../output/pbmc-purified/rds",
                         "fit-pbmc-purified-scd-ex-k=6.rds"))$fit
fit <- poisson2multinom(fit)

# Perform the DE analysis.
set.seed(1)
t0  <- proc.time()
de1 <- de_analysis(fit,counts,control = list(ns = 10000,nc = 20))
t1  <- proc.time()
cat(sprintf("Computation took %0.2f seconds.\n",(t1 - t0)["elapsed"]))

# Perform the DE analysis a second time to assess accuracy of the
# posterior calculations.
set.seed(2)
t0  <- proc.time()
de2 <- de_analysis(fit,counts,control = list(ns = 10000,nc = 20))
t1  <- proc.time()
cat(sprintf("Computation took %0.2f seconds.\n",(t1 - t0)["elapsed"]))

# Save the results.
save(list = c("genes","de1","de2"),
     file = "de-pbmc-purified.RData")
resaveRdaFiles("de-pbmc-purified.RData")
