# TO DO: Explain here what this script does, and how to use it.

# Load a couple packages.
library(Matrix)
library(fastTopics)

# Initialize the sequence of pseudorandom numbers.
set.seed(1)

# Load the previously prepared count data.
load("../data/droplet.RData")

# Select cells that have at least 10% membership to topic 6.
fit     <- readRDS("../output/droplet/rds/fit-droplet-scd-ex-k=7.rds")$fit
fit     <- poisson2multinom(fit)
rows    <- which(fit$L[,6] > 0.1)
samples <- samples[rows,]
counts  <- counts[rows,]

# FIT POISSON NMF
# ---------------
# Now we are ready to perform the main model-fitting step.
timing <- system.time({
  fit <- fit_poisson_nmf(counts,k = 5,numiter = 200,method = "scd",
                         innt.method = "random",
                         control = list(extrapolate = TRUE,numiter = 4,nc = 4))
})
cat(sprintf("Computation took %0.2f seconds.\n",timing["elapsed"]))

# SAVE RESULTS
# ------------
saveRDS(fit,file = "refit-droplet-scd-ex=k=7.rds")

