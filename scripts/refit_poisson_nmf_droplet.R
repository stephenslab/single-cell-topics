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

# Now we are ready to perform the main model-fitting step.
fits        <- vector("list",6)
names(fits) <- paste0("k",1:6)
for (k in 2:6)
  fits[[k]] <- fit_poisson_nmf(counts,k = k,numiter = 100,method = "scd",
                               init.method = "random",
                               control = list(extrapolate = TRUE,
                                              numiter = 4,nc = 4))
fits <- fits[-1]

plot_loglik_vs_rank(fits)

stop()

# Perform the DE analysis.
t0 <- proc.time()
de <- de_analysis(fit,counts,pseudocount = 0.1,control = list(ns = 1e4,nc = 4))
t1 <- proc.time()
timing <- t1 - t0
cat(sprintf("Computation took %0.2f seconds.\n",timing["elapsed"]))

# Perform a second DE analysis after merging topics 1 and 3.
t0 <- proc.time()
fit_merged <- merge_topics(fit,c("k1","k3"))
de_merged <- de_analysis(fit_merged,counts,pseudocount = 0.1,
                         control = list(ns = 1e4,nc = 4))
t1 <- proc.time()
timing <- t1 - t0
cat(sprintf("Computation took %0.2f seconds.\n",timing["elapsed"]))

# SAVE RESULTS
# ------------
saveRDS(list(fit = fit),file = "refit-droplet-scd-ex=k=7.rds")
