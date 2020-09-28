# Script for "re-fitting" a Poisson non-negative factorization to the
# droplet data after introducing a new topic for the T+N (tuft and
# neuroendocrine) cells based on the clustering of the topic
# proportions.
library(Matrix)
library(fastTopics)

# Load the previously prepared count data, the k = 7 Poisson NMF model
# fit, the 8 clusters identified in the clustering analysis.
load("../data/droplet.RData")
fit <- readRDS("../output/droplet/rds/fit-droplet-scd-ex-k=7.rds")$fit
samples <- readRDS("../output/droplet/clustering-droplet.rds")

# Add a new dimension to the Poisson NMF for the T+N (tuft and
# neuroendocrine) cluster. 
n    <- nrow(counts)
rows <- which(samples$cluster == "T+N")
f    <- matrix(pmax(0.001,colMeans(counts[rows,])))
l    <- matrix(1,n,1)
colnames(f) <- "k8"
colnames(l) <- "k8"
fit$F  <- cbind(fit$F,f)
fit$L  <- cbind(fit$L,l)
fit$Fy <- fit$F
fit$Ly <- fit$L
fit$Fn <- fit$F
fit$Ln <- fit$L
rm(n,rows,f,l)

# Begin by re-fitting the loadings only.
fit <- fit_poisson_nmf(counts,fit0 = fit,numiter = 20,
                        method = "scd",update.factors = NULL,
                        control = list(numiter = 4,nc = 8))

# Now we are ready to perform the main model re-fitting step.
timing <- system.time({
  fit <- fit_poisson_nmf(counts,fit0 = fit,numiter = 180,method = "scd",
                         control = list(extrapolate = TRUE,numiter = 4,nc = 8))
})
cat(sprintf("Computation took %0.2f seconds.\n",timing["elapsed"]))

# Save the results.
saveRDS(fit,file = "refit-droplet-scd-ex-k=8.rds")
