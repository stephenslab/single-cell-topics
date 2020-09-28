# TO DO: Revise this description.
#
# Script for "re-fitting" a Poisson non-negative factorization to the
# pulse-seq data after introducing a new topic for the T+N (tuft and
# neuroendocrine) cells based on the clustering of the topic
# proportions.
library(Matrix)
library(fastTopics)

# Load the previously prepared count data, the k = 7 Poisson NMF model
# fit, the 8 clusters identified in the clustering analysis.
load("../data/pulseseq.RData")
fit <- readRDS("../output/pulseseq/rds/fit-pulseseq-scd-ex-k=11.rds")$fit
samples <- readRDS("../output/pulseseq/clustering-pulseseq.rds")

# Add two new dimensions to the Poisson NMF for the P
# ("proliferatitng") and I ("ionocyte") clusters.
n     <- nrow(counts)
rows1 <- which(samples$cluster == "I")
rows2 <- which(samples$cluster == "P")
F     <- cbind(pmax(0.001,colMeans(counts[rows1,])),
               pmax(0.001,colMeans(counts[rows2,])))
L     <- matrix(1,n,2)
colnames(F) <- c("k12","k13")
colnames(L) <- c("k12","k13")
fit$F  <- cbind(fit$F,F)
fit$L  <- cbind(fit$L,L)
fit$Fy <- fit$F
fit$Ly <- fit$L
fit$Fn <- fit$F
fit$Ln <- fit$L
rm(n,rows1,rows2,F,L)

# Begin by re-fitting the loadings only.
fit <- fit_poisson_nmf(counts,fit0 = fit,numiter = 2,
                        method = "scd",update.factors = NULL,
                        control = list(numiter = 4,nc = 4))

# Now we are ready to perform the main model re-fitting step.
timing <- system.time({
  fit <- fit_poisson_nmf(counts,fit0 = fit,numiter = 10,method = "scd",
                         control = list(extrapolate = TRUE,numiter = 4,nc = 4))
})
cat(sprintf("Computation took %0.2f seconds.\n",timing["elapsed"]))

# Save the results.
saveRDS(fit,file = "refit-pulseseq-scd-ex-k=13.rds")
