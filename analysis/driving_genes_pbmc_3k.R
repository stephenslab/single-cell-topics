library(fastTopics)
library(ggplot2)
library(cowplot)
set.seed(1)

# Load the UMI counts, and the fitted multinomial topic model.
data(pbmc_facs)
counts <- pbmc_facs$counts
fit    <- pbmc_facs$fit

# Perform DE analysis using fitted multinomial topic model.
t0 <- proc.time()
de <- de_analysis(fit,counts)
t1 <- proc.time()
print(t1 - t0)

# Compute the KL-divergence based "distinctiveness" measure used in
# CountClust.
t0 <- proc.time()
e <- 1e-8
m <- nrow(fit$F)
k <- ncol(fit$F)
D <- matrix(0,m,k)
rownames(D) <- rownames(fit$F)
colnames(D) <- colnames(fit$F)
for (i in 1:m)
  for (j in 1:k) {
    kl <- rep(0,k)
    for (l in 1:k) {
      fj    <- fit$F[i,j] + e
      fl    <- fit$F[i,l] + e
      kl[l] <- fj*log(fj/fl) + fl - fj
    }
    D[i,j] <- min(kl[-j])
  }
t1 <- proc.time()
print(t1 - t0)
