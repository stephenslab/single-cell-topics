library(Matrix)
library(dplyr)
library(fastTopics)

# Load data and results.
load("../data/droplet.RData")
samples <- readRDS("../output/droplet/clustering-droplet.rds")
fit <- readRDS("../output/droplet/rds/fit-droplet-scd-ex-k=7.rds")$fit

totalvardist <- function (F) {
  n <- ncol(F)
  d <- matrix(0,n,n)
  rownames(d) <- colnames(F)
  colnames(d) <- colnames(F)
  for (i in 1:n)
    for (j in 1:n)
      d[i,j] <- sum(abs(F[,i] - F[,j]))/2
  return(d)
}

# Compute relative gene expression levels from Montoro et al (2018)
# clustering, and compute "distances" based on these relative
# expression levels.
fit_montoro <- init_poisson_nmf_from_clustering(counts,samples$tissue)
d_montoro   <- totalvardist(poisson2multinom(fit_montoro)$F)
print(d_montoro)

# Compute relative gene expression levels from clustering based on
# topic model, and compute "distances" based on these relative
# expression levels.
fit_cluster <- init_poisson_nmf_from_clustering(counts,samples$cluster)
d_cluster   <- totalvardist(poisson2multinom(fit_cluster)$F)
i <- c("B","C","Cil","G","T+N")
print(d_cluster[i,i])

# Compute relative gene expression levels from clustering based on
# topic model, and compute "distances" based on these relative
# expression levels.
fit_merge <- merge_topics(poisson2multinom(fit),c(5,7))
d_topics <- totalvardist(fit_merge$F)
i <- c("k1","k2","k4","k6","k5+k7")
print(d_topics[i,i])
