library(Matrix)
library(dplyr)
library(fastTopics)
library(ggplot2)
library(cowplot)

# Load the UMI count data, the K = 6 Poisson NMF model fit, and the
# clusters identified in the clustering analysis.
load("../data/pbmc_purified.RData")
fit <- readRDS(file.path("../output/pbmc-purified/rds",
                         "fit-pbmc-purified-scd-ex-k=6.rds"))$fit
samples <- readRDS("../output/pbmc-purified/clustering-pbmc-purified.rds")

# Create PCA plots.
cluster_colors <- c("dodgerblue","forestgreen","darkmagenta",
                    "red","gray","darkorange")
rows <- which(with(samples,is.element(cluster,c("NK","T","CD8+","U"))))
p1 <- pca_plot(poisson2multinom(fit),pcs = 1:2,fill = samples$cluster) +
  scale_fill_manual(values = cluster_colors) +
  labs(fill = "cluster")
p2 <- pca_plot(select(poisson2multinom(fit),loadings = rows),pcs = 1:2,
               fill = samples$cluster[rows,drop = FALSE]) +
  scale_fill_manual(values = cluster_colors,drop = FALSE) +
  labs(fill = "cluster")
p2 <- pca_hexbin_plot(poisson2multinom(fit))
p3 <- pca_hexbin_plot(poisson2multinom(fit))

# Create Structure plot.
