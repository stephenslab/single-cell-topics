# Here we perform PCA on the topic proportions to identify clusters in
# the mixture of FACS-purified PBMC data sets.
library(Matrix)
library(dplyr)
library(fastTopics)
library(ggplot2)
library(cowplot)
source("../code/plots.R")

# Load the count data.
load("../data/pbmc_purified.RData")

# Load the k = 6 Poisson NMF model fit.
fit <- readRDS(file.path("../output/pbmc-purified/rds",
                         "fit-pbmc-purified-scd-ex-k=6.rds"))$fit

# Define clusters.
pca <- prcomp(poisson2multinom(fit)$L)$x
pc1 <- pca[,1]
pc2 <- pca[,2]
x   <- rep("U",nrow(pca))
x[pc1 < -0.15 & pc2 > 0.15] <- "B"

# Project cells onto PCs 1 and 2.
colors <- c("dodgerblue","forestgreen","greenyellow","magenta",
            "firebrick","darkorange","gold","darkblue","darkmagenta","gray")
p1 <- pca_plot(poisson2multinom(fit),pcs = 1:2,fill = "none")
p2 <- pca_hexbin_plot(poisson2multinom(fit),pcs = 1:2)
p3 <- pca_plot(poisson2multinom(fit),fill = samples$celltype,pcs = 1:2) +
  scale_fill_manual(values = colors)
p4 <- pca_plot(poisson2multinom(fit),pcs = 1:2,fill = factor(x))
plot_grid(p1,p2,p3,p4)

# Compare clusters with FACS cell-types.
samples$cluster <- factor(x)
with(samples,table(celltype,cluster))
