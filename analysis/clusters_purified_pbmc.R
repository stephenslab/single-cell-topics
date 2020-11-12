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

# Define B-cell, CD14+ and CD34+ clusters.
pca <- prcomp(poisson2multinom(fit)$L)$x
n   <- nrow(pca)
x   <- rep("A",n)
pc1 <- pca[,1]
pc2 <- pca[,2]
pc3 <- pca[,3]
pc4 <- pca[,4]
pc5 <- pca[,5]
x[pc2 > 0.25] <- "B"
x[pc3 < -0.2 & pc4 < 0.2] <- "CD34+"
x[(pc4 + 0.1)^2 + (pc5 - 0.8)^2 < 0.07] <- "CD14+"

# TO DO: Explain here what these next lines of code do.
rows <- 1:94655
fit2 <- select(poisson2multinom(fit),loadings = rows)

# Project cells onto PCs 1 and 2.
colors <- c("dodgerblue","forestgreen","greenyellow","magenta",
            "firebrick","darkorange","gold","darkblue","darkmagenta","gray")
p1 <- pca_plot(fit2,pcs = 1:2,fill = "none")
p2 <- pca_hexbin_plot(fit2,pcs = 1:2,bins = 50) +
  scale_x_continuous(breaks = seq(-1,1,0.1)) +
  scale_y_continuous(breaks = seq(-1,1,0.1))
p3 <- pca_plot(fit2,fill = samples[rows,"celltype"],pcs = 1:2) +
  scale_fill_manual(values = colors)
p4 <- pca_plot(fit2,pcs = 1:2,fill = factor(x[rows]))
plot_grid(p1,p2,p3,p4)

# Compare our clusters with FACS cell-types.
samples$cluster <- factor(x)
with(samples,table(celltype,cluster))
