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
x   <- rep("U",n)
pc1 <- pca[,1]
pc2 <- pca[,2]
pc3 <- pca[,3]
pc4 <- pca[,4]
pc5 <- pca[,5]
x[pc2 > 0.25] <- "B"
x[pc3 < -0.2 & pc4 < 0.2] <- "CD34+"
x[(pc4 + 0.1)^2 + (pc5 - 0.8)^2 < 0.07] <- "CD14+"

# Define NK cells cluster.
rows <- which(x == "U")
n    <- length(rows)
fit2 <- select(poisson2multinom(fit),loadings = rows)
pca  <- prcomp(fit2$L)$x
y    <- rep("U",n)
pc1  <- pca[,1]
pc2  <- pca[,2]
y[pc1 < -0.3 & 1.1*pc1 < -pc2 - 0.57] <- "NK"
x[rows] <- y

# Define (less distinct) CD8+ cluster.
# TO DO.

# Plot the top two PCs of the remaining cells.
colors <- c("dodgerblue","forestgreen","greenyellow","magenta",
            "firebrick","darkorange","gold","darkblue","darkmagenta","gray")
rows <- which(x == "U")
fit2 <- select(poisson2multinom(fit),loadings = rows)
p1 <- pca_plot(fit2,fill = samples[rows,"celltype"],pcs = 1:2) +
  scale_fill_manual(values = colors)
p2 <- pca_plot(fit2,k = 1)
p3 <- pca_plot(fit2,k = 4)
p4 <- pca_plot(fit2,k = 6)
plot_grid(p1,p2,p3,p4)

# Compare our clusters with FACS cell-types.
samples$cluster <- factor(x)
with(samples,table(celltype,cluster))

# TO DO:
#
#  + Create PCA plots summarizing clustering.
#
#  + Create PCA plots showing FACS cell-types.
#
#  + Create Structure Plot summarizing clustering.
#
#  + Create scatterplots comparing likelihoods using a "hard"
#    clustering and using topic model.
#
