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

# Define (less distinct) CD8+ cluster, and label the remaining cells
# as T-cells.
rows <- which(x == "U")
n    <- length(rows)
fit2 <- select(poisson2multinom(fit),loadings = rows)
pca  <- prcomp(fit2$L)$x
y    <- rep("T",n)
pc1  <- pca[,1]
pc2  <- pca[,2]
y[pc1 < 0.25 & pc2 < -0.15] <- "CD8+"
x[rows] <- y

# Compare our clusters with FACS cell-types.
samples$cluster <- factor(x)
with(samples,table(celltype,cluster))

# Create PCA plots showing FACS cell-types.
colors <- c("dodgerblue","forestgreen","greenyellow","magenta",
            "firebrick","darkorange","gold","darkblue","darkmagenta","gray")
x     <- with(samples,cluster == "T" | cluster == "CD8+")
rows1 <- which(!x)
rows2 <- which(x)
p1 <- pca_plot(select(poisson2multinom(fit),loadings = rows1),
               fill = samples[rows1,"celltype"]) +
  scale_fill_manual(values = colors)
p2 <- pca_plot(select(poisson2multinom(fit),loadings = rows2),
               fill = samples[rows2,"celltype"]) +
  scale_fill_manual(values = colors)
plot_grid(p1,p2)

set.seed(5)
rows <- which(with(samples,
                   cluster == "B" |
                   cluster == "CD14+" |
                   cluster == "CD34+"))
rows <- sample(rows,2000)
fit2 <- select(poisson2multinom(fit),loadings = rows)
p3   <- pca_plot(fit2,fill = samples[rows,"cluster"])

# Run t-SNE, then plot the 2-d embedding.
library(Rtsne)
tsne <- Rtsne(fit2$L,dims = 2,pca = FALSE,normalize = FALSE,perplexity = 100,
              theta = 0.1,max_iter = 1000,eta = 200,verbose = TRUE)
tsne$x <- tsne$Y
colnames(tsne$x) <- c("tsne1","tsne2")
p4 <- pca_plot(fit2,out.pca = tsne,fill = samples[rows,"cluster"])

# Run UMAP, then plot the 2-d embedding.
library(uwot)
out.umap <- umap(fit2$L,n_neighbors = 30,metric = "euclidean",n_epochs = 1000,
                 min_dist = 0.1,scale = "none",learning_rate = 1,
                 verbose = TRUE)
out.umap <- list(x = out.umap)
colnames(out.umap$x) <- c("umap1","umap2")
p5 <- pca_plot(fit2,out.pca = out.umap,fill = samples[rows,"cluster"])
plot_grid(p3,p4,p5,nrow = 1)

# TO DO:
#
#  + Create Structure Plot summarizing clustering.
#
#  + Compare inter-cluster and inter-topic total variation distances.
# 
#  + Create scatterplots comparing likelihoods using a "hard"
#    clustering and using topic model.
#
