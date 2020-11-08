library(Matrix)
library(dplyr)
library(fastTopics)
library(Rtsne)
library(uwot)
library(ggplot2)
library(cowplot)
set.seed(1)

# Load data and results.
load("../data/droplet.RData")
fit <- readRDS("../output/droplet/rds/fit-droplet-scd-ex-k=7.rds")$fit
samples <- readRDS("../output/droplet/clustering-droplet.rds")
rows <- sort(c(sample(which(is.element(samples$cluster,
                                       c("B","C","B+C","H"))),1000),
               which(is.element(samples$cluster,c("Cil","T+N","U")))))
fit  <- poisson2multinom(fit)
fit  <- select(fit,loadings = rows)
fit  <- merge_topics(fit,c("k1","k3","k4","k5","k7"))
colnames(fit$L) <- c("k1","k2","k3")
p1 <- pca_plot(fit,k = 1) + guides(fill = "none")
p2 <- pca_plot(fit,k = 2) + guides(fill = "none")
p3 <- pca_plot(fit,k = 3) + guides(fill = "none")
plot_grid(p1,p2,p3,nrow = 1)

# Run t-SNE, then plot the 2-d embedding.
tsne <- Rtsne(fit$L,dims = 2,pca = FALSE,normalize = FALSE,perplexity = 30,
              theta = 0.1,max_iter = 1000,eta = 200,verbose = TRUE)
tsne$x <- tsne$Y
colnames(tsne$x) <- c("tsne1","tsne2")
p4 <- pca_plot(fit,out.pca = tsne,k = 1) + guides(fill = "none")
p5 <- pca_plot(fit,out.pca = tsne,k = 2) + guides(fill = "none")
p6 <- pca_plot(fit,out.pca = tsne,k = 3) + guides(fill = "none")

# Run UMAP, then plot the 2-d embedding.
out.umap <- umap(fit$L,n_neighbors = 30,metric = "euclidean",n_epochs = 1000,
                 min_dist = 0.1,scale = "none",learning_rate = 1,
                 verbose = TRUE)
out.umap <- list(x = out.umap)
colnames(out.umap$x) <- c("umap1","umap2")
p7 <- pca_plot(fit,out.pca = out.umap,k = 1) + guides(fill = "none")
p8 <- pca_plot(fit,out.pca = out.umap,k = 2) + guides(fill = "none")
p9 <- pca_plot(fit,out.pca = out.umap,k = 3) + guides(fill = "none")

# Arrange all nine plots.
print(plot_grid(p1,p2,p3,
                p4,p5,p6,
                p7,p8,p9))

