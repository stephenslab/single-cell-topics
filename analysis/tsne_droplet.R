library(Matrix)
library(dplyr)
library(fastTopics)
library(ggplot2)
library(cowplot)
set.seed(1)

# Run t-SNE on the droplet data.
load("../data/droplet.RData")
samples <- readRDS("../output/droplet/clustering-droplet.rds")
fit  <- readRDS("../output/droplet/rds/fit-droplet-scd-ex-k=7.rds")$fit
p1   <- pca_plot(poisson2multinom(fit),pcs = 1:2,fill = samples$tissue)
rows <- with(samples,cluster == "B" | cluster == "C" | cluster == "B+C")
rows <- sort(c(sample(which(rows),500),which(!rows)))
fit  <- select(poisson2multinom(fit),loadings = rows)
tsne <- tsne_from_topics(fit,perplexity = 30,pca = TRUE)
tsne <- list(x = tsne$Y)
p2   <- pca_plot(fit,out.pca = tsne,fill = samples[rows,"tissue"])
plot_grid(p1,p2)
