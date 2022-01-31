library(Matrix)
library(fastTopics)
library(ggplot2)
library(cowplot)
load("../data/pulseseq.RData")
fit <- readRDS("../output/pulseseq/rds/fit-pulseseq-scd-ex-k=9.rds")$fit
fit <- poisson2multinom(fit)
x <- as.character(samples$tissue)
x[x == "club (hillock-associated)"] <- "club"
x[x == "goblet.1" | x == "goblet.2" | x == "goblet.progenitor"] <- "goblet"
x[x == "tuft.1" | x == "tuft.2" | x == "tuft.progenitor"] <- "tuft"
samples$tissue <- factor(x)

# Define clusters from PCs of topic proportions.
set.seed(1)
p1 <- pca_plot(fit,pcs = 3:4,fill = samples$tissue)
p2 <- pca_hexbin_plot(fit,pcs = 3:4)
plot_grid(p1,p2)

set.seed(1)
pca <- prcomp(fit$L)$x
x   <- rep("A",nrow(pca))
pc3 <- pca[,3]
pc4 <- pca[,4]
x[1.5*pc3 - pc4 < -0.125] <- "B"
pca_plot(fit,pcs = 3:4,fill = factor(x))

set.seed(1)
rows <- which(x == "A")
fit2 <- select_loadings(fit,loadings = rows)
p1 <- pca_plot(fit2,pcs = 1:2,fill = samples[rows,"tissue"])
p2 <- pca_hexbin_plot(fit2,pcs = 1:2)
plot_grid(p1,p2)

set.seed(1)
rows <- which(x == "B")
fit2 <- select_loadings(fit,loadings = rows)
p1 <- pca_plot(fit2,pcs = 1:2,fill = samples[rows,"tissue"])
p2 <- pca_hexbin_plot(fit2,pcs = 1:2)
plot_grid(p1,p2)

samples$cluster <- factor(x)
