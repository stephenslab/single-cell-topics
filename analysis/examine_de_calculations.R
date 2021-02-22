library(Matrix)
library(fastTopics)
library(ggplot2)
library(cowplot)
set.seed(1)

# Load the count data, the K = 6 topic model fit, and the 7 clusters
# identified in the clustering analysis ("clusters_purified_pbmc").
load("../data/pbmc_purified.RData")
fit <- readRDS(file.path("../output/pbmc-purified/rds",
                         "fit-pbmc-purified-scd-ex-k=6.rds"))$fit
fit <- poisson2multinom(fit)
samples <- readRDS("../output/pbmc-purified/clustering-pbmc-purified.rds")

# Perform differential expression analysis using the multinomial topic
# model, after removing the dendritic cells cluster.
rows <- which(samples$cluster != "dendritic")
counts <- counts[rows,]
fit <- select_loadings(fit,loadings = rows)
de1 <- diff_count_analysis(fit,counts,pseudocount = 1,fit.method = "em")
j   <- which(is.element(genes$symbol,c("CD74","CD79A","CD79B","GNLY","IL32")))
fit2 <- fit
fit2$F <- fit$F[j,]
de2 <- diff_count_analysis(fit2,counts[,j],s = fit2$s,pseudocount = 1,
                           fit.method = "em")

plot(de1$Z[j,],de2$Z,pch = 20)
abline(a = 0,b = 1,col = "skyblue",lty = "dotted")

p1 <- volcano_plot(de1,"k3",genes$symbol)
