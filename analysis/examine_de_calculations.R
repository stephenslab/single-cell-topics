# TO DO: Explain here what this script is for, and how to use it.
library(Matrix)
library(fastTopics)
library(ggplot2)
library(cowplot)

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
fit_no_dendritic <- select_loadings(fit,loadings = rows)
de <- diff_count_analysis(fit_no_dendritic,counts[rows,],
                          shrink.method = "none")

i <- 5
x <- colSums(counts)
j <- which(x < 100)
pdat <- data.frame(count = x[j],
                   z     = abs(de$Z[j,i]),
                   beta  = de$beta[j,i])
ggplot(pdat,aes(x = count,y = z,fill = beta)) +
  geom_point(shape = 21,color = "white") +
  scale_fill_gradient2(low = "lightskyblue",mid = "gold",high = "orangered",
                       midpoint = 0) +
  labs(x = "total UMI count",y = "|z-score|",fill = "logFC") +
  theme_cowplot(font_size = 12) 
