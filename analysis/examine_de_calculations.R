# TO DO: Explain here what this script is for, and how to use it.
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

# Take a random subset of the cells.
n       <- nrow(counts)
rows    <- sample(n,9000)
samples <- samples[rows,]
counts  <- counts[rows,]
fit     <- select_loadings(fit,loadings = rows)

# Perform differential expression analysis using the multinomial topic
# model, after removing the dendritic cells cluster.
de <- diff_count_analysis(fit,counts,shrink.method = "none")

# Plot z-score vs. total UMI count.
i <- 5
x <- colSums(counts)
j <- which(x < 100)
pdat <- data.frame(count = x[j],
                   z     = abs(de$Z[j,i]),
                   beta  = de$beta[j,i])
p1 <- ggplot(pdat,aes(x = count,y = z,fill = beta)) +
  geom_point(shape = 21,color = "white") +
  scale_fill_gradient2(low = "lightskyblue",mid = "gold",high = "orangered",
                       midpoint = 0) +
  labs(x = "total UMI count",y = "|z-score|",fill = "logFC") +
  theme_cowplot(font_size = 12) 
print(p1)
