library(Matrix)
library(dplyr)
library(fastTopics)
library(ggplot2)
library(cowplot)
source("../code/plots.R")

# Load the data and results.
load("../data/pbmc_purified.RData")
load("../output/pbmc-purified/diff-count-pbmc-purified.RData")
fit <- readRDS(file.path("../output/pbmc-purified/rds",
                         "fit-pbmc-purified-scd-ex-k=6.rds"))$fit
samples <- readRDS("../output/pbmc-purified/clustering-pbmc-purified.rds")

# Volcano plot for B-cells.
p1 <- volcano_plot(diff_count_topics,k = 3,
                   labels = genes$symbol,
                   label_above_quantile = 0.995)

# Volcano plots for topics 1 and 6 in T-cells.
p2 <- volcano_plot(diff_count_t,k = 1,
                   labels = genes$symbol,
                   label_above_quantile = 0.995) +
  scale_fill_gradient2(low = "deepskyblue",mid = "gold",high = "orangered",
                       na.value = "gainsboro",midpoint = -2)
p3 <- volcano_plot(diff_count_t,k = 6,
                   labels = genes$symbol,
                   label_above_quantile = 0.99) +
  scale_fill_gradient2(low = "deepskyblue",mid = "gold",high = "orangered",
                       na.value = "gainsboro",midpoint = -2)
plot_grid(p2,p3)

loadings_plot(poisson2multinom(fit),x = samples$celltype,k = 1)
