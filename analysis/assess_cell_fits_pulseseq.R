library(Matrix)
library(dplyr)
library(fastTopics)
library(ggplot2)
library(cowplot)

colors <- c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#ffff33",
            "#a65628","#f781bf","#999999")

# Load data and results.
load("../data/pulseseq.RData")
samples <- readRDS("../output/pulseseq/clustering-pulseseq.rds")
fit <- readRDS("../output/pulseseq/rds/fit-pulseseq-scd-ex-k=11.rds")$fit
fit_cluster <- init_poisson_nmf_from_clustering(counts,samples$cluster)
samples$loglik1 <- loglik_poisson_nmf(counts,fit_cluster)
samples$loglik2 <- loglik_poisson_nmf(counts,fit)
p1 <- ggplot(subset(samples,loglik1 > -25000),
             aes(x = loglik1,y = loglik2,fill = cluster)) +
  geom_point(shape = 21,color = "white") +
  geom_abline(intercept = 0,slope = 1,linetype = "dotted",color = "gray") +
  scale_fill_manual(values = colors) +
  labs(x = "clusters",y = "topics") +
  theme_cowplot(12) 
