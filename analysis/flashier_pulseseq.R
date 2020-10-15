# Try running flashier on the topic proportions matrix estimated from
# the pulse-seq data.
library(Matrix)
library(dplyr)
library(fastTopics)
library(flashier)
library(ggplot2)
library(cowplot)
set.seed(6)

# Load the multinomial topic model fit and clustering.
samples <- readRDS("../output/pulseseq/clustering-pulseseq.rds")
fit <- readRDS("../output/pulseseq/rds/fit-pulseseq-scd-ex-k=11.rds")$fit
rows <-
  sort(c(sample(which(is.element(samples$cluster,c("B","C","P"))),2000),
         sample(which(samples$cluster == "Cil"),200),
         sample(which(samples$cluster == "T+N"),200),
         which(samples$cluster == "I")))
fit <- select(poisson2multinom(fit),loadings = rows)

# Fit flash model to topic proportions matrix.
k  <- 2
fl <- flash(fit$L,greedy.Kmax = k,backfit = TRUE,nullcheck = FALSE,
            prior.family = c(prior.unimodal(),prior.normal()),
            verbose.lvl = 3)
colnames(fl$loadings.pm[[1]]) <- paste0("d",1:k)

# Plot the flash loadings, and layer the clusters on top of these
# plots.
cluster_colors  <- c("royalblue",   # B
                     "forestgreen", # C
                     "peru",        # P
                     "firebrick",   # Cil
                     "darkorange",  # T+N
                     "darkmagenta", # I
                     "gainsboro")   # U
Y  <- fl$loadings.pm[[1]]
p1 <- pca_plot(fit,pcs = c("d1","d2"),out.pca = list(x = Y),
               fill = samples[rows,"cluster"]) +
  scale_fill_manual(values = cluster_colors)
