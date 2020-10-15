# Try running flashier on the topic proportions matrix estimated from
# the pulse-seq data.
library(Matrix)
library(fastTopics)
library(flashier)
library(ggplot2)
library(cowplot)
set.seed(1)

# Load the multinomial topic model fit and clustering.
samples <- readRDS("../output/pulseseq/clustering-pulseseq.rds")
fit <- readRDS("../output/pulseseq/rds/fit-pulseseq-scd-ex-k=11.rds")$fit

# Fit flash model to topic proportions matrix.
k  <- 4
n  <- nrow(fit$L)
m  <- ncol(fit$L)
fl <- flash.init(poisson2multinom(fit)$L)
fl <- flash.init.factors(fl,list(matrix(rnorm(n*k),n,k),
                                 matrix(rnorm(m*k),m,k)),
                         prior.family = c(prior.unimodal(),prior.normal()))
fl <- flash.backfit(fl,verbose.lvl = 3,maxiter = 50)
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
               fill = samples$cluster) +
  scale_fill_manual(values = cluster_colors)
