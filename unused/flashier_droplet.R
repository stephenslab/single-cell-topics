# Try running flashier on the topic proportions matrix estimated from
# the droplet data.
library(Matrix)
library(dplyr)
library(fastTopics)
library(flashier)
library(ggplot2)
library(cowplot)
set.seed(1)

# Load the multinomial topic model fit and clustering.
samples <- readRDS("../output/droplet/clustering-droplet.rds")
fit <- readRDS("../output/droplet/rds/fit-droplet-scd-ex-k=7.rds")$fit
fit <- poisson2multinom(fit)

# Fit flash model to topic proportions matrix.
k <- 2
fl <- flash(fit$L,greedy.Kmax = k,backfit = TRUE,nullcheck = FALSE,
            prior.family = c(prior.unimodal(),prior.normal()),
            verbose.lvl = 3)
colnames(fl$loadings.pm[[1]]) <- paste0("d",1:k)

# Plot the flash loadings, and layer the clusters on top of these
# plots.
cluster_colors <- c("royalblue",   # B
                    "forestgreen", # C
                    "slategray",   # B+C
                    "turquoise",   # H
                    "firebrick",   # Cil
                    "darkorange",  # T+N
                    "gold",        # G
                    "gainsboro")   # U
Y  <- fl$loadings.pm[[1]]
p1 <- pca_plot(fit,pcs = c("d1","d2"),out.pca = list(x = Y),
               fill = samples$cluster) +
  scale_fill_manual(values = cluster_colors)
p2 <- pca_hexbin_plot(fit,pcs = c(1,2),out.pca = list(x = Y))
print(plot_grid(p1,p2,rel_widths = c(1.1,1)))
