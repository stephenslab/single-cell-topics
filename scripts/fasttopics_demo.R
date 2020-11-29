# I will use this script to create the fastTopics single-cell RNA-seq
# vignette.
library(Matrix)
library(fastTopics)
library(ggplot2)
library(cowplot)
set.seed(1)

# Load the data.
load("pbmc_4k.RData")

# Remove genes with no expression.
j      <- which(colSums(counts > 0) > 0)
genes  <- genes[j,]
counts <- counts[,j]

# Fit the Poisson NMF model by running 100 EM updates.
fit0 <- fit_poisson_nmf(counts,k = 6,numiter = 100,method = "em",
                           control = list(numiter = 4,nc = 4))

# Refine the fit by running 100 extrapolated SCD updates.
fit <- fit_poisson_nmf(counts,fit0 = fit0,numiter = 100,method = "scd",
                       control = list(extrapolate = TRUE,numiter = 4,nc = 4))

# Plot improvement in solution over time.
p1 <- plot_progress_poisson_nmf(fit,x = "iter")
print(p1)

# Recover the multinomial topic model.
fit_multinom <- poisson2multinom(fit)

# Plot PCs of mixture proportions.
cluster_colors <- c("dodgerblue",  # B-cells
                    "forestgreen", # CD14+
                    "darkmagenta", # CD34+ 
                    "red",         # CD8+
                    "skyblue",     # dendritic
                    "gray",        # NK
                    "darkorange")  # T-cells
p2 <- pca_plot(fit_multinom,pcs = 1:2,fill = samples$cluster) +
  scale_fill_manual(values = cluster_colors)
p3 <- pca_plot(fit_multinom,pcs = 3:4,fill = samples$cluster) +
  scale_fill_manual(values = cluster_colors)
p4 <- pca_hexbin_plot(fit_multinom,pcs = 1:2)
p5 <- pca_hexbin_plot(fit_multinom,pcs = 3:4)
plot_grid(p2,p3,p4,p5)
