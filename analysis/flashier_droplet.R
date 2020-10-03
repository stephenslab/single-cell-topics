# Try running flashier on the topic proportions matrix estimated from
# the droplet data.
library(Matrix)
library(fastTopics)
library(flashier)
library(ggplot2)
library(cowplot)
set.seed(1)

# Script parameters.
flash_prior <- c(prior.point.normal(),prior.point.normal())

# Load the multinomial topic model fit and clustering.
samples <- readRDS("../output/droplet/clustering-droplet.rds")
fit <- readRDS("../output/droplet/rds/fit-droplet-scd-ex-k=7.rds")$fit

# Fit flash model.
X  <- poisson2multinom(fit)$L
fl <- flash(X,greedy.Kmax = 6,nullcheck = FALSE,prior.family = flash_prior)
colnames(fl$loadings.pm[[1]]) <- paste0("d",1:6)

# Visualize the flash results.
cluster_colors <- c("royalblue",   # B
                    "forestgreen", # C
                    "slategray",   # B+C
                    "turquoise",   # H
                    "firebrick",   # Cil
                    "darkorange",  # T+N
                    "gold",        # G
                    "gainsboro")   # U

labeled_pca_plot <- function (pca, pcs = 1:2, labels) {
  dat <- cbind(data.frame(label = factor(labels)),as.data.frame(pca))
  return(ggplot(dat,aes_string(x = pcs[1],y = pcs[2],fill = "label")) +
         geom_point(shape = 21,color = "white",size = 1.2,na.rm = TRUE) +
         scale_fill_manual(values = cluster_colors) +
         theme_cowplot(font_size = 10))
}

# Plot the flash loadings, and layer the clusters on top of these
# plots.
Y  <- fl$loadings.pm[[1]]
p1 <- labeled_pca_plot(Y,pcs = c("d1","d2"),samples$cluster)
p2 <- labeled_pca_plot(Y,pcs = c("d1","d4"),samples$cluster)
p3 <- labeled_pca_plot(Y,pcs = c("d1","d3"),samples$cluster)
p4 <- labeled_pca_plot(Y,pcs = c("d2","d5"),samples$cluster)
print(plot_grid(p1,p2,p3,p4))
