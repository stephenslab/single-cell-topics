library(Matrix)
library(fastTopics)
library(flashier)
library(ggplot2)
library(cowplot)
set.seed(1)

# Script parameters.
flash_prior <- c(prior.point.normal(),prior.point.normal())
# flash_prior <- c(prior.point.normal(),prior.nonnegative())
flash_prior <- c(prior.nonnegative(),prior.point.normal())
# flash_prior <- c(prior.nonnegative(),prior.nonnegative())

# Load the multinomial topic model fit.
samples <- readRDS("../output/droplet/clustering-droplet.rds")
fit <- readRDS("../output/droplet/rds/fit-droplet-scd-ex-k=7.rds")$fit

# Fit flash model.
X  <- poisson2multinom(fit)$L
fl <- flash(X,greedy.Kmax = 6,nullcheck = FALSE,prior.family = flash_prior)
colnames(fl$loadings.pm[[1]]) <- paste0("k",1:6)

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

Y  <- fl$loadings.pm[[1]]
print(range(Y))
p1 <- labeled_pca_plot(Y,pcs = c("k1","k2"),samples$cluster)
p2 <- labeled_pca_plot(Y,pcs = c("k3","k4"),samples$cluster)
p3 <- labeled_pca_plot(Y,pcs = c("k5","k6"),samples$cluster)
print(plot_grid(p1,p2,p3,nrow = 1))
ggsave("plot.png",dpi = 150)
