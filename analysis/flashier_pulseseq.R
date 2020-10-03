# Try running flashier on the topic proportions matrix estimated from
# the pulse-seq data.
library(Matrix)
library(fastTopics)
library(flashier)
library(ggplot2)
library(cowplot)
set.seed(1)

# Script parameters.
flash_prior <- c(prior.point.normal(),prior.point.normal())

# Load the multinomial topic model fit and clustering.
samples <- readRDS("../output/pulseseq/clustering-pulseseq.rds")
fit <- readRDS("../output/pulseseq/rds/fit-pulseseq-scd-ex-k=11.rds")$fit

# Fit flash model.
X  <- poisson2multinom(fit)$L
fl <- flash(X,greedy.Kmax = 10,nullcheck = FALSE,prior.family = flash_prior)
colnames(fl$loadings.pm[[1]]) <- paste0("d",1:10)

# Visualize the flash results.
cluster_colors  <- c("royalblue",   # B
                     "forestgreen", # C
                     "peru",        # P
                     "firebrick",   # Cil
                     "darkorange",  # T+N
                     "darkmagenta", # I
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
p2 <- labeled_pca_plot(Y,pcs = c("d3","d4"),samples$cluster)
p3 <- labeled_pca_plot(Y,pcs = c("d2","d5"),samples$cluster)
print(plot_grid(p1,p2,p3,nrow = 1))

pca_hexbin_plot <-
  function (pca, pcs = 1:2, n = 40, bins = c(0,1,10,100,1000,Inf),
            colors = c("gainsboro","lightskyblue","gold","orange","magenta")) {
  dat <- as.data.frame(pca)
  return(ggplot(dat,aes_string(x = pcs[1],y = pcs[2])) +
         stat_bin_hex(mapping = aes_q(fill = quote(cut(..count..,bins))),
                      bins = n) +
         scale_fill_manual(values = colors) +
         theme_cowplot(font_size = 10))
}

p4 <- pca_hexbin_plot(Y,pcs = c("d1","d2"))
p5 <- pca_hexbin_plot(Y,pcs = c("d3","d4"))
p6 <- pca_hexbin_plot(Y,pcs = c("d2","d5"))
print(plot_grid(p4,p5,p6,nrow = 1))
