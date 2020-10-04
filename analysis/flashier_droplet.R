# Try running flashier on the topic proportions matrix estimated from
# the droplet data.
library(Matrix)
library(fastTopics)
library(flashier)
library(ggplot2)
library(cowplot)
set.seed(1)

# These two functions are used to create the plots below.
labeled_pca_plot <- function (pca, pcs = 1:2, labels)
  ggplot(cbind(data.frame(label = factor(labels)),as.data.frame(pca)),
         aes_string(x = pcs[1],y = pcs[2],fill = "label")) +
         geom_point(shape = 21,color = "white",size = 1.2,na.rm = TRUE) +
         scale_fill_manual(values = cluster_colors) +
         theme_cowplot(font_size = 9)

pca_hexbin_plot <-
  function (pca, pcs = 1:2, n = 40, bins = c(0,1,10,100,1000,Inf),
            colors = c("gainsboro","lightskyblue","gold","orange","magenta"))
  ggplot(as.data.frame(pca),aes_string(x = pcs[1],y = pcs[2])) +
         stat_bin_hex(mapping = aes_q(fill = quote(cut(..count..,bins))),
                      bins = n) +
         scale_fill_manual(values = colors) +
         theme_cowplot(font_size = 8)

# Load the multinomial topic model fit and clustering.
samples <- readRDS("../output/droplet/clustering-droplet.rds")
fit <- readRDS("../output/droplet/rds/fit-droplet-scd-ex-k=7.rds")$fit

# Fit flash model to topic proportions matrix.
fl <- flash(poisson2multinom(fit)$L,greedy.Kmax = 6,nullcheck = FALSE,
            prior.family = c(prior.point.normal(),prior.point.normal()))
colnames(fl$loadings.pm[[1]]) <- paste0("d",1:6)

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
p1 <- labeled_pca_plot(Y,pcs = c("d1","d2"),samples$cluster)
p2 <- labeled_pca_plot(Y,pcs = c("d1","d4"),samples$cluster)
p3 <- labeled_pca_plot(Y,pcs = c("d1","d3"),samples$cluster)
p4 <- labeled_pca_plot(Y,pcs = c("d2","d5"),samples$cluster)
p5 <- pca_hexbin_plot(Y,pcs = c("d1","d2"))
p6 <- pca_hexbin_plot(Y,pcs = c("d1","d4"))
p7 <- pca_hexbin_plot(Y,pcs = c("d1","d3"))
p8 <- pca_hexbin_plot(Y,pcs = c("d2","d5"))
print(plot_grid(p1,p2,p3,p4,p5,p6,p7,p8,nrow = 2,ncol = 4))
