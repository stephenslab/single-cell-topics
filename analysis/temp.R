library(Matrix)
library(fastTopics)
library(dplyr)
library(ggplot2)
library(cowplot)
load("../data/pbmc_68k.RData")
fit <- readRDS("../output/pbmc-68k/rds/fit-pbmc-68k-scd-ex-k=12.rds")$fit
facs_colors <- c("forestgreen",
                 "dodgerblue",
                 "darkmagenta",
                 "firebrick",
                 "gray",
                 "tomato",
                 "yellow",
                 "magenta",
                 "darkorange",
                 "gold",
                 "darkblue",
                 "greenyellow")
pcs <- 5:6
p1 <- pca_plot(poisson2multinom(fit),pcs = pcs,fill = samples$celltype) +
  scale_fill_manual(values = facs_colors)
p2 <- pca_hexbin_plot(poisson2multinom(fit),pcs = pcs)
plot_grid(p1,p2,rel_widths = c(4,3))

qplot(samples$tsne1,samples$tsne2,color = samples$celltype,size = I(0.75)) +
  theme_cowplot(font_size = 10)
