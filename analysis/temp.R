# Quick inspection of the topic modeling results from the PBMC purified data.
library(Matrix)
library(fastTopics)
library(ggplot2)
library(cowplot)
load("../data/pbmc_purified.RData")
fit <- readRDS("../output/pbmc-purified/rds/fit-pbmc-purified-scd-ex-k=6.rds")$fit

# Loadings plot.
p1 <- loadings_plot(poisson2multinom(fit),samples$celltype)

# PCA plot
clrs <- c("dodgerblue",  # CD19+ B
          "forestgreen", # CD14+ Monocyte
          "darkmagenta", # CD34+
          "yellowgreen", # CD4+ T Helper2
          "gray",        # CD56+ NK
          "tomato",      # CD8+ Cytotoxic T
          "orange",      # CD4+/CD45RO+ Memory
          "magenta",     # CD8+/CD45RA+ Naive Cytotoxic
          "limegreen",   # CD4+/CD45RA+/CD25- Naive T
          "gold")        # CD4+/CD25 T Reg
fit2 <- poisson2multinom(fit)
pca  <- prcomp(fit2$L)
pdat <- cbind(samples,pca$x)
ggplot(pdat,aes(x = PC1,y = PC2,fill = celltype)) +
  geom_point(shape = 21,color = "white",size = 1.5) +
  scale_fill_manual(values = clrs) +
  theme_cowplot(font_size = 10)
