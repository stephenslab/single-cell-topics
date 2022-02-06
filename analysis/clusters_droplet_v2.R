library(Matrix)
library(fastTopics)
library(ggplot2)
library(cowplot)
load("../data/droplet.RData")
fit <- readRDS("../output/droplet/rds/fit-droplet-scd-ex-k=7.rds")$fit
fit <- poisson2multinom(fit)
pca <- prcomp(fit$L)$x
x   <- rep("rare",nrow(pca))
pc2 <- pca[,2]
pc6 <- pca[,6]
x[pc2 > -0.15] <- "common"
x[pc6 < -0.05] <- "rare"
samples$cluster <- factor(x)
pca_plot(fit,pcs = 1:2,fill = samples$cluster)
set.seed(3)
topic_colors <- c("gold","royalblue","salmon","turquoise","olivedrab",
                  "firebrick","forestgreen")
topics <- c(3,4,5,1,7,2,6)
rows <- sort(c(sample(which(samples$cluster == "common"),1200),
               which(samples$cluster == "rare")))
p <- structure_plot(select_loadings(fit,loadings = rows),
                    grouping = samples[rows,"cluster"],
                    topics = topics,colors = topic_colors,
                    perplexity = 70,n = Inf,gap = 20,
                    num_threads = 4,verbose = FALSE)
print(p)
