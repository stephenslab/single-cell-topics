library(Matrix)
library(fastTopics)
set.seed(1)

# Load the droplet data and K = 7 Poisson NMF model fit.
load("../data/droplet.RData")
fit <- readRDS("../output/droplet/rds/fit-droplet-scd-ex-k=7.rds")$fit
fit <- poisson2multinom(fit)

# Run k-means on the mixture proportions matrix.
out <- kmeans(fit$L,12)

# Create a Structure plot.
topic_colors <- c("gold","royalblue","salmon","turquoise","olivedrab",
                  "firebrick","forestgreen")
p1 <- structure_plot(fit,grouping = factor(out$cluster),topics = 1:7,
                     colors = topic_colors,gap = 10,num_threads = 4)
