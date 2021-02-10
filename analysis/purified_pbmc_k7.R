library(Matrix)
library(fastTopics)
library(ggplot2)
library(cowplot)

# Load the count data and the K = 7 Poisson NMF model fit.
load("../data/pbmc_purified.RData")
fit <- readRDS(file.path("../output/pbmc-purified/rds",
                         "fit-pbmc-purified-scd-ex-k=7.rds"))$fit
fit <- poisson2multinom(fit)

# Create a Structure plot.
topic_colors <- c("darkmagenta","gray","forestgreen","dodgerblue",
                  "darkorange","lightskyblue","gold")
celltype <- as.character(samples$celltype)
celltype[celltype == "CD4+ T Helper2" |
  celltype == "CD4+/CD45RO+ Memory" |
  celltype == "CD8+/CD45RA+ Naive Cytotoxic" |
  celltype == "CD4+/CD45RA+/CD25- Naive T" | 
  celltype == "CD4+/CD25 T Reg"] <- "T cell"
celltype <- factor(celltype)
set.seed(1)
rows <- sort(c(sample(which(celltype == "CD19+ B"),500),
               sample(which(celltype == "CD14+ Monocyte"),250),
               sample(which(celltype == "CD34+"),500),
               sample(which(celltype == "CD56+ NK"),400),
               sample(which(celltype == "CD8+ Cytotoxic T"),400),
               sample(which(celltype == "T cell"),1000)))
p1 <- structure_plot(select_loadings(fit,loadings = rows),
                     grouping = celltype[rows],topics = 1:7,
                     colors = topic_colors,gap = 50,num_threads = 4)
print(p1)

# Perform a differential expression (DE) analysis using the topic
# model.
res <- diff_count_analysis(fit,counts)

# plot the results of the DE analysis.
#
# topic meaning
#   1   CD34+ cells
#   2   NK cells
#   3   CD14+ cells
#   4   B cells
#   5
#   6
#   7
p2 <- volcano_plot(res,"k5",genes$symbol,label_above_quantile = 0.995,
                   subsample_below_quantile = 0.5,filter_low_counts = 5e-5)
print(p2)
