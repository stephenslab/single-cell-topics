# Load the packages.
library(Matrix)
library(fastTopics)
library(ggplot2)
library(cowplot)

# Set the seed.
set.seed(1)

# Load the data and results needed for this analysis.
load("../data/pbmc_purified.RData")
fit <- readRDS(file.path("../output/pbmc-purified/rds",
                         "fit-pbmc-purified-scd-ex-k=6.rds"))$fit
samples <- readRDS("../output/pbmc-purified/clustering-pbmc-purified.rds")

# Create a Structure plot in which the cells are grouped by FACS
# cell population.
topic_colors <- c("gold","forestgreen","dodgerblue","gray",
                  "darkmagenta","violet")
topics <- c(5,3,2,4,1,6)
celltype <- as.character(samples$celltype)
celltype[celltype == "CD4+ T Helper2" |
  celltype == "CD4+/CD45RO+ Memory" |
  celltype == "CD8+/CD45RA+ Naive Cytotoxic" |
  celltype == "CD4+/CD45RA+/CD25- Naive T" | 
  celltype == "CD4+/CD25 T Reg"] <- "T cell"
celltype <- factor(celltype)

# For each cell type, calculate the average observed expression (X)
# and the average predicted expression (Y).
m <- ncol(counts)
k <- nlevels(celltype)
X <- matrix(0,m,k)
Y <- matrix(0,m,k)
for (j in 1:k) {
  i     <- which(celltype == levels(celltype)[j])
  X[,j] <- colMeans(counts[i,])
  Y[,j] <- with(fit,tcrossprod(colMeans(L[i,]),F))
}
j <- 5
pdat <- data.frame(x = X[,j] + 0.001,y = Y[,j] + 0.001)
ggplot(pdat,aes(x = x,y = y)) +
  geom_point(color = "royalblue",size = 1) +
  geom_abline(slope = 1,intercept = 0,color = "black",linetype = "dotted") +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10") +
  labs(x = "observed",y = "predicted",title = levels(celltype)[j]) +
  theme_cowplot(font_size = 10)
