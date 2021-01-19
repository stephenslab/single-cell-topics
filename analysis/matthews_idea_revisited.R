# Load the packages.
library(Matrix)
library(fastTopics)
library(ggplot2)
library(ggrepel)
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
celltypes <- levels(celltype)
for (j in 1:k) {
  i     <- which(celltype == celltypes[j])
  X[,j] <- colMeans(counts[i,])
  Y[,j] <- with(fit,tcrossprod(colMeans(L[i,]),F))
}
j <- 5

# TO DO: Explain what this function does, and how to use it.
create_expression_scatterplot <- function (x, y, labels = NULL,
                                           title = NULL, e = 0.001,
                                           label_minbeta = 1,
                                           label_minx = 0.1) {
  x   <- x + e
  y   <- y + e
  dat <- data.frame(x = x,y = y,label = labels,stringsAsFactors = FALSE)
  r   <- log2(pmax(x/y,y/x))
  i   <- which(r < label_minbeta | x < label_minx)
  dat[i,"label"] <- ""
  
  return(ggplot(dat,aes_string(x = "x",y = "y",label = "label")) +
         geom_point(color = "dodgerblue",size = 1) +
         geom_abline(slope = 1,intercept = 0,color = "black",
                     linetype = "dotted") +
         geom_text_repel(color = "black",size = 2.2,fontface = "italic",
                         segment.color = "black",segment.size = 0.25,
                         max.overlaps = Inf,na.rm = TRUE) +
         scale_x_continuous(trans = "log10") +
         scale_y_continuous(trans = "log10") +
         labs(x = "observed",y = "predicted",title = title) +
         theme_cowplot(font_size = 10) +
         theme(plot.title = element_text(size = 10,face = "plain")))
}

# TO DO: Explain what these lines of code do.
p1 <- create_expression_scatterplot(X[,1],Y[,1],genes$symbol,celltypes[1],
                                    label_minbeta = 2,label_minx = 0.01)
p2 <- create_expression_scatterplot(X[,2],Y[,2],genes$symbol,celltypes[2],
                                    label_minbeta = 2,label_minx = 0.005)
p3 <- create_expression_scatterplot(X[,3],Y[,3],genes$symbol,celltypes[3],
                                    label_minbeta = 2,label_minx = 0.005)
p4 <- create_expression_scatterplot(X[,4],Y[,4],genes$symbol,celltypes[4],
                                    label_minbeta = 2,label_minx = 0.002)
p5 <- create_expression_scatterplot(X[,5],Y[,5],genes$symbol,celltypes[5],
                                    label_minbeta = 1.75,label_minx = 0.002)
p6 <- create_expression_scatterplot(X[,6],Y[,6],genes$symbol,celltypes[6],
                                    label_minbeta = 1.5,label_minx = 0.0015)
plot_grid(p1,p2,p3,p4,p5,p6,nrow = 2,ncol = 3)

