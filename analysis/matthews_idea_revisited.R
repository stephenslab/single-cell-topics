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

# Combine all the T cell subtypes, except for the CD8+ T cells, into a
# single label.
celltype <- as.character(samples$celltype)
celltype[celltype == "CD4+ T Helper2" |
  celltype == "CD4+/CD45RO+ Memory" |
  celltype == "CD8+/CD45RA+ Naive Cytotoxic" |
  celltype == "CD4+/CD45RA+/CD25- Naive T" | 
  celltype == "CD4+/CD25 T Reg"] <- "T cell"
celltype <- factor(celltype)

# Create a scatterplot comparing, for each gene, the sample mean
# expression against the mean predicted expression, where the mean is
# taken over all samples in the selected cell type (k). Genes are
# labeled if their mean expression is greater than obsmin and the
# distance from the diagonal in the scatterplot (on the base-2
# logarithmic scale) is greater than betamin.
create_expression_scatterplot <- function (counts, fit, celltype, k,
                                           labels = NULL, e = 0.001,
                                           betamin = 1, obsmin = 0.1) {

  # For each gene, compute the sample mean expression (x) and the
  # predicted mean expression (y).
  i   <- which(celltype == k)
  x   <- colMeans(counts[i,]) + e
  y   <- with(fit,drop(tcrossprod(colMeans(L[i,]),F))) + e
  dat <- data.frame(x = x,y = y,label = labels,stringsAsFactors = FALSE)

  # Label genes for which the mean expression is greater than
  # "obsmin", and the distance from the diagonal on the (base-2)
  # log-scale is greater than "betamin".
  r <- log2(pmax(x/y,y/x))
  i <- which(r < betamin | x < obsmin)
  dat[i,"label"] <- ""

  # Create the scatterplot.
  return(ggplot(dat,aes_string(x = "x",y = "y",label = "label")) +
         geom_point(color = "dodgerblue",size = 1) +
         geom_abline(slope = 1,intercept = 0,color = "black",
                     linetype = "dotted") +
         geom_text_repel(color = "black",size = 2.2,fontface = "italic",
                         segment.color = "black",segment.size = 0.25,
                         max.overlaps = Inf,na.rm = TRUE) +
         scale_x_continuous(trans = "log10",breaks = 10^seq(-4,4)) +
         scale_y_continuous(trans = "log10",breaks = 10^seq(-4,4)) +
         labs(x = "observed",y = "predicted",title = k) +
         theme_cowplot(font_size = 10) +
         theme(plot.title = element_text(size = 10,face = "plain")))
}

# Create 6 expression scatterplots, one for each of the FACS
# subpopulations.
p1 <- create_expression_scatterplot(counts,fit,celltype,"CD19+ B",
                                    genes$symbol,betamin = 2,obsmin = 0.005)
p2 <- create_expression_scatterplot(counts,fit,celltype,"CD14+ Monocyte",
                                    genes$symbol,betamin = 2,obsmin = 0.01)
p3 <- create_expression_scatterplot(counts,fit,celltype,"CD34+",
                                    genes$symbol,betamin = 2,obsmin = 0.005)
p4 <- create_expression_scatterplot(counts,fit,celltype,"CD56+ NK",
                                    genes$symbol,betamin = 2,obsmin = 0.005)
p5 <- create_expression_scatterplot(counts,fit,celltype,"CD8+ Cytotoxic T",
                                    genes$symbol,betamin = 1.75,obsmin = 0.005)
p6 <- create_expression_scatterplot(counts,fit,celltype,"T cell",
                                    genes$symbol,betamin = 1.5,obsmin = 0.005)
plot_grid(p1,p2,p3,p4,p5,p6,nrow = 2,ncol = 3)
