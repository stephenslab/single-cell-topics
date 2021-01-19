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
fit <- poisson2multinom(fit)

# Combine all the T cell subtypes, except for the CD8+ T cells, into a
# single label.
celltype <- as.character(samples$celltype)
celltype[celltype == "CD4+ T Helper2" |
  celltype == "CD4+/CD45RO+ Memory" |
  celltype == "CD8+/CD45RA+ Naive Cytotoxic" |
  celltype == "CD4+/CD45RA+/CD25- Naive T" | 
  celltype == "CD4+/CD25 T Reg"] <- "T cell"
celltype <- factor(celltype)

# Create a scatterplot comparing, for each gene, the
# maximum-likelihood estimate (MLE) for a simple multinomial model of
# expression against the mean multinomial probability under a topic
# model, in which the mean is taken over all samples in the selected
# cell type (k). Genes are labeled if their multinomial MLE is greater
# than mle.min and the distance from the diagonal in the scatterplot
# (on the base-2 logarithmic scale) is greater than beta.min.
create_expression_scatterplot <- function (counts, fit, celltype, k,
                                           labels = NULL, e = 1e-6,
                                           beta.min = 2, mle.min = 1e-6) {

  # For each gene, and for each cell, compute the multinomial
  # probability under the given multinomial topic model, then take the
  # average over all cells.
  i <- which(celltype == k)
  x <- drop(with(fit,tcrossprod(colMeans(L[i,]),F))) + e
  
  # Compute the maximum-likelihood estimate (MLE) for a basic
  # multinomial model of gene expression for each cell.
  y <- colSums(counts[i,])/sum(counts[i,]) + e

  # Create a data frame containing the data to be plotted.
  dat <- data.frame(x = x,y = y,label = labels,stringsAsFactors = FALSE)

  # Label genes for which the MLE is greater than mle.min, and the
  # distance from the diagonal on the (base-2) log-scale is greater
  # than beta.min.
  r <- log2(pmax(x/y,y/x))
  i <- which(r < beta.min | y < mle.min)
  dat[i,"label"] <- ""
  
  # Create the scatterplot.
  return(ggplot(dat,aes_string(x = "x",y = "y",label = "label")) +
         geom_point(color = "dodgerblue",size = 1) +
         geom_abline(slope = 1,intercept = 0,color = "black",
                     linetype = "dotted") +
         geom_text_repel(color = "black",size = 2.2,fontface = "italic",
                         segment.color = "black",segment.size = 0.25,
                         max.overlaps = Inf,na.rm = TRUE) +
         scale_x_continuous(trans = "log10",breaks = 10^seq(-8,0)) +
         scale_y_continuous(trans = "log10",breaks = 10^seq(-8,0)) +
         labs(x = "predicted",y = "observed",title = k) +
         theme_cowplot(font_size = 10) +
         theme(plot.title = element_text(size = 10,face = "plain")))
}

# Create six expression scatterplots, one for each of the FACS
# subpopulations.
p1 <- create_expression_scatterplot(counts,fit,celltype,"CD19+ B",
                                    genes$symbol,beta.min = 1.25,
                                    mle.min = 4e-6)
p2 <- create_expression_scatterplot(counts,fit,celltype,"CD14+ Monocyte",
                                    genes$symbol,beta.min = 1.5,
                                    mle.min = 1e-5)
p3 <- create_expression_scatterplot(counts,fit,celltype,"CD34+",
                                    genes$symbol,beta.min = 1.5,mle.min = 1e-5)
p4 <- create_expression_scatterplot(counts,fit,celltype,"CD56+ NK",
                                    genes$symbol,beta.min = 1.25,
                                    mle.min = 4e-6)
p5 <- create_expression_scatterplot(counts,fit,celltype,"CD8+ Cytotoxic T",
                                    genes$symbol,beta.min = 1.25,
                                    mle.min = 1e-5)
p6 <- create_expression_scatterplot(counts,fit,celltype,"T cell",
                                    genes$symbol,beta.min = 1.3,
                                    mle.min = 4e-6)
print(plot_grid(p1,p2,p3,p4,p5,p6,nrow = 2,ncol = 3))
