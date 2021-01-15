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
fit <- poisson2multinom(fit)
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
celltype <- factor(celltype,c("CD19+ B","CD14+ Monocyte","CD34+",
                              "CD56+ NK","CD8+ Cytotoxic T","T cell"))
rows <- sort(c(sample(which(celltype == "CD19+ B"),500),
               sample(which(celltype == "CD14+ Monocyte"),250),
               sample(which(celltype == "CD34+"),500),
               sample(which(celltype == "CD56+ NK"),400),
			   sample(which(celltype == "CD8+ Cytotoxic T"),400),
               sample(which(celltype == "T cell"),1000)))
p1 <- structure_plot(select_loadings(fit,loadings = rows),
                     grouping = celltype[rows],
                     topics = topics,colors = topic_colors[topics],
                     perplexity = c(70,30,30,30,30,70),n = Inf,gap = 30,
                     num_threads = 4,verbose = FALSE)

# Compute the likelihoods under this topic model.
loglik <- loglik_multinom_topic_model(counts,fit)

# Compute the likelihoods in which the underlying pattern of
# expression is estimated from all the cells in the FACS
# subpopulation.
fit_facs <- init_poisson_nmf_from_clustering(counts,samples$celltype)
fit_facs <- poisson2multinom(fit_facs)
loglik_facs <- loglik_multinom_topic_model(counts,fit_facs)

# The improvement in fit is greatest for cells in the CD34+ and T cell
# FACS subpopulations.
pdat <- data.frame(x = loglik_facs,y = loglik,celltype = samples$celltype)
p2 <- ggplot(pdat,aes(x = x,y = y)) +
  geom_point(shape = 21,color = "white",fill = "royalblue") +
  geom_abline(intercept = 0,slope = 1,linetype = "dotted") +
  facet_wrap(vars(celltype)) +
  scale_x_continuous(limits = c(-10000,0),breaks = seq(-10000,0,2500)) +
  scale_y_continuous(limits = c(-10000,0),breaks = seq(-10000,0,2500)) +
  labs(x = "cluster model",y = "topic model",fill = "FACS subpopulation") +
  theme_cowplot(font_size = 8)

# Next, let's examine the improvement in fit after splitting the
# a topic by FACS cell population.
# celltypes_to_split <- levels(celltype)
# k <- 1
celltypes_to_split <- c("CD56+ NK","CD8+ Cytotoxic T")
k <- 4
K <- ncol(fit$L)
n <- nrow(counts)
m <- length(celltypes_to_split) - 1
fit2 <- fit
L <- fit2$L
F <- fit2$F
L <- cbind(L,matrix(0,nrow(L),m))
F <- cbind(F,matrix(0,nrow(F),m))
colnames(L) <- paste0("k",1:ncol(L))
colnames(F) <- paste0("k",1:ncol(F))
for (j in 1:m) {
  i <- which(celltype == celltypes_to_split[j])
  F[,K+j]  <- F[,k]
  L[i,K+j] <- L[i,k]
  L[i,k]  <- 1e-15
}    
fit2$L  <- L
fit2$Ln <- L
fit2$Ly <- L
fit2$F  <- F
fit2$Fn <- F
fit2$Fy <- F
fit2 <- multinom2poisson(fit2)
fit2 <- fit_poisson_nmf(counts,fit0 = fit2,numiter = 20,update.loadings = NULL)

# Create a Structure plot in which the cells are grouped by FACS
# cell population.
fit2 <- poisson2multinom(fit2)
p3 <- structure_plot(select_loadings(fit2,loadings = rows),
                     grouping = celltype[rows],
                     perplexity = c(70,30,30,30,30,70),n = Inf,gap = 30,
                     num_threads = 4,verbose = FALSE)

# Compute the likelihoods under this new topic model.
loglik2 <- loglik_multinom_topic_model(counts,fit2)

# TO DO: Explain here what these lines of code do.
pdat <- data.frame(x = loglik,y = loglik2,celltype = samples$celltype)
p4 <- ggplot(pdat,aes(x = x,y = y)) +
  geom_point(shape = 21,color = "white",fill = "royalblue") +
  geom_abline(intercept = 0,slope = 1,linetype = "dotted") +
  facet_wrap(vars(celltype)) +
  scale_x_continuous(limits = c(-10000,0),breaks = seq(-10000,0,2500)) +
  scale_y_continuous(limits = c(-10000,0),breaks = seq(-10000,0,2500)) +
  labs(x = "topic model",y = "split topic model",fill = "FACS subpopulation") +
  theme_cowplot(font_size = 8)
