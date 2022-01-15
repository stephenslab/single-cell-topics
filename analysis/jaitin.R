# See stephenslab.github.io/count-clustering/project/src/
# jaitin_structure_genes.html
# Download data from
# github.com/jhsiao999/singleCellRNASeqMouseJaitinSpleen
# Download topic model fit from
# github.com/stephenslab/count-clustering/tree/master/project/rdas
library(Matrix)
library(Biobase)
library(fastTopics)
library(ggplot2)
library(cowplot)
library(CountClust)
set.seed(1)
load("MouseJaitinSpleen.rda")
meta_data <- pData(MouseJaitinSpleen)
gene_names <- rownames(fData(MouseJaitinSpleen))
counts <- exprs(MouseJaitinSpleen)
clust <- readRDS("MouseJaitinSpleen-topicFit.rds")$clust_7
ENSG_genes_index <- grep("ERCC",gene_names,invert = TRUE)
counts <- counts[ENSG_genes_index,]
filter_genes <- c("M34473","abParts","M13680","Tmsb4x",
                  "S100a4","B2m","Atpase6","Rpl23","Rps18",
                  "Rpl13","Rps19","H2-Ab1","Rplp1","Rpl4",
                  "Rps26","EF437368") 
fcounts <- counts[-match(filter_genes,rownames(counts)),]
sample_counts <- colSums(fcounts)
i <- which(meta_data$number_of_cells == 1 & 
           meta_data$group_name == "CD11c+" & 
           sample_counts > 600)
meta_data <- meta_data[i,]
counts <- t(counts)
counts <- as(counts,"dgCMatrix")
F <- clust$theta
L <- clust$omega
j <- match(rownames(F),colnames(counts))
counts <- counts[i,j]
j <- which(colSums(counts) > 0)
F <- F[j,]
counts <- counts[,j]
fit0 <- init_poisson_nmf(counts,F = F,L = L)
fit1 <- fit_poisson_nmf(counts,fit0 = fit0,numiter = 100,method = "em",
                        control = list(nc = 4,extrapolate = FALSE))
fit2 <- fit_poisson_nmf(counts,fit0 = fit0,numiter = 100,method = "scd",
                        control = list(nc = 4,extrapolate = TRUE))
p0 <- plot_progress(list(em = fit1,scd = fit2),add.point.every = 10)
set.seed(1)
p1 <- structure_plot(fit2,perplexity = 70)
pca_plot(poisson2multinom(fit2),k = 3,pcs = 2:3)
pca_plot(poisson2multinom(fit2),k = 4,pcs = 3:4)
pca_plot(poisson2multinom(fit2),k = 7,pcs = 4:5)
de <- de_analysis(fit2,counts,pseudocount = 0.1,
                  control = list(nc = 4,nc = 1e4))
volcano_plot(de,k = 1)
volcano_plot(de,k = 2)
volcano_plot(de,k = 3)
volcano_plot(de,k = 4) # cDC1 cells
volcano_plot(de,k = 5) 
volcano_plot(de,k = 7) # macrophages
