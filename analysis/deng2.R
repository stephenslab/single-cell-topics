library(Biobase)
library(fastTopics)
library(ggplot2)
library(cowplot)
library(CountClust)
set.seed(1)
load("Deng2014MouseEsc.rda")
load("deng_topic_fit.rda")
counts <- t(Biobase::exprs(Deng2014MouseESC))
F <- Topic_clus_list$clust_7$theta
L <- Topic_clus_list$clust_7$omega
i <- which(colSums(counts) > 0)
F <- F[i,]
counts <- counts[,i]
annotations <- factor(rownames(L),
                      c("zy","early2cell","mid2cell","late2cell","4cell",
                        "8cell","16cell","earlyblast","midblast","lateblast"))
levels(annotations)[1] <- "zygote"
fit0 <- init_poisson_nmf(counts,F = F,L = L)
fit1 <- fit_poisson_nmf(counts,fit0 = fit0,numiter = 40,method = "em",
                        control = list(nc = 4,extrapolate = FALSE))
set.seed(1)
fit2 <- fit_poisson_nmf(counts,k = 8,numiter = 100,method = "scd",
                        control = list(nc = 4,extrapolate = TRUE))
p0 <- plot_progress(list(em = fit1,scd = fit2),add.point.every = 10,e = 0.01)
p0 <- plot_progress(fit2,add.point.every = 5,e = 0.01)
stop()
set.seed(1)
f <- function (fit, ...)
  drop(pca_from_topics(fit,dims = 1,...))
p1 <- structure_plot(fit1,grouping = annotations,perplexity = 30,gap = 5,
                     embed_method = f)
p2 <- structure_plot(fit2,grouping = annotations,perplexity = 10,gap = 5,
                     embed_method = f)
plot_grid(p1,p2,nrow = 2,ncol = 1)
de <- de_analysis(fit2,counts,pseudocount = 0.1,
                  control = list(nc = 4,nc = 1e4))

# See Fiig. 7 of Guo et al (2010)
volcano_plot(de,k = 5,ymax = 500) # primitive endoderm (Cdx2, Gata3)
volcano_plot(de,k = 7,ymax = 500) # Tdgf1 = cripto, Sox2, Esrrb, Klf4 (inner cells)
volcano_plot(de,k = 8,ymax = 1000) # Id2, Mbnl3 (outer cells)
