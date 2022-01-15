# See stephenslab.github.io/count-clustering/project/src/
# deng_cluster_annotations.html
# Download .rda files from
# github.com/stephenslab/count-clustering/tree/master/project/rdas
library(Biobase)
library(fastTopics)
library(ggplot2)
library(cowplot)
library(CountClust)
set.seed(1)
load("Deng2014MouseEsc.rda")
load("deng_topic_fit.rda")
counts <- t(Biobase::exprs(Deng2014MouseESC))
F <- Topic_clus_list$clust_6$theta
L <- Topic_clus_list$clust_6$omega
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
fit2 <- fit_poisson_nmf(counts,fit0 = fit0,numiter = 40,method = "scd",
                        control = list(nc = 4,extrapolate = TRUE))
p1 <- plot_progress(list(em = fit1,scd = fit2),add.point.every = 5)
p2 <- structure_plot(fit1,grouping = annotations,perplexity = 30,gap = 5)
set.seed(1)
f <- function (fit, ...)
  drop(pca_from_topics(fit,dims = 1,...))
p3 <- structure_plot(fit2,grouping = annotations,perplexity = 10,gap = 5,
                     embed_method = f)
plot_grid(p2,p3,nrow = 2,ncol = 1)
de <- de_analysis(fit2,counts,pseudocount = 0.1,
                  control = list(nc = 4,nc = 1e4))
top_features <- ExtractTopFeatures(F,top_features = 100,
                                   method = "poisson",options = "min")
kl <- fastTopics:::min_kl_poisson(F)
k <- 1
plot(kl[top_features$indices[k,],k],top_features$scores[k,],
     pch = 20,log = "xy")
dat <- de
dat$lfsr[,1] <- 0
dat$lfsr[top_features$indices[k,],k] <- 0.01
volcano_plot(dat,k = k,ymax = 200)

