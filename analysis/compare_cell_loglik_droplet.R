library(Matrix)
library(dplyr)
library(fastTopics)
library(ggplot2)
library(cowplot)

colors = c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2",
           "#D55E00","#CC79A7")

# Load data and results.
load("../data/droplet.RData")
samples <- readRDS("../output/droplet/clustering-droplet.rds")
fit <- readRDS("../output/droplet/rds/fit-droplet-scd-ex-k=7.rds")$fit
fit_merge   <- merge_topics(poisson2multinom(fit),c("k5","k7"))
fit_cluster <- init_poisson_nmf_from_clustering(counts,samples$cluster)
fit_montoro <- init_poisson_nmf_from_clustering(counts,samples$tissue)
fit_cluster <- poisson2multinom(fit_cluster)
fit_montoro <- poisson2multinom(fit_montoro)

# Compute the likelihood of each cell under the two hard clusterings.
samples$loglik1 <- loglik_multinom_topic_model(counts,fit_montoro)
samples$loglik2 <- loglik_multinom_topic_model(counts,fit_cluster)

# TO DO: Explain here what this function does, and how to use it.
create_loglik_scatterplot <- function (clusters, k, color = "black",
                                       minloglik = -Inf) {
  dat <- samples
  dat$cluster <- clusters
  dat    <- subset(dat,loglik1 > minloglik & loglik2 > minloglik)
  limits <- range(c(dat$loglik1,dat$loglik2))
  dat    <- subset(dat,cluster == k)
  return(ggplot(dat,aes(x = loglik1,y = loglik2)) +
         geom_point(shape = 21,color = "white",fill = color) +
         geom_abline(intercept = 0,slope = 1,linetype = "dotted",
                     color = "black") +
         xlim(limits) +
         ylim(limits) +
         labs(x = "clusters",y = "topics",
              title = sprintf("%s (n = %d)",k,nrow(dat))) +
         theme_cowplot(font_size = 10) +
         theme(plot.title = element_text(size = 10,face = "plain")))
}

minloglik <- -20000
p1 <- create_loglik_scatterplot(samples$tissue,"Basal","royalblue",-20000)
p2 <- create_loglik_scatterplot(samples$tissue,"Club","forestgreen",-20000)
p3 <- create_loglik_scatterplot(samples$tissue,"Ciliated","firebrick",-20000)
p4 <- create_loglik_scatterplot(samples$tissue,"Neuroendocrine","darkorange",
                                -20000)
p5 <- create_loglik_scatterplot(samples$tissue,"Tuft","dodgerblue",-20000)
p6 <- create_loglik_scatterplot(samples$tissue,"Goblet","gold",-20000)
plot_grid(p1,p2,p3,
          p4,p5,p6,
          nrow = 2,ncol = 3)

# Compute the likelihood of each cell under the topic model and the
# hard clustering.
samples$loglik1 <- loglik_multinom_topic_model(counts,fit_cluster)
samples$loglik2 <- loglik_multinom_topic_model(counts,poisson2multinom(fit))

p7  <- create_loglik_scatterplot(samples$cluster,"B","royalblue",minloglik)
p8  <- create_loglik_scatterplot(samples$cluster,"C","forestgreen",minloglik)
p9  <- create_loglik_scatterplot(samples$cluster,"Cil","firebrick",minloglik)
p10 <- create_loglik_scatterplot(samples$cluster,"H","turquoise",minloglik)
p11 <- create_loglik_scatterplot(samples$cluster,"T+N","darkorange",minloglik)
p12 <- create_loglik_scatterplot(samples$cluster,"G","gold",minloglik)
plot_grid(p7,p8,p9,
          p10,p11,p12,
          nrow = 2,ncol = 3)

# Compute the likelihood of each cell under the Montoro et al (2018)
# hard clustering.
samples$loglik1 <- loglik_multinom_topic_model(counts,fit_montoro)
p13 <- create_loglik_scatterplot(samples$tissue,"Neuroendocrine","darkorange",
                                 minloglik)
p14 <- create_loglik_scatterplot(samples$tissue,"Tuft","dodgerblue",minloglik)
p15 <- create_loglik_scatterplot(samples$tissue,"Ionocyte","darkmagenta",
                                 minloglik)
plot_grid(p13,p14,p15,
          nrow = 1,ncol = 3)
