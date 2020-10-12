library(Matrix)
library(dplyr)
library(fastTopics)
library(ggplot2)
library(cowplot)

# Load data and results.
load("../data/pulseseq.RData")
samples <- readRDS("../output/pulseseq/clustering-pulseseq.rds")
fit <- readRDS("../output/pulseseq/rds/fit-pulseseq-scd-ex-k=11.rds")$fit

# Compute the likelihood of each cell under the topic model and the
# two hard clusterings.
fit <- merge_topics(poisson2multinom(fit),c("k4","k5","k6","k8","k10"))
fit <- merge_topics(fit,c("k1","k3","k9"))
fit_cluster    <- init_poisson_nmf_from_clustering(counts,samples$cluster)
fit_montoro    <- init_poisson_nmf_from_clustering(counts,samples$tissue)
fit_cluster    <- poisson2multinom(fit_cluster)
fit_montoro    <- poisson2multinom(fit_montoro)
loglik_topics  <- loglik_multinom_topic_model(counts,fit)
loglik_cluster <- loglik_multinom_topic_model(counts,fit_cluster)
loglik_montoro <- loglik_multinom_topic_model(counts,fit_montoro)

create_loglik_scatterplot <- function (dat, k, color = "black",
                                       minloglik = -Inf) {
  dat    <- subset(dat,loglik1 > minloglik & loglik2 > minloglik)
  limits <- range(c(dat$loglik1,dat$loglik2))
  dat    <- subset(dat,cluster == k)
  return(ggplot(dat,aes(x = loglik1,y = loglik2)) +
         geom_point(shape = 21,color = "white",fill = color) +
         geom_abline(intercept = 0,slope = 1,linetype = "dotted",
                     color = "black") +
         xlim(limits) +
         ylim(limits) +
         labs(title = sprintf("%s (n = %d)",k,nrow(dat))) +
         theme_cowplot(font_size = 9) +
         theme(plot.title = element_text(size = 9,face = "plain")))
}

# Compare the two hard clusterings.
minloglik <- -30000
dat <- data.frame(cluster = samples$tissue,
                  loglik1 = loglik_montoro,
                  loglik2 = loglik_cluster)
p1 <- create_loglik_scatterplot(dat,"basal","royalblue",minloglik)
p2 <- create_loglik_scatterplot(dat,"club","forestgreen",minloglik)
p3 <- create_loglik_scatterplot(dat,"ciliated","firebrick",minloglik)
p4 <- create_loglik_scatterplot(dat,"neuroendocrine","darkorange",minloglik)
p5 <- create_loglik_scatterplot(dat,"tuft","dodgerblue",minloglik)
p6 <- create_loglik_scatterplot(dat,"ionocyte","darkmagenta",minloglik)
p7 <- create_loglik_scatterplot(dat,"proliferating","peru",minloglik)
p8 <- create_loglik_scatterplot(dat,"goblet","gold",-Inf)
plot_grid(p1 + labs(x = "Montoro et al cluster",y = "our clusters"),
          p2 + labs(x = "Montoro et al cluster",y = "our clusters"),
          p3 + labs(x = "Montoro et al cluster",y = "our clusters"),
          p4 + labs(x = "Montoro et al cluster",y = "our clusters"),
          p5 + labs(x = "Montoro et al cluster",y = "our clusters"),
          p6 + labs(x = "Montoro et al cluster",y = "our clusters"),
          p7 + labs(x = "Montoro et al cluster",y = "our clusters"),
          p8 + labs(x = "Montoro et al cluster",y = "our clusters"),
          nrow = 3,ncol = 3)

# Compare the hard clustering against the topics.
minloglik <- -25000
dat <- data.frame(cluster = samples$cluster,
                  loglik1 = loglik_cluster,
                  loglik2 = loglik_topics)
p9  <- create_loglik_scatterplot(dat,"B","royalblue",minloglik)
p10 <- create_loglik_scatterplot(dat,"C","forestgreen",minloglik)
p11 <- create_loglik_scatterplot(dat,"Cil","firebrick",minloglik)
p12 <- create_loglik_scatterplot(dat,"T+N","darkorange",minloglik)
p13 <- create_loglik_scatterplot(dat,"I","darkmagenta",minloglik)
p14 <- create_loglik_scatterplot(dat,"P","peru",minloglik)
plot_grid(p9,p10,p11,p12,p13,p14,nrow = 2,ncol = 3)

# Compare the Montoro et al (2018) hard clustering against the topics.
minloglik <- -30000
dat <- data.frame(cluster = samples$tissue,
                  loglik1 = loglik_montoro,
                  loglik2 = loglik_topics)
p15 <- create_loglik_scatterplot(dat,"basal","royalblue",minloglik)
p16 <- create_loglik_scatterplot(dat,"club","forestgreen",minloglik)
p17 <- create_loglik_scatterplot(dat,"ciliated","firebrick",minloglik)
p18 <- create_loglik_scatterplot(dat,"neuroendocrine","darkorange",minloglik)
p19 <- create_loglik_scatterplot(dat,"tuft","dodgerblue",minloglik)
p20 <- create_loglik_scatterplot(dat,"ionocyte","darkmagenta",minloglik)
p21 <- create_loglik_scatterplot(dat,"proliferating","peru",minloglik)
p22 <- create_loglik_scatterplot(dat,"goblet","gold",-Inf)
plot_grid(p15 + labs(x = "Montoro et al cluster",y = "topics"),
          p16 + labs(x = "Montoro et al cluster",y = "topics"),
          p17 + labs(x = "Montoro et al cluster",y = "topics"),
          p18 + labs(x = "Montoro et al cluster",y = "topics"),
          p19 + labs(x = "Montoro et al cluster",y = "topics"),
          p20 + labs(x = "Montoro et al cluster",y = "topics"),
          p21 + labs(x = "Montoro et al cluster",y = "topics"),
          p22 + labs(x = "Montoro et al cluster",y = "topics"),
          nrow = 3,ncol = 3)
