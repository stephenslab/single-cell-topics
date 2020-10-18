library(Matrix)
library(dplyr)
library(fastTopics)
library(ggplot2)
library(cowplot)

# Load data and results.
load("../data/droplet.RData")
samples <- readRDS("../output/droplet/clustering-droplet.rds")
fit <- readRDS("../output/droplet/rds/fit-droplet-scd-ex-k=7.rds")$fit

# Compute the likelihood of each cell under the topic model and the
# two hard clusterings.
fit_merge      <- merge_topics(poisson2multinom(fit),c("k5","k7"))
fit_cluster    <- init_poisson_nmf_from_clustering(counts,samples$cluster)
fit_montoro    <- init_poisson_nmf_from_clustering(counts,samples$tissue)
fit_cluster    <- poisson2multinom(fit_cluster)
fit_montoro    <- poisson2multinom(fit_montoro)
loglik_topics  <- loglik_multinom_topic_model(counts,fit_merge)
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
minloglik <- -20000
dat <- data.frame(cluster = samples$tissue,
                  loglik1 = loglik_montoro,
                  loglik2 = loglik_cluster)
p1 <- create_loglik_scatterplot(dat,"Basal","royalblue",minloglik)
p2 <- create_loglik_scatterplot(dat,"Club","forestgreen",minloglik)
p3 <- create_loglik_scatterplot(dat,"Ciliated","firebrick",minloglik)
p4 <- create_loglik_scatterplot(dat,"Neuroendocrine","darkorange",minloglik)
p5 <- create_loglik_scatterplot(dat,"Tuft","dodgerblue",minloglik)
p6 <- create_loglik_scatterplot(dat,"Goblet","gold",minloglik)
plot_grid(p1 + labs(x = "Montoro et al cluster",y = "our clusters"),
          p2 + labs(x = "Montoro et al cluster",y = "our clusters"),
          p3 + labs(x = "Montoro et al cluster",y = "our clusters"),
          p4 + labs(x = "Montoro et al cluster",y = "our clusters"),
          p5 + labs(x = "Montoro et al cluster",y = "our clusters"),
          p6 + labs(x = "Montoro et al cluster",y = "our clusters"),
          nrow = 2,ncol = 3)

# Compare the hard clustering against the topics.
dat <- data.frame(cluster = samples$cluster,
                  loglik1 = loglik_cluster,
                  loglik2 = loglik_topics)
p7  <- create_loglik_scatterplot(dat,"B","royalblue",minloglik)
p8  <- create_loglik_scatterplot(dat,"C","forestgreen",minloglik)
p9  <- create_loglik_scatterplot(dat,"Cil","firebrick",minloglik)
p10 <- create_loglik_scatterplot(dat,"H","turquoise",minloglik)
p11 <- create_loglik_scatterplot(dat,"T+N","darkorange",minloglik)
p12 <- create_loglik_scatterplot(dat,"G","gold",minloglik)
plot_grid(p7 + labs(x = "our cluster",y = "topics"),
          p8 + labs(x = "our cluster",y = "topics"),
          p9 + labs(x = "our cluster",y = "topics"),
          p10 + labs(x = "our cluster",y = "topics"),
          p11 + labs(x = "our cluster",y = "topics"),
          p12 + labs(x = "our cluster",y = "topics"),
          nrow = 2,ncol = 3)

# Compare the Montoro et al (2018) hard clustering against the topics.
dat <- data.frame(cluster = samples$tissue,
                  loglik1 = loglik_montoro,
                  loglik2 = loglik_topics)
p13 <- create_loglik_scatterplot(dat,"Basal","royalblue",minloglik)
p14 <- create_loglik_scatterplot(dat,"Club","forestgreen",minloglik)
p15 <- create_loglik_scatterplot(dat,"Ciliated","firebrick",minloglik)
p16 <- create_loglik_scatterplot(dat,"Neuroendocrine","darkorange",minloglik)
p17 <- create_loglik_scatterplot(dat,"Tuft","dodgerblue",minloglik)
p18 <- create_loglik_scatterplot(dat,"Ionocyte","darkmagenta",minloglik)
p19 <- create_loglik_scatterplot(dat,"Goblet","gold",minloglik)
plot_grid(p13 + labs(x = "Montoro et al cluster",y = "topics"),
          p14 + labs(x = "Montoro et al cluster",y = "topics"),
          p15 + labs(x = "Montoro et al cluster",y = "topics"),
          p16 + labs(x = "Montoro et al cluster",y = "topics"),
          p17 + labs(x = "Montoro et al cluster",y = "topics"),
          p18 + labs(x = "Montoro et al cluster",y = "topics"),
          p19 + labs(x = "Montoro et al cluster",y = "topics"),
          nrow = 3,ncol = 3)
