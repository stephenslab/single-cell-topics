library(fastTopics)
library(ggplot2)
library(cowplot)
k <- c(10,11,12,13,15,20,25,30)
loglik <- c(-40022356.71064168,
            -39759467.82240736,
            -39536064.99617583,
            -39341529.73828159,
            -38910982.63832612,
            -38423380.53641042,
            -38110854.54406747,
            -37847817.65065855)
plot(k,loglik - min(loglik),pch = 20)
lines(k,loglik - min(loglik))

samples <- readRDS("../output/droplet/clustering-droplet.rds")
fit <- readRDS("../output/droplet/rds/fit-droplet-scd-ex-k=15.rds")$fit
p1 <- pca_plot(poisson2multinom(fit),pcs = 1:2,fill = samples$tissue)
p2 <- pca_plot(poisson2multinom(fit),pcs = 3:4,fill = samples$tissue)

