library(fastTopics)
library(dplyr)
library(ggplot2)
library(cowplot)
k <- c(10,11,12,13,15,20,25,30)
loglik <- c(-399413776.0089223,
            -397218007.3648255,
            -395641654.2734525,
            -393746499.9804388,
            -391327750.4356079,
            -386673304.3196605,
            -383575254.6774732,
            -381414168.2207881)
plot(k,loglik - min(loglik),pch = 20)
lines(k,loglik - min(loglik))

samples <- readRDS("../output/pulseseq/clustering-pulseseq.rds")
fit <- readRDS("../output/pulseseq/rds/fit-pulseseq-scd-ex-k=25.rds")$fit
p1 <- pca_plot(poisson2multinom(fit),pcs = 5:6,fill = samples$tissue)
p2 <- pca_plot(poisson2multinom(fit),pcs = 5:6,k = 1)
p3 <- pca_plot(poisson2multinom(fit),pcs = 5:6,k = 21)

rows <- which(samples$cluster == "I")
fit2 <- select(poisson2multinom(fit),loadings = rows)
fit2 <- merge_topics(fit2,k = c("k3","k4","k5","k7","k8","k9","k10","k11",
                                "k14","k16","k17","k18","k19","k20","k23",
                                "k24","k25"))
p5 <- structure_plot(fit2)
