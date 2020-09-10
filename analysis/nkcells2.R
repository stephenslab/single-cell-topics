library(Matrix)
library(fastTopics)
load("../data/pbmc_purified.RData")
fit <-
  readRDS("../output/pbmc-purified/rds/fit-pbmc-purified-scd-ex-k=6.rds")$fit
fit    <- poisson2multinom(fit)
k      <- 6
j      <- c(2824,20864)
fit$F  <- fit$F[j,]
s      <- rowSums(counts)
counts <- as.matrix(counts[,j])
out    <- diff_count_analysis(fit,counts,s = s)

