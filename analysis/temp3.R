library(Matrix)
library(fastTopics)
fit1 <- readRDS("../output/pbmc-68k/rds/fit-pbmc-68k-em-k=12.rds")$fit
fit2 <- readRDS("../output/pbmc-68k/rds/fit-pbmc-68k-scd-ex-k=12.rds")$fit
fit1 <- poisson2multinom(fit1)
fit2 <- poisson2multinom(fit2)
k <- 7
plot(fit1$L[,k],fit2$L[,k],pch = 20)
abline(a = 0,b = 1,col = "skyblue",lty = "dotted")
