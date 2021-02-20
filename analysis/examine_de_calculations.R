# TO DO: Explain here what this script is for, and how to use it.
library(Matrix)
# library(fastTopics)
library(ggplot2)
library(cowplot)
library(progress)
set.seed(1)

# Load the count data, the K = 6 topic model fit, and the 7 clusters
# identified in the clustering analysis ("clusters_purified_pbmc").
load("../data/pbmc_purified.RData")
fit <- readRDS(file.path("../output/pbmc-purified/rds",
                         "fit-pbmc-purified-scd-ex-k=6.rds"))$fit
fit <- poisson2multinom(fit)
samples <- readRDS("../output/pbmc-purified/clustering-pbmc-purified.rds")

# Perform differential expression analysis using the multinomial topic
# model, after removing the dendritic cells cluster.
rows <- sample(which(samples$cluster != "dendritic"),9000)
counts <- counts[rows,]
fit <- select_loadings(fit,loadings = rows)
de1 <- diff_count_analysis(fit,counts,pseudocount = 1,fit.method = "em")

# Implement the differential expression analysis using glm with family
# = poisson(link = "identity"). Note that the parameterization is
# slightly different: b0 = f0 and b = f1 - f0.
fit_poisson_glm <- function (x, s, q, pseudocount = 0.01) {
  x <- c(x,pseudocount,pseudocount)
  s <- c(s,1,1)
  q <- c(q,0,1)
  dat <- data.frame(x = x,b0 = s,b = s*q)
  res <- suppressWarnings(glm(x ~ b0 + b - 1,
                              family = poisson(link = "identity"),
                              data = dat,start = c(0.5,0),maxit = 100))
  out <- summary(res)$coefficients
  colnames(out) <- c("coef","se","z","pval")
  b0   <- out["b0","coef"]
  b    <- out["b","coef"]
  u    <- b0 + q*b
  se   <- sqrt(diag(solve(rbind(c(sum(x/u^2),sum(x*q/u^2)),
                                c(sum(x*q/u^2),sum(x*(q/u)^2))))))
  z    <- coef(res)/se
  pval <- 2*pnorm(-abs(z))
  # out[,"z"] <- z
  return(out)
}
fit_univar_poisson_models_glm <- function (X, q, s, pseudocount = 0.01)  {
  m <- ncol(X)
  out <- matrix(0,m,5)
  rownames(out) <- colnames(X)
  colnames(out) <- c("b0","b","se","z","pval")
  pb <- progress_bar$new(total = m)
  for (i in 1:m) {
    pb$tick()
    res <- fit_poisson_glm(X[,i],s,q,pseudocount)
    out[i,] <- c(res["b0","coef"],res["b",c("coef","se","z","pval")])
  }
  return(out)
}
i <- 3
j <- sample(which(de1$beta[,i] > 0),100)
s <- rowSums(counts)
de2 <- fit_univar_poisson_models_glm(counts[,j],fit$L[,i],s)
fit2 <- fit
fit2$F <- fit$F[j,]
de3 <- diff_count_analysis(fit2,counts[,j],s = s,shrink.method = "none",
                           pseudocount = 1,fit.method = "glm")
pdat <- data.frame(x = de3$Z[j,i],y = de2[,"z"],lfc = de1$beta[j,i])
p1 <- ggplot(pdat,aes(x = x,y = y,fill = lfc)) +
  geom_point(color = "white",shape = 21,size = 1.5) +
  scale_fill_gradient2(low = "darkblue",mid = "lightskyblue",
                       high = "orangered") +
  geom_abline(intercept = 0,slope = 1,color = "black",linetype = "dotted") +
  theme_cowplot(font_size = 12)

p2 <- volcano_plot(de1,"k3",genes$symbol,
                   label_above_quantile = 0.999,
                   subsample_below_quantile = 0.5)

