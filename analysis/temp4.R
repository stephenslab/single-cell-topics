# Look for genes in the dendritic of CD8+ cells that are poorly
# captured by the topic model.
rows <- which(samples$cluster == "CD8+")
fit1 <- select_loadings(fit,loadings = rows)
fit1 <- multinom2poisson(fit1)
fit1 <- list(F = fit1$L,L = fit1$F)
class(fit1) <- c("poisson_nmf_fit","list")
f0 <- colMeans(counts[rows,])
f1 <- rowMeans(with(fit1,tcrossprod(L,F)))
ll <- loglik_poisson_nmf(t(counts[rows,]),fit1)
j <- which(ll < (-4000))
plot(log(f0),log(f1),pch = 20)
