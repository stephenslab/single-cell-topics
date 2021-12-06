# library(Matrix)
# load("../data/pbmc_purified.RData")
l11 <- rep(0,13)
ll2 <- rep(0,13)
ll3 <- rep(0,13)
for (i in 2:13) {
  fit    <- fits[[paste0("fit-pbmc-purified-scd-ex-k=",i)]]
  ll1[i] <- max(fit$progress$loglik)
  ll2[i] <- sum(loglik_poisson_nmf(counts,fit))
  ll3[i] <- sum(loglik_multinom_topic_model(counts,fit))
}
