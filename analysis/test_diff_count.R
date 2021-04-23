library(Matrix)
library(ggplot2)
library(cowplot)
# devtools::load_all("~/git/fastTopics")
set.seed(1)
load("../data/pbmc_purified.RData")
fit     <- readRDS(file.path("../output/pbmc-purified/rds",
                             "fit-pbmc-purified-scd-ex-k=6.rds"))$fit
fit     <- poisson2multinom(fit)
fit     <- merge_topics(fit,c(1:2,4:6))
i       <- sample(94655,1000)
j       <- c(1386,2824,6542,11281,21434,18349,20547)
samples <- samples[i,]
genes   <- genes[j,]
counts  <- counts[i,j]
F       <- fit$F[j,]
s       <- fit$s[i]
L       <- s*fit$L[i,]

# Add "pseudocounts" to the data.
dat    <- add_pseudocounts(counts,L,0.01)
counts <- dat$X
L      <- dat$L

# Fit a Poisson model for each gene.
F <- fit_poisson_models(counts,L,method = "scd",numiter = 200)

# Compute the log-fold change statistics.
# mu <- colSums(counts)/sum(s)
n  <- ncol(counts)
se <- matrix(0,n,2)
for (i in 1:n)
  se[i,] <- diag(compute_poisson_covariance(counts[,i],L,F[i,]))
rownames(F) <- genes$symbol
colnames(F) <- c("k1","k2")
rownames(se) <- genes$symbol
colnames(se) <- c("k1","k2")
print(log(F))
print(se)

# Plot the likelihood surface.
j   <- 7
dat <- expand.grid(t1 = seq(-5.85,-5.7,0.002),t2 = seq(-20,-14,0.02))
dat$lik <- 0
n <- nrow(dat)
x <- counts[,j]
loglik_poisson <- function (x, y, e = 1e-15)
  return(sum(x*log(y + e) - y))
for (i in 1:n) {
  f <- exp(c(dat[i,"t1"],dat[i,"t2"]))
  u <- drop(L %*% f)
  dat[i,"lik"] <- loglik_poisson(x,u)
}
dat$lik <- exp(dat$lik - max(dat$lik))
ggplot(dat,aes(x = t1,y = t2,z = lik)) +
  geom_contour(color = "black",bins = 20) +
  geom_point(data = as.data.frame(t(log(F[j,]))),
             mapping = aes(x = k1,y = k2),
             color = "royalblue",shape = 4,
             inherit.aes = FALSE) +
  theme_cowplot(font_size = 10)
dat2 <- subset(dat,t1 == -5.786)
plot(dat2$t2,dat2$lik,type = "l",col = "royalblue",lwd = 2)
points(-17.68,1,col = "red",pch = 4,cex = 0.75)

# Compute Monte Carlo estimates of the 90% HPD intervals.
samples <- simulate_posterior_poisson(counts[,j],L,F[j,],ns = 1000,s = 0.3)
print(hpd(log(samples[,1]),0.95))
print(hpd(log(samples[,2]),0.95))
