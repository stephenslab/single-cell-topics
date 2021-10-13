# These are some preliminary experiments using simulated data to
# verify implementation of the DE analysis methods in fastTopics.
library(Matrix)
library(MCMCpack)
library(fastTopics)
library(ggplot2)
library(cowplot)
source("../code/de_eval_functions.R")
set.seed(1)

# Simulate data.
set.seed(1)
n  <- 200
m  <- 10000
k  <- 2
p  <- 0.5
s  <- rep(1,n)
se <- 1
L  <- fastTopics:::generate_mixture_proportions(n,k)
F  <- matrix(0,m,k)
for (j in 1:m) {
  y <- rnorm(1,-4,2)
  e <- rnorm(1,0,se)
  u <- runif(1)
  w <- runif(1)
  if (u > p)
    F[j,] <- 2^y
  else if (w < 0.5)
    F[j,] <- 2^(y + c(0,e))
  else
    F[j,] <- 2^c(y + c(e,0))
}
X <- fastTopics:::generate_poisson_nmf_counts(F,s*L)
rownames(X) <- paste0("c",1:n)
colnames(X) <- paste0("g",1:m)
rownames(L) <- paste0("c",1:n)
rownames(F) <- paste0("g",1:m)
colnames(L) <- paste0("k",1:k)
colnames(F) <- paste0("k",1:k)

# Fit a multinomial topic model using the ground-truth topic
# proportions.
fit0 <- init_poisson_nmf(X,F = F,L = s*L)
fit <- fit_poisson_nmf(X,fit0 = fit0,numiter = 40,method = "scd",
                       update.loadings = NULL)
fit <- poisson2multinom(fit)

# Compute the KL-divergence based "distinctiveness" measure used in
# CountClust.
# D <- min_kl_poisson(fit$F)

# Perform DE analysis without adaptive shrinkage.
de.noshrink <- de_analysis(fit,X,shrink.method = "none")

# Perform DE analysis with adaptive shrinkage.
de <- de_analysis(fit,X,shrink.method = "ash")

# Compare distributions of z-scores.
lfc  <- log2(F[,1]/F[,2])
pdat <- data.frame(noshrink = clamp(de.noshrink$z[,1],-4,+4),
                   shrink   = clamp(de$z[,1],-4,+4),
                   de       = factor(abs(F[,1] - F[,2]) > 1e-8))
p1 <- ggplot(pdat,aes(x = noshrink,color = de,fill = de)) +
  geom_histogram(bins = 64) +
  scale_color_manual(values = c("darkorange","darkblue")) +
  scale_fill_manual(values = c("darkorange","darkblue")) +
  theme_cowplot()
p2 <- ggplot(pdat,aes(x = shrink,color = de,fill = de)) +
  geom_histogram(bins = 64) +
  scale_color_manual(values = c("darkorange","darkblue")) +
  scale_fill_manual(values = c("darkorange","darkblue")) +
  theme_cowplot()
print(plot_grid(p1,p2,nrow = 2,ncol = 1))

stop()

lfc_true <- log2(F[,1]/F[,2])
i <- which(de$lfsr[,1] < 0.1)
pdat <- data.frame(true = lfc_true,est = de$est[,1],lf0 = log10(de$f0))[i,]
ggplot(pdat,aes(x = true,y = est,fill = lf0)) +
  geom_point(shape = 21,color = "white") +
  geom_abline(intercept = 0,slope = 1,color = "magenta",linetype = "dotted") +
  scale_fill_gradient2(low = "deepskyblue",mid = "gold",high = "orangered",
                       na.value = "gainsboro",midpoint = -3) +
  theme_cowplot()

pdat <- data.frame(true = lfc_true,
                   postmean = de$postmean[,1],
                   lf0 = log10(de$f0))[i,]
ggplot(pdat,aes(x = true,y = postmean,fill = lf0)) +
  geom_point(shape = 21,color = "white") +
  geom_abline(intercept = 0,slope = 1,color = "magenta",linetype = "dotted") +
  scale_fill_gradient2(low = "deepskyblue",mid = "gold",high = "orangered",
                       na.value = "gainsboro",midpoint = -3) +
  theme_cowplot()

plot(lfc_true[i],de$postmean[i,1],pch = 20)
abline(a = 0,b = 1,col = "cyan",lty = "dotted")

pdat <- data.frame(z  = de$z[,1],
                   de = factor(abs(F[,1] - F[,2]) > 1e-15))
pdat <- transform(pdat,
                  z = sign(z) * pmin(4,abs(z)))
ggplot(pdat,aes(x = z,color = de,fill = de)) +
  geom_histogram(bins = 64) +
  scale_color_manual(values = c("darkorange","darkblue")) +
  scale_fill_manual(values = c("darkorange","darkblue")) +
  theme_cowplot()

pdat <- data.frame(z  = de.noshrink$z[,1],
                   de = factor(abs(F[,1] - F[,2]) > 1e-15))
pdat <- transform(pdat,
                  z = sign(z) * pmin(4,abs(z)))
ggplot(pdat,aes(x = z,color = de,fill = de)) +
  geom_histogram(bins = 64) +
  scale_color_manual(values = c("darkorange","darkblue")) +
  scale_fill_manual(values = c("darkorange","darkblue")) +
  theme_cowplot()

i    <- which(de$est[,1] >= 0)
pdat <- data.frame(d  = D[i,1],
                   de = factor(abs(F[i,1] - F[i,2]) > 1e-15))
pdat <- transform(pdat,
                  d = pmin(0.001,abs(d)))
ggplot(pdat,aes(x = d,color = de,fill = de)) +
  geom_histogram(bins = 64) +
  scale_color_manual(values = c("darkorange","darkblue")) +
  scale_fill_manual(values = c("darkorange","darkblue")) +
  theme_cowplot()

i    <- which(de$est[,1] >= 0)
pdat <- data.frame(lfsr = de$lfsr[i,1],
                   de   = factor(abs(F[i,1] - F[i,2]) > 1e-15))
ggplot(pdat,aes(x = lfsr,color = de,fill = de)) +
  geom_histogram(bins = 64) +
  scale_color_manual(values = c("darkorange","darkblue")) +
  scale_fill_manual(values = c("darkorange","darkblue")) +
  theme_cowplot()

pdat <- data.frame(pval = 10^(-de$lpval[,1]),
                   de   = factor(abs(F[,1] - F[,2]) > 1e-15))
ggplot(pdat,aes(x = pval,color = de,fill = de)) +
  geom_histogram(bins = 64) +
  scale_color_manual(values = c("darkorange","darkblue")) +
  scale_fill_manual(values = c("darkorange","darkblue")) +
  theme_cowplot()

stop()

# Perform DE analysis using fitted multinomial topic model, first
# without adaptive shrinkage.
de1 <- de_analysis(fit,X,shrink.method = "none")


# First check that the adaptive shrinkage is doing what it is supposed
# to do.
pdat <- data.frame(x   = as.vector(de1$est),
                   y   = as.vector(de2$est),
                   lf0 = log10(de1$f0))
ggplot(pdat,aes(x = x,y = y,fill = lf0)) +
  geom_point(shape = 21,size = 2.5,color = "white") +
  geom_abline(intercept = 0,slope = 1,color = "magenta",linetype = "dotted") +
  scale_fill_gradient2(low = "deepskyblue",mid = "gold",high = "orangered",
                       na.value = "gainsboro",midpoint = -3) +
  labs(x = "LFC estimate before shrinkage",
       y = "LFC estimate after shrinkage",
       fill = "log10 mean") +
  theme_cowplot()

# Compare the topic 1 LFC estimates against the true values; that is,
# the values used to simulate the data.
fit_true <- poisson2multinom(fit0)
pdat <- data.frame(x   = log2(fit_true$F[,1]/fit_true$F[,2]),
                   y   = de2$z[,1],
                   lf0 = log10(de1$f0))
ggplot(pdat,aes(x = x,y = y,fill = lf0)) +
  geom_point(shape = 21,size = 2.5,color = "white") +
    scale_fill_gradient2(low = "deepskyblue",mid = "gold",high = "orangered",
                         na.value = "gainsboro",midpoint = -3) +
  theme_cowplot()

