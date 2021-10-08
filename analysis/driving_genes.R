# These are some preliminary experiments using simulated data to
# verify implementation of the DE analysis methods in fastTopics.
library(Matrix)
library(ggplot2)
library(cowplot)
set.seed(1)

# Simulate data.
set.seed(1)
n  <- 800
m  <- 10000
k  <- 2
p  <- 0.5
s  <- rep(1,n)
se <- 1
L  <- generate_mixture_proportions(n,k)
F  <- matrix(0,m,k)
for (j in 1:m) {
  y     <- rnorm(1,-4,2)
  F[j,] <- 2^c(y,y + (runif(1) < p) * rnorm(1,0,se))
}
X <- generate_poisson_nmf_counts(F,s*L)

# Fit a multinomial topic model using the ground-truth topic
# proportions.
fit0 <- init_poisson_nmf(X,F = F,L = L)
fit <- fit_poisson_nmf(X,fit0 = fit0,numiter = 40,method = "scd",
                       update.loadings = NULL)
fit <- poisson2multinom(fit)

# Compute the KL-divergence based "distinctiveness" measure used in
# CountClust.
D <- min_kl_poisson(fit$F)

# Perform DE analysis with adaptive shrinkage.
de <- de_analysis(fit,X,shrink.method = "ash")

lfc_true <- log2(F[,1]/F[,2])
plot(lfc_true,pmax(-2,pmin(2,de$est[,1])),pch = 20)

pdat <- data.frame(z  = de$z[,1],
                   de = factor(abs(F[,1] - F[,2]) > 1e-15))
pdat <- transform(pdat,
                  z = sign(z) * pmin(4,abs(z)))
ggplot(pdat,aes(x = z,color = de,fill = de)) +
  geom_histogram(bins = 64) +
  scale_color_manual(values = c("darkorange","darkblue")) +
  scale_fill_manual(values = c("darkorange","darkblue")) +
  theme_cowplot()

pdat <- data.frame(d  = D[,1],
                   de = factor(abs(F[,1] - F[,2]) > 1e-15))
pdat <- transform(pdat,
                  d = pmin(0.0001,abs(d)))
ggplot(pdat,aes(x = d,color = de,fill = de)) +
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

