library(Matrix)
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
