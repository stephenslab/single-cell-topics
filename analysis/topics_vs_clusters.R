library(Matrix)
library(Ternary)
library(ggplot2)
library(cowplot)
library(Rtsne)
library(fastTopics)

# SCRIPT PARAMETERS
# -----------------
n <- 400
m <- 40
s <- 1000
k <- 3

# SIMULATE DATA
# -------------
set.seed(1)

# Generate the m x k matrix of gene frequencies.
F <- matrix(runif(m*k),m,k)

# Generate the topic proportions.
L <- matrix(0,n,k)
for (i in 1:n) {
  x <- rep(0,k)
  if (runif(1) < 0.9) {
    j     <- sample(k,1)
    x[j]  <- 1
    x     <- x + abs(rnorm(k,0,0.15))
  } else
    x <- rep(1/k,k) + abs(rnorm(k,0,0.15))
  L[i,] <- x/sum(x)
}
colnames(L) <- paste0("k",1:k)

# Generate the counts.
X <- matrix(0,n,m)
P <- tcrossprod(L,F)
for (i in 1:n)
  X[i,] <- rmultinom(1,s,P[i,])

fit <- init_poisson_nmf(X,L = L,init.method = "random")
fit <- fit_poisson_nmf(X,fit0 = fit,numiter = 20,update.loadings = NULL)
fit <- fit_poisson_nmf(X,fit0 = fit,numiter = 100,
                       control = list(extrapolate = TRUE))
Lsim <- L
L    <- fit$L

# Plot the topic proportions.
pdat <- as.data.frame(L)
par(mar = c(0,0,0,0))
TernaryPlot(alab = "k1",blab = "k2",clab = "k3",
            grid.col = "skyblue",grid.minor.lines = 0)
TernaryPoints(pdat,pch = 21,cex = 0.75)

# Create a t-SNE plot.
tsne <- Rtsne(X,2,pca = FALSE,normalize = FALSE,theta = 0.1,
              perplexity = 100,max_iter = 1000)
colnames(tsne$Y) <- c("d1","d2")
pdat <- cbind(as.data.frame(tsne$Y),
              data.frame(k = factor(apply(L,1,which.max))))
p1 <- ggplot(pdat,aes(x = d1,y = d2,fill = k)) +
  geom_point(shape = 21,color = "white",size = 1.5) +
  scale_fill_manual(values = c("dodgerblue","darkorange","darkblue")) +
  labs(x = "tsne 1",y = "tsne 2") +
  theme_cowplot(font_size = 10) 

# t-SNE from loadings.
tsne <- Rtsne(L,2,pca = FALSE,normalize = FALSE,theta = 0.1,
              perplexity = 100,max_iter = 1000)
colnames(tsne$Y) <- c("d1","d2")
pdat[,1:2] <- tsne$Y
p2 <- ggplot(pdat,aes(x = d1,y = d2,fill = k)) +
  geom_point(shape = 21,color = "white",size = 1.5) +
  scale_fill_manual(values = c("dodgerblue","darkorange","darkblue")) +
  labs(x = "tsne 1",y = "tsne 2") +
  theme_cowplot(font_size = 10) 

# Implements the "Laplacian eigenmaps" method, in which samples are
# "projected" nonlinearly onto an k-dimension embedding by computing the
# largest k + 1 eigenvalues of the (normalized) graph Laplacian. Input
# W must be an n x n matrix of class 'dgCMatrix' (see function
# 'sparseMatrix'), where n is the number of samples.
calc.eigenmap <- function (W, k, tol = 1e-6, maxiter = 1e4) {

  # Compute the largest k + 1 eigenvalues (and corresponding
  # eigenvectors) of the normalized Laplacian matrix.
  out <- eigs(calc.laplacian(W),k + 1,which = "LM",
              opts = list(tol = tol,maxitr = maxiter,ncv = max(20,4*k)))

  # Get the k largest eigenvalues (except for the largest eigenvalue,
  # which is typically 1, or very close to 1).
  v        <- out$values[-1]
  names(v) <- paste0("d",1:k)

  # Add row and column names to the eigenvectors corresponding to the
  # k largest eigenvalues.
  R           <- out$vectors[,-1]
  rownames(R) <- rownames(W)
  colnames(R) <- paste0("d",1:k)

  # Output (1) the eigenvalues, (2) a data frame containing the
  # corresponding eigenvectors, and (3) the number of iterations it
  # took to reach convergence.
  return(list(niter = out$niter,v = v,R = as.data.frame(R)))
}

# If W21 = W11, the return value is the normalized Laplacian of the n
# x n sparse, symmetric affinity matrix W11. In general, if W21 is an
# m x n matrix, the return value is the m x n matrix of "normalized"
# affinities L(i,j) between all rows i of W21 and all rows/columns j
# of W11.
calc.laplacian <- function (W11, W21 = W11) {
  D1 <- .symDiagonal(nrow(W11),1/sqrt(rowSums(W11)))
  D2 <- .symDiagonal(nrow(W21),1/sqrt(rowSums(W21)))
  return(D2 %*% W21 %*% D1)
}

# Laplacian eigenmaps from loadings.
library(rARPACK)
out <- calc.eigenmap(tcrossprod(L),2)
plot(out$R[,1],out$R[,2],pch = 20,col = pdat$k)
ans <- kmeans(out$R,4,iter.max = 100)
points(ans$centers[,1],ans$centers[,2],pch = 4,col = "dodgerblue")
