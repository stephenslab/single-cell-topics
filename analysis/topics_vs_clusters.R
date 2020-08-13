library(Ternary)
library(ggplot2)
library(cowplot)
library(Rtsne)

# SCRIPT PARAMETERS
# -----------------
n <- 400
m <- 20
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
  if (runif(1) < 0.5) {
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

# Plot the topic proportions.
pdat <- as.data.frame(L)
par(mar = c(0,0,0,0))
TernaryPlot(alab = "k1",blab = "k2",clab = "k3",
            grid.col = "skyblue",grid.minor.lines = 0)
TernaryPoints(pdat,pch = 21,cex = 0.75)

# Create a t-SNE plot
tsne <- Rtsne(X,2,pca = FALSE,normalize = FALSE,theta = 0.1,
              perplexity = 100,max_iter = 1000)
colnames(tsne$Y) <- c("d1","d2")
pdat <- cbind(as.data.frame(tsne$Y),L)
p1 <- ggplot(pdat,aes(x = d1,y = d2,fill = k1)) +
  geom_point(shape = 21,color = "white",size = 1.5) +
  theme_cowplot(font_size = 10)
print(p1)

# TO DO:
#
#  + Compute t-SNE from loadings.
#
#  + Run Laplacian eigenmaps + k-means on loadings.
# 
