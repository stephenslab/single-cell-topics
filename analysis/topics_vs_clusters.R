suppressMessages(library(ggtern))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))

# SCRIPT PARAMETERS
# -----------------
n <- 200
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
    x     <- x + abs(rnorm(k,0,0.125))
  } else
    x <- rep(1/k,k) + abs(rnorm(k,0,0.2))
  L[i,] <- x/sum(x)
}
    
# Generate the counts.
X <- matrix(0,n,m)
P <- tcrossprod(L,F)
for (i in 1:n)
  X[i,] <- rmultinom(1,s,P[i,])

# Plot the topic proportions.
pdat <- as.data.frame(L)
names(pdat) <- paste0("k",1:k)
p1 <- ggtern(pdat,aes(x = k1,y = k2,z = k3)) +
  geom_point(shape = 21,color = "white",fill = "dodgerblue",size = 1.5) +
  theme_classic(base_size = 10) +
  theme_showarrows() +
  theme(tern.panel.mask.show = FALSE)
print(p1)
