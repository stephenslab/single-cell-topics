suppressMessages(library(MCMCpack))
suppressMessages(library(ggtern))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))

# SCRIPT PARAMETERS
# -----------------
n <- 100
m <- 20
s <- 1000
k <- 3

# SIMULATE DATA
# -------------
set.seed(1)

# Generate the m x k matrix of gene frequencies.
F <- matrix(runif(m*k),m,k)

# Generate the topic proportions.
L <- rdirichlet(n,rep(1,k))
    
# Generate the counts.
X <- matrix(0,n,m)
P <- tcrossprod(L,F)
for (i in 1:n)
  X[i,] <- rmultinom(1,s,P[i,])

# Plot the topic proportions.
pdat <- as.data.frame(L)
names(pdat) <- paste0("k",1:k)
p1 <- ggtern(pdat,aes(x = k1,y = k2,z = k3)) +
  geom_point(color = "dodgerblue",size = 1.5) +
  theme_classic(base_size = 10) +
  theme_showarrows()

