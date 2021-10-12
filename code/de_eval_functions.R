# TO DO: Explain what this function does, and how to use it.
simulate_twotopic_umi_data <- function (m = 1e4, s = 10^rnorm(200,-1,0.2),
                                        p = 0.5, a = c(1,1)) {

  # Get the number of cells to simulate.
  n <- length(s)
    
  # For each sample (row), generate the two topic proportions
  # according to the following procedure: (1) Choose uniformly at
  # random the configuration of the nonzero proportions---both topics,
  # the first topic only, or the second topic only; (2) If both topics
  # have nonzero proportions, generate the proportions from the
  # Dirichlet distribution with shape parameter a.
  L           <- matrix(1,n,2)
  k1          <- sample(3,n,replace = TRUE)
  L[k1 == 3,] <- rdirichlet(sum(k1 == 3),a)

  # Generate the relative expression levels.
  #
  # TO DO: Add details describing how this simulation is done.
  # 
  F <- matrix(0,m,2)
  for (i in 1:m) {
    y <- runif(1,-10,3)
    e <- pmax(-5,pmin(5,0.7 * rt(1,df = 3)))
    u <- runif(1)
    w <- runif(1)
    if (u > p)
      F[i,] <- 2^y
    else if (w < 0.5)
      F[i,] <- 2^(y + c(0,e))
    else
      F[i,] <- 2^c(y + c(e,0))
  }

  # Generate the counts.
  X <- matrix(as.double(rpois(n*m,tcrossprod(s*L,F))),n,m)

  # Return the counts matrix (X), and the parameters of the Poisson
  # NMF model used to simulate the data (F, L).
  rownames(X) <- paste0("c",1:n)
  colnames(X) <- paste0("g",1:m)
  rownames(L) <- paste0("c",1:n)
  colnames(L) <- paste0("k",1:2)
  rownames(F) <- paste0("g",1:m)
  colnames(F) <- paste0("k",1:2)
  return(list(X = X,F = F,L = L))
}
