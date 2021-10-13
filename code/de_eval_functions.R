# Set all entries of x less than a to a, and set alll entries of x
# greater than b to b.
clamp <- function (x, a, b)
  pmax(a,pmin(b,x))

# Simulate counts from a Poisson NMF model intended to "mimic" UMI
# count data from a single-cell RNA sequencing experiment. The
# simulation parameters were chosen based on examining expression
# levels (mean UMI counts) in B cells vs. other cells in the
# "purified" PBMC data.
simulate_twotopic_umi_data <- function (m = 1e4, s = 10^rnorm(200,0.2,0.2),
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

  # Generate the expression rates. The process for generating the
  # expression rates F[i,j] is as follows: with probability 0.5, the
  # rates F[i,1] and F[i,2] are the same, F[i,1] = F[i,2] = 2^y;
  # otherwise, the rates are different, with (multiplicative, base-2)
  # difference given by 0.7*e (with lower and upper limits of -10 and
  # +10 on the differences). Here, y is uniform on [-10,+3] and e is
  # t-distributed with 3 degrees of freedom.
  F <- matrix(0,m,2)
  for (i in 1:m) {
    y <- rnorm(1,-4,2)
    e <- clamp(0.7*rt(1,df = 3),-10,10)
    u <- runif(1)
    w <- runif(1)
    if (u > p)
      F[i,] <- 2^y
    else if (w < 0.5)
      F[i,] <- 2^(y + c(0,e))
    else
      F[i,] <- 2^c(y + c(e,0))
  }

  # Generate the counts from a Poisson NMF model with factors matrix F
  # and loadings matrix s*L.
  X <- matrix(as.double(rpois(n*m,tcrossprod(s*L,F))),n,m)

  # Return the counts matrix (X), the multinomial sample sizes (s),
  # and the parameters of the Poisson NMF model used to simulate the
  # data (F and L).
  names(s)    <- paste0("c",1:n)
  rownames(X) <- paste0("c",1:n)
  colnames(X) <- paste0("g",1:m)
  rownames(L) <- paste0("c",1:n)
  colnames(L) <- paste0("k",1:2)
  rownames(F) <- paste0("g",1:m)
  colnames(F) <- paste0("k",1:2)
  return(list(X = X,s = s,F = F,L = L))
}
