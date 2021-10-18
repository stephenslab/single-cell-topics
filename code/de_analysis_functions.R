# Set all entries of x less than a to a, and set alll entries of x
# greater than b to b.
clamp <- function (x, a, b)
  pmax(a,pmin(b,x))

# Simulate counts from a K=2 Poisson NMF model so that the count data
# look roughly like UMI counts from a single-cell RNA sequencing
# experiment. Input argument m is the number of genes to simulate;
# input s is the vector of scaling factors (one scaling factor for
# each cell); and p is the probability of a expression rates being
# different in the two topics.
simulate_twotopic_umi_data <- function (m = 10000, s = 10^rnorm(200,0,0.2),
                                        p = 0.5) {

  # Get the number of cells to simulate.
  n <- length(s)
    
  # For each sample (row), generate the two topic proportions.
  L <- fastTopics:::generate_mixture_proportions(n,2)

  # Generate the expression rates. The process for generating the
  # expression rates F[i,j] is as follows: with probability 0.5, the
  # rates F[i,1] and F[i,2] are the same, F[i,1] = F[i,2] = 2^y;
  # otherwise, the rates are different, with (multiplicative, base-2)
  # difference e. Here, y is normal with mean -4 and s.d., and e is
  # normal with mean zero and s.d. 1.
  #
  # There might very well be more efficient ways to implement this
  # step, but this current implementation in nonetheless fast enough
  # for my purposes.
  F <- matrix(0,m,2)
  for (i in 1:m) {
    y <- rnorm(1,-4,2)
    e <- rnorm(1,0,1)
    u <- runif(1)
    w <- runif(1)
    if (u > p)
      f <- 2^y
    else if (w < 0.5)
      f <- 2^(y + c(0,e))
    else
      f <- 2^c(y + c(e,0))
    F[i,] <- f
  }

  # Generate the counts from a Poisson NMF model with factors matrix F
  # and loadings matrix s*L.
  X <- fastTopics:::generate_poisson_nmf_counts(F,s*L)

  # Return the counts matrix (X), the scaling factors (s), and the
  # parameters of the Poisson NMF model used to simulate the data (F
  # and L).
  names(s)    <- paste0("c",1:n)
  rownames(X) <- paste0("c",1:n)
  colnames(X) <- paste0("g",1:m)
  rownames(L) <- paste0("c",1:n)
  colnames(L) <- paste0("k",1:2)
  rownames(F) <- paste0("g",1:m)
  colnames(F) <- paste0("k",1:2)
  return(list(X = X,s = s,F = F,L = L))
}

# TO DO: Explain here what this function does, and how to use it.
simulate_manytopic_umi_data <- function (m = 10000, s = 10^rnorm(1000,0,0.2),
                                         k = 6, p = 0.5) {

  # Get the number of cells to simulate.
  n <- length(s)
    
  # For each sample (row), generate the topic proportions.
  L <- fastTopics:::generate_mixture_proportions(n,k)

  # Generate the expression rates. The process for generating the
  # expression rates F[i,j] is as follows: with probability 0.5, the
  # rates F[i,j], for j = 1,...,k, are the same; otherwise, all the
  # rates are the same except for a randomly chosen j, with
  # (multiplicative, base-2) difference e. Here, y is normal with mean
  # -4 and s.d., and e is normal with mean zero and s.d. 1.
  #
  # There might very well be more efficient ways to implement this
  # step, but this current implementation in nonetheless fast enough
  # for my purposes.
  F <- matrix(0,m,k)
  for (i in 1:m) {
    y <- rnorm(1,-4,2)
    e <- rnorm(1,0,1)
    u <- runif(1)
    j <- sample(k,1)
    f <- rep(2^y,k)
    if (u < p)
      f[j] <- 2^(y + e)
    F[i,] <- f
  }

  # Generate the counts from a Poisson NMF model with factors matrix F
  # and loadings matrix s*L.
  X <- fastTopics:::generate_poisson_nmf_counts(F,s*L)

  # Return the counts matrix (X), the scaling factors (s), and the
  # parameters of the Poisson NMF model used to simulate the data (F
  # and L).
  names(s)    <- paste0("c",1:n)
  rownames(X) <- paste0("c",1:n)
  colnames(X) <- paste0("g",1:m)
  rownames(L) <- paste0("c",1:n)
  colnames(L) <- paste0("k",1:k)
  rownames(F) <- paste0("g",1:m)
  colnames(F) <- paste0("k",1:k)
  return(list(X = X,s = s,F = F,L = L))
}

# Create a data frame used to plot a power vs. FDR curve. Input
# argument "true" is a logical vector in which true[i] is TRUE if and
# only if the ith element is a true discovery; input argument "pval"
# may be either a vector of p-values of the same length as "true", or
# an equivalent ranking in which the discoveries with the strongest
# support are ordered first in an ordering obtained from order(pval).
# The output is a data frame with three columns: "power", "fdr" (false
# discovery rate") and "t", the latter being the threshold used to
# calculate power and FDR. Note power = TP/(TP + TN) and fdr = FP/(TP
# + FP), where TP is the number of true positives, FP is the number of
# false positives, and TN is the number of true negatives.
create_fdr_vs_power_curve <- function (pval, true) {
  pval[is.na(pval)] <- max(pval,na.rm = TRUE)
  t   <- sort(unique(pval))
  n   <- length(t)
  out <- data.frame(t = t,power = 0,fdr = 0)
  for (i in 1:n) {
    pos <- pval <= t[i]
    tp  <- sum(pos & true)
    fp  <- sum(pos & !true)
    fn  <- sum(!pos & true)
    out[i,"power"] <- tp/(tp + fn)
    out[i,"fdr"]   <- fp/(tp + fp)
  }
  return(out)
}
