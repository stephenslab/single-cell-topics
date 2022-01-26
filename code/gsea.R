# Recover a list of gene sets from n x m adjacency matrix A, in which
# n is the number of genes and m is the number of gene sets, and
# A[i,j] is greater than zero if and only if gene i is included in
# gene set j. This function can be used to prepare the "pathways"
# input to fgsea from a collection of gene sets encoded as a matrix.
# Note that the rows and columns of the matrix, A, should be named.
matrix2list <- function (A) {
  n          <- ncol(A)
  out        <- vector("list",n)
  names(out) <- colnames(A)
  genes      <- rownames(A)
  for (i in 1:n)
    out[[i]] <- which(A[,i] > 0)
  return(out)
}
