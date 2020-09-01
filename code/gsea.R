# This function aligns the gene-set data (gene_info, gene_sets) with
# the gene-wise statistics (genes, diff_count_res) by their Ensembl
# ids to prepare these data for a gene-set enrichment analysis. It is
# assumed that gene_info and genes each have an "ensembl" column
# containing the Ensembl gene ids.
align_gene_data_by_ensembl <- function (gene_info, gene_sets, 
                                        genes, diff_count_res) {
  ids       <- intersect(gene_info$ensembl,genes$ensembl)
  i         <- match(ids,gene_info$ensembl)
  j         <- match(ids,genes$ensembl)
  gene_info <- gene_info[i,]
  gene_sets <- gene_sets[i,]
  genes     <- genes[j,]
  diff_count_res$colmeans <- diff_count_res$colmeans[j]
  diff_count_res$F0       <- diff_count_res$F0[j,]
  diff_count_res$F1       <- diff_count_res$F1[j,]
  diff_count_res$beta     <- diff_count_res$beta[j,]
  diff_count_res$se       <- diff_count_res$se[j,]
  diff_count_res$Z        <- diff_count_res$Z[j,]
  diff_count_res$pval     <- diff_count_res$pval[j,]
  return(list(gene_info = gene_info,gene_sets = gene_sets,
              genes = genes,diff_count_res = diff_count_res))
}

# Recover a list of gene sets from n x m adjacency matrix A, in which
# n is the number of genes and m is the number of gene sets, and
# A[i,j] is greater than zero if and only if gene i is included in
# gene set j. This function can be used to prepare the "pathways"
# input to fgsea from a collection of gene sets encoded as a matrix.
# Note that the rows and columns of the matrix, A, should be named.
matrix2pathways <- function (A) {
  n          <- ncol(A)
  out        <- vector("list",n)
  names(out) <- colnames(A)
  genes      <- rownames(A)
  for (i in 1:n)
    out[[i]] <- genes[A[,i] > 0]
  return(out)
}

# Perform gene-set enrichment analysis using the fast method described
# in Korotkevich et al (2016). This method uses permutation testing
# to account for correlations between genes. The p-values are based on
# Kolmogorov-Smirnov test statistics which are "normalized" to account
# for differences in gene-set sizes; the normalized test statistics
# are provided in the "NES" (normalized enrichment score) column of
# the output. See also Subramanian et al (2005) for background on the
# method.
#
# Input argument "gene_sets" should be an n x m matrix, in which n is
# the number of genes, and m is the number of gene sets, and entry
# (i,j) is greater than zero if and only if gene i is included in gene
# set j. Input argument "gene_scores" should be a numeric vector of
# length n containing gene-wise statistics, such as z-scores. The
# elements of gene_scores should be named.
#
# Input argument "eps" controls the accuracy of the small p-values.
# Here I set it to be much lower than the default setting since some
# of the gene-set enrichment p-values can be quite small.
perform_gsea <- function (gene_sets, gene_scores, eps = 1e-32,
                          nproc = 1, ...) {

  # Convert the gene-sets adjacency matrix into the fgsea gene-sets
  # format.
  pathways <- matrix2pathways(gene_sets)

  # Perform the gene-set enrichment analysis using fgsea.
  out <- suppressWarnings(fgsea(pathways,gene_scores,eps = eps,
                                nproc = nproc,...))
  class(out) <- "data.frame"

  # Post-process the fgsea output.
  # rownames(out) <- out$pathway
  # out <- out[c("pval","log2err","ES","NES")]
  # out <- out[colnames(gene_sets),]
  # out[is.na(out$log2err),] <- NA
  return(out)
}

# Perform fgsea gene-set enrichment analysis once per column of the
# gene_scores matrix.
perform_gsea_all_topics <- function (gene_sets, gene_scores, eps = 1e-32,
                                     nproc = 1, ...) {

  # Get the number of gene sets (n) and the number of topics (k).
  n <- ncol(gene_sets)
  k <- ncol(gene_scores)

  # Initialize the outputs.
  out <- list(pval    = matrix(0,n,k),
              log2err = matrix(0,n,k),
              ES      = matrix(0,n,k),
              NES     = matrix(0,n,k))
  rownames(out$pval)    <- colnames(gene_sets)
  rownames(out$log2err) <- colnames(gene_sets)
  rownames(out$ES)      <- colnames(gene_sets)
  rownames(out$NES)     <- colnames(gene_sets)
  colnames(out$pval)    <- colnames(gene_scores)
  colnames(out$log2err) <- colnames(gene_scores)
  colnames(out$ES)      <- colnames(gene_scores)
  colnames(out$NES)     <- colnames(gene_scores)
  
  # Run the gene-set enrichment analysis for each topic.
  for (i in 1:k) {
    z               <- gene_scores[,i]
    names(z)        <- rownames(gene_scores)
    ans             <- perform_gsea(gene_sets,z,eps,nproc,...)
    out$pval[,i]    <- ans$pval
    out$log2err[,i] <- ans$log2err
    out$ES[,i]      <- ans$ES
    out$NES[,i]     <- ans$NES
  }

  # Output the gene-set enrichment results.
  return(out)
}
