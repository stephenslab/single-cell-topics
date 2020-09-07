# Select the differential expression analysis results for the
# specified genes only (these are the names of rows in the
# diff_count_analysis output matrices).
select_diff_count_res <- function (diff_count_res, genes) {
  diff_count_res$colmeans <- diff_count_res$colmeans[genes]
  diff_count_res$F0       <- diff_count_res$F0[genes,]
  diff_count_res$F1       <- diff_count_res$F1[genes,]
  diff_count_res$pval     <- diff_count_res$pval[genes,]
  diff_count_res$beta     <- diff_count_res$beta[genes,]
  diff_count_res$se       <- diff_count_res$se[genes,]
  diff_count_res$Z        <- diff_count_res$Z[genes,]
  return(diff_count_res)
}
