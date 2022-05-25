# Generate a data frame summarizing the results of the susie gene set
# enrichment analysis.
compile_gsea_table <- function (s, X, y, gene_set_info, ntop = 10) {

  # Initialize the output.
  out <- NULL
    
  # Get the credible sets.
  cs <- s$sets$cs

  # Add labels to some susie outputs.
  n                 <- length(s$lbf)
  names(s$lbf)      <- paste0("L",1:n)
  rownames(s$alpha) <- paste0("L",1:n)
  rownames(s$mu)    <- paste0("L",1:n)
  
  # Repeat for each CS.
  for (i in names(cs)) {
      
    # Get information about the variables (i.e., the gene sets)
    # included in the CS.
    j <- cs[[i]]
    n <- length(j)
    x <- data.frame(CS = rep(i,n),
                    lbf = rep(s$lbf[i],n),
                    stringsAsFactors = FALSE)
    x <- cbind(x,
               data.frame(pip  = s$alpha[i,j],
                          coef = s$mu[i,j]),
               gene_set_info[j,])

    # For each gene set, extract the "top genes" (genes in the gene
    # set with the Y's that are the largest in magnitude).
    top_genes_all <- vector("character",n)
    for (k in 1:n) {
      genes <- which(X[,j[k]] > 0)
      genes <- genes[order(abs(y[genes]),decreasing = TRUE)]
      genes <- genes[seq(1,ntop)]
      top_genes <- paste0(names(genes)," (",round(y[genes],digits = 3),")")
      top_genes <- paste(top_genes,collapse = " ")
      top_genes_all[k] <- top_genes
    }
    x <- cbind(x,data.frame(top_genes = top_genes_all,
                            stringsAsFactors = FALSE))
    out <- rbind(out,x)
  }

  rownames(out) <- NULL
  return(out)
}
