gsea_plot_shapes <- c(23,24,21)
gsea_plot_colors <- c("gold",         # biocyc
                      "yellowgreen",  # C1
                      "dodgerblue",   # C2
                      "olivedrab",    # C3
                      "firebrick",    # C4
                      "darkorange",   # C5
                      "magenta",      # C6
                      "darkblue",     # C7
                      "gold",         # H
                      "tomato",       # humancyc
                      "olivedrab",    # inoh
                      "darkblue",     # kegg
                      "royalblue",    # netpath
                      "darkorange",   # panther
                      "yellowgreen",  # pathbank
                      "magenta",      # pid
                      "dodgerblue")   # reactome

# This function aligns the gene-set data (gene_sets) with the
# gene-wise statistics (diff_count_res) by their Ensembl ids to
# prepare these data for a gene-set enrichment analysis. It is assumed
# that the row names of gene_sets and the row names of diff_count_res
# matrices give the unique gene symbols or ids (e.g., Ensembl ids).
align_gene_data <- function (gene_sets, diff_count_res) {
  x   <- rownames(gene_sets)
  y   <- rownames(diff_count_res$Z)
  ids <- intersect(x,y)
  i   <- match(ids,x)
  j   <- match(ids,y)
  gene_sets               <- gene_sets[i,]
  diff_count_res$colmeans <- diff_count_res$colmeans[j]
  diff_count_res$F0       <- diff_count_res$F0[j,]
  diff_count_res$F1       <- diff_count_res$F1[j,]
  diff_count_res$beta     <- diff_count_res$beta[j,]
  diff_count_res$se       <- diff_count_res$se[j,]
  diff_count_res$Z        <- diff_count_res$Z[j,]
  diff_count_res$pval     <- diff_count_res$pval[j,]
  return(list(gene_sets = gene_sets,diff_count_res = diff_count_res))
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
# set j. Input argument z should be a numeric vector of length n
# containing gene-wise statistics, such as z-scores. The elements of
# z should be named.
#
# Input argument "eps" controls the accuracy of the small p-values.
# Here I set it to be much lower than the default setting since some
# of the gene-set enrichment p-values can be quite small.
perform_gsea <- function (gene_sets, z, eps = 1e-32, nproc = 1, ...) {

  # Convert the gene-sets adjacency matrix into the fgsea gene-sets
  # format.
  pathways <- matrix2pathways(gene_sets)

  # Perform the gene-set enrichment analysis using fgsea.
  out <- suppressWarnings(fgsea(pathways,z,eps = eps,nproc = nproc,...))
  class(out) <- "data.frame"

  # Post-process the fgsea output.
  rownames(out) <- out$pathway
  out <- out[c("pval","log2err","ES","NES")]
  out <- out[colnames(gene_sets),]
  out[is.na(out$log2err),] <- NA
  return(out)
}

# Perform fgsea gene-set enrichment analysis once per column of the
# gene_scores matrix.
perform_gsea_all_topics <- function (gene_sets, diff_count_res,
                                     eps = 1e-32, nproc = 1, ...) {

  # Get the number of gene sets (n) and the number of topics (k).
  Z <- diff_count_res$Z
  n <- ncol(gene_sets)
  k <- ncol(Z)

  # Initialize the outputs.
  out <- list(pval    = matrix(0,n,k),
              log2err = matrix(0,n,k),
              ES      = matrix(0,n,k),
              NES     = matrix(0,n,k))
  rownames(out$pval)    <- colnames(gene_sets)
  rownames(out$log2err) <- colnames(gene_sets)
  rownames(out$ES)      <- colnames(gene_sets)
  rownames(out$NES)     <- colnames(gene_sets)
  colnames(out$pval)    <- colnames(Z)
  colnames(out$log2err) <- colnames(Z)
  colnames(out$ES)      <- colnames(Z)
  colnames(out$NES)     <- colnames(Z)
  
  # Run the gene-set enrichment analysis for each topic.
  cat("k = ")
  for (i in 1:k) {
    if (i == 1)
      cat("1")
    else
      cat(",",i)
    z               <- Z[,i]
    names(z)        <- rownames(Z)
    ans             <- perform_gsea(gene_sets,z,eps,nproc,...)
    out$pval[,i]    <- ans$pval
    out$log2err[,i] <- ans$log2err
    out$ES[,i]      <- ans$ES
    out$NES[,i]     <- ans$NES
  }
  cat("\n")
  
  # Output the results of the gene-set enrichment analysis.
  return(out)
}

# Compile the data frame used for create_gsea_plotly.
compile_data_for_gsea_plot <- function (gene_set_info, gsea_res, k,
                                        max_name_len = 44) {

  # Compute the "signed p-values".
  P <- with(gsea_res,-sign(ES) * log10(pval))

  # Compute the "most extreme" signed p-value among other topics.
  n  <- nrow(P)
  p0 <- rep(0,n)
  for (i in 1:n) {
    if (sign(gsea_res$ES[i,k]) > 0)
      p0[i] <- max(P[i,-k],na.rm = TRUE)
    else
      p0[i] <- min(P[i,-k],na.rm = TRUE)
  }
  
  # Compile the data for plotting.
  dat <- data.frame(p1         = P[,k],
                    p0         = p0,
                    name       = substr(gene_set_info$name,1,max_name_len),
                    id         = gene_set_info$id,
                    database   = gene_set_info$database,
                    collection = with(gene_set_info,
                                      ifelse(database == "MSigDB",
                                             as.character(category_code),
                                             as.character(data_source))),
                    stringsAsFactors = FALSE)
  
  # Subsample the gene sets with p-values close to zero because those
  # gene sets are not particularly interesting, and there are many of
  # them, which slows down the plotting.
  i   <- which(with(dat,!(abs(p0) < 3 & abs(p1) < 3)))
  j   <- which(with(dat,abs(p0) < 3 & abs(p1) < 3))
  n   <- length(j)
  j   <- sample(j,ceiling(n/10))
  dat <- dat[c(i,j),]
  
  # Set any missing signed p-values to zero.
  dat[is.na(dat$p1),"p1"] <- 0
  dat[is.na(dat$p0),"p0"] <- 0
  return(dat)
}

# Create an interactive scatterplot using plotly to explore the
# gene-set enrichment results for a given topic k. Input argument
# "gene_set_info" is a data frame containing information about the
# gene sets; "gsea_res" is an output from function "perform_gsea"; and
# "label_gene_sets" is a vector of gene set ids to be labeled in the
# plot.
create_gsea_plotly <- function (gene_set_info, gsea_res, k, margin = 1,
                                height = 675, width = 800,
                                title = paste("topic",k)) {

  # Compile the data for plotting.
  dat <- compile_data_for_gsea_plot(gene_set_info,gsea_res,k)
  
  # Create the plotly plot.
  pmin <- min(c(dat$p0,dat$p1))
  pmax <- max(c(dat$p0,dat$p1))
  p <- plot_ly(data = dat,x = ~p0,y = ~p1,symbol = ~database,
               color = ~collection,
               text = ~sprintf("%s\nid: %s\np-value: %0.2e",name,id,10^(-p1)),
               type = "scatter",mode = "markers",hoverinfo = "text",
               symbols = gsea_plot_shapes,colors = gsea_plot_colors,
               marker = list(size = 8,line = list(color = "white",width = 1)),
               height = height,width = width)
  p <- add_trace(p,data = data.frame(x = c(pmin,pmax),y = c(pmin,pmax)),
                 x = ~x,y = ~y,mode = "lines",type = "scatter",
                 inherit = FALSE,showlegend = FALSE,
                 line = list(color = "lightgray",dash = "dash",size = 0.3))
  p <- layout(p,
              legend = list(font = list(size = 10)),
              xaxis = list(title="most extreme signed p-value in other topics",
                           showgrid = FALSE,showline = FALSE,zeroline = FALSE),
              yaxis = list(title = "signed p-value in topic",
                           showgrid = FALSE,showline = FALSE,zeroline = FALSE),
              hoverlabel = list(bgcolor = "white",bordercolor = "black",
                                font = list(color = "black",family = "arial",
                                            size = 12)),
              font = list(family = "arial",size = 12),
              title = title)
  return(p)
}

