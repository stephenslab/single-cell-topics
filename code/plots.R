# This function transforms the data to a standard normal distribution
# via "quantile normalization". It randomly assigns rank to ties.
qt_random_tie <- function (x) {
  y = rank(x,ties.method = "random")
  return(qqnorm(y,plot.it = FALSE)$x)
}

# Given a matrix F of probability vectors (each column should sum to
# 1), compute the total variation distance between each pair of
# columns. The output is an n x n matrix, where n is the number of
# columns in F.
totalvardist <- function (F) {
  n <- ncol(F)
  d <- matrix(0,n,n)
  rownames(d) <- colnames(F)
  colnames(d) <- colnames(F)
  for (i in 1:n)
    for (j in 1:n)
      d[i,j] <- sum(abs(F[,i] - F[,j]))/2
  return(d)
}

# This is a refinement of the volcano plot implemented in the
# fastTopics package. Here, an additional set of genes is highlighted
# by labeling the points with darker text.
volcano_plot_with_highlighted_genes <- function (diff_count_res, k, 
                                                 genes, ...) {
  dat <- data.frame(beta  = diff_count_res$beta[genes,k],
                    y     = abs(diff_count_res$Z[genes,k]),
                    label = genes) 
  rows <- match(genes,rownames(diff_count_res$beta))
  rownames(diff_count_res$beta)[rows] <- ""
  rownames(diff_count_res$Z)[rows] <- ""
  return(volcano_plot(diff_count_res,k = k,
                      ggplot_call = ggplot_call_for_volcano_plot,...) +
         geom_text_repel(data = dat,
                         mapping = aes(x = beta,y = y,label = label),
                         inherit.aes = FALSE,color = "black",size = 2.25,
                         fontface = "italic",segment.color = "black",
                         segment.size = 0.25,max.overlaps = Inf,
                         na.rm = TRUE))
}

# This is used by volcano_plot_with_highlighted_genes to create the
# volcano plot.
ggplot_call_for_volcano_plot <- function (dat, y.label, topic.label)
  ggplot(dat,aes_string(x = "beta",y = "y",fill = "mean",label = "label")) +
    geom_point(color = "white",stroke = 0.3,shape = 21,na.rm = TRUE) +
    scale_y_continuous(trans = "sqrt",
      breaks = c(0,1,2,5,10,20,50,100,200,500,1e3,2e3,5e3,1e4,2e4,5e4)) +
    scale_fill_gradient2(low = "deepskyblue",mid = "gold",high = "orangered",
                         midpoint = mean(range(dat$mean))) +
    geom_text_repel(color = "gray",size = 2.25,fontface = "italic",
                    segment.color = "gray",segment.size = 0.25,
                    max.overlaps = 10,na.rm = TRUE) +
    labs(x = "log-fold change (\u03b2)",y = y.label,fill = "log10 mean") +
    theme_cowplot(font_size = 9) +
    theme(plot.title = element_text(size = 9,face = "plain"))

# TO DO: Explain here what this function does, and how to use it.
loglik_scatterplot <- function (x, y, cluster, k, color = "black",
                                minloglik = -Inf, xlab = "loglik1",
                                ylab = "loglik2") {
  dat <- data.frame(x = x,y = y,cluster = cluster)
  dat <- subset(dat,x > minloglik & y > minloglik)
  xy  <- with(dat,range(c(x,y)))
  dat <- subset(dat,cluster == k)
  return(ggplot(dat,aes_string(x = "x",y = "y")) +
         geom_point(shape = 21,color = "white",fill = color) +
         geom_abline(intercept = 0,slope = 1,linetype = "dotted",
                     color = "black") +
         xlim(xy) +
         ylim(xy) +
         labs(x = xlab,y = ylab,title = sprintf("%s (n = %d)",k,nrow(dat))) +
         theme_cowplot(font_size = 9) +
         theme(plot.title = element_text(size = 9,face = "plain")))
}

# Create a scatterplot comparing two sets of log-fold change
# statistics generated from two different differential expression
# analyses of the same data. Only log-fold change statistics with
# z-scores greater than zmin are shown. Points with the largest
# z-scores (in magnitude) are labeled (controlled by the
# "label_above_score" argument). An additional set of genes can be
# highlighted with the "genes" argument.
beta_scatterplot <- function (res1, res2, k1, k2, genes = NULL,
                              label_above_score = Inf, zmin = 10,
                              betamax = 10,xlab = "beta1",ylab = "beta2") {
  z1   <- res1$Z[,k1]
  z2   <- res2$Z[,k2]
  pdat <- data.frame(mean = cut(res1$colmeans,c(0,0.01,0.1,1,10,Inf)),
                     b1   = res1$beta[,k1],
                     b2   = res2$beta[,k2],
                     gene = rownames(res1$Z),
                     stringsAsFactors = FALSE)
  pdat <- transform(pdat,
                    b1 = sign(b1) * pmin(abs(b1),betamax),
                    b2 = sign(b2) * pmin(abs(b2),betamax))
  rows <- which(!(is.element(rownames(res1$Z),genes) |
                  abs(z1) > label_above_score |
                  abs(z2) > label_above_score))
  pdat[rows,"gene"] <- ""
  rows <- which(abs(z1) > zmin | abs(z2) > zmin)
  pdat <- pdat[rows,]
  return(ggplot(pdat,aes_string(x = "b1",y = "b2",fill = "mean",
                                label = "gene")) +
         geom_point(shape = 21,size = 1.5,color = "white") +
         geom_text_repel(color = "black",size = 2.25,fontface = "italic",
                         segment.color = "black",segment.size = 0.25,
                         max.overlaps = Inf,na.rm = TRUE) +
         scale_fill_manual(values = c("skyblue","cornflowerblue","orange",
                                      "tomato","firebrick")) +
         geom_abline(slope = 1,intercept = 0,color = "black",
                     linetype = "dotted") +
         xlim(range(c(pdat$b1,pdat$b2))) + 
         ylim(range(c(pdat$b1,pdat$b2))) +
         labs(x = xlab,y = ylab,title = "log-fold change (\u03b2)") +
         theme_cowplot(font_size = 10) +
         theme(plot.title = element_text(size = 10,face = "plain")))
}

# Create a scatterplot comparing two sets of z-scores generated from
# two different differential expression analyses of the same data.
# Points with the largest z-scores (in magnitude) are labeled
# (controlled by the "label_above_score" argument). The "zmax" argument
# is useful for creating a nice scatterplot when a small number of the
# z-scores are much larger than the others. An additional set of genes
# can be highlighted with the "genes" argument.
zscores_scatterplot <- function (res1, res2, k1, k2, genes = NULL,
                                 label_above_score = Inf, zmax = Inf,
                                 xlab = "z1", ylab = "z2") {
  z1   <- pmin(res1$Z[,k1],zmax)
  z2   <- pmin(res2$Z[,k2],zmax)
  pdat <- data.frame(mean = cut(res1$colmeans,c(0,0.01,0.1,1,10,Inf)),
                     z1   = z1,
                     z2   = z2,
                     gene = rownames(res1$Z),
                     stringsAsFactors = FALSE)
  rows <- which(!(is.element(rownames(res1$Z),genes) |
                  abs(z1) > label_above_score |
                  abs(z2) > label_above_score))
  pdat[rows,"gene"] <- ""
  pdat <- transform(pdat,
                    z1 = sign(z1)*sqrt(abs(z1)),
                    z2 = sign(z2)*sqrt(abs(z2)))     
  return(ggplot(pdat,aes_string(x = "z1",y = "z2",fill = "mean",
                                label = "gene")) +
         geom_point(shape = 21,size = 1.5,color = "white") +
         geom_text_repel(color = "black",size = 2.25,fontface = "italic",
                         segment.color = "black",segment.size = 0.25,
                         max.overlaps = Inf,na.rm = TRUE) +
         scale_fill_manual(values = c("skyblue","cornflowerblue","orange",
                                      "tomato","firebrick")) +
         geom_abline(slope = 1,intercept = 0,color = "black",
                     linetype = "dotted") +
         xlim(range(c(pdat$z1,pdat$z2))) +
         ylim(range(c(pdat$z1,pdat$z2))) +
         labs(x = xlab,y = ylab,title = "sqrt(z-score)") +
         theme_cowplot(font_size = 10) +
         theme(plot.title = element_text(size = 10,face = "plain")))
}

# Create a basic scatterplot showing the topic proportions projected
# onto two principal components (PCs), and the colour of the points is
# varied according to to a simple "score" for cell-cycle genes.
cellcycle_pca_plot <- function (fit, counts, pcs = 1:2,
                                cell_cycle_genes = c("Cdk1","Ube2c","Top2a"),
                                font_size = 10) {
  X <- counts[,cell_cycle_genes]
  Y <- apply(X,2,qt_random_tie)
  rownames(Y) <- rownames(X)
  score <- rowSums(Y)
  return(suppressMessages(pca_plot(fit,pcs = pcs,fill = score) +
           scale_fill_gradientn(colors = c("darkblue","royalblue",
                                           "lightskyblue","darkorange",
                                           "firebrick"))))
}
