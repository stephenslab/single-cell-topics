# Create a basic scatterplot showing the topic proportions projected
# onto two principal components (PCs).
basic_pca_plot <- function (fit, pcs = 1:2) {
  if (inherits(fit,"poisson_nmf_fit"))
    fit <- poisson2multinom(fit)
  dat <- as.data.frame(prcomp(fit$L)$x)
  if (is.numeric(pcs))
    pcs <- names(dat)[pcs]
  return(ggplot(dat,aes_string(x = pcs[1],y = pcs[2])) +
         geom_point(shape = 21,color = "white",fill = "black",size = 1.25) +
         theme_cowplot(font_size = 10))
}

# Create a basic scatterplot showing the topic proportions projected
# onto two principal components (PCs), and the colour of the points is
# varied according to a factor ("labels").
labeled_pca_plot <-
  function (fit, pcs = 1:2, labels,
            colors = c("firebrick","dodgerblue","forestgreen","darkmagenta",
                       "darkorange","gold","darkblue","peru","greenyellow"),
            legend_label = "label", font_size = 10) {
  if (inherits(fit,"poisson_nmf_fit"))
    fit <- poisson2multinom(fit)
  dat <- as.data.frame(prcomp(fit$L)$x)
  if (is.numeric(pcs))
    pcs <- names(dat)[pcs]
  dat <- cbind(data.frame(label = factor(labels)),dat)
  return(ggplot(dat,aes_string(x = pcs[1],y = pcs[2],fill = "label")) +
         geom_point(shape = 21,color = "white",size = 1.2,na.rm = TRUE) +
         scale_fill_manual(values = colors) +
         labs(fill = legend_label) +
         theme_cowplot(font_size = font_size))
}

# Create a "hexbin plot" showing the density of the data points
# (specifically, the topic proportions) as they are projected onto two
# principal components (PCs).
pca_hexbin_plot <-
  function (fit, pcs = 1:2, n = 40, bins = c(0,1,10,100,1000,Inf),
            colors = c("gainsboro","lightskyblue","gold","orange","magenta")) {
  if (inherits(fit,"poisson_nmf_fit"))
    fit <- poisson2multinom(fit)
  dat <- as.data.frame(prcomp(fit$L)$x)
  if (is.numeric(pcs))
    pcs <- names(dat)[pcs]
  return(ggplot(dat,aes_string(x = pcs[1],y = pcs[2])) +
         stat_bin_hex(mapping = aes_q(fill = quote(cut(..count..,bins))),
                      bins = n) +
         scale_fill_manual(values = colors) +
         labs(fill = "count") +
         theme_cowplot(font_size = 10))
}

# TO DO: Explain here what this function does, and how to use it.
volcano_plot_with_highlighted_genes <- function (diff_count_res, k, genes,
                                                 ...) {
  dat <- data.frame(beta  = diff_count_res$beta[genes,k],
                    y     = diff_count_res$Z[genes,k],
                    label = genes)
  return(volcano_plot(diff_count_res,k = k,
                      ggplot_call = ggplot_call_for_volcano_plot,...) +
         geom_text_repel(data = dat,
                         mapping = aes(x = beta,y = y,label = label),
                         inherit.aes = FALSE,color = "black",size = 2.25,
                         fontface = "italic",segment.color = "black",
                         segment.size = 0.25,max.overlaps = Inf,
                         na.rm = TRUE))
}

# TO DO: Explain here what this function is for, and how to use it.
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
    theme_cowplot(font_size = 9)

# Create a scatterplot comparing two sets of log-fold change
# statistics generated from two different differential expression
# analyses of the same data. Argument "betamin" is useful for not
# including very large, negative log-fold change ("beta") estimates.
logfoldchange_scatterplot <- function (x, y, colmeans, betamin = -20) {
  i        <- which(x >= betamin & y >= betamin)
  x        <- x[i]
  y        <- y[i]
  colmeans <- colmeans[i]
  pdat     <- data.frame(mean = cut(colmeans,c(0,0.01,0.1,1,10,Inf)),
                         x = x,y = y,stringsAsFactors = FALSE)
  return(ggplot(pdat,aes_string(x = "x",y = "y",fill = "mean")) +
         geom_point(shape = 21,size = 1.5,color = "white") +
         scale_fill_manual(values = c("skyblue","cornflowerblue","orange",
                                      "tomato","firebrick")) +
         geom_abline(slope = 1,intercept = 0,color = "black",
                     linetype = "dotted") +
         xlim(range(c(x,y))) + 
         ylim(range(c(x,y))) +
         theme_cowplot(font_size = 10) +
         theme(plot.title = element_text(size = 10,face = "plain")))
}

# Create a scatterplot comparing two sets of z-scores generated from
# two different differential expression analyses of the same data.
# Points with the largest z-scores (in magnitude) are labeled
# (according ot the "label_above_score" argument). The "zmax" argument
# is useful when a small number of the z-scores are much larger than
# the others.
zscores_scatterplot <- function (z1, z2, colmeans, genes,
                                 label_above_score = 100,
                                 zmax = 1000) {
  z1   <- pmin(z1,zmax)
  z2   <- pmin(z2,zmax)
  pdat <- data.frame(mean = cut(colmeans,c(0,0.01,0.1,1,10,Inf)),
                     z1 = z1,z2 = z2,gene = genes,stringsAsFactors = FALSE)
  rows <- which(abs(z1) < label_above_score &
                abs(z2) < label_above_score)
  pdat[rows,"gene"] <- ""
  return(ggplot(pdat,aes_string(x = "z1",y = "z2",fill = "mean",
                                label = "gene")) +
         geom_point(shape = 21,size = 1.5,color = "white") +
         geom_text_repel(color = "black",size = 2.25,fontface = "italic",
                         segment.color = "black",segment.size = 0.25,
                         na.rm = TRUE) +
         scale_fill_manual(values = c("skyblue","cornflowerblue","orange",
                                      "tomato","firebrick")) +
         geom_abline(slope = 1,intercept = 0,color = "black",
                     linetype = "dotted") +
         xlim(range(c(z1,z2))) +
         ylim(range(c(z1,z2))) +
         theme_cowplot(font_size = 10) +
         theme(plot.title = element_text(size = 10,face = "plain")))
}
