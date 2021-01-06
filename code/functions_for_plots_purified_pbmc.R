# Automates a few adjustments to the volcano_plot function from
# fastTopics to create the volcano plots in the plots_purified_pbmc
# analysis.
volcano_plot_better <- function (diff_count_res, k, labels,
                                 label_above_quantile = 0.995,
                                 title = "", show_legend = FALSE) {
  i <- which(diff_count_res$beta >= 9.9)
  n <- length(i)
  diff_count_res$beta[i] <- diff_count_res$beta[i] + runif(n)
  out <- volcano_plot(diff_count_res,k,genes$symbol,betamax = 11,
                      label_above_quantile = label_above_quantile) +
    theme(plot.title = element_text(size = 9,face = "plain")) +
    ggtitle(title)
  if (!show_legend)
    out <- out + guides(fill = "none")
  return(out)
}

# Create a scatterplot comparing two sets of log-fold change
# statistics (beta). Estimates of beta are not used if both z-scores
# are less than zmin in magnitude. Data points in which the z-scores
# from the second set of results (res2) exceeding the specified
# quantile (label_above_quantile) are labeled.
#
# Additional details: beta estimates less than -10 are set to -10;
# beta estimates 10 or greater are randomly sampled to be between 10
# and 11 to better visualize the many data points with large log-fold
# change values.
lfc_scatterplot <- function (res1, res2, k1, k2, labels, zmin = 6,
                             label_above_quantile = 0.995,
                             xlab = "", ylab = "") {
  z2  <- abs(res2$Z[,k2])
  z0  <- quantile(z2,label_above_quantile)
  dat <- data.frame(label = labels,
                    mean  = res1$colmeans,
                    x     = res1$beta[,k1],
                    y     = res2$beta[,k2],
                    stringsAsFactors = FALSE)
  dat$mean <- cut(dat$mean,c(0,0.001,0.01,0.1,1,Inf))
  dat$x    <- pmax(dat$x,-10)
  dat$y    <- pmax(dat$y,-10)
  i        <- which(dat$x >= 9.9)
  dat$x[i] <- dat$x[i] + runif(length(i))
  i        <- which(dat$y >= 9.9)
  dat$y[i] <- dat$y[i] + runif(length(i))
  dat$label[abs(z2) < z0] <- ""
  rows <- which(abs(res1$Z[,k1]) > zmin |
                abs(res2$Z[,k2]) > zmin)
  dat <- dat[rows,]
  return(ggplot(dat,aes_string(x = "x",y = "y",fill = "mean",
                               label = "label")) +
         geom_point(shape = 21,color = "white") +
         geom_abline(slope = 1,intercept = 0,linetype = "dotted") +
         geom_text_repel(color = "black",size = 2,fontface = "italic",
                         segment.color = "black",segment.size = 0.25,
                         max.overlaps = Inf,na.rm = TRUE) +
         scale_x_continuous(limits = c(-10,11),breaks = seq(-10,10,5)) +
         scale_y_continuous(limits = c(-10,11),breaks = seq(-10,10,5)) +
         scale_fill_manual(values = c("skyblue","cornflowerblue","orange",
                                      "tomato","firebrick")) +
         labs(x = xlab,y = ylab) +
         theme_cowplot(font_size = 9))
}
