# Create a scatterplot comparing two sets of log-fold change (LFC)
# statistics. LFC estimates are not plotted if both z-scores are less
# than zmin in magnitude. Data points in which the LFC estimates and
# z-scores from the second set of results (res2) exceeding
# label_above_lfc and label_above_quantile, respectively, are labeled.
lfc_scatterplot <- function (res1, res2, k1, k2, labels, zmin = 2.5,
                             betamax = 15, label_above_lfc = 0,
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
  dat$x    <- pmin(pmax(dat$x,-betamax),betamax)
  dat$y    <- pmin(pmax(dat$y,-betamax),betamax)
  dat$label[abs(z2) < z0 | dat$y < label_above_lfc] <- ""
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
         scale_x_continuous(limits = c(-betamax,betamax + 1),
                            breaks = seq(-20,20,5)) +
         scale_y_continuous(limits = c(-betamax,betamax + 1)
                            ,breaks = seq(-20,20,5)) +
         scale_fill_manual(values = c("skyblue","cornflowerblue","orange",
                                      "tomato","firebrick")) +
         labs(x = xlab,y = ylab) +
         theme_cowplot(font_size = 9))
}
