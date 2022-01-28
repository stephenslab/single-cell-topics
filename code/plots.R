# A simple scatterplot to compare likelihoods x and y computed from
# models.
loglik_scatterplot <- function (x, y, cluster, k, minloglik = -Inf,
                                xlab = "loglik1", ylab = "loglik2") {
  dat <- data.frame(x = x,y = y,cluster = cluster)
  dat <- subset(dat,x > minloglik & y > minloglik)
  xy  <- with(dat,range(c(x,y)))
  dat <- subset(dat,cluster == k)
  return(ggplot(dat,aes_string(x = "x",y = "y")) +
         geom_point(shape = 21,color = "white",fill = "dodgerblue") +
         geom_abline(intercept = 0,slope = 1,linetype = "dotted",
                     color = "black") +
         xlim(xy) +
         ylim(xy) +
         labs(x = xlab,y = ylab,title = sprintf("%s (n = %d)",k,nrow(dat))) +
         theme_cowplot(font_size = 9) +
         theme(plot.title = element_text(size = 9,face = "plain")))
}

