# TO DO: Explain here what this function does, and how to use it.
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
                         box.padding = 0.1,point.padding = 0.1,
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
