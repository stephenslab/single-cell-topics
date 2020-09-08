# TO DO: Explain here what this function does, and how to use it.
diff_count_scatterplot <- function (z1, z2, colmeans, genes,
                                    label_above_score = 200) {
  pdat <- data.frame(mean = cut(colmeans,c(0,0.01,0.1,1,10,Inf)),
                     z1   = z1,
                     z2   = z2,
                     gene = genes,
                     stringsAsFactors = FALSE)
  rows <- which(abs(z1) < label_above_score & abs(z2) < label_above_score)
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
