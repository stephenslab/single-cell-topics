p1 <- pca_plot_with_counts(fit2,counts[rows,"Cdk1"],5:6,log = TRUE) +
  labs(title = "Cdk1")
p2 <- pca_plot_with_counts(fit2,counts[rows,"Ube2c"],5:6,log = TRUE) +
  labs(title = "Ube2c")
p3 <- pca_plot_with_counts(fit2,counts[rows,"Top2a"],5:6,log = TRUE) +
  labs(title = "Top2a")
plot_grid(p1,p2,p3,nrow = 1)

qt_random_tie <- function (x) {
  y = rank(x,ties.method = "random")
  return(qqnorm(y,plot.it = FALSE)$x)
}
X <- counts[rows,c("Cdk1","Ube2c","Top2a")]
Y <- apply(X,2,qt_random_tie)
rownames(Y) <- rownames(X)           
score <- rowSums(Y)

pdat <- as.data.frame(prcomp(fit2$L)$x)
pdat <- cbind(pdat,data.frame(score = score))
p4 <- ggplot(pdat,aes(x = PC5,y = PC6,fill = score)) +
  geom_point(shape = 21,color = "white",size = 1.25) +
  scale_fill_gradientn(colors = c("darkblue","royalblue",
                       "lightskyblue","darkorange","firebrick")) +
  theme_cowplot(font_size = 10)
