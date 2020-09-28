p1 <- pca_plot_with_counts(fit2,counts[rows,"Cdk1"],5:6,log = TRUE) +
  labs(title = "Cdk1")
p2 <- pca_plot_with_counts(fit2,counts[rows,"Ube2c"],5:6,log = TRUE) +
  labs(title = "Ube2c")
p3 <- pca_plot_with_counts(fit2,counts[rows,"Top2a"],5:6,log = TRUE) +
  labs(title = "Top2a")
plot_grid(p1,p2,p3,nrow = 1)
