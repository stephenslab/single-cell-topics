rows <- which(x == "D")
fit2 <- select(poisson2multinom(fit),loadings = rows)
plot_grid(labeled_pca_plot(fit2,1:2,samples[rows,"tissue"]),
          labeled_pca_plot(fit2,3:4,samples[rows,"tissue"]),
          labeled_pca_plot(fit2,5:6,samples[rows,"tissue"]),
          labeled_pca_plot(fit2,7:8,samples[rows,"tissue"]),
          labeled_pca_plot(fit2,9:10,samples[rows,"tissue"]))
p13  <- pca_hexbin_plot(fit2,1:2) + guides(fill = "none")
