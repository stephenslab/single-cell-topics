rows <- which(is.element(x,c("C","H","B+C")))
fit2 <- select(poisson2multinom(fit),loadings = rows)
res  <- diff_count_analysis(fit2,counts[rows,])
pca_plot(fit2,k = c(4,5,7))
pca_plot_with_counts(fit2,counts[rows,"Muc5b"],pcs = 1:2,log = TRUE)
res$colmeans <- res$colmeans + 1e-4
plot_grid(volcano_plot(res,k = 5)
          volcano_plot(res,k = 7))
