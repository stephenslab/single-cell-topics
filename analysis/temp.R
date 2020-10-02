Y <- matrix(0,ncol(counts),nlevels(samples$tissue))
colnames(Y) <- levels(samples$tissue)
for (i in levels(samples$tissue)) {
  rows  <- which(samples$tissue == i)
  Y[,i] <- colMeans(counts[rows,])
}
R <- cor(Y)

# -------

plot_grid(labeled_pca_plot(fit,1:2,samples$tissue) +
          xlim(-0.1,0.3) + ylim(-0.8,-0.4),
          pca_plot_with_counts(fit,counts[,"Chga"],1:2,log = TRUE) +
          xlim(-0.1,0.3) + ylim(-0.8,-0.4),
          pca_plot(poisson2multinom(fit),k = 6) +
          xlim(-0.1,0.3) + ylim(-0.8,-0.4),
          nrow = 1)

# -------

labeled_pca_plot(fit,5:6,samples$tissue)

# -------

rows <- which(with(samples_droplet,
                   tissue == "Tuft" | tissue == "Neuroendocrine"))
fit2 <- init_poisson_nmf_from_clustering(counts_droplet[rows,],
          factor(samples_droplet[rows,"tissue"]))
res1 <- diff_count_analysis(fit2,counts_droplet[rows,])

rows <- which(samples_droplet$cluster == "T+N")
fit2 <- select(poisson2multinom(fit_droplet),loadings = rows)
res2 <- diff_count_analysis(fit2,counts_droplet[rows,])

volcano_plot(res1,k = "Tuft")
volcano_plot(res2,k = 2)
beta_scatterplot(res2,res2,"k2","k6",tuft_neuroendocrine_genes,zmin = 1)
