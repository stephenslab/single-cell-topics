library(uwot)
set.seed(1)
out <- umap(poisson2multinom(fit_droplet)$L,n_neighbors = 10,n_epochs = 1000,
            min_dist = 0.1,scale = "none",learning_rate = 1,verbose = TRUE)
colnames(out) <- c("d1","d2")
pca_plot(fit_droplet,out.pca = list(x = out),fill = samples_droplet$cluster)
