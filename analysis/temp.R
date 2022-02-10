library(cowplot)
clusters <- as.character(samples$tissue)
clusters[clusters == "Basal" |
         clusters == "Club" |
         clusters == "Goblet"] <- "Other"
clusters <- factor(clusters)
pca_plot(fit,pcs = 1:2,fill = clusters)
set.seed(1)
p1 <- structure_plot(fit,
                     topics = c(4,5,2,3,1),
                     colors = c("firebrick","darkorange","salmon","gainsboro",
                                "skyblue"),
                     perplexity = 70)

# ----

volcano_plot(de_merged,k = "k1+k3",ymax = 200) # Ciliated
volcano_plot(de_merged,k = "k2",ymax = 100) # Neuroendocrine + ionocyte
volcano_plot(de_merged,k = "k5",ymax = 100) # Tuft
plot_grid(volcano_plot(de,k = "k1",ymax = 60),
          volcano_plot(de,k = "k3",ymax = 60))


