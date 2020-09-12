# NOTE: Run plots_pbmc analysis first before running this code.
library(ggrepel)
load("../data/pbmc_purified.RData")
counts_purified <- counts
genes_purified <- genes
rm(samples,genes,counts)

# Perform differential expression analysis with the topic model
diff_count_purified <- diff_count_analysis(fit_purified,counts_purified)

# Perform differential expression analysis using the clusters
fit_clusters_purified <-
  init_poisson_nmf_from_clustering(counts_purified,samples_purified$cluster)
diff_count_clusters_purified <- diff_count_analysis(fit_clusters_purified,
                                                    counts_purified)

# Volcano plot for B-cells.
p1 <- volcano_plot(diff_count_clusters_purified,k = "B",
                   labels = genes_purified$symbol,
                   label_above_quantile = 0.998)
print(p1)

# Volcano plot for natural killer (NK) cells.
p2 <- volcano_plot(diff_count_clusters_purified,k = "NK",
                   labels = genes_purified$symbol,
                   label_above_quantile = 0.995)
print(p2)

# Volcano plot for T-cells.
p3 <- volcano_plot(diff_count_clusters_purified,k = "T",
                   labels = genes_purified$symbol,
                   label_above_quantile = 0.995)
print(p3)

# Volcano plot for CD34+ cells.
genes <- genes_purified$symbol
i     <- which(diff_count_clusters_purified$beta[,"CD34+"] < 7.5 &
               diff_count_clusters_purified$Z[,"CD34+"] < 80)
genes[i] <- ""
p4 <- volcano_plot(diff_count_clusters_purified,k = "CD34+",
                   labels = genes,label_above_quantile = 0)
print(p4)

# Volcano plot for CD14+ monocytes.
p5 <- volcano_plot(diff_count_clusters_purified,k = "CD14+",
                   labels = genes_purified$symbol,
                   label_above_quantile = 0.995)

# Volcano plot for CD8+ cytotoxic cells.
p6 <- volcano_plot(diff_count_clusters_purified,k = "CD8+",
                   labels = genes_purified$symbol,
                   label_above_quantile = 0.998)
print(p6)

# Volcano plots for topics 1 and 6 in T-cells.
rows <- which(samples_purified$cluster == "T")
fit  <- select(poisson2multinom(fit_purified),loadings = rows)
diff_count_purified_T <- diff_count_analysis(fit,counts_purified[rows,])
diff_count_purified_T$colmeans <- pmax(diff_count_purified_T$colmeans,0.001)
p7a <- volcano_plot(diff_count_purified_T,k = 1,
                    labels = genes_purified$symbol,
                    label_above_quantile = 0.995)
p7b <- volcano_plot(diff_count_purified_T,k = 6,
                    labels = genes_purified$symbol,
                    label_above_quantile = 0.995)
plot_grid(p7a,p7b)
