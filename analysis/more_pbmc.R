# NOTE: Run plots_pbmc analysis first before running this code.
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
