# NOTE: Run plots_pbmc analysis first before running this code.
load("../data/pbmc_68k.RData")
counts_68k <- counts
genes_68k <- genes
rm(samples,genes,counts)

# Perform differential expression analysis with the topic model.
diff_count_68k <- diff_count_analysis(fit_68k,counts_68k)

# Perform differential expression analysis using the clusters.
fit_clusters_68k <-
  init_poisson_nmf_from_clustering(counts_68k,samples_68k$cluster)
diff_count_clusters_68k <- diff_count_analysis(fit_clusters_68k,counts_68k)

# Volcano plot for B-cells.
p1 <- volcano_plot(diff_count_clusters_68k,k = "A2",labels = genes_68k$symbol,
                   label_above_quantile = 0.995)
print(p1)
