# This is an exploratory analysis of the topic modeling results
# regarding B-cells in the PBMC data sets.
#
# TO DO:
#
#   + Compare z-scores from topic (k=3) and cluster (B) in purified data.
#
#   + Compare z-scores from topic (k=5) and cluster (A1) in 68k data.
#
#   + In both data sets, examine in greater detail relationship
#     between topic proportions and expression level of key genes, in
#     both purified and 68k data.
#

# Perform differential expression analysis on the clusters identified
# in the purified data.
load("../data/pbmc_purified.RData")
rm(samples,genes)
fit_clusters_purified <-
  init_poisson_nmf_from_clustering(counts,samples_purified$cluster)
diff_count_clusters_purified <- diff_count_analysis(fit_clusters_purified,
                                                    counts)
diff_count_clusters_purified <-
  select_diff_count_res(diff_count_clusters_purified,genes_purified)

# Perform differential expression analysis on the clusters identified
# in the 68k data.
load("../data/pbmc_68k.RData")
rm(samples,genes)
fit_clusters_68k <- 
  init_poisson_nmf_from_clustering(counts,samples_68k$cluster)
diff_count_clusters_68k <- diff_count_analysis(fit_clusters_68k,counts)
diff_count_clusters_68k <- select_diff_count_res(diff_count_clusters_68k,
                                                 genes_68k)

# Compare z-scores from topic 3 and cluster B in the purified data.
pdat <- data.frame(x     = diff_count_res_clusters$Z[,"B"],
                   y     = diff_count_res$Z[,"k3"],
                   label = genes$symbol,
                   stringsAsFactors = FALSE)
pdat[with(pdat,abs(x) < 150 & abs(y) < 150),"label"] <- ""
p <- ggplot(pdat,aes(x = x,y = y,label = label)) +
  geom_point(shape = 21,fill = "dodgerblue",color = "white") +
  geom_text_repel(color = "black",size = 2.25,fontface = "italic",
                  box.padding = 0.1,point.padding = 0.1,
                  segment.color = "black",segment.size = 0.25,
                  na.rm = TRUE) +
  geom_abline(slope = 1,intercept = 0,color = "lightskyblue",
              linetype = "dotted") + 
  theme_cowplot(font_size = 10)
