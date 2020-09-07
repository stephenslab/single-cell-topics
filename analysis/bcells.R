# This is an exploratory analysis of the topic modeling results
# regarding B-cells in the PBMC data sets.
#
# TO DO:
#
#   + Compare z-scores from topic (k=5) and cluster (A2) in 68k data.
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
  select_diff_count_res(diff_count_clusters_purified,genes_purified$ensembl)

# Perform differential expression analysis on the clusters identified
# in the 68k data.
load("../data/pbmc_68k.RData")
rm(samples,genes)
fit_clusters_68k <- 
  init_poisson_nmf_from_clustering(counts,samples_68k$cluster)
diff_count_clusters_68k <- diff_count_analysis(fit_clusters_68k,counts)
diff_count_clusters_68k <- select_diff_count_res(diff_count_clusters_68k,
                                                 genes_68k$ensembl)

# Compare z-scores from topic 3 and cluster B in the purified data.
pdat <- data.frame(x     = diff_count_clusters_purified$Z[,"B"],
                   y     = diff_count_purified$Z[,"k3"],
                   label = genes_purified$symbol,
                   stringsAsFactors = FALSE)
pdat[with(pdat,abs(x) < 100 & abs(y) < 100),"label"] <- ""
p1 <- ggplot(pdat,aes(x = x,y = y,label = label)) +
  geom_point(shape = 21,fill = "dodgerblue",color = "white") +
  geom_text_repel(color = "black",size = 2.25,fontface = "italic",
                  box.padding = 0.1,point.padding = 0.1,
                  segment.color = "black",segment.size = 0.25,
                  na.rm = TRUE) +
  geom_abline(slope = 1,intercept = 0,color = "lightskyblue",
              linetype = "dotted") + 
  theme_cowplot(font_size = 10)

# Compare z-scores from topic 5 and cluster A2 in the 68k data.
pdat <- data.frame(x     = diff_count_clusters_68k$Z[,"A2"],
                   y     = diff_count_68k$Z[,"k5"],
                   label = genes_68k$symbol,
                   stringsAsFactors = FALSE)
pdat[with(pdat,abs(x) < 100 & abs(y) < 100),"label"] <- ""
p2 <- ggplot(pdat,aes(x = x,y = y,label = label)) +
  geom_point(shape = 21,fill = "dodgerblue",color = "white") +
  geom_text_repel(color = "black",size = 2.25,fontface = "italic",
                  box.padding = 0.1,point.padding = 0.1,
                  segment.color = "black",segment.size = 0.25,
                  na.rm = TRUE) +
  geom_abline(slope = 1,intercept = 0,color = "lightskyblue",
              linetype = "dotted") + 
  theme_cowplot(font_size = 10)

# Examine in greater detail relationship between topic 3 proportions
# and expression level of CD79A in purified data.
load("../data/pbmc_purified.RData")
counts <- counts[,genes_purified$ensembl]
rm(samples,genes)
i    <- which(genes_purified$symbol == "CD79A")
fit2 <- poisson2multinom(fit_purified)
table(cut(fit2$L[samples_purified$cluster == "B","k3"],seq(0,1,0.1)))
table(cut(fit2$L[samples_purified$cluster != "B","k3"],seq(0,1,0.1)))
pdat <- data.frame(x = cut(fit2$L[,"k3"],
                           c(-1,seq(0.1,1,0.1))),
                   y = counts[,i])
p3 <- ggplot(pdat,aes(x = x,y = y)) +
  geom_boxplot(width = 0.25,size = 0.4,outlier.shape = NA) +
  theme_cowplot(font_size = 10)

# Examine in greater detail relationship between topic 5 proportions
# and expression level of CD79A in 68k data.
load("../data/pbmc_68k.RData")
counts <- counts[,genes_68k$ensembl]
rm(samples,genes)
i    <- which(genes_68k$symbol == "CD79A")
fit2 <- poisson2multinom(fit_68k)
table(cut(fit2$L[,"k5"],seq(0,1,0.1)))
table(cut(fit2$L[samples_68k$cluster == "A2","k5"],seq(0,1,0.1)))
table(cut(fit2$L[samples_68k$cluster != "A2","k5"],seq(0,1,0.1)))
pdat <- data.frame(x = cut(fit2$L[,"k5"],
                           c(-1,seq(0.1,1,0.1))),
                   y = counts[,i])
p4 <- ggplot(pdat,aes(x = x,y = y)) +
  geom_boxplot(width = 0.25,size = 0.4,outlier.shape = NA) +
  theme_cowplot(font_size = 10)

