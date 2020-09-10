# NKG7 in purified data.
j  <- which(diff_count_purified$colmeans > 0.01  &
            diff_count_clusters_purified$beta[,"A1"] > -20)
p5 <- diff_count_scatterplot(diff_count_clusters_purified$beta[j,"A1"],
                             diff_count_purified$beta[j,4],
			     diff_count_purified$colmeans[j],
			     genes_purified$symbol[j],
			     label_above_score = 10) +
  labs(x = "cluster A1",y = "topic 4",title = "log-fold change (beta)")
print(p5)
