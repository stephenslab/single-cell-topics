fit_merged <- merge_topics(fit,c("k3","k5","k10"))
fit_merged <- merge_topics(fit_merged,c("k2","k7"))
fit_merged <- merge_topics(fit_merged,c("k9","k11"))
fit_merged <- merge_topics(fit_merged,c("k6","k8"))
rows <- which(samples$cluster == "CD8+" |
              samples$cluster == "TA" |
              samples$cluster == "TB" |
              samples$cluster == "NK")
fit_t <- select_loadings(fit,loadings = rows)
out1 <- diff_count_analysis(fit_merged,counts)
out2 <- diff_count_analysis(fit,counts)
out3 <- diff_count_clusters(samples$cluster,counts)
out4 <- diff_count_analysis(fit_t,counts[rows,])
volcano_plot(out3,k = "CD8+",labels = genes$symbol,
             label_above_quantile = 0.998)
