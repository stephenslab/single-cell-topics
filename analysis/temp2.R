volcano_plot(diff_count_facs,k = "T cell",labels = genes$symbol,
             label_above_quantile = 0.995)
volcano_plot(diff_count_topics,"k3",genes$symbol,
             label_above_quantile = 0.998) +
  scale_x_continuous(expand = expansion(mult = 0.5))
