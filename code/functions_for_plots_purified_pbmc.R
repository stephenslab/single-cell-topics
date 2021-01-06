# Automates a few adjustments to the volcano_plot function from
# fastTopics to create the volcano plots in the plots_purified_pbmc
# analysis.
volcano_plot_better <- function (diff_count_res, k, labels,
                                 label_above_quantile = 0.995,
                                 title = "", show_legend = FALSE) {
  i <- which(diff_count_res$beta >= 9.9)
  n <- length(i)
  diff_count_res$beta[i] <- diff_count_res$beta[i] + runif(n)
  out <- volcano_plot(diff_count_res,k,genes$symbol,betamax = 11,
                      label_above_quantile = label_above_quantile) +
    theme(plot.title = element_text(size = 9,face = "plain")) +
    ggtitle(title)
  if (!show_legend)
    out <- out + guides(fill = "none")
  return(out)
}
