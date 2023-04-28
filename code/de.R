# Generate a data frame containing the results of a GoM DE analysis.
compile_de_table <- function (de, topics = colnames(de$postmean)) {
  out <- NULL
  for (k in topics) {
    x <- data.frame(topic   = k,
                    gene    = rownames(de$postmean),
                    lfc     = de$postmean[,k],
                    z       = de$z[,k],
                    lfsr    = de$lfsr[,k],
                    stringsAsFactors = FALSE)
    out <- rbind(out,x)
  }
  return(transform(out,topic = factor(topic)))
}

