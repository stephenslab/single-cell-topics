# This function imports and prepares the "pulse-seq" data from the
# Montoro et al (2018) paper for analysis in R. The return value is
# list with two list elements: (1) a data frame containing sample
# attributes (specifically, sample ids, barcodes, inferred tissue
# labels, mouse ids, lineage labels, and profiling time-point), and
# (2) an n x p matrix containing gene expression data (read counts),
# where n is the number of samples, and p is the number of genes.
read_pulseseq_data <- function (file) {
  counts <- t(readRDS(file))

  # Extract some information from the sample ids: barcode, tissue,
  # mouse, lineage and time-point.
  ids <- rownames(counts)
  ids <- strsplit(ids,"_")
  samples <-
    data.frame(id      = sapply(ids,function (x) paste(x[1:4],collapse = "-")),
               barcode = sapply(ids,function (x) x[1]),
               tissue  = factor(sapply(ids,function (x) tolower(x[5]))),
               mouse   = factor(sapply(ids,function (x) x[3])),
               lineage = factor(sapply(ids,function (x) x[4])),
               tp = factor(sapply(ids,function(x) substr(x[2],3,length(x)))),
               stringsAsFactors = FALSE)
  rownames(counts) <- samples$id

  # Return a data frame containing the sample attributes (samples),
  # and a matrix containing the gene expression data (counts).
  return(list(samples = samples,counts = counts))
}
