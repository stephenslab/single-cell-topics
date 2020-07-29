# This function imports and prepares the "droplet" data from the
# Montoro et al (2018) paper for analysis in R. The return value is
# list with two list elements: (1) a data frame containing sample
# attributes (specifically, mouse ids, barcodes and tissue labels),
# and (2) an n x m matrix containing gene expression data (UMI
# counts), where n is the number of samples, and m is the number of
# genes.
read_droplet_data <- function (file) {
    
  # Load the gene expression data (i.e., the UMI counts) as an n x m
  # double-precision matrix, where n is the number of samples, and m
  # is the number of genes for which we have expression data (UMI
  # counts).
  suppressMessages(
    counts <- read_delim(file,delim = "\t",progress = FALSE,
                         col_types = cols(.default = col_double(),
                                          gene = col_character())))
  class(counts)    <- "data.frame"
  genes            <- counts$gene
  rownames(counts) <- genes
  counts           <- t(as.matrix(counts[-1]))

  # Extract the barcodes, mouse ids and tissue labels, and create a
  # new data frame containing this info.
  ids       <- strsplit(rownames(counts),"_")
  mouse.ids <- factor(sapply(ids,function (x) x[1]))
  barcodes  <- sapply(ids,function (x) x[2])
  tissues   <- factor(sapply(ids,function (x) x[3]))
  samples   <- data.frame(mouse.id = mouse.ids,
                          barcode  = barcodes,
                          tissue   = tissues,
                          stringsAsFactors = FALSE)

  # Set the row labels to be the combination of the mouse ids and
  # barcodes.
  rownames(counts) <- paste(mouse.ids,barcodes,sep = "_")
  
  # Return a data frame containing the sample attributes (samples),
  # and a matrix containing the gene expression data (counts).
  return(list(samples = samples,counts = counts))
}
