# This file implements functions to load and prepare the 68k and
# FACS-purified PBMC data sets described in the paper, "Massively
# parallel digital transcriptional profiling of single cells" (Zheng
# et al, Nature Communications 8, 2017, doi:10.1038/ncomms14049).

# Read the barcodes from the tab-delimited text file, and return a
# data frame with a single column containing the barcodes.
read_barcodes <- function (filename) {
  out <- suppressMessages(read_tsv(filename,col_names = FALSE))
  class(out) <- "data.frame"
  names(out) <- "barcode"
  return(out)
}

# Read the 68k PBMC barcodes and annotations from the tab-delimited
# text file. The return value is a data frame with four columns: the
# first two columns give the 2-d t-SNE projection; column 3 gives the
# barcode; and the last column gives the estimated "cell-type".
read_barcodes_with_annotations <- function (filename) {
  out <- read_tsv("../data/pbmc_68k/68k_pbmc_barcodes_annotation.tsv.gz")
  class(out) <- "data.frame"
  names(out) <- c("tsne1","tsne2","barcode","celltype")
  return(transform(out,celltype = factor(celltype)))
}
  
# Read the gene information from the tab-delimited text file, and
# return a data frame with two columns containing the Ensembl ids and
# gene symbols.
read_genes <- function (filename) {
  out <- suppressMessages(read_tsv(filename,col_names = FALSE))
  class(out) <- "data.frame"
  names(out) <- c("ensembl","symbol")
  return(out)
}

# Read the count data from a file in "Matrix Market" (.mtx) format,
# and return a n x m sparse matrix containing the counts, where n is
# the number of samples and m is the number of genes. The row and
# column names are set to the barcodes and Ensembl gene ids,
# respectively.
create_counts_matrix <- function (filename, sample.names = NULL,
                                  gene.names = NULL) {
  out <- suppressMessages(read_delim(filename,delim = " ",comment = "%",
                                     col_names = FALSE,progress = FALSE))
  class(out) <- "data.frame"
  names(out) <- c("j","i","x")
  n   <- out[1,2]
  m   <- out[1,1]
  out <- out[-1,]
  return(sparseMatrix(out$i,out$j,x = out$x,dims = c(n,m),
                      dimnames = list(sample = sample.names,
                                      gene = gene.names)))
}

# Read and combine multiple 10x Genomics gene expression data sets,
# each corresponding to "cell type". These data sets are stored in
# "datadir". The "datasets" argument should be a named list in which
# each list element is a subdirectory of "datadir"; the names of the
# elements give the "cell types", or the labels assigned to the
# samples from each of the data sets.
#
# Each of the subdirectories should contain three files:
# barcodes.tsv.gz, a tab-delimited text file containing the sample
# information (read and processed by read_barcodes); genes.tsv.gz,
# another tab-delimited text file containing the gene information
# (read and processed by read_genes); and matrix.mtx,gz, a "Matrix
# Market" file containing count data (read into a sparse matrix by
# create_counts_matrix).
#
# Each data set should have exactly the same gene information (each
# gene.tsv.gz should be the same). Further, the UMI barcodes in each
# data set should be unique.
#
read_purified_pbmc_data <- function (datadir, datasets) {

  # Get the number of data sets to combine.
  n         <- length(datasets)
  celltypes <- names(datasets)
  samples   <- NULL
  counts    <- NULL

  # Repeat for each data set.
  cat("Importing data from these files:\n")
  for (i in 1:n) {

    # Read the barcodes and gene information from the tab-delimited
    # text files, and read the count data from the "Matrix Market"
    # (.mtx) file.
    dataset      <- datasets[[i]]
    genes.file   <- file.path(datadir,dataset,"genes.tsv.gz")
    samples.file <- file.path(datadir,dataset,"barcodes.tsv.gz")
    counts.file  <- file.path(datadir,dataset,"matrix.mtx.gz")
    cat(genes.file,"\n")
    genes <- read_genes(genes.file)
    cat(samples.file,"\n")
    samples1 <- read_barcodes(samples.file)
    cat(counts.file,"\n")
    counts1 <- create_counts_matrix(counts.file,
                 sample.names = paste(samples1$barcode,dataset,sep="-"),
                 gene.names = genes$ensembl)

    # Add the sample information to the larger data frame.
    samples1 <- data.frame(barcode  = samples1$barcode,
                           dataset  = dataset,
                           celltype = celltypes[i],
                           stringsAsFactors = FALSE)
    samples <- rbind(samples,samples1)

    # Add the counts to the larger matrix.
    counts <- rbind(counts,counts1)
  }

  # Return a list containing: (1) "samples", the data frame with
  # information about the samples; (2) "genes", the data frame with
  # information about the genes; and (3) "counts", the n x m sparse
  # matrix of UMI counts.
  samples <- transform(samples,
                       dataset  = factor(dataset,datasets),
                       celltype = factor(celltype,celltypes))
  return(list(samples = samples,
              genes   = genes,
              counts  = counts))
}
