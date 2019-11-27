# Notes on Zheng et al PBMC data:
#
#  + Downloaded "Fresh 68k PBMCs (Donor A)" data set from:
#    https://support.10xgenomics.com/single-cell-gene-expression/datasets
#
#  + Associated publication: Zheng et al (2017). Massively parallel
#    digital transcriptional profiling of single cells. Nature
#    Communications 8, 14049. doi:10.1038/ncomms14049
#
#  + Specifically, I downloaded the "Gene/cell matrix (filtered)"
#    tar.gz file. Then I moved the files to
#    data/fresh_68k_pbmc_donor_a_filtered and compressed them with gzip.
#
#  + A summary of the data is given here:
#    http://cf.10xgenomics.com/samples/cell-exp/1.1.0/fresh_68k_pbmc_donor_a/fresh_68k_pbmc_donor_a_web_summary.html
#   
#  + These are the numbers given from that webpage:
#    n = 68,579
#    mean reads per cell = 20,491
#    median read per cell = 1,292
#    median genes per cell = 525
#    number of reads = 1,405,303,609
#
library(Matrix)
library(readr)
library(rsvd)

# SCRIPT PARAMETERS
# -----------------
data.dir      <- "../data/fresh_68k_pbmc_donor_a_filtered"
barcodes.file <- "barcodes.tsv.gz"
genes.file    <- "genes.tsv.gz"
counts.file   <- "matrix.mtx.gz"

# LOAD THE DATA
# -------------
barcodes.file <- file.path(data.dir,barcodes.file)
genes.file    <- file.path(data.dir,genes.file)
counts.file   <- file.path(data.dir,counts.file)
barcodes      <- read_tsv(barcodes.file,col_names = FALSE)[[1]]
genes         <- read_tsv(genes.file,col_names = FALSE)
counts        <- read_delim(counts.file,delim = " ",comment = "%",
                            col_names = FALSE)
class(genes)  <- "data.frame"
class(counts) <- "data.frame"
names(genes)  <- c("ensembl","symbol")
names(counts) <- c("j","i","x")
n      <- counts[1,2]
m      <- counts[1,1]
counts <- counts[-1,]
counts <- sparseMatrix(i = counts$i,j = counts$j,x = counts$x,dims = c(n,m),
                       dimnames = list(sample = barcodes,gene = genes$ensembl))

# EXAMINE DATA
# ------------
