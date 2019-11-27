# Notes on Zheng et al PBMC data:
#
#  + Downloaded "Fresh 68k PBMCs (Donor A)" data set from:
#    https://support.10xgenomics.com/single-cell-gene-expression/datasets
#
#  + See also: https://community.10xgenomics.com/t5/Data-Sharing/10x-Single-Cell-3-Paper-Zheng-et-al-2016-Datasets/td-p/231
#
#  + Associated publication: Zheng et al (2017). Massively parallel
#    digital transcriptional profiling of single cells. Nature
#    Communications 8, 14049. doi:10.1038/ncomms14049
#
#  + Specifically, I downloaded the "Gene/cell matrix (filtered)"
#    tar.gz file. Then I moved the files to
#    data/fresh_68k_pbmc_donor_a_filtered and compressed them with gzip.
#
#  + Additional information about the samples was downloaded from here:
#    https://github.com/10XGenomics/single-cell-3prime-paper
#    Specifically, file 68k_pbmc_barcodes_annotation.tsv was downloaded.
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
library(ggplot2)
library(cowplot)

# SCRIPT PARAMETERS
# -----------------
data.dir      <- "../data/fresh_68k_pbmc_donor_a_filtered"
samples.file  <- "68k_pbmc_barcodes_annotation.tsv.gz"
genes.file    <- "genes.tsv.gz"
counts.file   <- "matrix.mtx.gz"
out.file      <- "68k_pbmc.RData"

# SET UP ENVIRONMENT
# ------------------
set.seed(1)

# LOAD THE DATA
# -------------
samples.file   <- file.path(data.dir,samples.file)
genes.file     <- file.path(data.dir,genes.file)
counts.file    <- file.path(data.dir,counts.file)
samples        <- read_tsv(samples.file)
genes          <- read_tsv(genes.file,col_names = FALSE)
counts         <- read_delim(counts.file,delim = " ",comment = "%",
                            col_names = FALSE,progress = FALSE)
class(samples) <- "data.frame"
class(genes)   <- "data.frame"
class(counts)  <- "data.frame"
names(samples) <- c("tsne1","tsne2","barcode","celltype")
names(genes)   <- c("ensembl","symbol")
names(counts)  <- c("j","i","x")
samples        <- transform(samples,celltype = factor(celltype))
n      <- counts[1,2]
m      <- counts[1,1]
counts <- counts[-1,]
counts <- sparseMatrix(i = counts$i,j = counts$j,x = counts$x,dims = c(n,m),
                       dimnames = list(sample = samples$barcode,
                                       gene = genes$ensembl))
rm(n,m)

# SUMMARIZE DATA
# --------------
cat(sprintf("Number of samples:     %d\n",nrow(counts)))
cat(sprintf("Number of genes:       %d\n",ncol(counts)))
cat(sprintf("Total number of reads: %d\n",nnzero(counts)))

# EXAMINE DATA
# ------------
ggplot(data.frame(x = log10(1 + colSums(counts > 0))),aes(x)) +
  geom_density(fill = "darkorange",color = "darkorange") +
  scale_x_continuous(breaks = 0:6) +
  labs(x = "log10(1 + num. nonzero read counts)",
       y = "proportion of genes") +
  theme_cowplot(font_size = 12)

# Remove genes that are expressed in fewer than 4 cells.
j      <- which(colSums(counts > 0) >= 4)
genes  <- genes[j,]
counts <- counts[,j]

# SAVE PROCESSED DATA
# -------------------
save(list = c("samples","genes","counts"),file = out.file)
