# TO DO: Explain here what this package is for, and how to use it.

# Load a few packages.
library(Matrix)
library(Seurat)

# Initialize the sequence of pseudorandom numbers.
set.seed(1)

# Load the count data.
load("pbmc_purified_for_de.RData")

# Remove genes that are expressed in fewer than 10 cells. It is
# doubtful that we will be able to obtain accurate estimates of
# differential expression for these genes anyhow.
j      <- which(colSums(counts > 0) >= 10)
genes  <- genes[j,]
counts <- counts[,j]

i <- "CD19+ B"

counts <- t(counts)
names(celltype) <- colnames(counts)
pbmc <- CreateSeuratObject(counts = counts)
Idents(pbmc) <- celltype
de <- FindMarkers(pbmc,ident.1 = "CD19+ B",ident.2 = NULL,test.use = "MAST")
