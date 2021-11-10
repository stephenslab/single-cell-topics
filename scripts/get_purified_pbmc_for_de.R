# TO DO: Explain here what this script is for, and how to use it.

# Load a couple packages.
library(tools)
library(Matrix)

# Initialize the sequence of pseudorandom numbers.
set.seed(1)

# Load the count data.
load("../data/pbmc_purified.RData")

# Create a single label for all the T cells.
celltype <- as.character(samples$celltype)
celltype[celltype == "CD4+ T Helper2" |
  celltype == "CD4+/CD45RO+ Memory" |
  celltype == "CD8+/CD45RA+ Naive Cytotoxic" |
  celltype == "CD4+/CD45RA+/CD25- Naive T" | 
  celltype == "CD4+/CD25 T Reg" |
  celltype == "CD8+ Cytotoxic T"] <- "T cell"
celltype <- factor(celltype)

# Subsample the cell types so that there are at most 4,000 of each
# cell type.
i <- c(which(celltype == "CD14+ Monocyte"),
       sample(which(celltype == "CD19+ B"),4000),
       sample(which(celltype == "CD34+"),4000),
       sample(which(celltype == "CD56+ NK"),4000),
       sample(which(celltype == "T cell"),4000))
i <- sort(i)
celltype <- celltype[i]
samples  <- samples[i,]
counts   <- counts[i,]

# Remove genes that are expressed in fewer than 10 cells. It is
# doubtful that we will be able to obtain accurate estimates of
# differential expression for these genes anyhow.
j      <- which(colSums(counts > 0) >= 10)
genes  <- genes[j,]
counts <- counts[,j]

# Save the new data set.
save(list = c("genes","celltype","samples","counts"),
     file = "pbmc_purified_for_de.RData")
resaveRdaFiles("pbmc_purified_for_de.RData")
