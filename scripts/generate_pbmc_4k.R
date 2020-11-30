# Generate the single-cell RNA-seq data set used for the fastTopics
# demo.
library(tools)
library(Matrix)
library(fastTopics)
set.seed(1)

# Load the UMI count data.
load("../data/pbmc_purified.RData")

# Load the clustering.
samples <- readRDS("../output/pbmc-purified/clustering-pbmc-purified.rds")

# Select a random subset of the cells.
x       <- samples$cluster
i       <- c(sample(which(x == "T"),1212),
             sample(which(x != "T" & x != "dendritic"),2562))
i       <- sample(i)
samples <- samples[i,]
counts  <- counts[i,]

# Save the data.
# samples <- samples[c("barcode","dataset","celltype")]
save(list = c("samples","genes","counts"),
     file = "pbmc_4k.RData")
resaveRdaFiles("pbmc_4k.RData")
