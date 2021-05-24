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
x <- samples$cluster
i <- c(sample(which(x == "T"),1212),
       sample(which(x != "T" & x != "dendritic"),2562))
j <- c(sample(setdiff(which(x == "T"),i),321),
       sample(setdiff(which(x != "T" & x != "dendritic"),i),679))
j <- sample(j)
samples_test <- samples[j,]
counts_test  <- counts[j,]
samples      <- samples[i,]
counts       <- counts[i,]

# Remove genes with no expression in the training set.
j           <- which(colSums(counts) > 0)
genes       <- genes[i,]
counts      <- counts[,j]
counts_test <- counts_test[,j]

# Save the data.
save(list = c("samples","samples_test","genes","counts","counts_test"),
     file = "pbmc_4k.RData")
resaveRdaFiles("pbmc_4k.RData")
