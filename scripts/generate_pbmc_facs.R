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
subpop <- as.character(samples$celltype)
subpop[subpop == "CD19+ B"] <- "B cell"
subpop[subpop == "CD14+ Monocyte"] <- "CD14+"
subpop[subpop == "CD34+ Monocyte"] <- "CD34+"
subpop[subpop == "CD56+ NK"] <- "NK cell"
subpop[subpop == "CD4+ T Helper2" |
       subpop == "CD8+ Cytotoxic T" |
       subpop == "CD8+/CD45RA+" |
       subpop == "CD4+/CD45RO+ Memory" |
       subpop == "CD8+/CD45RA+ Naive Cytotoxic" |
       subpop == "CD4+/CD45RA+/CD25- Naive T" |
       subpop == "CD4+/CD25 T Reg" |
       subpop == "Naive Cytotoxic"] <- "T cell"
samples$subpop <- factor(subpop)

# Select a random subset of the cells.
x <- samples$cluster
i <- c(sample(which(x == "T"),1212),
       sample(which(x != "T" & x != "dendritic"),2562))
i <- sample(i)
j <- c(sample(setdiff(which(x == "T"),i),321),
       sample(setdiff(which(x != "T" & x != "dendritic"),i),679))
j <- sample(j)
samples_test <- samples[j,]
counts_test  <- counts[j,]
samples      <- samples[i,]
counts       <- counts[i,]

# Remove genes with no expression in the training set.
j           <- which(colSums(counts) > 0)
genes       <- genes[j,]
counts      <- counts[,j]
counts_test <- counts_test[,j]

# Save the data.
samples <- samples[,c("barcode","celltype","subpop")]
samples_test <- samples_test[,c("barcode","celltype","subpop")]
rownames(samples) <- NULL
rownames(samples_test) <- NULL
pbmc_facs <- list(samples      = samples,
                  samples_test = samples_test,
                  genes        = genes,
                  counts       = counts,
                  counts_test  = counts_test)
save(list = "pbmc_facs",file = "pbmc_facs.RData",version = 2)
resaveRdaFiles("pbmc_facs.RData")
