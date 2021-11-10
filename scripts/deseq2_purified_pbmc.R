# TO DO: Explain here what this script is for, and how to use it.
# sinteractive -p mstephens --account=pi-mstephens -c 4 --mem=24G \
#   --time=24:00:00
# module load R/4.1.0

# Load a few packages.
library(Matrix)
library(scran)
library(DESeq2)

# Initialize the sequence of pseudorandom numbers.
set.seed(1)

# Load the count data.
load("../data/pbmc_purified.RData")

i <- "CD19+ B"

# Prepare the UMI count data for analysis with DESeq2.
t0 <- proc.time()
coldata <- data.frame(celltype = factor(celltype == i))
levels(coldata$celltype) <- 1:2
counts <- t(counts)
deseq <- DESeqDataSetFromMatrix(counts,coldata,~celltype)
sizeFactors(deseq) <- calculateSumFactors(counts)

# Now we perform the DE analysis using DESeq2, using the settings
# recommended for single-cell RNA-seq data (see the main DESeq2
# vignette for details). To replicate the fastTopics analysis as
# closely as possible, we also shrink the LFC estimates using adaptive
# shrinkage.
deseq <- DESeq(deseq,test = "LRT",reduced=~1,useT = TRUE,minmu = 1e-6,
               minReplicatesForReplace = Inf)
deseq <- lfcShrink(deseq,coef = "cluster_2_vs_1",type = "ashr")
t1 <- proc.time()
cat(sprintf("Computation took %0.2f seconds.\n",(t1 - t0)["elapsed"]))

# Save the results.
save(list = c("genes","deseq"),
     file = "deseq2-pbmc-purified.RData")
resaveRdaFiles("deseq2-pbmc-purified.RData")
