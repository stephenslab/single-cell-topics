# TO DO: Explain here what this script is for, and how to use it.
#
# sinteractive -p mstephens --account=pi-mstephens -c 4 --mem=264G \
#   --time=30:00:00
# module load R/4.1.0
# export MEM_CHECK_INTERVAL=0.01
# python3 monitor_memory.py Rscript deseq2_purified_pbmc.R
# DESeq + lfcShrink took 92,235 s.
# max rss_memory: 262 GB

# Load a few packages.
library(Matrix)
library(scran)
library(DESeq2)

# Initialize the sequence of pseudorandom numbers.
set.seed(1)

# Load the count data.
load("../data/pbmc_purified.RData")

# Remove genes that are expressed in fewer than 10 cells. It is
# doubtful that we will be able to obtain accurate estimates of
# differential expression for these genes anyhow.
j      <- which(colSums(counts > 0) >= 10)
genes  <- genes[j,]
counts <- counts[,j]

i <- "CD19+ B"

# Prepare the UMI count data for analysis with DESeq2.
celltype <- samples$celltype
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
t0 <- proc.time()
deseq <- DESeq(deseq,test = "LRT",reduced = ~1,useT = TRUE,minmu = 1e-6,
               minReplicatesForReplace = Inf)
t1 <- proc.time()
timing <- t1 - t0
cat(sprintf("DESeq + lfcShrink took %0.2f seconds.\n",timing["elapsed"]))

# Save the results.
save(list = c("genes","deseq"),
     file = "deseq2-pbmc-purified.RData")
