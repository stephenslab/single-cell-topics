# This is the script that was used to perform differential expression
# (DE) analysis for a single (FACS) cell-type in the mixture of
# purified PBMC data, using DESeq2. The DESeq2 analysis took roughly
# 30 h and required a lot of memory, roughly 262 GB. These were the
# steps taken to allocate computing resources on midway2 and run the
# DESeq2 analysis while also assessing memory usage:
#
#   sinteractive -c 4 --mem=264G --time=48:00:00
#   module load R/4.1.0
#   export MEM_CHECK_INTERVAL=0.01
#   python3 monitor_memory.py Rscript deseq2_purified_pbmc.R
#

# Load a few packages.
library(Matrix)
library(scran)
library(DESeq2)
i <- "CD19+ B" # CD19+ B, CD56+ NK, T cell, CD14+ Monocyte, CD34+
outfile <- "deseq2-pbmc-purified-bcells.RData"
# deseq2-pbmc-purified-bcells.RData
# deseq2-pbmc-purified-nkcells.RData
# deseq2-pbmc-purified-tcells.RData
# deseq2-pbmc-purified-cd14+.RData
# deseq2-pbmc-purified-cd34+.RData
print(i)
print(outfile)

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

# Prepare the UMI count data for analysis with DESeq2.
celltype <- samples$celltype
if (i == "T cell") {
  t_cell_subtypes <- c("CD4+ T Helper2",
                       "CD8+ Cytotoxic T",
                       "CD4+/CD45RO+ Memory",
                       "CD8+/CD45RA+ Naive Cytotoxic",
                       "CD4+/CD45RA+/CD25- Naive T",
                       "CD4+/CD25 T Reg")
  coldata <- data.frame(celltype=factor(is.element(celltype,t_cell_subtypes)))
} else {
  coldata <- data.frame(celltype = factor(celltype == i))
}
levels(coldata$celltype) <- 1:2
print(summary(coldata$celltype))
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
deseq <- lfcShrink(deseq,coef = "celltype_2_vs_1",type = "ashr")
t1 <- proc.time()
timing <- t1 - t0
cat(sprintf("DESeq + lfcShrink took %0.2f seconds.\n",timing["elapsed"]))

# Save the results.
save(list = c("genes","deseq"),
     file = outfile)

