# A short script used to perform the gene set enrichment analysis
# using the DE analysis results on the purified PBMC data, with k = 6
# topics. These were the steps taken to load R and allocate computing
# resources for this analysis:
#
#   sinteractive -p broadwl -c 20 --mem=36G --time=60:00:00
#   module load R/3.5.1
#
library(Matrix)
library(tools)
library(susieR)
library(pathways)
set.seed(1)

# Load the gene sets.
data(gene_sets_human)
X <- gene_sets_human$gene_sets

# Load the results of the DE analysis.
load("../output/pbmc-purified/de-pbmc-purified-seed=1.RData")
Y <- de$postmean

# Remove gene sets in several MSigDB collections that clearly aren't
# relevant.
i <- which(!is.element(gene_sets_human$gene_set_info$database,
                       c("MSigDB-ARCHIVED","MSigDB-C1","MSigDB-C3",
                         "MSigDB-C4","MSigDB-C6")))
X <- X[,i]

# Align the gene-set data with the gene-wise statistics.
ids <- intersect(rownames(X),rownames(Y))
X <- X[ids,]
Y <- Y[ids,]

# Next, remove gene sets with fewer than 10 genes and with more than
# 400 genes. Gene sets with a large number of genes are less likely to
# be interesting, and slow down the enrichment analysis, so they are
# removed.
i <- which(colSums(X) >= 10 & colSums(X) <= 400)
X <- X[,i]

# Perform a gene set enrichment analysis using susieR.
k <- ncol(Y)
gsea <- vector("list",k)
names(gsea) <- colnames(Y)
t0 <- proc.time()
for (i in 1:k)
  gsea[[i]] <- susie(X,Y[,i],L = 10,intercept = TRUE,standardize = FALSE,
                     estimate_residual_variance = TRUE,refine = FALSE,
                     compute_univariate_zscore = FALSE,verbose = TRUE)
t1 <- proc.time()
timing <- t1 - t0
cat(sprintf("Computation took %0.2f seconds.\n",timing["elapsed"]))

# Save the results to file.
save(list = "gsea",file = "gsea-pbmc-purified.RData")
resaveRdaFiles("gsea-pbmc-purified.RData")
