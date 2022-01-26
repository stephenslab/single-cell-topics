# TO DO: Explain here what this script does, and how to use it.
library(Matrix)
library(fgsea)
library(iDEA)
source("../code/gsea.R")
set.seed(1)

# Script parameters
k <- "k3"

# Load the gene sets.
load("../data/gene_sets_human.RData")
gene_sets <- gene_sets_human$gene_sets

# *** TESTING ***
gene_sets <- as.matrix(gene_sets[,1:20])

# Load the results of the DE analysis, and create the summary data
# object, a data frame with two rows: (1) the LFC estimates, and (2)
# the variances of the LFC estimates.
load("../output/pbmc-purified/de-pbmc-purified-seed=1.RData")
zscores <- de$z[,k]
sdat    <- with(de,data.frame(log2FoldChange = postmean[,k],
                                   lfcSE2         = (postmean[,k]/z[,k])^2))
sdat <- subset(sdat,!is.na(lfcSE2))

# Align the gene-set data with the gene-wise statistics.
ids <- intersect(rownames(gene_sets),rownames(sdat))
i   <- match(ids,rownames(gene_sets))
j   <- match(ids,rownames(sdat))
gene_sets <- gene_sets[i,]
sdat      <- sdat[j,]
zscores   <- zscores[j]

# Next, remove gene sets with fewer than 4 genes, and with more than
# 400 genes. Gene sets with a large number of genes are less likely to
# be interesting, and slow down the enrichment analysis, so they are
# removed.
i <- which(colSums(gene_sets) >= 4 & colSums(gene_sets) <= 400)
gene_sets <- gene_sets[,i]

# Convert the sparse matrix representation of the gene sets to a list.
gene_sets <- matrix2list(gene_sets)

# Perform a gene set enrichment analysis using fgsea.
out <- fgsea(lapply(gene_sets,function (x) names(x)),zscores,
             eps = 1e-32,nproc = 2)
class(out)    <- "data.frame"
rownames(out) <- names(gene_sets)
out <- out[c("pval","padj","log2err","ES","NES","size")]
out <- out[names(gene_sets),]
out[is.na(out$log2err),] <- NA
out.fgsea <- out

# Perform a gene set enrichment analysis using iDEA.
idea <- new(Class      = "iDEA",
            gene_id    = rownames(sdat),
            annot_id   = names(gene_sets),
            summary    = sdat,
            annotation = gene_sets,
            num_gene   = nrow(sdat),
            num_core   = 2,
            project    = "idea")
idea <- iDEA.fit(idea,min_degene = 4,em_iter = 20,mcmc_iter = 1000, 
	         fit.tol = 1e-6,modelVariant = TRUE,verbose = TRUE)
idea <- iDEA.louis(idea)

# Save the results to file.
# TO DO.
