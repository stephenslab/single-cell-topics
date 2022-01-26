# TO DO: Explain here what this script does, and how to use it.
#
# NOTES:
#
#   - Why does iDEA not represent gene sets as a dgCMatrix?
#
library(Matrix)
library(iDEA)

# Script parameters
k <- "k3"

# Load the gene sets.
load("../data/gene_sets_human.RData")
gene_sets <- gene_sets_human$gene_sets

# *** TESTING ***
gene_sets <- as.matrix(gene_sets[,1:20])

# Load the results of the DE analysis.
load("../output/pbmc-purified/de-pbmc-purified-seed=1.RData")
summary_data <- with(de,data.frame(log2FoldChange = postmean[,k],
                                   lfcSE2         = postmean[,k]/z[,k]))

# align the gene-set data with the gene-wise statistics.
x   <- rownames(gene_sets)
y   <- rownames(summary_data)
ids <- intersect(x,y)
i   <- match(ids,x)
j   <- match(ids,y)
gene_sets    <- gene_sets[i,]
summary_data <- summary_data[j,]

# Run the gene-set enrichment analysis.
idea <- CreateiDEAObject(summary_data,gene_sets,max_var_beta = Inf,
                         min_precent_annot = 0,num_core = 2)
idea <- iDEA.fit(idea,min_degene = 4,em_iter = 20,mcmc_iter = 1000, 
	         fit.tol = 1e-6,modelVariant = TRUE,verbose = TRUE)
