# A short script used to perform the gene set enrichment analysis
# using the DE analysis results on the purified PBMC data, with k = 6
# topics. These were the steps taken to load R and allocate computing
# resources for this analysis:
#
#   sinteractive -p broadwl -c 20 --mem=36G --time=60:00:00
#   module load R/4.1.0
#
library(Matrix)
library(fgsea)
library(iDEA)
library(tools)
source("../code/gsea.R")

# Script parameters.
k       <- "k3"
outfile <- "gsea-pbmc-purified-k=3.RData"
print(k)
print(outfile)
set.seed(1)

# Load the gene sets.
load("../data/gene_sets_human.RData")
gene_sets <- gene_sets_human$gene_sets

# Load the results of the DE analysis, and create the summary data
# object, a data frame with two rows: (1) the LFC estimates, and (2)
# the variances of the LFC estimates.
load("../output/pbmc-purified/de-pbmc-purified-seed=1.RData")
zscores <- de$z[,k]
sdat    <- with(de,data.frame(log2FoldChange = postmean[,k],
                              lfcSE2         = (postmean[,k]/z[,k])^2))
i       <- which(!(is.na(zscores) | is.na(sdat$lfcSE2)))
zscores <- zscores[i]
sdat    <- sdat[i,]

# Align the gene-set data with the gene-wise statistics.
ids <- intersect(rownames(gene_sets),rownames(sdat))
i   <- match(ids,rownames(gene_sets))
j   <- match(ids,rownames(sdat))
gene_sets <- gene_sets[i,]
sdat      <- sdat[j,]
zscores   <- zscores[j]

# Remove gene sets in several collections that aren't relevant.
i <- which(!is.element(gene_sets_human$gene_set_info$database,
                       c("MSigDB-ARCHIVED","MSigDB-C1","MSigDB-C3",
                         "MSigDB-C4","MSigDB-C6")))
gene_sets <- gene_sets[,i]

# Next, remove gene sets with fewer than 10 genes, and with more than
# 400 genes. Gene sets with a large number of genes are less likely to
# be interesting, and slow down the enrichment analysis, so they are
# removed.
i <- which(colSums(gene_sets) >= 10 & colSums(gene_sets) <= 400)
gene_sets <- gene_sets[,i]

# Convert the sparse matrix representation of the gene sets to a list.
gene_sets <- matrix2list(gene_sets)

# Perform a gene set enrichment analysis using fgsea.
t0 <- proc.time()
out <- fgsea(lapply(gene_sets,function (x) names(x)),zscores,
             eps = 1e-32,nproc = 20)
class(out)    <- "data.frame"
rownames(out) <- names(gene_sets)
out <- out[c("pval","padj","log2err","ES","NES","size")]
out <- out[names(gene_sets),]
out[is.na(out$log2err),] <- NA
out.fgsea <- out
t1 <- proc.time()
timing <- t1 - t0
cat(sprintf("Computation took %0.2f seconds.\n",timing["elapsed"]))

stop()

# Perform a gene set enrichment analysis using iDEA.
t0 <- proc.time()
idea <- new(Class      = "iDEA",
            gene_id    = rownames(sdat),
            annot_id   = names(gene_sets),
            summary    = sdat,
            annotation = gene_sets,
            num_gene   = nrow(sdat),
            num_core   = 20,
            project    = "idea")
# Error in names(res_idea) <- object@annot_id :
#   attempt to set an attribute on NULL
set.seed(1)
idea1 <- iDEA.fit(idea,min_degene = 4,em_iter = 10,mcmc_iter = 100, 
 	          fit.tol = 1e-6,modelVariant = TRUE,verbose = TRUE)
set.seed(1)
idea2 <- iDEA.fit(idea,min_degene = 4,em_iter = 10,mcmc_iter = 100, 
 	          fit.tol = 1e-6,modelVariant = TRUE,verbose = TRUE)
idea <- iDEA.louis(idea)
t1 <- proc.time()
timing <- t1 - t0
cat(sprintf("Computation took %0.2f seconds.\n",timing["elapsed"]))

# Save the results to file.
idea_gsea <- idea@gsea
save(list = c("out.fgsea","idea_gsea"),file = outfile)
resaveRdaFiles(outfile)
