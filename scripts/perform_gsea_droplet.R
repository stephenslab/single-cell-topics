# A short script used to perform the gene set enrichment analysis
# using the GoM DE analysis results on the droplet data.
library(Matrix)
library(tools)
library(susieR)
library(pathways)

# Initialize the sequence of pseudorandom numbers.
set.seed(1)

# Load the gene sets.
data(gene_sets_mouse)
X <- gene_sets_mouse$gene_sets

# Load the results of the topic-model-based DE analysis.
load("../output/droplet/de-droplet.RData")

# Compile the DE results for all topics into a single matrix.
Y <- cbind(de_merged$postmean,de$postmean[,c("k5","k7")])

# Remove gene sets that aren't relevant.
dat <- gene_sets_mouse$gene_set_info
i <- which(with(dat,
                grepl("BioSystems-",database,fixed = TRUE) |
                (database == "MSigDB-C2" &
                 grepl("CP",sub_category_code,fixed = TRUE)) |
                (database == "MSigDB-C5" &
                 grepl("GO",sub_category_code,fixed = TRUE))))
X <- X[,i]

# Align the gene-set data with the gene-wise statistics.
rownames(X) <- gene_sets_mouse$gene_info$Symbol
genes <- intersect(rownames(X),rownames(Y))
X <- X[genes,]
Y <- Y[genes,]

# Next, remove gene sets with fewer than 10 genes and with more than
# 400 genes. Gene sets with a large number of genes are less likely to
# be interesting, and slow down the enrichment analysis, so they are
# removed.
i <- which(colSums(X) >= 10 & colSums(X) <= 400)
X <- X[,i]

# Perform a gene set enrichment analysis using susieR.
topics <- colnames(Y)
gsea <- vector("list",ncol(Y))
names(gsea) <- topics
t0 <- proc.time()
for (i in topics) {
  cat("topic",i,"\n")
  out <- susie(X,Y[,i],L = 10,intercept = TRUE,standardize = FALSE,
               estimate_residual_variance = TRUE,refine = FALSE,
               compute_univariate_zscore = FALSE,verbose = TRUE,
               min_abs_corr = 0)
  gsea[[i]] <- out[c("KL","lbf","sigma2","V","elbo","sets","pip","alpha","mu")]
}
t1 <- proc.time()
timing <- t1 - t0
cat(sprintf("Computation took %0.2f seconds.\n",timing["elapsed"]))

# Save the results to file.
save(list = c("X","Y","gsea"),file = "gsea-droplet-curated-only.RData")
resaveRdaFiles("gsea-droplet-curated-only.RData")
