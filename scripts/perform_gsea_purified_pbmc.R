# A short script used to perform the gene set enrichment analysis
# using the GoM DE analysis results on the purified PBMC data, with 
# k = 6 topics.
library(Matrix)
library(tools)
library(susieR)
library(pathways)

# Initialize the sequence of pseudorandom numbers.
set.seed(1)

# Load the gene sets.
data(gene_sets_human)
X <- gene_sets_human$gene_sets

# Load the results of the topic-model-based DE analysis.
load("../output/pbmc-purified/de-pbmc-purified-seed=1.RData")
de1 <- de
load("../output/pbmc-purified/de-pbmc-purified-seed=2.RData")
de2 <- de
rm(de)

# When the two z-scores disagree, use the one that is nearer to zero.
de             <- de1[c("postmean","z")]
i              <- which(abs(de2$z) < abs(de1$z))
de$postmean[i] <- de2$postmean[i]
de$z[i]        <- de2$z[i]
Y              <- de$postmean

# Remove gene sets that aren't relevant.
dat <- gene_sets_human$gene_set_info
i <- which(with(dat,
                grepl("BioSystems-",database,fixed = TRUE) |
                grepl("PC-",database,fixed = TRUE) |
                (database == "MSigDB-C2" &
                 grepl("CP",sub_category_code,fixed = TRUE)) |
                (database == "MSigDB-C5" &
                 grepl("GO",sub_category_code,fixed = TRUE))))
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
save(list = c("X","Y","gsea"),
     file = "gsea-pbmc-purified-curated-only.RData")
resaveRdaFiles("gsea-pbmc-purified-curated-only.RData")
