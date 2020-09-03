# TO DO: Explain here what this script does, and how to use it. Also,
# give Slurm settings I used to run script, and combinations of script
# parameters settings I used.
library(tools)
library(Matrix)
library(fastTopics)
library(fgsea)
source("../code/gsea.R")

# SCRIPT PARAMETERS
# -----------------
# TO DO: Explain here what these script parameters are for.
genesetfile  <- "../data/gene_sets_human.RData"
countsfile   <- "../data/pbmc_purified.RData"
modelfitfile <- "../output/pbmc-purified/rds/fit-pbmc-purified-scd-ex-k=6.rds"
outfile      <- "postfit-pbmc-purified-scd-ex-k=6.rds"

# LOAD DATA
# ---------
# Load the gene-set data.
cat(sprintf("Loading gene-set data from %s.\n",genesetfile))
load(genesetfile)
rownames(gene_sets) <- gene_info$Ensembl
cat(sprintf("Loaded data for %d gene sets.\n",nrow(gene_set_info)))

# Load the previously prepared count data.
cat(sprintf("Loading data from %s.\n",countsfile))
load(countsfile)
cat(sprintf("Loaded %d x %d counts matrix.\n",nrow(counts),ncol(counts)))

# LOAD MODEL FIT
# --------------
cat(sprintf("Loading Poisson NMF model fit from %s\n",modelfitfile))
fit <- readRDS(modelfitfile)$fit

# COMPUTE Z-SCORES
# ----------------
# Perform differential expression analysis with the topic model.
cat("Computing differential expression statistics from topic model.\n")
timing <- system.time(diff_count_res <- diff_count_analysis(fit,counts))
cat(sprintf("Computation took %0.2f seconds.\n",timing["elapsed"]))

# PREPARE DATA FOR GSEA
# ---------------------
# Prepare the gene-set data and gene-wise statistics for the gene-set
# enrichment analysis. First, align the gene-set data with the
# gene-wise statistics.
out            <- align_gene_data(gene_sets,diff_count_res)
gene_sets      <- out$gene_sets
diff_count_res <- out$diff_count_res
rm(out)

# Next, remove gene sets with fewer than 4 genes, and with more than
# 400 genes. Gene sets with a large number of genes are less likely to
# be interesting, and slow down the enrichment analysis, so they are
# removed.
i <- which(colSums(gene_sets) >= 4 & colSums(gene_sets) <= 400)
gene_set_info <- gene_set_info[i,]
gene_sets     <- gene_sets[,i]
rm(i)

# PERFORM GSEA
# ------------
# For each topic, perform a gene-set enrichment analysis using fgsea.
cat("Performing gene-set enrichment analysis.\n")
rownames(gene_sets) <- rownames(diff_count_res$Z)
out <- perform_gsea_all_topics(gene_sets,diff_count_res$Z,nproc = 8)

# SAVE RESULTS
# ------------
cat("Saving results.\n")
# TO DO.
# resaveRdaFiles
