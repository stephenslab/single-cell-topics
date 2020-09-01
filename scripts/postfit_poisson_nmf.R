# TO DO: Explain here what this script does, and how to use it. Also,
# give Slurm settings I used to run script, and combinations of script
# parameters settings I used.

# SCRIPT PARAMETERS
# -----------------
# TO DO: Explain here what these script parameters are for.
genesetfile  <- "../data/gene_sets_human.RData"
countsfile   <- "../data/pbmc_purified.RData"
modelfitfile <- "../output/pbmc-purified/rds/fit-pbmc-purified-scd-ex-k=6.rds"

# Load a few packages.
library(Matrix)
library(fastTopics)

# LOAD DATA
# ---------
# Load the gene-set data.
cat(sprintf("Loading gene-set data from %s.\n",genesetfile))
load(genesetfile)

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
timing <- system.time(diff_count_res <- diff_count_analysis(fit,counts))
cat(sprintf("Computation took %0.2f seconds.\n",timing["elapsed"]))
