# Based on the topic modeling results in the droplet data, with k = 7
# topics, one of the topics, topic 6, appears to be capturing rare,
# specialized epithelial cells, including ionocyte, tuft and
# neuroendocrine cells. Therefore, here we reanalyze the cells with
# membership to topic 6.
#
# These were the steps taken to load R and allocate computing
# resources for this analysis:
#
#   sinteractive -p broadwl -c 8 --mem=16G --time=4:00:00
#   module load R/3.5.1
#

# Load a couple packages.
library(Matrix)
library(fastTopics)

# Initialize the sequence of pseudorandom numbers.
set.seed(1)

# Load the previously prepared count data.
load("../data/droplet.RData")

# Select cells that have at least 10% membership in topic 6.
fit     <- readRDS("../output/droplet/rds/fit-droplet-scd-ex-k=7.rds")$fit
fit     <- poisson2multinom(fit)
rows    <- which(fit$L[,6] > 0.1)
samples <- samples[rows,]
counts  <- counts[rows,]

# Remove genes that are not expressed in any of the cells.
j      <- which(colSums(counts > 0) >= 1)
counts <- counts[,j]
print(dim(counts))

# Now we are ready to perform the main model-fitting step. We repeat
# the model fitting for k = 2 through 6.
fits        <- vector("list",6)
names(fits) <- paste0("k",1:6)
for (k in 2:6)
  fits[[k]] <- fit_poisson_nmf(counts,k = k,numiter = 200,method = "scd",
                               init.method = "random",
                               control = list(extrapolate = TRUE,
                                              numiter = 4,nc = 8))
fits <- fits[-1]

# Perform a DE analysis using the topic model with k = 5 topics.
fit <- poisson2multinom(fits$k5)
t0 <- proc.time()
de <- de_analysis(fit,counts,pseudocount = 0.1,control = list(ns = 1e5,nc = 8))
t1 <- proc.time()
timing <- t1 - t0
cat(sprintf("Computation took %0.2f seconds.\n",timing["elapsed"]))

# Perform a second DE analysis with the k = 5 topic model after
# merging topics 3 and 4.
t0 <- proc.time()
fit_merged <- merge_topics(fit,c("k3","k4"))
de_merged <- de_analysis(fit_merged,counts,pseudocount = 0.1,
                         control = list(ns = 1e5,nc = 8))
t1 <- proc.time()
timing <- t1 - t0
cat(sprintf("Computation took %0.2f seconds.\n",timing["elapsed"]))

# SAVE RESULTS
# ------------
saveRDS(list(fits = fits,de = de,de_merged = de_merged),
        file = "refit-droplet.rds")
