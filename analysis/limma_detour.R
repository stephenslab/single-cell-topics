library(Matrix)
library(fastTopics)
library(limma)
library(edgeR)
set.seed(1)

# Load the count data and K = 6 topic model fit.
load("../data/pbmc_purified.RData")
fit <- readRDS(file.path("../output/pbmc-purified/rds",
                         "fit-pbmc-purified-scd-ex-k=6.rds"))$fit
fit <- poisson2multinom(fit)

# Take a random subset of the cells.
n       <- nrow(counts)
rows    <- sample(n,4000)
samples <- samples[rows,]
counts  <- counts[rows,]
cols    <- which(colSums(counts > 0) >= 1)
counts  <- counts[,cols]
fit     <- select_loadings(fit,loadings = rows)

# Perform differential expression analysis using fastTopics.
celltype <- factor(samples$celltype == "CD19+ B")
levels(celltype) <- c("A","B")
res <- diff_count_clusters(celltype,counts)

# Perform differential expression analysis using limma.
dge    <- DGEList(t(counts))
logCPM <- cpm(dge,log = TRUE,prior.count = 3)
design <- model.matrix(~celltype)
eb     <- lmFit(logCPM,design)
eb     <- eBayes(eb,trend = TRUE)
topTable(eb,number = 20)

# Perform differential expression analysis using limma/voom.
design <- model.matrix(~0 + celltype)
v <- voom(t(counts),design)
v <- lmFit(v,design)
v <- eBayes(v)
plot(v$coefficients[,"celltypeCD19+ B"],res$beta[,"CD19+ B"],
     pch = 20,ylim = c(-10,10))
