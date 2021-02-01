library(Matrix)
library(fastTopics)
library(limma)
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
#
# IDEA: Perform differential expression analysis for B-cells
# vs. non-B-cells only.
#
celltype <- as.character(samples$celltype)
celltype[celltype == "CD4+/CD45RA+/CD25- Naive T" |
         celltype == "CD4+/CD45RO+ Memory" |
         celltype == "CD8+/CD45RA+ Naive Cytotoxic" |
         celltype == "CD4+ T Helper2" |
         celltype == "CD4+/CD25 T Reg"] <- "T cell"
celltype <- factor(celltype)
res <- diff_count_clusters(celltype,counts)

# Perform differential expression analysis using limma.
# TO DO.

# Perform differential expression analysis using limma/voom.
design <- model.matrix(~celltype)
colnames(design) <- letters[1:6]
v <- voom(t(counts),design)
v <- lmFit(v,design)
contrast.matrix <- makeContrasts("b - (a + c + d + e + f)",
                                 levels = design)
v <- contrasts.fit(v,contrast.matrix)
v <- eBayes(v)
topTable(v,number = 20)
plot(v$coefficients[,"celltypeCD19+ B"],res$beta[,"CD19+ B"],
     pch = 20,ylim = c(-10,10))
