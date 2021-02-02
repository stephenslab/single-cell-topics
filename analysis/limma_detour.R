library(Matrix)
library(fastTopics)
library(limma)
library(edgeR)
set.seed(1)

# Load the count data.
load("../data/pbmc_purified.RData")

# Take a random subset of the cells.
n       <- nrow(counts)
rows    <- sample(n,4000)
samples <- samples[rows,]
counts  <- counts[rows,]

# Perform differential expression analysis using fastTopics.
celltype <- factor(samples$celltype == "CD19+ B")
levels(celltype) <- c("A","B")
res <- diff_count_clusters(celltype,counts)

out <- ash(res$beta[,"B"],res$se[,"B"],mixcompdist = "normal")
plot(quantile(abs(res$Z[,"B"]),seq(0,1,0.01)),
     quantile(abs(out$result$PosteriorMean/out$result$PosteriorSD),
              seq(0,1,0.01)),
     pch = 20,xlim = c(-10,10),ylim = c(-10,10),
     xlab = "raw",ylab = "ashr")
abline(a = 0,b = 1,col = "skyblue",lty = "dotted")

n   <- 4000 - 2
out <- squeezeVar(res$se[,"B"]^2,df = n,robust = FALSE)
plot(res$se[,"B"],sqrt(out$var.post),pch = 20,log = "xy")

# Perform differential expression analysis using limma/voom.
design <- model.matrix(~0 + celltype)
contrast.matrix <- makeContrasts(AvsB = celltypeB - celltypeA,
                                 levels = colnames(design))
v <- voom(t(counts),design,plot = TRUE)
v <- lmFit(v,design)
v <- contrasts.fit(v,contrasts = contrast.matrix)
v <- eBayes(v)
plot(v$t,res$Z[,"B"],pch = 20,xlim = c(-25,25),ylim = c(-25,25))
abline(a = 0,b = 1,col = "skyblue",lty = "dotted")
