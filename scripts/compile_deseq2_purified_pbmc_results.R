# TO DO: Explain here what this script is for, and how to use it.
library(tools)
library(DESeq2)
celltypes <- c("CD19+ B","CD56+ NK","T cell","CD14+ Monocyte","CD34+")
datfiles <- c("deseq2-pbmc-purified-bcells.RData",
              "deseq2-pbmc-purified-nkcells.RData",
              "deseq2-pbmc-purified-tcells.RData",
              "deseq2-pbmc-purified-cd14+.RData",
              "deseq2-pbmc-purified-cd34+.RData")
n <- length(celltypes)
res <- vector("list",n)
names(res) <- celltypes
for (i in 1:n) {
  load(datfiles[i])
  res[[i]] <- deseq
}
deseq <- res
save(list = c("genes","deseq"),file = "deseq2-pbmc-purified.RData")
resaveRdaFiles("deseq2-pbmc-purified.RData")
