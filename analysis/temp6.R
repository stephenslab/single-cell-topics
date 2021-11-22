library(DESeq2)
library(fastTopics)
library(ggplot2)
library(cowplot)
clamp <- function (x, a, b)
  pmax(a,pmin(b,x)) 
load("deseq2-pbmc-purified-bcells.RData")
load("de-pbmc-purified-seed=1.RData")
de1 <- de
load("de-pbmc-purified-seed=2.RData")
de2 <- de
k   <- 3
z   <- de1$z
i   <- which(abs(de2$z) < abs(de1$z))
z[i] <- de2$z[i]
z   <- z[abs(de1$est[,k]) > 0.25 & de1$f0 > 1e-5,]
i   <- intersect(rownames(deseq),rownames(z))
deseq$z <- with(deseq,log2FoldChange/lfcSE)
rownames(genes) <- genes$ensembl
pdat <- data.frame(deseq      = deseq[i,"z"],
                   fasttopics = z[i,k],
                   deseq.est  = deseq[i,"log2FoldChange"],
                   ft.est     = de1$postmean[i,k],
                   f0         = log10(de1$f0[i]),
                   postmean   = de1$postmean[i,k],
                   gene       = genes[i,"symbol"])
pdat <- transform(pdat,
                  deseq      = clamp(deseq,-200,+200),
                  fasttopics = clamp(fasttopics,-200,+200))
ggplot(pdat,aes(x = deseq,y = fasttopics,fill = f0)) +
  geom_point(shape = 21,color = "white") +
  geom_abline(intercept = 0,slope = 1,color = "magenta",linetype = "dotted") +
  scale_fill_gradient2(low = "deepskyblue",mid = "gold",high = "orangered",
                       na.value = "gainsboro",midpoint = -3) +
  theme_cowplot()
pdat2 <- subset(pdat,abs(fasttopics) > 4)
ggplot(pdat2,aes(x = deseq.est,y = ft.est,fill = f0)) +
  geom_point(shape = 21,color = "white") +
  geom_abline(intercept = 0,slope = 1,color = "magenta",linetype = "dotted") +
  scale_fill_gradient2(low = "deepskyblue",mid = "gold",high = "orangered",
                       na.value = "gainsboro",midpoint = -3) +
  theme_cowplot()
hist(clamp(deseq$z,-20,+20),n = 64)
hist(clamp(z,-20,+20),n = 64)
pdat2 <- subset(pdat,fasttopics > 4 & postmean > 4)
pdat2[order(pdat2$postmean,decreasing = TRUE),]
pdat2 <- subset(pdat,fasttopics > 4 & deseq < 4 & postmean > 1)
pdat2[order(pdat2$postmean,decreasing = TRUE),]
pdat2 <- subset(pdat,deseq > 50 & fasttopics < 4)
plot(de1$postmean,de2$postmean,pch = 20)
abline(a = 0,b = 1,col = "dodgerblue",lty = "dotted")
